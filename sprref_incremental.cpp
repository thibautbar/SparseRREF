/*
    sprref_incremental.cpp

    Incremental sparse RREF over GF(p) with eager backward substitution.

    State model:
      - basis: pivot column -> normalised row (off-pivot entries) + rhs
      - col_to_pivots: column -> set of pivot columns whose row touches it
        (reverse index, used to find rows needing eager backsub when a new
         pivot column appears)
      - pivot_keys: optional Laporta ordering for pivot selection
      - masters: columns that may never be pivoted

    On every successful insert the basis is left in true RREF (each pivot
    column appears with coefficient 0 in all other pivot rows). This matches
    SpotlightSolverSparseGF semantics and makes is_solved a near-O(1) lookup.

    Concurrency: a handle is NOT thread-safe. Use one handle per
    mp.Pool worker.
*/

#include "sprref_incremental.h"

#include <flint/ulong_extras.h>

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

using Row = std::unordered_map<uint32_t, uint64_t>;

struct PivotRow {
    Row coeffs;     /* off-pivot columns only (the pivot column itself is implicit, with value 1) */
    uint64_t rhs;
};

inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t p) {
    /* Use 128-bit intermediate to avoid overflow for primes up to 2^63. */
    return (uint64_t)((__uint128_t)a * b % p);
}

inline uint64_t submod(uint64_t a, uint64_t b, uint64_t p) {
    return (a >= b) ? (a - b) : (a + p - b);
}

} /* namespace */

struct sprref_inc {
    uint64_t prime = 0;
    uint32_t nvars = 0;

    std::unordered_map<uint32_t, PivotRow> basis;
    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> col_to_pivots;
    std::unordered_map<uint32_t, uint64_t> pivot_keys;
    std::unordered_set<uint32_t> masters;
};

/* =========================================================================
   Internal helpers
   ========================================================================= */

namespace {

/* col_to_pivots[c] = { pivot rows whose stored .coeffs (off-pivot part) contains c }.
   Used by eager backsub to find rows that need their entry at a new pivot column
   zeroed out. The pivot column itself is NOT inserted into col_to_pivots[pivot]
   (its presence in basis is checked directly via h->basis.count). */
void index_row(sprref_inc_t* h, uint32_t pivot, const Row& row) {
    for (auto& [c, v] : row) {
        (void)v;
        h->col_to_pivots[c].insert(pivot);
    }
}

void deindex_row(sprref_inc_t* h, uint32_t pivot, const Row& row) {
    for (auto& [c, v] : row) {
        (void)v;
        auto it = h->col_to_pivots.find(c);
        if (it == h->col_to_pivots.end()) continue;
        it->second.erase(pivot);
        if (it->second.empty()) h->col_to_pivots.erase(it);
    }
}

/* dst -= factor * src ; arithmetic mod p. */
void saxpy(Row& dst_row, uint64_t& dst_rhs,
           uint64_t factor,
           const Row& src_row, uint64_t src_rhs,
           uint64_t p) {
    if (factor == 0) return;
    for (auto& [c, v] : src_row) {
        uint64_t prod = mulmod(factor, v, p);
        if (prod == 0) continue;
        auto it = dst_row.find(c);
        if (it == dst_row.end()) {
            dst_row.emplace(c, p - prod); /* 0 - prod */
        } else {
            uint64_t nv = submod(it->second, prod, p);
            if (nv == 0) dst_row.erase(it);
            else it->second = nv;
        }
    }
    if (src_rhs != 0) {
        uint64_t prod = mulmod(factor, src_rhs, p);
        dst_rhs = submod(dst_rhs, prod, p);
    }
}

/* Forward elimination: while `row` contains a pivot column, saxpy against
   basis. Loops until fixed point (each saxpy may introduce columns from
   src that are themselves pivots — though in true RREF this never happens,
   we don't rely on that invariant here so the impl is robust to lazy state). */
void forward_eliminate(sprref_inc_t* h, Row& row, uint64_t& rhs) {
    const uint64_t p = h->prime;
    while (true) {
        uint32_t pivot = 0;
        bool found = false;
        for (auto& [c, v] : row) {
            (void)v;
            if (h->basis.count(c)) {
                pivot = c;
                found = true;
                break;
            }
        }
        if (!found) return;

        uint64_t factor = row[pivot];
        /* The pivot is normalised to 1, so eliminating it means subtracting
           factor * (basis_row + e_pivot). The pivot column itself is implicit
           in the stored basis row (coeff 1). */
        row.erase(pivot);  /* equivalent to subtracting factor*1 from row[pivot] */
        const PivotRow& src = h->basis.at(pivot);
        saxpy(row, rhs, factor, src.coeffs, src.rhs, p);
    }
}

bool choose_pivot(const sprref_inc_t* h,
                  const Row& row,
                  uint32_t preferred_pivot,
                  uint32_t& out_pivot) {
    if (row.empty()) return false;

    /* Spotlight semantics: preferred_pivot wins iff present and non-master. */
    if (preferred_pivot != SPRREF_INC_NO_PREF) {
        auto it = row.find(preferred_pivot);
        if (it != row.end() && !h->masters.count(preferred_pivot)) {
            out_pivot = preferred_pivot;
            return true;
        }
    }

    bool found = false;
    uint32_t best_col = 0;
    uint64_t best_key = 0;
    auto key_of = [&](uint32_t c) -> uint64_t {
        auto it = h->pivot_keys.find(c);
        return it != h->pivot_keys.end() ? it->second : (uint64_t)c;
    };

    for (auto& [c, v] : row) {
        (void)v;
        if (h->masters.count(c)) continue;
        uint64_t k = key_of(c);
        /* Tie-break by column index for determinism. */
        if (!found || k < best_key || (k == best_key && c < best_col)) {
            found = true;
            best_key = k;
            best_col = c;
        }
    }
    if (found) out_pivot = best_col;
    return found;
}

/* Multiply row+rhs by inverse of row[pivot], then erase the pivot entry
   (it's stored implicitly as coefficient 1). */
void normalize_around_pivot(Row& row, uint64_t& rhs, uint32_t pivot, uint64_t p) {
    auto it = row.find(pivot);
    if (it == row.end()) return; /* shouldn't happen: caller chose this pivot */
    uint64_t pv = it->second;
    if (pv == 1) {
        row.erase(it);
        return;
    }
    uint64_t inv = n_invmod(pv, p);
    /* In-place scale (skip pivot — we'll erase it). */
    for (auto& kv : row) {
        if (kv.first == pivot) continue;
        kv.second = mulmod(kv.second, inv, p);
    }
    rhs = mulmod(rhs, inv, p);
    row.erase(pivot);
}

} /* anonymous namespace */

/* =========================================================================
   API
   ========================================================================= */

extern "C" {

const char* sprref_inc_version(void) {
    return "sprref_incremental v0.3.0 (delta pivot keys)";
}

sprref_inc_t* sprref_inc_init(uint64_t field_order, int n_threads) {
    if (field_order < 3) return nullptr;
    if (n_threads < 1) n_threads = 1;
    auto* h = new (std::nothrow) sprref_inc_t();
    if (!h) return nullptr;
    h->prime = field_order;
    h->nvars = 0;
    (void)n_threads; /* M2 is single-threaded; M5 may parallelise saxpy. */
    return h;
}

void sprref_inc_free(sprref_inc_t* h) {
    delete h;
}

void sprref_inc_resize_nvars(sprref_inc_t* h, uint32_t new_nvars) {
    if (!h) return;
    if (new_nvars > h->nvars) h->nvars = new_nvars;
}

void sprref_inc_set_masters(sprref_inc_t* h, const uint32_t* master_cols, size_t n) {
    if (!h) return;
    h->masters.clear();
    h->masters.reserve(n);
    for (size_t i = 0; i < n; ++i) h->masters.insert(master_cols[i]);
}

void sprref_inc_set_pivot_order(sprref_inc_t* h,
                                const uint32_t* cols,
                                const uint64_t* keys,
                                size_t n) {
    if (!h) return;
    h->pivot_keys.clear();
    h->pivot_keys.reserve(n);
    for (size_t i = 0; i < n; ++i) h->pivot_keys[cols[i]] = keys[i];
}

void sprref_inc_add_pivot_keys(sprref_inc_t* h,
                               const uint32_t* cols,
                               const uint64_t* keys,
                               size_t n) {
    if (!h) return;
    h->pivot_keys.reserve(h->pivot_keys.size() + n);
    for (size_t i = 0; i < n; ++i) h->pivot_keys[cols[i]] = keys[i];
}

int sprref_inc_insert(sprref_inc_t* h,
                      const uint32_t* cols,
                      const uint64_t* vals,
                      size_t nnz,
                      uint64_t rhs,
                      uint32_t preferred_pivot) {
    if (!h) return SPRREF_INC_DEPENDENT;
    const uint64_t p = h->prime;

    Row row;
    row.reserve(nnz);
    for (size_t i = 0; i < nnz; ++i) {
        uint64_t v = vals[i] % p;
        if (v != 0) {
            uint32_t c = cols[i];
            /* Combine duplicates if any (mod p sum). */
            auto it = row.find(c);
            if (it == row.end()) row.emplace(c, v);
            else {
                uint64_t sum = it->second + v;
                if (sum >= p) sum -= p;
                if (sum == 0) row.erase(it);
                else it->second = sum;
            }
            if (c >= h->nvars) h->nvars = c + 1;
        }
    }
    rhs %= p;

    forward_eliminate(h, row, rhs);

    if (row.empty()) {
        return rhs == 0 ? SPRREF_INC_DEPENDENT : SPRREF_INC_INCONSISTENT;
    }

    uint32_t pivot_col = 0;
    if (!choose_pivot(h, row, preferred_pivot, pivot_col)) {
        /* All remaining entries are masters: no pivotable column. */
        return SPRREF_INC_DEPENDENT;
    }

    normalize_around_pivot(row, rhs, pivot_col, p);

    /* Eager backward substitution: zero out pivot_col in every existing
       basis row that touches it. We use the reverse index, then re-index
       each modified row since saxpy may add/remove columns. */
    auto rev_it = h->col_to_pivots.find(pivot_col);
    if (rev_it != h->col_to_pivots.end()) {
        /* Snapshot: saxpy will mutate col_to_pivots, so iterate over a copy. */
        std::vector<uint32_t> affected(rev_it->second.begin(), rev_it->second.end());
        for (uint32_t other_pivot : affected) {
            if (other_pivot == pivot_col) continue; /* not in basis yet, but defensive */
            auto bit = h->basis.find(other_pivot);
            if (bit == h->basis.end()) continue;
            PivotRow& other = bit->second;
            auto cit = other.coeffs.find(pivot_col);
            if (cit == other.coeffs.end() || cit->second == 0) continue;
            uint64_t factor = cit->second;
            /* Deindex, saxpy against (row, rhs) [which has implicit pivot 1 — but the
               new pivot row's stored coeffs do NOT include pivot_col since we erased it,
               and the implicit "1 at pivot_col" cancels exactly with `factor` at that column].
               After saxpy + erasing pivot_col entry, the row no longer touches pivot_col. */
            deindex_row(h, other_pivot, other.coeffs);
            other.coeffs.erase(cit); /* implicit "1" cancels factor exactly */
            saxpy(other.coeffs, other.rhs, factor, row, rhs, p);
            index_row(h, other_pivot, other.coeffs);
        }
    }

    /* Insert new pivot. */
    h->basis.emplace(pivot_col, PivotRow{std::move(row), rhs});
    index_row(h, pivot_col, h->basis[pivot_col].coeffs);

    if (pivot_col >= h->nvars) h->nvars = pivot_col + 1;

    return SPRREF_INC_INDEPENDENT;
}

int sprref_inc_is_solved(sprref_inc_t* h,
                         uint32_t var_idx,
                         const uint32_t* acceptable_free,
                         size_t n_free,
                         uint32_t** out_cols,
                         uint64_t** out_vals,
                         size_t* out_nnz,
                         uint64_t* out_rhs) {
    if (!h) return SPRREF_INC_NOT_PIVOT;

    auto it = h->basis.find(var_idx);
    if (it == h->basis.end()) return SPRREF_INC_NOT_PIVOT;

    const PivotRow& pr = it->second;

    /* Decide solved-ness: every off-pivot column in the stored row must
       either be a master or be in acceptable_free. (Because the basis is
       maintained in true RREF, no off-pivot column can itself be a pivot —
       forward elim + eager backsub guarantee this.) */
    std::unordered_set<uint32_t> free_set;
    if (acceptable_free && n_free > 0) {
        free_set.reserve(n_free);
        for (size_t i = 0; i < n_free; ++i) free_set.insert(acceptable_free[i]);
    }

    bool has_free = false;
    for (auto& [c, v] : pr.coeffs) {
        (void)v;
        if (h->masters.count(c)) continue;
        if (!free_set.empty() && free_set.count(c)) continue;
        has_free = true;
        break;
    }

    /* Still allocate and emit the expression even when has_free is true:
       the caller (Python wrapper) may want the partial expression for
       intermediate logs / pretty-printing. The bool return signals
       solved-ness; the buffers carry the coefficients regardless. */
    const size_t nnz = pr.coeffs.size();
    if (out_nnz) *out_nnz = nnz;
    if (out_rhs) *out_rhs = pr.rhs;

    if (out_cols && out_vals) {
        if (nnz == 0) {
            *out_cols = nullptr;
            *out_vals = nullptr;
        } else {
            auto* cols_buf = (uint32_t*)std::malloc(sizeof(uint32_t) * nnz);
            auto* vals_buf = (uint64_t*)std::malloc(sizeof(uint64_t) * nnz);
            if (!cols_buf || !vals_buf) {
                std::free(cols_buf);
                std::free(vals_buf);
                *out_cols = nullptr;
                *out_vals = nullptr;
                if (out_nnz) *out_nnz = 0;
                return SPRREF_INC_NOT_PIVOT; /* allocation failure: degrade safely */
            }
            size_t i = 0;
            for (auto& [c, v] : pr.coeffs) {
                cols_buf[i] = c;
                vals_buf[i] = v;
                ++i;
            }
            *out_cols = cols_buf;
            *out_vals = vals_buf;
        }
    }

    return has_free ? SPRREF_INC_HAS_FREE : SPRREF_INC_SOLVED;
}

void sprref_inc_buffer_free(uint32_t* cols, uint64_t* vals) {
    std::free(cols);
    std::free(vals);
}

size_t sprref_inc_rank(const sprref_inc_t* h) {
    if (!h) return 0;
    return h->basis.size();
}

uint32_t sprref_inc_nvars(const sprref_inc_t* h) {
    if (!h) return 0;
    return h->nvars;
}

} /* extern "C" */
