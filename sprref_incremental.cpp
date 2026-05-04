/*
    sprref_incremental.cpp

    M1 scaffolding for the incremental sparse RREF C-API. This compiles into
    a shared library (libsparse_rref_inc.{dylib,so}) but the algorithmic core
    is stubbed out and will be implemented in M2/M3.

    The state layout is the minimum the M2/M3 work will need:
      - field setup (prime / Flint field_t)
      - per-handle thread pool (kept as ptr for ABI compat across translation
        units; constructed in init, destroyed in free)
      - pivot ordering and master set
      - the RREF basis: a map from pivot column -> normalised sparse row
      - reverse index col -> set of pivot-rows containing that column (for
        cheap forward elimination of new rows)
      - last_insert_changed_basis flag (controls is_solved cache invalidation)

    Concurrent calls on the same handle are NOT supported. Callers (Python
    wrappers, mp.Pool workers) should each own their own handle.
*/

#include "sprref_incremental.h"

#include "sparse_mat.h"
#include "sparse_rref.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using SparseRREF::field_t;
using SparseRREF::RING;

namespace {

using Row = std::unordered_map<uint32_t, uint64_t>; /* col -> coeff (mod p) */

struct PivotRow {
    /* Normalised so pivot coefficient == 1. */
    Row coeffs;     /* off-pivot columns only (pivot column omitted) */
    uint64_t rhs;
};

} /* namespace */

struct sprref_inc {
    uint64_t prime;
    field_t F;
    int n_threads;

    uint32_t nvars = 0;

    /* pivot column -> normalised row */
    std::unordered_map<uint32_t, PivotRow> basis;

    /* reverse index: column -> set of pivot columns whose row touches that column */
    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> col_to_pivots;

    /* pivot ordering: smaller key preferred. Columns not listed fall back to
       sort by numeric column id. */
    std::unordered_map<uint32_t, uint64_t> pivot_keys;

    std::unordered_set<uint32_t> masters;

    bool dirty_since_last_solve = false;
};

/* =========================================================================
   Helpers (M1: minimal — full impl in M2/M3)
   ========================================================================= */

namespace {

inline uint64_t mod_reduce(uint64_t x, uint64_t p) {
    return x % p;
}

/* Stub: forward eliminate `row` against existing basis pivots. No-op for M1. */
void forward_eliminate(sprref_inc_t* /*h*/, Row& /*row*/, uint64_t& /*rhs*/) {
    /* TODO(M2): for each pivot column c present in row, saxpy row -= row[c] * basis[c]. */
}

/* Stub: choose pivot column from a reduced row. Returns false if row empty.
   M1 picks the lowest-key non-master column; full Laporta logic comes in M2. */
bool choose_pivot(const sprref_inc_t* h, const Row& row, uint32_t& out_pivot) {
    if (row.empty()) return false;

    bool found = false;
    uint32_t best_col = 0;
    uint64_t best_key = 0;
    auto key_of = [&](uint32_t c) -> uint64_t {
        auto it = h->pivot_keys.find(c);
        return it != h->pivot_keys.end() ? it->second : (uint64_t)c;
    };

    for (auto& kv : row) {
        uint32_t c = kv.first;
        if (h->masters.count(c)) continue;
        uint64_t k = key_of(c);
        if (!found || k < best_key) {
            found = true;
            best_key = k;
            best_col = c;
        }
    }
    if (found) out_pivot = best_col;
    return found;
}

} /* namespace */

/* =========================================================================
   API
   ========================================================================= */

extern "C" {

const char* sprref_inc_version(void) {
    return "sprref_incremental v0.0.1 (M1 scaffolding)";
}

sprref_inc_t* sprref_inc_init(uint64_t field_order, int n_threads) {
    if (field_order < 3) return nullptr;
    if (n_threads < 1) n_threads = 1;

    auto* h = new (std::nothrow) sprref_inc_t();
    if (!h) return nullptr;

    h->prime = field_order;
    h->n_threads = n_threads;
    /* Flint Fp field setup. */
    h->F = field_t(RING::FIELD_Fp, field_order);
    h->nvars = 0;
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

int sprref_inc_insert(sprref_inc_t* h,
                      const uint32_t* cols,
                      const uint64_t* vals,
                      size_t nnz,
                      uint64_t rhs) {
    if (!h) return SPRREF_INC_DEPENDENT;
    const uint64_t p = h->prime;

    Row row;
    row.reserve(nnz);
    for (size_t i = 0; i < nnz; ++i) {
        uint64_t v = mod_reduce(vals[i], p);
        if (v != 0) row[cols[i]] = v;
    }
    rhs = mod_reduce(rhs, p);

    forward_eliminate(h, row, rhs);

    if (row.empty()) {
        if (rhs == 0) return SPRREF_INC_DEPENDENT;
        return SPRREF_INC_INCONSISTENT;
    }

    uint32_t pivot_col = 0;
    if (!choose_pivot(h, row, pivot_col)) {
        /* All remaining columns are masters: treat as dependent (cannot pivot). */
        return SPRREF_INC_DEPENDENT;
    }

    /* M1: skip normalisation and basis insertion to keep this a pure stub.
       M2 will normalise (multiply by inverse of row[pivot_col]) and store. */
    (void)pivot_col;
    h->dirty_since_last_solve = true;
    return SPRREF_INC_INDEPENDENT;
}

int sprref_inc_is_solved(sprref_inc_t* h,
                         uint32_t var_idx,
                         const uint32_t* /*acceptable_free*/,
                         size_t /*n_free*/,
                         uint32_t** /*out_cols*/,
                         uint64_t** /*out_vals*/,
                         size_t* /*out_nnz*/,
                         uint64_t* /*out_rhs*/) {
    if (!h) return SPRREF_INC_NOT_SOLVED;
    /* M1 stub: nothing is ever solved. M3 will add lazy backsub + extraction. */
    (void)var_idx;
    return SPRREF_INC_NOT_SOLVED;
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
