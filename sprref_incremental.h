/*
    sprref_incremental.h

    C-callable API for incremental sparse RREF over GF(p).

    Unlike the one-shot batch RREF in main.cpp, this API exposes a persistent
    handle that supports row-by-row insertion (with on-the-fly forward
    elimination) and lazy is_solved queries with backward substitution.

    Designed to be loaded as a shared library by Python (ctypes/cffi) for the
    SparseRREFBatchSolver "incremental" mode in project-feynman.

    All operations are stateless w.r.t. globals: state lives entirely in the
    opaque sprref_inc_t* handle. Multiple independent handles can coexist
    (e.g. one per worker process).

    M1 milestone: only the API surface is defined. Implementation is stubbed.
    M2/M3 will fill in forward elimination and lazy backsub.
*/

#ifndef SPRREF_INCREMENTAL_H
#define SPRREF_INCREMENTAL_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) || defined(__CYGWIN__)
  #define SPRREF_API __declspec(dllexport)
#else
  #define SPRREF_API __attribute__((visibility("default")))
#endif

typedef struct sprref_inc sprref_inc_t;

/* Status codes returned by sprref_inc_insert. */
#define SPRREF_INC_INDEPENDENT 0
#define SPRREF_INC_DEPENDENT   1
#define SPRREF_INC_INCONSISTENT 2

/* Returns for sprref_inc_is_solved (three-way: distinguishes "var_idx is
   not a pivot at all" from "is a pivot but has unaccounted free columns"). */
#define SPRREF_INC_NOT_PIVOT  0  /* var_idx has no pivot row in the basis */
#define SPRREF_INC_HAS_FREE   1  /* is a pivot, but some non-master/non-free col remains */
#define SPRREF_INC_SOLVED     2  /* is a pivot AND fully reduced (all off-pivot cols accounted for) */

/* Extra status code for sprref_inc_insert_forward. */
#define SPRREF_INC_NEEDS_PIVOT 3  /* row reduced to non-empty form; caller must pick a pivot column */

/* Lifecycle ----------------------------------------------------------------- */

/*
   Create a new incremental solver over GF(field_order).
   field_order must be an odd prime (> 2).
   n_threads: thread pool size for internal parallelism (>= 1).
   Returns NULL on failure (e.g. invalid prime).
*/
SPRREF_API sprref_inc_t* sprref_inc_init(uint64_t field_order, int n_threads);

SPRREF_API void sprref_inc_free(sprref_inc_t* h);

SPRREF_API const char* sprref_inc_version(void);

/* Configuration ------------------------------------------------------------- */

SPRREF_API void sprref_inc_resize_nvars(sprref_inc_t* h, uint32_t new_nvars);

/* Replaces the master-column set. Master columns are never chosen as pivots. */
SPRREF_API void sprref_inc_set_masters(sprref_inc_t* h,
                                       const uint32_t* master_cols,
                                       size_t n);

/*
   Replaces the pivot ordering. cols[i] gets sort-key keys[i].
   Smaller keys are preferred when picking a pivot from a row's nonzero set.
   Columns not listed here use the fallback (numeric column index).
*/
SPRREF_API void sprref_inc_set_pivot_order(sprref_inc_t* h,
                                           const uint32_t* cols,
                                           const uint64_t* keys,
                                           size_t n);

/*
   Upsert version of set_pivot_order: adds (cols[i] -> keys[i]) entries
   without clearing existing ones. Existing keys for the same column are
   overwritten. This is the hot path for monotonically growing pivot_order
   maps (e.g. as new integrals are inserted into the solver) — callers can
   push only the delta on each step instead of the full map.
*/
SPRREF_API void sprref_inc_add_pivot_keys(sprref_inc_t* h,
                                          const uint32_t* cols,
                                          const uint64_t* keys,
                                          size_t n);

/*
   Two-component variant of add_pivot_keys: each column gets a *pair* of
   uint64 keys, compared lexicographically (key0 first, key1 as tiebreaker).
   This is necessary when the natural pivot priority is a tuple that doesn't
   fit in 64 bits (e.g. a Laporta key whose first four scalar components fit
   in key0 but whose tail components need additional bits in key1).

   Calling this with a column that already has a single-uint64 key (set via
   set_pivot_order or add_pivot_keys) overwrites both halves; mixing the two
   forms on different columns is allowed (a column with no key1 is treated
   as key1==0 for comparison).
*/
SPRREF_API void sprref_inc_add_pivot_keys2(sprref_inc_t* h,
                                           const uint32_t* cols,
                                           const uint64_t* key0s,
                                           const uint64_t* key1s,
                                           size_t n);

/* Two-phase insert ---------------------------------------------------------- */

/*
   Phase 1: forward eliminate the incoming row against existing pivots.
   Stashes the reduced row in the handle's pending slot.

   Returns one of:
     SPRREF_INC_DEPENDENT    - row reduced to 0 = 0; nothing pending. out_* untouched.
     SPRREF_INC_INCONSISTENT - row reduced to 0 = nonzero; nothing pending. out_* untouched.
     SPRREF_INC_NEEDS_PIVOT  - row pivotable. out_candidate_cols / out_n_candidates are
                               populated with the column indices of the reduced row's
                               non-zero entries. The caller must pick one of them
                               (typically using its own pivot ordering, applying
                               master / acceptable_free filtering Python-side) and
                               then call sprref_inc_commit_pivot.
                               out_candidate_cols is heap-allocated; free via
                               sprref_inc_buffer_free_u32.

   Designed to keep the actual pivot ordering decision in the calling
   language (Python in our case), which matches SpotlightSolverSparseGF
   semantics bit-for-bit without needing to translate Python comparison
   keys into the C side.

   At most one pending row may exist on a handle at a time. Calling this
   while a pending row exists overwrites the pending slot.
*/
SPRREF_API int sprref_inc_insert_forward(
    sprref_inc_t* h,
    const uint32_t* cols,
    const uint64_t* vals,
    size_t nnz,
    uint64_t rhs,
    uint32_t** out_candidate_cols,
    size_t* out_n_candidates);

/*
   Phase 2 (success): commit the pending row using pivot_col as the pivot.
   Performs normalisation (so the pivot coefficient becomes 1) and eager
   backward substitution against existing basis rows (zeroing out pivot_col
   in every row that touches it). After this the basis is in true RREF.

   pivot_col MUST be one of the candidate columns returned by the matching
   sprref_inc_insert_forward call; passing a column not in the pending row
   is undefined behaviour.
*/
SPRREF_API void sprref_inc_commit_pivot(sprref_inc_t* h, uint32_t pivot_col);

/*
   Phase 2 (cancel): discard the pending row without modifying the basis.
   Use when the caller decides the row can't be pivoted (e.g. all candidate
   columns are masters): treat it as DEPENDENT instead.
*/
SPRREF_API void sprref_inc_abort_pending(sprref_inc_t* h);

/* Free a uint32 buffer allocated by sprref_inc_insert_forward. */
SPRREF_API void sprref_inc_buffer_free_u32(uint32_t* ptr);

/* Single-call insert (legacy / direct C-API users) -------------------------- */

/* Sentinel for "no preferred pivot" in sprref_inc_insert. */
#define SPRREF_INC_NO_PREF UINT32_MAX

/*
   Insert one sparse augmented row [a | b] where a has nnz entries
   (cols[i], vals[i]) and b is rhs. All values are reduced modulo field_order
   on entry (pass any nonneg representative; will be canonicalised).

   preferred_pivot: if not SPRREF_INC_NO_PREF and that column is present in
     the (forward-reduced) row and is not a master, it is used as the pivot.
     Otherwise the lowest-pivot-key non-master column wins (with column index
     as fallback tie-breaker).

   The row is reduced against existing pivots (forward elimination). If the
   reduced row is empty:
     - rhs == 0  -> SPRREF_INC_DEPENDENT
     - rhs != 0  -> SPRREF_INC_INCONSISTENT
   Otherwise a new pivot is chosen, the row is normalised so the pivot
   coefficient is 1, eager backward substitution removes the new pivot column
   from all existing basis rows, and SPRREF_INC_INDEPENDENT is returned.

   Post-condition (on independent return): the basis is in true reduced row
   echelon form — each pivot column has 0 in all other pivot rows.
*/
SPRREF_API int sprref_inc_insert(sprref_inc_t* h,
                                 const uint32_t* cols,
                                 const uint64_t* vals,
                                 size_t nnz,
                                 uint64_t rhs,
                                 uint32_t preferred_pivot);

/* Query --------------------------------------------------------------------- */

/*
   Check if var_idx is solved. acceptable_free (length n_free) is an optional
   set of column indices that may appear free in the expression without
   counting as "unsolved". Pass NULL/0 to require strict full solving (only
   masters allowed as free).

   Output buffers (regardless of solved-ness, when var_idx is a pivot):
     - *out_cols, *out_vals: heap-allocated arrays of length *out_nnz holding
       the (col, coeff) pairs of the stored basis row's off-pivot part.
       Caller must free via sprref_inc_buffer_free. May be NULL when nnz==0.
     - *out_rhs: constant term of the stored row.

   Semantics: the basis row encodes x_{var_idx} + sum_j coeff_j * x_j = rhs,
   so the caller derives x_{var_idx} = rhs - sum_j coeff_j * x_j. Coefficients
   are returned as STORED (positive mod p), no negation applied — the Python
   wrapper applies (-coeff) % p to match the existing solver-result convention.

   If var_idx is not a pivot at all, returns SPRREF_INC_NOT_SOLVED and leaves
   out_* untouched. If var_idx is a pivot but has unaccounted free columns,
   the buffers are still populated (useful for intermediate logs) but the
   return is SPRREF_INC_NOT_SOLVED.
*/
SPRREF_API int sprref_inc_is_solved(sprref_inc_t* h,
                                    uint32_t var_idx,
                                    const uint32_t* acceptable_free,
                                    size_t n_free,
                                    uint32_t** out_cols,
                                    uint64_t** out_vals,
                                    size_t* out_nnz,
                                    uint64_t* out_rhs);

/*
   Free buffers allocated by sprref_inc_is_solved.
   Both pointers must come from the same is_solved call (or be NULL).
*/
SPRREF_API void sprref_inc_buffer_free(uint32_t* cols, uint64_t* vals);

/* Diagnostics --------------------------------------------------------------- */

SPRREF_API size_t sprref_inc_rank(const sprref_inc_t* h);
SPRREF_API uint32_t sprref_inc_nvars(const sprref_inc_t* h);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SPRREF_INCREMENTAL_H */
