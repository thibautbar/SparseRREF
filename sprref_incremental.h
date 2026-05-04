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

/* Boolean returns for sprref_inc_is_solved. */
#define SPRREF_INC_NOT_SOLVED 0
#define SPRREF_INC_SOLVED     1

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

/* Insert -------------------------------------------------------------------- */

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

   On SPRREF_INC_SOLVED, the expression for var_idx is written to out_*:
     - *out_cols, *out_vals are heap-allocated arrays of length *out_nnz.
       Caller must free them via sprref_inc_buffer_free.
     - *out_rhs is the constant term.
   On SPRREF_INC_NOT_SOLVED, out_* are left untouched (and the caller should
   not read them).

   The returned expression is the *negated* off-pivot part: i.e. for the
   normalised pivot row  x_v = rhs - sum_j coeff_j * x_j , out_cols/vals
   contain (j, coeff_j) pairs and out_rhs = rhs.
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
