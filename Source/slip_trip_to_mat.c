# include "SLIP_LU_internal.h"
/*
 * Purpose: This function converts triplet matrix into compressed column
 * matrix A
 */

#define SLIP_FREE_WORKSPACE              \
    SLIP_FREE(w);

SLIP_info slip_trip_to_mat
(
    SLIP_sparse *A,     //matrix stored in ccf that will take B
    int32_t *I,         // Row indices.
    int32_t *J,         // Column indices
    mpz_t *x,           // Values in the matrix
    int32_t n,          // Dimension of the matrix
    int32_t nz          // Number of nonzeros in the matrix
)
{
    if (!I || !J || !x || !A)
    {
        return SLIP_INCORRECT_INPUT;
    }
    int32_t k, p;
    SLIP_info ok;
    int32_t* w = (int32_t*) SLIP_calloc(n, sizeof(int32_t));
    if (!w)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    SLIP_CHECK(slip_sparse_alloc(A, n, n, nz));

    // Column pointers
    for (k = 0; k < nz; k++)
    {
        w[J[k]]++;
    }
    // Column sums for A, w = A->p[0..n-1]
    slip_cumsum(A->p, w, A->n);
    for (k = 0; k < nz; k++)
    {
        p = w[J[k]]++;
        // Place values of i
        A->i[p] = I[k];
        // Place values of x
        if (A->i[p] < 0)
        {
            SLIP_FREE_WORKSPACE;
            return SLIP_INCORRECT_INPUT;
        }
        SLIP_CHECK(slip_mpz_set(A->x[p], x[k]));
    }
    // Number of nonzeros in A
    A->nz = A->nzmax;
    // Delete w
    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}
#undef SLIP_FREE_WORKSPACE