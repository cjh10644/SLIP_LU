# include "SLIP_LU_internal.h"


/* 
 * Purpose: This function populates the SLIP_sparse A by the ccf stored vectors I,
 * p, and x
 */
SLIP_info slip_mpz_populate_mat
(
    SLIP_sparse* A,   // matrix to be populated
    int32_t* I,       // row indices
    int32_t* p,       // column pointers
    mpz_t* x,         // set of values in A
    int32_t n,        // size of the matrix A
    int32_t nz        // number of nonzeros in the matrix A
)
{
    if (!I || !p || n <= 0 || nz < 0 || !x || !A ) 
    {
        return SLIP_INCORRECT_INPUT;
    }
    SLIP_info ok = slip_sparse_alloc(A, n, n, nz);// Allocate the matrix A
    if (ok != SLIP_OK) 
    {
        return ok;
    }
    A->nz = nz;
    for (int32_t k = 0; k <= n; k++)              // Set A->p
    {
        A->p[k] = p[k];
    }
    for (int32_t k = 0; k < nz; k++)              // Set A->i and A->x
    {
        A->i[k] = I[k];
        if (A->i[k] < 0)
        {
            return SLIP_INCORRECT_INPUT;
        }
        ok = slip_mpz_set(A->x[k],x[k]);
        if (ok != SLIP_OK) {return ok;}
    }
    return ok;
}