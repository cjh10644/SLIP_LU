#include "SLIP_LU_internal.h"
#include <assert.h>

SLIP_info slip_update_LU
(
    SLIP_sparse *L,     // partially built L matrix upto column k-1
    SLIP_sparse *U,     // partially built U matrix upto column k-1
    mpz_t  *rhos,       // sequence of pivots
    size_t *rhos_bs,    // bit size of the sequence of pivots
    int32_t *rc,
    slip_column *col,   // the selected candidate column
    int32_t rpiv,       // the row index for the new pivot
    int32_t k           // the column index that col will be in final LU
)
{
    SLIP_info ok = SLIP_OK;
    int32_t p, i, ci, hp;
    int32_t lnz = L->p[k], unz = U->p[k], nz = col->nz, new_unz = col->unz;
    int32_t *h = col->h;
    mpz_t *Ax = col->x;
    size_t size, *Axb = col->bs;

    //--------------------------------------------------------------------------
    // Reallocate memory if necessary
    //--------------------------------------------------------------------------
    if (lnz+nz-new_unz > L->nzmax)
    {
        // Set L->nz = lnz
        L->nz = lnz;
        SLIP_CHECK(slip_sparse_realloc(L));
    }
    if (unz+new_unz+1 > U->nzmax)
    {
        // Set U->nz = unz
        U->nz = unz;
        SLIP_CHECK(slip_sparse_realloc(U));
    }

    //--------------------------------------------------------------------------
    // update the k-th column of L, U, rhos and rhos_bs
    //--------------------------------------------------------------------------
    // Iterate accross the nonzeros in col
    for (p = 0; p < nz; p++)
    {            
        i = col->i[p];
        if (p < new_unz)           // entries belongs to U
        {
            // Place the i location of the U->nz nonzero
            U->i[unz] = i;
            // Place the x value of the U->nz nonzero
            SLIP_CHECK(SLIP_mpz_init_set(U->x[unz], Ax[p]));
            // Increment U->nz
            unz++;
        }
        else if (Axb[p] != 0)      // entries belongs to L
        {
            // remove this nz from row count
            rc[i] --;
            assert(rc[i] >= 0);
            // finish history update
            hp = h[p];
            if (hp < k-1)
            {
                // x[j] = x[j] * rho[k-1]
                SLIP_CHECK(SLIP_mpz_mul(Ax[p], Ax[p], rhos[k-1]));
                if (hp > -1)
                {
                    // x[j] = x[j] / rho[h[j]]
                    SLIP_CHECK(SLIP_mpz_divexact(Ax[p], Ax[p], rhos[hp]));
                }
            }
            // Place the i location of the L->nz nonzero
            L->i[lnz] = i;
            // Place the x value of the L->nz nonzero
            SLIP_CHECK(SLIP_mpz_init_set(L->x[lnz], Ax[p]));
            // Increment L->nz
            lnz++;

            if (i == rpiv)          // pivot element
            {
                SLIP_CHECK(SLIP_mpz_sizeinbase(&size, Ax[p], 2));
                // set rhos[k] = L(k,k)
                SLIP_CHECK(SLIP_mpz_set(rhos[k], Ax[p]));
            //SLIP_gmp_printf("rhos[%d]= %Zd %d\n",k,rhos[k],col->bs[p]);
                rhos_bs[k] = size;

                // add the pivot element to U as well
                // Place the i location of the U->nz nonzero
                U->i[unz] = i;
                // Place the x value of the U->nz nonzero
                SLIP_CHECK(SLIP_mpz_init_set(U->x[unz], Ax[p]));
                // Increment U->nz
                unz++;
                if(abs(col->bs[p]-size)>1)
                {
                    printf("[est real diff ] = [%zu %zu %u ]<----------------------------------------pivot\n",col->bs[p], size, abs(col->bs[p]-size));
                    return SLIP_INCORRECT;
                }
            }
            else
            {
                SLIP_CHECK(SLIP_mpz_sizeinbase(&size, col->x[p], 2));
                if(abs(col->bs[p]-size)>1)
                {
                    printf("[est real diff ] = [%zu %zu %u ]\n",col->bs[p], size, abs(col->bs[p]-size));
                    return SLIP_INCORRECT;
                }
            }
        }
    }

    L->p[k+1] = lnz;
    U->p[k+1] = unz;
    return SLIP_OK;
}
