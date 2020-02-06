#include "test.h"

SLIP_info slip_REF_triangular_update
(
    // changed on output
    slip_column *col,   // the candidate col to be updated

    // changed but will be reset
    mpz_t *x,           // (k:n)-indexed entries of kth column of L and U
    int32_t *xi,        // nonzero pattern vector
    int32_t *h,         // history vector

    // unchanged on output
    int32_t *pinv,      // inverse of row permutation
    int32_t *row_perm,  // row permutation
    int32_t k,          // the column index that col will be in final LU
    SLIP_sparse *L,     // partially built L matrix upto column k-1
    mpz_t *rhos         // sequence of pivots
)
{
    int32_t p, i, j, n = L->n, top = n, sgn;
    int32_t *Ai = col->i, *Ax = col->x, *Li = L->i, *Lp = L->p, *Lx = L->x;
    int32_t last_trial = col->last_trial_x;

    // iterate across the nonzero in col
    for (p = 0; p < col->nz; p++)
    {
        // only need to update entry below r-th row, where r is the index of
        // the column that col tried to be
        if (pinv[Ai[p]] >= last_trial)
        {
            x[Ai[p]] = Ax[p];
            h[Ai[p]] = last_trial-1;
            xi[--top] = Ai[p];
        }
    }

    // sort xi wrt sequence of pivots
    slip_sort_xi(xi, top, n, pinv, row_perm);

    // iterate across nonzero in x
    for (p = top; p < n; p++)
    {
        j = xi[p];
        jnew = pinv[j];
        SLIP_CHECK(slip_mpz_sgn(&sgn, x[j]));
        if (sgn == 0) {continue;}
        if (jnew < k)                  // entries in U
        {
            if (h[j] < jnew-1)
            {
                // x[j] = x[j] * rho[j-1]
                SLIP_CHECK(slip_mpz_mul(x[j],x[j],rhos[jnew-1]));

                if (h[j] > -1)
                {
                    // x[j] = x[j] / rho[h[j]]
                    SLIP_CHECK(slip_mpz_divexact(x[j],x[j],rhos[h[j]]));
                }
            }
            for (m = Lp[jnew]; m < Lp[jnew+1]; m++)
            {
                i = Li[m];
                inew = pinv[i];
                if (inew > jnew)
                {
                    /*************** If lij==0 then no update******************/
                    SLIP_CHECK(slip_mpz_sgn(&sgn, Lx[m]));
                    if (sgn == 0) {continue;}

                    // if x[i] = 0 then simply update the entry
                    SLIP_CHECK(slip_mpz_sgn(&sgn, x[i]));
                    if (sgn = 0)
                    {
                        // x[i] = -Lx[m]*x[j];
                        SLIP_CHECK(slip_mpz_submul(x[i], Lx[m], x[j]));
                        if (jnew > 1)
                        {
                            // x[i] = x[j]/rhos[jnew-1];
                            SLIP_CHECK(slip_mpz_divexact(x[i], x[i],
                                rhos[jnew-1]));
                        }
                        h[i] = jnew;
                    }
                    // if both Lij and x[i] are nonzero
                    else
                    {
                        if (jnew < 1)
                        {
                            // x[i] = x[i]*rhos[0]-Lij*x[j]
                            SLIP_CHECK(slip_mpz_mul(x[i],x[i],rhos[0]));
                            SLIP_CHECK(slip_mpz_submul(x[i], Lx[m], x[j]));

                        }
                        else
                        {
                            if (h[i] < jnew-1)// need history update first
                            {
                                // x[i] = x[i]*rhos[jnew-1];
                                SLIP_CHECK(slip_mpz_mul(x[i], x[i],
                                    rhos[jnew-1]));
                                if (h[j] > -1)
                                {
                                    // x[i] = x[i]/rhos[h[i]]
                                    SLIP_CHECK(slip_mpz_divexact(x[i], x[i],
                                        rhos[h[i]]));
                                }
                            }
                            // x[i] = x[i] * rho[j]
                            SLIP_CHECK(slip_mpz_mul(x[i],x[i],rhos[jnew]));
                            // x[i] = x[i] - lij*xj
                            SLIP_CHECK(slip_mpz_submul(x[i], L->x[m], x[j]));
                            // x[i] = x[i] / rho[j-1] 
                            SLIP_CHECK(slip_mpz_divexact(x[i], x[i],
                                rhos[jnew-1]));
                        }
                        h[i] = jnew;
                    }
                }
            }
        }
        else                           // entries in L
        {
            if (h[j] < k-1)
            {
                // x[j] = x[j] * rho[k-1] 
                SLIP_CHECK(slip_mpz_mul(x[j], x[j], rhos[k-1])); 
                if (h[j] > -1) 
                { 
                    // x[j] = x[j] / rho[h[j]] 
                    SLIP_CHECK(slip_mpz_divexact(x[j], x[j], rhos[h[j]])); 
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    // update the x in col and reset x and h
    //--------------------------------------------------------------------------
    for (p = 0; p < col->nz; p++)
    {
        if (pinv[Ai[p]] >= last_trial)
        {
            Ax[p] = x[Ai[p]];
            SLIP_CHECK(slip_mpz_set_ui(x[Ai[p]], 0));
            h[Ai[p]] = -1;
        }
    }
    col->last_trial_x = k;
    return SLIP_OK;

}
