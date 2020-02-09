#include "SLIP_LU_internal.h"

SLIP_info slip_REF_triangular_update
(
    // changed on output
    slip_column *col,   // the candidate col to be updated
    int32_t *x,         // array of indices of nonzeros in col, whose nonzero
                        // pattern is found using xi. garbage on input,
                        // modified on output
    int32_t *xi,        // nonzero pattern vector, first nonzero can be found at
                        // n - col->nz. garbage on input, modified on output
    int32_t *h,         // history vector, garbage on input, modified on output

    // unchanged on output
    int32_t *pinv,      // inverse of row permutation
    int32_t *row_perm,  // row permutation
    int32_t k,          // the column index that col will be in final LU
    SLIP_sparse *L,     // partially built L matrix upto column k-1
    mpz_t *rhos         // sequence of pivots
)
{
    SLIP_info ok = SLIP_OK;
    int32_t p, m, i, j, n = L->n, top = n, sgn, jnew, inew, ci, cj;
    int32_t *Ai = col->i, *Li = L->i, *Lp = L->p;
    mpz_t *Ax = col->x, *Lx = L->x;
    int32_t last_trial = col->last_trial_x;

    // if col has been updated to be k-th pivot column, then nothing
    // needs to be done. This happens when this column has multiple minimums
    // of bit size estimation
    if (last_trial == k) {return SLIP_OK;}

    // iterate across the nonzero mpz in col
    for (ci = 0; ci < col->nz_mpz; ci++)
    {
        i = Ai[ci];
        x[i] = ci; // only get the index of this entry in col
        h[i] = last_trial-1;
        xi[--top] = i;
    }
    // then add other reachable nodes to the nonzero pattern list
    // and set x[i] = 0
    for (ci = col->nz_mpz; ci < col->nz; ci++)
    {
        i = Ai[ci];
        x[i] = ci; // only get the index of this entry in col
        SLIP_CHECK(SLIP_mpz_set_ui(Ax[ci], 0));
        h[i] = last_trial-1;
        xi[--top] = i;
    }
    col->nz_mpz = col->nz;

    // sort xi wrt sequence of pivots
    slip_sort_xi(xi, top, n, pinv, row_perm);

    // iterate across nonzero in x
    for (p = top; p < n; p++)
    {
        j = xi[p];
        jnew = pinv[j];
        cj = x[j];                     // index in col

        // only need to update entry L(r:n-1,r:n-1), where r is the index of
        // the column that col tried to be
        if (jnew < last_trial) {continue;}

        SLIP_CHECK(SLIP_mpz_sgn(&sgn, Ax[cj]));
        if (sgn == 0) {continue;}
        if (jnew < k)                  // entries in U
        {
            if (h[j] < jnew-1)
            {
                // x[j] = x[j] * rho[j-1]
                SLIP_CHECK(SLIP_mpz_mul(Ax[cj],Ax[cj],rhos[jnew-1]));

                if (h[j] > -1)
                {
                    // x[j] = x[j] / rho[h[j]]
                    SLIP_CHECK(SLIP_mpz_divexact(Ax[cj], Ax[cj], rhos[h[j]]));
                }
            }
            for (m = Lp[jnew]; m < Lp[jnew+1]; m++)
            {
                i = Li[m];
                inew = pinv[i];
                ci = x[i];             // index in col
                if (inew > jnew)
                {
                    /*************** If lij==0 then no update******************/
                    SLIP_CHECK(SLIP_mpz_sgn(&sgn, Lx[m]));
                    if (sgn == 0) {continue;}

                    // if x[i] = 0 then simply update the entry
                    SLIP_CHECK(SLIP_mpz_sgn(&sgn, Ax[ci]));
                    if (sgn == 0)
                    {
                        // x[i] = -Lx[m]*x[j];
                        SLIP_CHECK(SLIP_mpz_submul(Ax[ci], Lx[m], Ax[cj]));
                        if (jnew >= 1)
                        {
                            // x[i] = x[j]/rhos[jnew-1];
                            SLIP_CHECK(SLIP_mpz_divexact(Ax[ci],
                                Ax[ci], rhos[jnew-1]));
                        }
                    }
                    // if both Lij and x[i] are nonzero
                    else
                    {
                        if (jnew < 1)
                        {
                            // x[i] = x[i]*rhos[0]-Lij*x[j]
                            SLIP_CHECK(SLIP_mpz_mul(Ax[ci], Ax[ci], rhos[0]));
                            SLIP_CHECK(SLIP_mpz_submul(Ax[ci], Lx[m], Ax[cj]));

                        }
                        else
                        {
                            if (h[i] < jnew-1)// need history update first
                            {
                                // x[i] = x[i]*rhos[jnew-1];
                                SLIP_CHECK(SLIP_mpz_mul(Ax[ci], Ax[ci],
                                    rhos[jnew-1]));
                                if (h[i] > -1)
                                {
                                    // x[i] = x[i]/rhos[h[i]]
                                    SLIP_CHECK(SLIP_mpz_divexact(Ax[ci],
                                        Ax[ci], rhos[h[i]]));
                                }
                            }
                            // x[i] = x[i] * rho[j]
                            SLIP_CHECK(SLIP_mpz_mul(Ax[ci], Ax[ci],rhos[jnew]));
                            // x[i] = x[i] - lij*xj
                            SLIP_CHECK(SLIP_mpz_submul(Ax[ci], Lx[m],
                                Ax[cj]));
                            // x[i] = x[i] / rho[j-1] 
                            SLIP_CHECK(SLIP_mpz_divexact(Ax[ci], Ax[ci],
                                rhos[jnew-1]));
                        }
                    }
                    // update h for all inew>jnew
                    h[i] = jnew;
                }
            }
        }
        else                           // entries in L
        {
            if (h[j] < k-1)
            {
                // x[j] = x[j] * rho[k-1] 
                SLIP_CHECK(SLIP_mpz_mul(Ax[cj], Ax[cj], rhos[k-1])); 
                if (h[j] > -1) 
                { 
                    // x[j] = x[j] / rho[h[j]] 
                    SLIP_CHECK(SLIP_mpz_divexact(Ax[cj], Ax[cj], rhos[h[j]])); 
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    // update column info
    //--------------------------------------------------------------------------
    col->last_trial_x = k;

    return SLIP_OK;

}
