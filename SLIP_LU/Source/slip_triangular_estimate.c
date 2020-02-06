#include "test.h"

SLIP_info slip_triangular_estimate
(
    // changed on output
    slip_column *col,   // the candidate col whose bit size need to be estimated

    // changed but will be reset
    int32_t *x,         // estimation result
    int32_t *xi,        // nonzero pattern vector
    int32_t *h,         // history vector

    // unchanged on output
    int32_t *pinv,      // inverse of row permutation
    int32_t *row_perm,  // row permutation
    int32_t k,          // the column index that col tend to be in final LU
    SLIP_sparse *L,     // partially built L matrix upto column k-1
    int32_t *Lb,        // the bit size of the known L->x
    int32_t *rhos       // bit size of the sequence of pivots
)
{
    int32_t p, i, j, n = L->n, top = n;
    int32_t *Ai = col->i, *Axb = col->bs, *Li = L->i, *Lp = L->p;
    int32_t last_trial = col->last_trial_bs;

    // iterate across the nonzero in col
    for (p = 0; p < col->nz; p++)
    {
        // only need to update entry below r-th row, where r is the index of
        // the column that col tried to be
        if (pinv[Ai[p]] >= last_trial)
        {
            x[Ai[p]] = Axb[p];
            h[Ai[p]] = last_trial-1;

            // find the reachable of this node
            if(!SLIP_MARKED(L->p, Ai[p]))
            {
                slip_dfs(&top, Ai[p], L, xi, xi+n, pinv);
            }
        }
    }

    // Restore L
    for ( p = top; p < n; p++)
    {   
        SLIP_MARK(L->p, xi[p]);
    }

    // sort xi wrt sequence of pivots
    slip_sort_xi(xi, top, n, pinv, row_perm);

    // iterate across nonzeros in x
    for (p = top; p < n; p++)
    {
        j = xi[p];
        jnew = pinv[j];
        if (x[j] == 0) {continue;}     // this woulnd not occur mostly
        if (jnew < k)                  // entries in U
        {
            if (h[j] < jnew-1)
            {
                if (h[j] > -1)
                {
                    // x[j] = x[j]*rhos[jnew-1]/rhos[h[j]];
                    x[j] = x[j]+rhos[jnew-1]-rhos[h[j]];
                }
                else
                {
                    // x[j] = x[j]*rhos[jnew-1];
                    x[j] = x[j]+rhos[jnew-1];
                }
            }
            for (m = Lp[jnew]; m < Lp[jnew+1]; m++)
            {
                i = Li[m];
                inew = pinv[i];
                if (inew > jnew)
                {
                    if (Lb[m] == 0) {continue;}
                    // if x[i] = 0 then simply update the entry
                    if (x[i] = 0)
                    {
                        if (jnew < 1)
                        {
                            // x[i] = -Lb[m]*x[j];
                            x[i] = Lb[m]+x[j];
                        }
                        else
                        {
                            // x[i] = -Lb[m]*x[j]/rhos[jnew-1];
                            x[i] = Lb[m]+x[j]-rhos[jnew-1];
                        }
                        h[i] = jnew;
                    }
                    // if both Lij and x[i] are nonzero
                    else
                    {
                        if (jnew < 1)
                        {
                            // x[i] = x[i]*rhos[0]-Lb[m]*x[j]
                            x[i] = max(x[i]+rhos[0], Lb[m]+x[j])+1;
                        }
                        else
                        {
                            if (h[i] < jnew-1)// need history update first
                            {
                                if (h[j] > -1)
                                {
                                    // x[i] = x[i]*rhos[jnew-1]/rhos[h[i]];
                                    // x[i] = x[i]*rhos[jnew]/rhos[jnew-1]
                                    //        -Lb[m]*x[j]/rhos[jnew-1];
                                    x[i] = max(x[i]+rhos[jnew]-rhos[h[i]],
                                               Lb[m]+x[j]-rhos[jnew-1])+1;
                                }
                                else
                                {
                                    // x[i] = x[i]*rhos[jnew-1];
                                    // x[i] = x[i]*rhos[jnew]/rhos[jnew-1]
                                    //        -Lb[m]*x[j]/rhos[jnew-1];
                                    x[i] = max(x[i]+rhos[jnew],
                                               Lb[m]+x[j]-rhos[jnew-1])+1;
                                }
                            }
                            else// no need for history update
                            {
                                // x[i] = x[i]*rhos[jnew]/rhos[jnew-1]
                                //        -Lb[m]*x[j]/rhos[jnew-1];
                                x[i] = max(x[i]+rhos[jnew]-rhos[jnew-1],
                                           Lb[m]+x[j]-rhos[jnew-1])+1;
                            }
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
                // x[j] = x[j]*rhos[k-1]
                x[j] = x[j]+rhos[k-1];
                if (h[j] > -1)
                {
                    // x[j] = x[j]/rhos[h[j]];
                    x[j] = x[j]-rhos[h[j]];
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // update the bit size estimation in col and reset x and h
    //--------------------------------------------------------------------------
    // firstly update the entries that have been in col
    for (p = 0; p < col->nz; p++)
    {
        if (x[Ai[p]] != 0)
        {
            Axb[p] = x[Ai[p]];
            x[Ai[p]] = 0;
            h[Ai[p]] = -1;
        }
    }

    // then add entries that were newly found
    for (p = top; p < n; p++)
    {
        if (x[xi[p]] != 0)
        {
            Ai[col->nz] = xi[p];
            Axb[col->nz] = x[xi[p]];
            x[xi[p]] = 0;
            h[xi[p]] = -1;
            (col->nz)++;
        }
    }
    col->last_trial_bs = k;
    return SLIP_OK;
}
