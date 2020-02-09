#include "SLIP_LU_internal.h"

SLIP_info slip_triangular_estimate
(
    // changed on output
    slip_column *col,   // the candidate col whose bit size need to be estimated
    int32_t *x,         // array of indices of nonzeros in col, whose nonzero 
                        // pattern is found using xi. defaulted -1 for each
                        // entry. garbage on input, modified on output
    int32_t *xi,        // nonzero pattern vector, first nonzero can be found at
                        // n - col->nz. garbage on input, modified on output
    int32_t *h,         // history vector, garbage on input, modified on output

    // unchanged on output
    int32_t *pinv,      // inverse of row permutation
    int32_t *row_perm,  // row permutation
    int32_t k,          // the column index that col tend to be in final LU
    SLIP_sparse *L,     // partially built L matrix upto column k-1
    size_t *Lb,         // the bit size of the known L->x
    size_t *rhos        // bit size of the sequence of pivots
)
{
    int32_t p, m, i, j, n = L->n, top = n, jnew, inew, ci, cj;
    int32_t *Ai = col->i, *Li = L->i, *Lp = L->p;
    size_t *Axb = col->bs;
    int32_t last_trial = col->last_trial_bs;

    // iterate across the nonzero in col
    for (ci = 0; ci < col->nz; ci++)
    {
        i = Ai[ci];
        x[i] = ci; // only get the index of this entry in col
        h[i] = last_trial-1;

        // find the reachable of this node
        if(!SLIP_MARKED(Lp, i))
        {
            slip_dfs(&top, i, L, xi, xi+n, pinv);
        }
    }

    // Iterate all across the nonzero pattern and restore L and reset x
    for ( p = top; p < n; p++)
    {
        i = xi[p];                     // row index
        SLIP_MARK(Lp, i);

        // find new reachable nodes and add them to col, and initialize
        // corresponding x and history value
        ci = x[i];                     // index in col
        if (ci == -1 ||                // if x[i] has not been used, or
            ci >= col->nz ||           // x[i] points to non-existing entry
                                       // in col, or
           (ci != -1 && ci < col->nz &&// x[i] points to existing entry
            i != Ai[ci]))              // incorrectly, then it is new reachable
        {
            Axb[col->nz] = 0;
            Ai [col->nz] = i;
            x[i] = col->nz;
            h[i] = last_trial-1;       // this is only for consistency, the h
                                       // value for new reachable doesn't matter
            col->nz ++;
        }
    }

    // sort xi wrt sequence of pivots
    slip_sort_xi(xi, top, n, pinv, row_perm);

    // iterate across nonzeros in x
    for (p = top; p < n; p++)
    {
        j = xi[p];
        jnew = pinv[j];
        cj = x[j];                     // index in col

        // only need to update entry L(r:n-1,r:n-1), where r is the index of
        // the column that col tried to be
        if (jnew < last_trial) {continue;}

        if (Axb[cj] == 0) {continue;}  // this woulnd not occur mostly
        if (jnew < k)                  // entries in U
        {
            if (h[j] < jnew-1)
            {
                if (h[j] > -1)
                {
                    // x[j] = x[j]*rhos[jnew-1]/rhos[h[j]];
                    Axb[cj] = Axb[cj]+rhos[jnew-1]-rhos[h[j]];
                }
                else
                {
                    // x[j] = x[j]*rhos[jnew-1];
                    Axb[cj] = Axb[cj]+rhos[jnew-1];
                }
            }
            for (m = Lp[jnew]; m < Lp[jnew+1]; m++)
            {
                i = Li[m];
                inew = pinv[i];
                ci = x[i];             // index in col
                if (inew > jnew)
                {
                    if (Lb[m] == 0) {continue;}
                    // if x[i] = 0 then simply update the entry
                    if (Axb[ci] == 0)
                    {
                        if (jnew < 1)
                        {
                            // x[i] = -Lb[m]*x[j];
                            Axb[ci] = Lb[m]+Axb[cj];
                        }
                        else
                        {
                            // x[i] = -Lb[m]*x[j]/rhos[jnew-1];
                            Axb[ci] = Lb[m]+Axb[cj]-rhos[jnew-1];
                        }
                    }
                    // if both Lij and x[i] are nonzero
                    else
                    {
                        if (jnew < 1)
                        {
                            // x[i] = x[i]*rhos[0]-Lb[m]*x[j]
                            Axb[ci] = SLIP_MAX(Axb[ci]+rhos[0],Lb[m]+Axb[cj])+1;
                        }
                        else
                        {
                            if (h[i] < jnew-1)// need history update first
                            {
                                if (h[i] > -1)
                                {
                                    // x[i] = x[i]*rhos[jnew-1]/rhos[h[i]];
                                    // x[i] = x[i]*rhos[jnew]/rhos[jnew-1]
                                    //        -Lb[m]*x[j]/rhos[jnew-1];
                                    Axb[ci] = SLIP_MAX(
                                               Axb[ci]+rhos[jnew]-rhos[h[i]],
                                               Lb[m]+Axb[cj]-rhos[jnew-1])+1;
                                }
                                else
                                {
                                    // x[i] = x[i]*rhos[jnew-1];
                                    // x[i] = x[i]*rhos[jnew]/rhos[jnew-1]
                                    //        -Lb[m]*x[j]/rhos[jnew-1];
                                    Axb[ci] = SLIP_MAX(Axb[ci]+rhos[jnew],
                                               Lb[m]+Axb[cj]-rhos[jnew-1])+1;
                                }
                            }
                            else// no need for history update
                            {
                                // x[i] = x[i]*rhos[jnew]/rhos[jnew-1]
                                //        -Lb[m]*x[j]/rhos[jnew-1];
                                Axb[ci] = SLIP_MAX(
                                           Axb[ci]+rhos[jnew]-rhos[jnew-1],
                                           Lb[m]+Axb[cj]-rhos[jnew-1])+1;
                            }
                        }
                    }
                    // update h for all inew > jnew
                    h[i] = jnew;
                }
            }
        }
        else                           // entries in L
        {
            if (h[j] < k-1)
            {
                // x[j] = x[j]*rhos[k-1]
                Axb[cj] = Axb[cj]+rhos[k-1];
                if (h[j] > -1)
                {
                    // x[j] = x[j]/rhos[h[j]];
                    Axb[cj] = Axb[cj]-rhos[h[j]];
                }
            }
        }
        if(Axb[cj] <0){printf("x[%d] is negative!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",j);}
    }

    //--------------------------------------------------------------------------
    // update column info
    //--------------------------------------------------------------------------

    col->last_trial_bs = k;

    return SLIP_OK;
}
