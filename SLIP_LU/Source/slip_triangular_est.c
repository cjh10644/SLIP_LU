#include "SLIP_LU_internal.h"
#include <assert.h>

SLIP_info slip_triangular_est
(
    // changed on output
    slip_column *col,   // the candidate col whose bit size need to be estimated
    int32_t *x,         // array of indices of nonzeros in col, whose nonzero 
                        // pattern is found using xi. defaulted -1 for each
                        // entry. garbage on input and output
    int32_t *xi,        // nonzero pattern vector, first nonzero can be found at
                        // n - col->nz. garbage on input and output
    int32_t *rc,

    // unchanged on output
    int32_t *pinv,      // inverse of row permutation
    int32_t *row_perm,  // row permutation
    int32_t k,          // the column index that col tend to be in final LU
    SLIP_sparse *L,     // partially built L matrix upto column k-1
    mpz_t  *rhos,       // sequence of pivots
    size_t *rhos_bs,    // bit size of the sequence of pivots
    int32_t bound       // worst case bit length for each mpz entry
)
{
    SLIP_info ok = SLIP_OK;
    int32_t p, m, i, j, n = L->n, top = n, jnew, inew, ci, cj, sgn;
    int32_t nz = col->nz, unz, unz_new = 0, max_mpz = col->max_mpz,
            real_lnz = 0;
    int32_t *Ai = col->i, *h = col->h, *Li = L->i, *Lp = L->p;
    size_t *Axb = col->bs, size;
    mpz_t *Ax = col->x, *Lx = L->x;
    int32_t last_trial;

    if (col->unz == -1)                // this col is newly added to candidates
    {
        last_trial = 0;
        unz = 0;
    }
    else
    {
        last_trial = k-1;              // col is updated in every iteration
        unz = col->unz;
    }
    //--------------------------------------------------------------------------
    // iterate across the all nonzero in col that would be in L
    //--------------------------------------------------------------------------
    for (ci = unz; ci < col->nz; ci++)
    {
        i = Ai[ci];
        x[i] = ci; // only get the index of this entry in col

        // find the reachable of this node
        if(!SLIP_MARKED(Lp, i))
        {
            slip_dfs(&top, i, L, xi, xi+n, pinv);
        }
    }

    //--------------------------------------------------------------------------
    // Iterate all across the nonzero pattern and restore L and reset x
    //--------------------------------------------------------------------------
    for ( p = top; p < n; p++)
    {
        i = xi[p];                     // row index

        // restore L
        SLIP_MARK(Lp, i);

        //----------------------------------------------------------------------
        // find new reachable nodes and add them to col, and initialize
        // corresponding x and history value
        //----------------------------------------------------------------------
        ci = x[i];                     // index in col
        if (ci == -1 ||                // if x[i] has not been used, or
            ci >= nz ||                // x[i] points to non-existing entry
                                       // in col, or
           (ci != -1 && ci < nz &&     // x[i] points to existing entry
            i != Ai[ci]))              // incorrectly, then it is new reachable
        {
            Axb[nz] = 0;
            Ai [nz] = i;
            h  [nz] =-1;               // this is only for consistency, the h
                                       // value for new reachable doesn't matter
            x[i] = nz;

            // max_mpz always >= nz
            if (max_mpz == nz)         // allocate memory for newly found entry
            {                          // if needed
                // allocate worse case bit size for each entry to
                // avoid possible reallocation
                if (SLIP_mpz_init2(Ax[max_mpz], bound) != SLIP_OK)
                {
                    SLIP_MPZ_SET_NULL(Ax[max_mpz]);
                    // update the number of mpz allocated to help delete
                    // mpz array
                    col->max_mpz = max_mpz;
                    return SLIP_OUT_OF_MEMORY;
                }
                max_mpz ++;
            }
            else
            {
                SLIP_CHECK(SLIP_mpz_set_ui(Ax[nz], 0));
            }
            nz ++;
        }

        //----------------------------------------------------------------------
        // move all entries that would be in U to the top of xi
        //----------------------------------------------------------------------
        if (pinv[i] < k)
        {
            m = top+unz_new;           // index in xi for newly found nz in U
            xi[p] = xi[m];
            xi[m] = i;
            unz_new++;
        }
    }
    col->max_mpz = max_mpz;

    //--------------------------------------------------------------------------
    // sort the first unz_new entries if there are more than 1 newly found nz
    // in U, this happens for most newly added col
    //--------------------------------------------------------------------------
    if (unz_new > 1)
    {
        // sort xi wrt sequence of pivots
        slip_sort_xi(xi, top, top+unz_new, pinv, row_perm);
    }
    //--------------------------------------------------------------------------
    // iterate across nonzeros in x
    //--------------------------------------------------------------------------
    for (p = top; p < n; p++)
    {
        j = xi[p];
        jnew = pinv[j];
        cj = x[j];                     // index in col

        if (jnew < k)                  // entries in U
        {
            //------------------------------------------------------------------
            // get rid of the zero entry from the list
            //------------------------------------------------------------------
            SLIP_CHECK(SLIP_mpz_sgn(&sgn, Ax[cj]));
            if (sgn == 0)
            {
                // move the zero entry to the last of the list
                nz --;
                if (cj != nz)
                {
                    // change the index pointer first
                    x[Ai[nz]] = cj;

                    // copy all values from index of nz to cj
                    SLIP_mpz_swap(Ax[cj], Ax[nz]);
                    Ai[cj] = Ai[nz];
                    Axb[cj] = Axb[nz];
                    h[cj] = h[nz];
                }
                continue;
            }
            //------------------------------------------------------------------
            // finish history update if needed
            //------------------------------------------------------------------
            h[cj] = SLIP_UNFLIP(h[cj]);
            if (h[cj] < jnew-1)
            {
                // x[j] = x[j] * rho[j-1]
                SLIP_CHECK(SLIP_mpz_mul(Ax[cj],Ax[cj],rhos[jnew-1]));

                if (h[cj] > -1)
                {
                    // x[j] = x[j] / rho[h[j]]
                    SLIP_CHECK(SLIP_mpz_divexact(Ax[cj], Ax[cj], rhos[h[cj]]));
                }
            }

            //------------------------------------------------------------------
            // perform IPGE
            //------------------------------------------------------------------
            for (m = Lp[jnew]; m < Lp[jnew+1]; m++)
            {
                i = Li[m];
                inew = pinv[i];
                ci = x[i];             // index in col
                if (inew > jnew)
                {
                    //----------------------------------------------------------
                    // If lij==0 then no update
                    //----------------------------------------------------------
                    SLIP_CHECK(SLIP_mpz_sgn(&sgn, Lx[m]));
                    if (sgn == 0) {continue;}

                    //----------------------------------------------------------
                    // if x[i] = 0 then simply update the entry
                    //----------------------------------------------------------
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
                        if (h[ci] == -1 && inew >= k)
                        {
                            rc[i] ++;
                            assert(rc[i] <= 10);
                        }
                    }
                    //----------------------------------------------------------
                    // if both Lij and x[i] are nonzero
                    //----------------------------------------------------------
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
                            h[ci] = SLIP_UNFLIP(h[ci]);
                            if (h[ci] < jnew-1)// need history update first
                            {
                                // x[i] = x[i]*rhos[jnew-1];
                                SLIP_CHECK(SLIP_mpz_mul(Ax[ci], Ax[ci],
                                    rhos[jnew-1]));
                                if (h[ci] > -1)
                                {
                                    // x[i] = x[i]/rhos[h[i]]
                                    SLIP_CHECK(SLIP_mpz_divexact(Ax[ci],
                                        Ax[ci], rhos[h[ci]]));
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

                    //----------------------------------------------------------
                    // update h for all inew > jnew, but flip to mark as updated
                    //----------------------------------------------------------
                    h[ci] = SLIP_FLIP(jnew);
                }
            }

            //------------------------------------------------------------------
            // swap entries in col with indices of cj and col->unz
            //------------------------------------------------------------------
            if (cj != unz)
            {
                // change the index pointer first
                x[Ai[unz]] = cj;

                // both row indices and mpz values need to be swapped
                SLIP_mpz_swap(Ax[cj], Ax[unz]);
                Ai[cj] = Ai[unz];
                Ai[unz] = j;

                // entries in U do not care about bit size nor h
                Axb[cj] = Axb[unz];
                h[cj] = h[unz];
            }
            unz ++;
        }

        //----------------------------------------------------------------------
        // only need to update bit size of entries that would be in L
        //----------------------------------------------------------------------
        else
        {
            if (h[cj] < -1)            // entries has just been updated
            {
                SLIP_MARK(h, cj);
                SLIP_CHECK(SLIP_mpz_sgn(&sgn, Ax[cj]));
                if (sgn == 0)
                {
                    Axb[cj] = 0;
                    h[cj] = -1;        // h won't affect explicit 0 entries
                    rc[j] --;
                    assert(rc[j] >= 0);
                    continue;
                }
                SLIP_CHECK(SLIP_mpz_sizeinbase(&size, Ax[cj], 2));
                if (h[cj] < k-1)
                {
                    // x[j] = x[j]*rhos[k-1]
                    Axb[cj] = size+rhos_bs[k-1];
                    if (h[cj] > -1)
                    {
                        // x[j] = x[j]/rhos[h[j]];
                        Axb[cj] = Axb[cj]-rhos_bs[h[cj]];
                    }
                }
                else
                {
                    Axb[cj] = size;
                }
            }
            else
            {
                if (Axb[cj] == 0)
                {
                    continue;
                }
                if (last_trial == 0)   // this col is newly added
                {
                    Axb[cj] = Axb[cj]+rhos_bs[k-1];
                }
                else
                {
                    Axb[cj] = Axb[cj]+rhos_bs[k-1]-rhos_bs[k-2];
                }
            }
            real_lnz ++;
        }
    }
    //--------------------------------------------------------------------------
    // update column info
    //--------------------------------------------------------------------------
    col->nz = nz;
    col->real_lnz = real_lnz;
    col->unz = unz;

    return SLIP_OK;
}
