//TODO consider test case with explicit 0
#define CAND_SIZE 10
#define SLIP_FREE_WORKSPACE         \
{                              \
    slip_delete_cand_columns(&cands);\
    SLIP_FREE(priority);       \
    SLIP_FREE(rc);       \
    SLIP_FREE(row_perm);       \
    SLIP_FREE(xi);             \
    SLIP_FREE(x);          \
    SLIP_FREE(rhos_bs);         \
    SLIP_MPFR_CLEAR(temp);      \
}

#include"SLIP_LU_internal.h"
#include <assert.h>

SLIP_info SLIP_LU_analyze_and_factorize1
(
    SLIP_sparse *L,
    SLIP_sparse *U,
    SLIP_sparse *A,
    SLIP_LU_analysis *S,
    mpz_t *rhos,
    int32_t *pinv,
    SLIP_options *option
)
{
    if (L == NULL || U == NULL || rhos == NULL || pinv == NULL ||
        A == NULL || S == NULL || S->q == NULL || option == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // initialize variables and allocate memory
    //--------------------------------------------------------------------------
    double total_time = 0, est_time = 0, init_time = 0;
    // infomation of the A matrix
    int32_t n = A->n, nz = A->nz;

    // number of nonzeros in L and U
    int32_t lnz = 0, unz = 0;

    // the nth minimum digit size, which is used to determine the column to
    // perform REF triangular update
    int32_t nthMin = 0;

    // whether or not the pivot column is found
    bool foundPivotColumn = false;

    // number of columns that have been taken as pivots or in the candidate list
    int32_t usedcol = 0;

    // other miscellaneous variables
    SLIP_info ok = SLIP_OK;
    int32_t sgn, i, j, k, ci;
    size_t size;
    slip_column *col;
    int32_t sigma_index;
    mpfr_t temp; SLIP_MPFR_SET_NULL(temp);

    // vectors to be used, corresponding usage can be found below
    int32_t *xi = NULL, *row_perm = NULL, *col_perm = S->q, *x = NULL,
            *priority = NULL, *rc = NULL;
    size_t *rhos_bs = NULL;
    slip_cand_columns *cands;

    // candidate columns for next pivot column
    int32_t ncand = SLIP_MIN(CAND_SIZE, n), nthcand;
    cands = slip_initialize_cand_columns(ncand, n);
    priority = (int32_t *) SLIP_malloc(ncand*sizeof(int32_t));
    rc = (int32_t *) SLIP_malloc(n*sizeof(int32_t));

    // row permutaion, i.e., inew = row_perm[i] is
    // the row index in A that will be permuted to L(i,:)
    row_perm = (int32_t *) SLIP_malloc(n*sizeof(int32_t));
    // nonzero pattern vector
    xi = (int32_t *) SLIP_malloc(2*n*sizeof(int32_t));
    // used in slip_triangular_estimate and slip_REF_triangular_update as the
    // pointers to corresponding entries in slip_column struct. That is, x[i]
    // indicates nonzero at i-th row, while its value is the x[i]-th entry in
    // slip_column struct
    x = (int32_t *) SLIP_malloc(n*sizeof(int32_t));

    // bit size of pivot elements
    rhos_bs = (size_t *) SLIP_malloc(n*sizeof(size_t));

    if (!cands || !row_perm || !xi || !rhos_bs || !x ||!priority)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    // initialize each entry of x_bs and h as required by
    // slip_triangular_estimate and slip_REF_triangular_update
    for (i = 0; i < n; i++)
    {
        x[i] = -1;
        rc[i] = 0;
    }

    // Assume the number of nonzeros in L or U is 10 times of that in A
    // only allocate memory for L->x and U->x but not initializing each entry
    SLIP_CHECK(slip_sparse_alloc2(L, n, n, 10*nz));
    SLIP_CHECK(slip_sparse_alloc2(U, n, n, 10*nz));

    //--------------------------------------------------------------------------
    // This section of the code computes a bound for the worst case bit-length
    // of each entry in the matrix. This bound is used to allocate the size of
    // each mpz number in the candidate columns (as slip_column struct, col->x
    // to be more specifically). As a result of this allocation, computing the
    // values in L and U via repeated triangular solves will not require
    // intermediate memory reallocations from the GMP library.
    //
    // This bound is based on a relaxation of sparse Hadamard's bound
    //--------------------------------------------------------------------------
    SLIP_CHECK(SLIP_mpfr_init2(temp, 256));

    // sigma_index is the index of the largest initial entry in A. First we
    // initialize sigma_index to be the index of the first nonzero in A
    sigma_index = 0;

    // Iterate throughout A and set sigma_index = index of max (A)
    for (i = 1; i < A->nz; i++)
    {
        if(mpz_cmpabs(A->x[sigma_index],A->x[i]) < 0)
        {
            sigma_index = i;
        }
    }

    // gamma is the number of nonzeros in the most dense column of A. First, we
    // initialize gamma to be the number of nonzeros in A(:,1).
    int32_t gamma = A->p[1];

    // Iterate throughout A and obtain gamma as the most dense column of A
    for (i = 1; i<n; i++)
    {
        if( gamma < A->p[i+1] - A->p[i])
        {
            gamma = A->p[i+1]-A->p[i];
        }
    }

    // temp = |sigma|
    SLIP_CHECK(SLIP_mpfr_set_z(temp,A->x[sigma_index],option->SLIP_MPFR_ROUND));
    SLIP_CHECK(SLIP_mpfr_abs(temp,temp,option->SLIP_MPFR_ROUND));

    //--------------------------------------------------------------------------
    // The bound is given as: gamma*log2(sigma sqrt(gamma))
    //--------------------------------------------------------------------------
    // temp = sigma*sqrt(gamma)
    SLIP_CHECK(SLIP_mpfr_mul_d(temp, temp, (double) sqrt(gamma),
        option->SLIP_MPFR_ROUND));
    // temp = log2(temp)
    SLIP_CHECK(SLIP_mpfr_log2(temp, temp, option->SLIP_MPFR_ROUND));
    // inner2 = temp
    double inner2;
    SLIP_CHECK(SLIP_mpfr_get_d(&inner2, temp, option->SLIP_MPFR_ROUND));
    // Free cache from log2. Even though mpfr_free_cache is called in
    // SLIP_LU_final(), it has to be called here to prevent memory leak in
    // some rare situations.
    SLIP_mpfr_free_cache();
    // bound = gamma * inner2+1. We add 1 to inner2 because log2(1) = 0
    int32_t bound = ceil(gamma*(inner2+1));
    // Ensure bound is at least 64 bit. In some rare cases the bound is very
    // small so we default to 64 bits.
    if (bound < 64) {bound = 64;}

    //--------------------------------------------------------------------------
    // find initial ordering using COLAMD and initialize row_perm and pinv
    //--------------------------------------------------------------------------
    // Declared as per COLAMD documentation
    int32_t Alen = 2*nz + 6 *(n+1) + 6*(n+1) + n;
    // allocate new vector Ai = A->i which will be modified in colamd
    int32_t* Ai = (int32_t*) SLIP_malloc(Alen* sizeof(int32_t));
    if (!Ai)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    // Initialize S->q as per COLAMD documentation
    for (i = 0; i < n; i++)
    {
        col_perm[i] = A->p[i];
        // initialize row_perm and pinv
        row_perm[i] = i;
        pinv[i] = i;
    }
    col_perm[n] = A->p[n];

    // Initialize Ai per COLAMD documentation
    for (i = 0; i < nz; i++)
    {
        Ai[i] = A->i[i];
    }
    int32_t stats [COLAMD_STATS];
    colamd(n, n, Alen, Ai, col_perm, (double *) NULL, stats);

#if 0
    printf("init col_perm: ");
    for (int32_t ii =0; ii<n; ii++)
    {
        printf("%d ", col_perm[ii]);
    }
    printf("\n");
#endif
    SLIP_FREE(Ai);

    // initialize L and U
    L->p[0] = 0;
    U->p[0] = 0;

    //--------------------------------------------------------------------------
    // try to find the k-th pivot column
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        //printf("***********************iteration %d(%d)**********************\n\n\n\n",k,n);

        //----------------------------------------------------------------------
        // Update the bit size estimation for candidate columns
        //----------------------------------------------------------------------
        // reset the values related to the minimum bit size
        int32_t mincol = -1;
        int32_t minrow = -1;
        int32_t cur_nz;//min_nz = 1;
        double min_w = SIZE_MAX, w;

        for (nthcand = 0; nthcand < ncand; nthcand++)
        {
            col = cands->columns[nthcand];
clock_t tic = clock();
            // if the nth candidate has not initialized or has been taken
            if (cands->col_index[nthcand] == -1)
            {
                // grap a new column from remaining list if there is any
                if (usedcol < n)
                {
                    j = col_perm[usedcol];
                    cands->col_index[nthcand] = j;
                    col->nz = 0;
                    col->unz = -1;
                    priority[nthcand] = 1;
                    int32_t cnz = 0;
                    for (i = A->p[j]; i < A->p[j+1]; i++)
                    {
                        // init col->x[cnz] if more mpz entries needed than
                        // previous allocated number
                        if (cnz >= col->max_mpz)
                        {
                            // allocate worse case bit size for each entry to
                            // avoid possible reallocation
                            if (SLIP_mpz_init2(col->x[cnz], bound) != SLIP_OK)
                            {
                                SLIP_MPZ_SET_NULL(col->x[cnz]);
                                col->max_mpz = cnz;
                                SLIP_FREE_WORKSPACE;
                                return SLIP_OUT_OF_MEMORY;
                            }
                        }

                        // set the value for col->x[cnz]=A[i][j]
                        SLIP_CHECK (SLIP_mpz_set(col->x[cnz], A->x[i]));

                        // get the bit size of entry A[i][j]
                        SLIP_mpz_sizeinbase(&size, A->x[i], 2);
                        col->bs[cnz] = size;

                        // set the row index for col
                        col->i[cnz] = A->i[i];

                        // set history vector
                        col->h[cnz] = -1;
                        cnz++;

                        // increase number of nonzero in current row
                        if (size != 0)
                        {
                            rc[A->i[i]] ++;
                        }
                    }
                    col->nz = cnz;
                    col->real_lnz = cnz;

                    // if new candidate column has less nz than previous
                    // allocated size (this would happen since the previous
                    // candidate column that was just selected for k-1 pivot
                    // column was only pretended to delete, the initialized mpz
                    // entries still exist)
                    if (cnz > col->max_mpz)
                    {
                        col->max_mpz = cnz;
                    }
                    usedcol++;
                }

                // otherwise, go for next candidate column
                else
                {
                    continue;
                }

            }
            else
            {
                priority[nthcand] ++;
            }

                clock_t toc = clock();
            init_time = init_time+(double)(toc-tic)/CLOCKS_PER_SEC;
            //printf("-------------------%d-th cand: A(:,%d) %d-----------------\n",nthcand,cands->col_index[nthcand], priority[nthcand]);

            //------------------------------------------------------------------
            // estimate number of digits for column cand
            //------------------------------------------------------------------
            tic = clock();
            // Since the bit size is exact as candidate for first col,
            // the estimation process is skipped.
            if (k != 0)
            {
                SLIP_CHECK(slip_triangular_est(col, x, xi, rc, pinv,
                    row_perm, k, L, rhos, rhos_bs, bound));
            }
            else
            {
                col->unz = 0;
            }
            toc = clock();
            est_time = est_time+(double)(toc-tic)/CLOCKS_PER_SEC;

            //------------------------------------------------------------------
            // check if the matrix is singular
            //------------------------------------------------------------------
            cur_nz = col->real_lnz;
            if (cur_nz == 0)
            {
                SLIP_FREE_WORKSPACE;
                return SLIP_SINGULAR;
            }
 
            // find the minimum from all entries that would be in L
            for (ci = col->unz; ci < col->nz; ci++)
            {
                if (rc[col->i[ci]] > ncand || rc[col->i[ci]] < 0)
                {
                    printf("rc too large!\n");
                    return SLIP_INCORRECT;
                }
                //w = col->bs[ci];
                //w = cur_nz*col->bs[ci];
                w = cur_nz*col->bs[ci]/((double)priority[nthcand]);
                //w = cur_nz*col->bs[ci]*rc[col->i[ci]];
                //w = cur_nz*col->bs[ci]*rc[col->i[ci]]/((double)priority[nthcand]);
                //printf("%lf(%d)>%lf ", w,ci,min_w);
                if (w > 0 && (w < min_w || mincol == -1))
                {
                    min_w = w;
                    mincol = nthcand;
                    minrow = ci;
                }
            }
        }

        //----------------------------------------------------------------------
        // add the selected candidate col to L and U
        //----------------------------------------------------------------------
        col = cands->columns[mincol];
        int32_t tmp, rpiv = col->i[minrow],
                cpiv = cands->col_index[mincol];
//printf("%s line %d using %d-th cand A(%d,%d) as pivot col\n",__FILE__,__LINE__,mincol,rpiv,cpiv);
//fprintf(stderr,"%s line %d using %d-th cand A(%d,%d) as pivot col\n",__FILE__,__LINE__,mincol,rpiv,cpiv);

        // use the column that has the minimum of the estimated digit size
        // to update L and U
        clock_t tic  = clock();
        SLIP_CHECK(slip_update_LU(L, U, rhos, rhos_bs, rc, col, rpiv, k));
        clock_t toc = clock();
        total_time = total_time+(double)(toc-tic)/CLOCKS_PER_SEC;
//printf("time used in REF update: %lf\n", total_time);

        //------------------------------------------------------------------
        // update row_perm, col_perm and pinv
        //------------------------------------------------------------------
        if (row_perm[k] != rpiv)
        {
            tmp = row_perm[k];
            row_perm[k] = rpiv;
            row_perm[pinv[rpiv]] = tmp;

            pinv[tmp] = pinv[rpiv];
            pinv[rpiv] = k;
        }

        if (col_perm[k] != cpiv)
        {
            tmp = col_perm[k];
            col_perm[k] = cpiv;
            for (i = k+1; i < k+ncand; i++)
            {
                if (col_perm[i] == cpiv)
                {
                    col_perm[i] = tmp;
                }
            }
        }

#if 0
        printf("row_perm: ");
        for (int32_t ii =0; ii<=k; ii++)
        {
            printf("%d(%d) ", row_perm[ii],ii);
        }
        printf("\n");
        printf("pinv: ");
        for (int32_t ii =0; ii<=k; ii++)
        {
            printf("%d(%d) ", pinv[ii],ii);
        }
        printf("\n");
        printf("col_perm: ");
        for (int32_t ii =0; ii<=k; ii++)
        {
            printf("%d(%d) ", col_perm[ii],ii);
        }
        printf("\n");
#endif

#if 0
        // keep the candidates in order to make sure the new candidate
        // appears at the end
        slip_column *pending_col = cands->columns[mincol];
        for (nthcand = mincol; nthcand < ncand-1; nthcand++)
        {
            if (cands->col_index[nthcand+1] != -1)
            {
                cands->columns[nthcand] = cands->columns[nthcand+1];
                cands->col_index[nthcand] = cands->col_index[nthcand+1];
                priority[nthcand] = priority[nthcand+1];
            }
            else
            {
                break;
            }
        }
        // pretend to delete the candidate column that has been used
        cands->col_index[nthcand] = -1;
        cands->columns[nthcand] = pending_col;
#endif
        // pretend to delete the candidate column that has been used
        cands->col_index[mincol] = -1;
    }
    SLIP_FREE_WORKSPACE;

    // Finalize L and U
    L->nz = L->p[n];
    U->nz = U->p[n];

    // shrink L and U
    slip_sparse_collapse(L);
    slip_sparse_collapse(U);

    //--------------------------------------------------------------------------
    // finalize the row indices in L and U
    //--------------------------------------------------------------------------

    // Permute entries in L
    for (i = 0; i < L->nz; i++)
    {
        L->i[i] = pinv[L->i[i]];
    }
    // Permute entries in U
    for (i = 0; i < U->nz; i++)
    {
        U->i[i] = pinv[U->i[i]];
    }
/*printf("time used in est: %lf\n", est_time);
printf("time used in init: %lf\n", init_time);
printf("time used in REF update: %lf\n", total_time);
*/    return SLIP_OK;
}
#undef CAND_SIZE
