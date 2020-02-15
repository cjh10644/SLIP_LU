//TODO
#define CAND_SIZE 10
#define SLIP_FREE_WORKSPACE         \
{                              \
    slip_delete_cand_columns(&cands);\
    SLIP_FREE(index);          \
    SLIP_FREE(bs);             \
    SLIP_FREE(range_stack);    \
    SLIP_FREE(row_perm);       \
    SLIP_FREE(xi);             \
    SLIP_FREE(x);          \
    SLIP_FREE(h);              \
    SLIP_FREE(Lb);             \
    SLIP_FREE(rhos_bs);         \
    SLIP_MPFR_CLEAR(temp);      \
}

#include"SLIP_LU_internal.h"

SLIP_info SLIP_LU_analyze_and_factorize
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
        A == NULL || S == NULL ||option == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // initialize variables and allocate memory
    //--------------------------------------------------------------------------
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
    int32_t sgn, i, j, k, jnew, loc;
    size_t size;
    slip_column *col;
    int32_t sigma_index;
    mpfr_t temp; SLIP_MPFR_SET_NULL(temp);

    // vectors to be used, corresponding usage can be found below
    int32_t *range_stack = NULL, *index = NULL, *xi = NULL, *h = NULL,
            *row_perm = NULL, *col_perm = S->q, *x = NULL;
    size_t *Lb = NULL, *bs = NULL, *rhos_bs = NULL;
    slip_cand_columns *cands;

    // candidate columns for next pivot column
    int32_t ncand = SLIP_MIN(CAND_SIZE, n);
    cands = slip_initialize_cand_columns(ncand, n);

    // bs is an array of estimated bit size extracted from candidate columns.
    // the location of bs[i] in cands can be obtained from index[i].
    // Specifically, index[i]/n gives the column index of bs[i] in cands and
    // index[i]%n gives the row index
    int32_t bs_max_size = n, bs_nz = 0;
    index = (int32_t *) SLIP_malloc(bs_max_size*sizeof(int32_t));
    bs    = (size_t *)  SLIP_malloc(bs_max_size*sizeof(size_t));

    // range_stack={...,r1,r0} gives the range of indices of bs to be sorted
    // when calling slip_quicksort, specifically bs[r0:r1]. Upon returning from
    // slip_quicksort, range_stack will be updated to {...,r1',r0'}, then
    // bs[i], i = r0:r0'-1, are sorted in ascending order
    // To get bs fully sorted, initially range_stack = {bs_nz-1, 0}, and
    // and iterate until range_stack becomes size of 1
    int32_t stack_max_size = n, stack_nz = 0;
    range_stack = (int32_t *) SLIP_malloc(stack_max_size*sizeof(int32_t));

    // row permutaion, i.e., inew = row_perm[i] is
    // the row index in A that will be permuted to L(i,:)
    row_perm = (int32_t *) SLIP_malloc(n*sizeof(int32_t));
    // inverse of row permutaion, i.e., i = pinv[inew] is
    // the row index in L from A(inew,:)
    // pinv     = (int32_t *) SLIP_malloc(n*sizeof(int32_t));
    // column permutation, which should be size of n+1 to be used in COLAMD
    // col_perm = (int32_t *) SLIP_malloc((n+1)*sizeof(int32_t));
    // nonzero pattern vector
    xi = (int32_t *) SLIP_malloc(2*n*sizeof(int32_t));
    // used in slip_triangular_estimate and slip_REF_triangular_update as the
    // pointers to corresponding entries in slip_column struct. That is, x[i]
    // indicates nonzero at i-th row, while its value is the x[i]-th entry in
    // slip_column struct
    x = (int32_t *) SLIP_malloc(n*sizeof(int32_t));
    // history vector
    h = (int32_t *) SLIP_malloc(n*sizeof(int32_t));

    // only allocate memory for pivot elements here without initializing each
    // entry. Entries are initialized and set using mpz_init_set when the pivot
    // is found and calculated. Using calloc instead of malloc is to avoid error
    // caused by clearing uninitialized mpz_t value with SLIP_MPZ_CLEAR
    // rhos    = (mpz_t*)   SLIP_calloc(n, sizeof(mpz_t));
    // bit size of pivot elements
    rhos_bs = (size_t *) SLIP_malloc(n*sizeof(size_t));

    if (!cands || !index || !bs || !range_stack || !row_perm || !pinv ||
        !col_perm || !xi || !h || !rhos_bs )
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }

    // initialize each entry of x_bs and h as required by
    // slip_triangular_estimate and slip_REF_triangular_update
    for (i = 0; i < n; i++)
    {
        x[i] = -1;
        h[i] = -1;
    }

    // Assume the number of nonzeros in L or U is 10 times of that in A
    // only allocate memory for L->x and U->x but not initializing each entry
    SLIP_CHECK(slip_sparse_alloc2(L, n, n, 10*nz));
    SLIP_CHECK(slip_sparse_alloc2(U, n, n, 10*nz));

    // bit size of L->x
    Lb = (size_t *) SLIP_malloc(10*nz* sizeof(size_t));

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

    //--------------------------------------------------------------------------
    // try to find the k-th pivot column
    //--------------------------------------------------------------------------
    for (k = 0; k < n; k++)
    {
        //printf("***********************iteration %d(%d)**********************\n\n\n\n",k,n);
        // Column pointers for column k of L and U
        L->p[k] = lnz;
        U->p[k] = unz;

        //----------------------------------------------------------------------
        // Reallocate memory if necessary
        //----------------------------------------------------------------------
        if (lnz + n > L->nzmax)
        {
            // Set L->nz = lnz
            L->nz = lnz;
            Lb = (size_t *) SLIP_realloc(Lb, L->nzmax*sizeof(size_t),
                2*(L->nzmax)*sizeof(size_t));
            if (!Lb)
            {
                SLIP_FREE_WORKSPACE;
                return SLIP_OUT_OF_MEMORY;
            }
            SLIP_CHECK(slip_sparse_realloc(L));
        }
        if (unz + n > U->nzmax)
        {
            // Set U->nz = unz
            U->nz = unz;
            SLIP_CHECK(slip_sparse_realloc(U));
        }

        //----------------------------------------------------------------------
        // Update the bit size estimation for candidate columns
        //----------------------------------------------------------------------
        // reset bs_nz
        bs_nz = 0;

        for (int32_t nthcand = 0; nthcand < ncand; nthcand++)
        {
            col = cands->columns[nthcand];

            // if the nth candidate has not initialized or has been taken
            if (cands->col_index[nthcand] == -1)
            {
                // grap a new column from remaining list if there is any
                if (usedcol < n)
                {
                    j = col_perm[usedcol];
                    cands->col_index[nthcand] = j;
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
                        cnz++;
                    }
                    col->nz = cnz;
                    col->nz_mpz = cnz;

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
        //printf("-------------------%d-th cand: A(:,%d)-----------------\n",nthcand,cands->col_index[nthcand]);

            // Since the bit size is exact as candidate for first col,
            // the estimation process is skipped.
            if (k != 0)
            {
                // estimate number of digits for column cand
                SLIP_CHECK(slip_triangular_estimate(col, x, xi, h, pinv,
                    row_perm, k, L, Lb, rhos_bs));

                // allocate memory for newly found reachable entries if there
                // is not enough allocated mpz_t values in col->x
                if (col->max_mpz < col->nz)
                {
                    for (i = col->max_mpz; i < col->nz; i++)
                    {
                        // allocate worse case bit size for each entry to
                        // avoid possible reallocation
                        if (SLIP_mpz_init2(col->x[i], bound) != SLIP_OK)
                        {
                            SLIP_MPZ_SET_NULL(col->x[i]);
                            // update the number of mpz allocated to help delete
                            // mpz array
                            col->max_mpz = i;
                            SLIP_FREE_WORKSPACE;
                            return SLIP_OUT_OF_MEMORY;
                        }
                    }
                    // update the max number of mpz allocated in col
                    col->max_mpz = col->nz;
                }
            }

            // check if there are too much nonzero added in, reallocate for
            // index array if necessary
            if (bs_nz+col->nz > bs_max_size)
            {
                index = (int32_t *) SLIP_realloc(index,
                      bs_max_size*sizeof(int32_t),
                    2*bs_max_size*sizeof(int32_t));
                bs    = (size_t *) SLIP_realloc(bs,
                      bs_max_size*sizeof(size_t),
                    2*bs_max_size*sizeof(size_t));

                if (index == NULL || bs == NULL)
                {
                    SLIP_FREE_WORKSPACE;
                    return SLIP_OUT_OF_MEMORY;
                }
                bs_max_size *= 2;
            }
            slip_get_new_bs(bs, index, &bs_nz, col, pinv, nthcand*n);
        }

        //----------------------------------------------------------------------
        // try the column that has the n-th minimum of the estimated digit size
        //----------------------------------------------------------------------
        nthMin = 0;
        range_stack[0] = bs_nz-1;
        range_stack[1] = 0;
        stack_nz = 2;
        // k-th pivot column has not been found
        foundPivotColumn = false;

        while (!foundPivotColumn)
        {
            // try the column with smallest number of digits
            // slip_quicksort only has the first range_stack[stack_nz-1] entries
            // in array well sorted
            if (nthMin >= range_stack[stack_nz-1])
            {
                SLIP_CHECK(slip_quicksort(bs, index, &range_stack, &stack_nz,
                    &stack_max_size));
            }

            // get the column index in the candidate list
            int32_t mincol = index[nthMin]/n;
            int32_t minrow = index[nthMin]%n;
            col = cands->columns[mincol];

            // only need to update for column index > 0
            if (k != 0)// TODO is this necessray?
            {
                SLIP_CHECK(slip_REF_triangular_update(col, x, xi, h, pinv,
                    row_perm, k, L, rhos));
            }

            // check if the selected pivot is zero
            SLIP_CHECK (SLIP_mpz_sgn (&sgn, col->x[minrow]));

            if (sgn != 0)
            {
                // iterate thru all entries below k to make sure that we
                // selected the smallest nonzero
                int32_t min_index = minrow;
                for (i = 0; i < col->nz_mpz; i++)
                {
                    if(pinv[col->i[i]] >= k &&
                       mpz_cmpabs(col->x[min_index],col->x[i]) > 0)
                    {
                        SLIP_CHECK (SLIP_mpz_sgn (&sgn, col->x[i]));
                        if (sgn != 0) {min_index = i;}
                    }
                }
                /*if (min_index != minrow)
                {
                    SLIP_gmp_printf("new minium found! %d~?(%Zd) index=%d\n",
                        col->i[min_index], col->x[min_index],min_index);
                    SLIP_gmp_printf("                  %d~?(%Zd) index=%d\n",
                        col->i[minrow], col->x[minrow],minrow);
                }*/
                // TODO
                //minrow = min_index;

                foundPivotColumn = true;
                //--------------------------------------------------------------
                // update row_perm, col_perm and pinv
                //--------------------------------------------------------------
                int32_t tmp, rpiv = col->i[minrow],
                        cpiv = cands->col_index[mincol];
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
                for (int32_t ii =0; ii<n; ii++)
                {
                    printf("%d ", row_perm[ii]);
                }
                printf("\n");
                printf("pinv: ");
                for (int32_t ii =0; ii<n; ii++)
                {
                    printf("%d ", pinv[ii]);
                }
                printf("\n");
                printf("col_perm: ");
                for (int32_t ii =0; ii<n; ii++)
                {
                    printf("%d ", col_perm[ii]);
                }
                printf("\n");
#endif

                //--------------------------------------------------------------
                // update the k-th column of L, Lb and U rhos
                //--------------------------------------------------------------
                // Iterate accross the nonzeros in col->x
                for (i = n-col->nz_mpz; i < n; i++)//TODO FIXME
                {
                    // j is the index of entry in col, while jnew = col->i[j]
                    // indicates this entry is at jnew-th row
                    if (k != 0)
                    {
                        // if k!=0, nonzero pattern is initialized and ordered
                        // based on row permutation, so we iterate thru nonzero
                        // in the known row permutation order using xi
                        jnew = xi[i];           // row index
                        j = x[jnew];            // index in col
                    }
                    else
                    {
                        // otherwise, nonzero pattern is not initialized, so
                        // we simply iterate in the same order as they show in
                        // the slip_column struct
                        j = i-n+col->nz_mpz;   // index in col
                        jnew = col->i[j];      // row index
                    }
                    // Location of x[j] in final matrix
                    loc = pinv[jnew];

                    //----------------------------------------------------------
                    // loc <= k are rows above k, thus go to U
                    //----------------------------------------------------------
                    if (loc <= k)
                    {
                        // Place the i location of the U->nz nonzero
                        U->i[unz] = jnew;
                        // Place the x value of the U->nz nonzero
                        SLIP_CHECK(SLIP_mpz_init_set(U->x[unz], col->x[j]));
                        // Increment U->nz
                        unz++;
                    }

                    //----------------------------------------------------------
                    // loc >= k are rows below k, thus go to L
                    //----------------------------------------------------------
                    if (loc >= k)
                    {
                        // Place the i location of the L->nz nonzero
                        L->i[lnz] = jnew;
                        // Place the x value of the L->nz nonzero
                        SLIP_CHECK(SLIP_mpz_init_set(L->x[lnz], col->x[j]));
                        // update Lb
                        SLIP_CHECK(SLIP_mpz_sizeinbase(&size, col->x[j], 2));
                        Lb[lnz] = size;
                        // Increment L->nz
                        lnz++;
                        if (loc == k)
                        {
                            //SLIP_gmp_printf("rhos[%d]=%Zd,rpiv=%d(%d)(%d)\n",k,col->x[j], rpiv,col->i[minrow],col->i[j]);
                            // set rhos[k] = L(k,k)
                            SLIP_CHECK(SLIP_mpz_set(rhos[k],col->x[j]));
                            rhos_bs[k] = size;
                        }
                    }
                }

                //--------------------------------------------------------------
                // pretend to delete the candidate column that has been used
                //--------------------------------------------------------------
                cands->col_index[mincol] = -1;
                col->nz = 0;
                col->nz_mpz = 0;
                col->last_trial_x = 0;
                col->last_trial_bs = 0;
            }

            // use the (n+1)-th minimum
            nthMin++;

            // if all columns are tested but none can be used,
            // then the matrix is singular
            if (nthMin >= bs_nz && !foundPivotColumn)
            {
                printf("%s line %d\n",__FILE__,__LINE__);

                for (int32_t nthcand = 0; nthcand < ncand; nthcand++)
                {
                    col = cands->columns[nthcand];

                    if (cands->col_index[nthcand] != -1)
                    {
                        printf("k=%d,last trial=%d, A(:,%d)\n",k,col->last_trial_x,cands->col_index[nthcand]);
                        for (i = 0; i < col->nz_mpz; i++)
                        {
                            if(pinv[col->i[i]] >= k)
                            SLIP_gmp_printf("%Zd ",col->x[i]);
                            else
                            SLIP_gmp_printf("(%d)%Zd ",pinv[col->i[i]],col->x[i]);
                        }
                        printf("\n");
                    }
                }

                SLIP_FREE_WORKSPACE;
                return SLIP_SINGULAR;
            }
        }
    }

    SLIP_FREE_WORKSPACE;
    // Finalize L and U
    L->nz = lnz;
    U->nz = unz;
    L->p[n] = lnz;
    U->p[n] = unz;

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

    return SLIP_OK;
}
#undef CAND_SIZE
