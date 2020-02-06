#include"test.h"

typedef struct
{
    int32_t last_trial_bs;// column index that bs last updated, defaulted 0
    int32_t last_trial_x; // column index that x last updated, defaulted 0
    int32_t nz;           // number of nonzero in this column
    //int32_t nzmax;        // allocated size for the vectors
    int32_t *i;           // row indices
    size_t *bs;           // estimated bit size in i-th entry
    mpz_t *x;             // previously calculated to be the (last_trial)-th
                          // column of LU
} slip_column;

typedef struct
{
    slip_column **columns;  // columns of the matrix
    int32_t n;             // number of columns
} slip_columns_of_M;

#define FREE_WORKSPACE         \
{                              \
    slip_delete_columnsofM(&M);\
    SLIP_FREE(cands);          \
    SLIP_FREE(index);          \
    SLIP_FREE(bs);             \
    SLIP_FREE(range_stack);    \
    SLIP_FREE(row_perm);       \
    SLIP_FREE(col_perm);       \
    SLIP_FREE(xi);             \
    SLIP_FREE(x_int);          \
    SLIP_delete_mpz_array(&x_mpz, n);\
    SLIP_FREE(h);              \
    SLIP_FREE(Lb);             \
    SLIP_delete_mpz_array(&rhos, n);\
    SLIP_FREE(rhos_b);         \
}

SLIP_info SLIP_LU_analyze_and_factorize
(
    SLIP_sparse *L,
    SLIP_sparse *U,
    SLIP_sparse *A
)
{
    if (L == NULL || U == NULL || A == NULL)
    {
        return SLIP_INCORRECT_INPUT;
    }

    // initialize variables
    int32_t n = A->n, nz = A->nz;
    int32_t lnz = 0, unz = 0, bs_max_size = n, stack_max_size = n, sgn;
    int32_t nthMin = 0;
    bool foundColumnk = false;
    int32_t *cands = NULL, *range_stack = NULL, *index = NULL,
            *row_perm = NULL, *col_perm = NULL, *xi = NULL, *x_int = NULL,
            h = NULL, rhos_b = NULL;
    mpz_t *x_mpz = NULL, *rhos = NULL;
    int32_t ncand = 0, bs_size = 0;
    size_t size, *Lb = NULL, *bs = NULL;
    slip_column *col;

    // rebuild matrix A in column form
    slip_columns_of_M *M = slip_initialize_columnsofM(A);
    // the indices of candidate columns for next pivot column
    cands = (int32_t *) SLIP_malloc(ncand*sizeof(int32_t));

    // bs is an anrray of estimated bit size extracted from candidate columns.
    // the location of bs[i] in M can be obtained from index[i]. Specifically,
    // index[i]/n gives the column index of bs[i] in M and index[i]%n gives
    // the row index
    index = (int32_t *) SLIP_malloc(bs_max_size*sizeof(int32_t));
    bs    = (size_t *) SLIP_malloc(bs_max_size*sizeof(size_t));

    // range_stack={...,r2,r1} gives the range of indices of bs to be sorted
    // when calling slip_quicksort, specifically bs[r1:r2]. Upon returning from
    // slip_quicksort, range_stack will be updated to {...,r2',r1'}, then
    // bs[i], i = r1:r1'-1, are sorted in ascending order
    // To get bs fully sorted, initially range_stack = {bs_size-1, 0}, and
    // and iterate until range_stack becomes size of 1
    range_stack = (int32_t *) SLIP_malloc(stack_max_size*sizeof(int32_t));

    // row permutaion, i.e., inew = row_perm[i] is
    // the column index in L from A(:,i)
    row_perm = (int32_t *) SLIP_malloc(n*sizeof(int32_t));
    // column permutation
    col_perm = (int32_t *) SLIP_malloc(n*sizeof(int32_t));
    // nonzero pattern vector
    xi = (int32_t *) SLIP_malloc(2*n*sizeof(int32_t));
    // only used in slip_triangular_estimate as the bit size of current column
    x_int = (int32_t *) SLIP_malloc(n*sizeof(int32_t));
    // history vector
    h = (int32_t *) SLIP_malloc(n*sizeof(int32_t));
    // bit size of pivot elements
    rhos_b = (int32_t *) SLIP_malloc(n*sizeof(int32_t));

    L = slip_sparse_alloc(L, n, n, 10*nz);
    U = slip_sparse_alloc(U, n, n, 10*nz);
    // bit size of L->x
    Lb = (size_t *) SLIP_malloc(10*nz* sizeof(size_t));

    // initialize x_mpz refer to SLIP_LU_factorize.c
    x_mpz = slip_create_mpz_array2(n,bound);
    // pivot elements
    rhos = SLIP_create_mpz_array(n);

    if (NULL)
    {
        return SLIP_OUT_OF_MEMORY;
    }

    // try to find the k-th column
    for (int32_t k = 0; k < A->n; k++) // TODO k= 0:n-2?
    {
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
            SLIP_CHECK(slip_sparse_realloc(L));
        }
        if (unz + n > U->nzmax)
        {
            // Set U->nz = unz
            U->nz = unz;
            SLIP_CHECK(slip_sparse_realloc(U));
        }

        // find candidate columns using COLAMD
        SLIP_COLAMD(cands, &ncand, A);

        // reset bs_size
        bs_size = 0;

        for (int32_t nthcand = 0; nthcand < ncand; nthcand++)
        {
            col = M->columns[cands[nthcand]];

            // estimate number of digits for column cand
            if (k != 0) // the bit size is exact as candidate for first col
            {
                SLIP_CHECK(slip_triangular_estimate(col, x_int, xi, h, pinv,
                    row_perm, k, L, Lb, rhos));
            }

            // check if there are too much nonzero added in, reallocate for
            // index array if necessary
            if (bs_size+col->nz > bs_max_size)
            {
                index = (int32_t *) SLIP_realloc(index, bs_max_size,
                    2*bs_max_size);
                bs    = (int32_t *) SLIP_realloc(bs,    bs_max_size,
                    2*bs_max_size);

                if (index == NULL || bs == NULL)
                {
                    return SLIP_OUT_OF_MEMORY;
                }
                bs_max_size *= 2;
            }
            slip_get_new_bs(bs, index, &bs_size, col, pinv, cands[nthcand]);
        }

        range_stack[0] = bs_size-1;
        range_stack[1] = 0;
        n_stack = 2;
        // k-th pivot column has not been found
        foundColumnk = false;
        // try the column that has the n-th minimum of the estimated digit size
        nthMin = 0;

        while (!foundColumnk)
        {
            // try the column with smallest number of digits
            // slip_quicksort only has the first range_stack[n_stack-1] entries
            // in array well sorted
            if (nthMin >= range_stack[n_stack-1])
            {
                SLIP_CHECK(slip_quicksort(bs, index, &range_stack, n_stack,
                    &stack_max_size));
            }

            // get the column index in the original A
            q[k] = ordered_index[nthMin]/n;
            row_perm[k] = ordered_index[nthMin]%n;
            col = M->columns[q[k]];
            
            // only need to update for column index > 0
            if (k != 0)
            {
                SLIP_CHECK(slip_REF_triangular_update(col, x_mpz, xi, h, pinv,
                    row_perm, k, L, rhos));
            }

            // check if the selected pivot is zero
            SLIP_CHECK (slip_mpz_sgn (&sgn, col.x[row_perm[k]]));

            // TODO iterate thru all entries to make sure that
            //      we selected the smallest
            if (sgn != 0)
            {
                foundColumnk = true;
                // update the k-th column of L, Lb and U rhos
                // delete column[i]
                //--------------------------------------------------------------
                // Iterate accross the nonzeros in x
                //--------------------------------------------------------------
                for (j = 0; j < col->nz; j++)
                {
                    jnew = col->i[j];
                    // Location of x[j] in final matrix
                    loc = pinv[jnew];

                    //----------------------------------------------------------
                    // loc <= k are rows above k, thus go to U
                    //----------------------------------------------------------
                    if (loc <= k)
                    {
                        // Place the i location of the U->nz nonzero
                        U->i[unz] = jnew;
                        SLIP_CHECK(slip_mpz_sizeinbase(&size, col->x[jnew], 2));
                        // GMP manual: Allocated size should be size+2
                        SLIP_CHECK(slip_mpz_init2(U->x[unz], size+2));
                        // Place the x value of the U->nz nonzero
                        SLIP_CHECK(slip_mpz_set(U->x[unz], col->x[jnew]));
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
                        SLIP_CHECK(slip_mpz_sizeinbase(&size, col->x[jnew], 2));
                        // GMP manual: Allocated size should be size+2
                        SLIP_CHECK(slip_mpz_init2(L->x[lnz], size+2));
                        // Place the x value of the L->nz nonzero
                        SLIP_CHECK(slip_mpz_set(L->x[lnz], col->x[jnew]));
                        // update Lb
                        Lb[lnz] = size;
                        // Increment L->nz
                        lnz++;
                        if (loc == k)
                        {
                            // GMP manual: Allocated size should be size+2
                            SLIP_CHECK(slip_mpz_init2(rhos[k], size+2));
                            // Place the x value of the L->nz nonzero
                            SLIP_CHECK(slip_mpz_set(rhos[k], col->x[jnew]));
                        }
                    }
                }

            }

            // use the (n+1)-th minimum
            nthMin++;

            // if all columns are tested but none can be used,
            // then the matrix is singular
            if (nthMin >= bs_size)
            {
                FREE_WORKSPACE;
                return SLIP_SINGULAR;
            }
        }
    }

    FREE_WORKSPACE;
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
