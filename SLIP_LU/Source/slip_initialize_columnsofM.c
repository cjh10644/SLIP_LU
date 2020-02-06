#include "test.h"

slip_columns_of_M slip_initialize_columnsofM
(
    SLIP_sparse *A
)
{
    if (A == NULL || A->n == 0 || A->m == 0)
    {
        return NULL;
    }

    slip_columns_of_M *M;
    M = (slip_columns_of_M*) SLIP_malloc(sizeof(slip_columns_of_M));

    if (M == NULL) {return M;}

    M->n = A->n;
    M->columns = (slip_column**) SLIP_malloc((M->n) * sizeof(slip_column*));
    if (!M->columns) {return NULL;}

    // set each of the column pointers to be NULL
    for (int32_t i = 0; i < M->n; i++)
    {
        M->columns[i] = NULL;
    }

    // allocate memeory for and initialize each column
    for (int32_t i = 0; i < M->n; i++)
    {
        // initialize the column size as twice size of nnz(A[:,i])
        int32_t col_nz = (A->p[i+1]-A->p[i])*2;
        M->columns[i] = slip_initialize_column(col_nz);

        // if there is no enough memory, free the whole matrix
        if (M->columns[i] == NULL)
        {
            slip_delete_columnsofM(&M);
            return NULL;
        }
    }

    // fill the bs array in column i with the bit size of each entry
    // and update the row index of each entry
    size_t size;
    for (int32_t i = 0; i < M->n; i++)
    {
        int32_t nz = 0;
        for (int32_t j = A->p[i]; j < A->p[i+1]; j++)
        {
            // get the bit size of entry A[j][i]
            slip_mpz_sizeinbase(&size, A->x[j], 2);
            M->columns[i]->i[nz] = A->i[j];
            M->columns[i]->bs[nz] = size;
            nz++;
        }
        M->columns[i]->nz = nz;
    }
    return M;
}

