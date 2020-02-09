#include "SLIP_LU_internal.h"

slip_cand_columns *slip_initialize_cand_columns
(
    int32_t ncand,      // number of candidates
    int32_t nrow        // number of rows in each column
)
{
    if (ncand <= 0 || nrow <= 0)    {   return NULL;    }

    slip_cand_columns *M;
    M = (slip_cand_columns*) SLIP_malloc(sizeof(slip_cand_columns));

    if (M == NULL) {return M;}

    M->n = ncand;
    M->nrow = nrow;
    M->col_index = (int32_t*) SLIP_malloc(ncand * sizeof(int32_t));
    if (M->col_index == NULL)
    {
        SLIP_FREE(M);
        return NULL;
    }

    M->columns = (slip_column**) SLIP_malloc(ncand * sizeof(slip_column*));
    if (!M->columns)
    {
        SLIP_FREE(M->col_index);
        SLIP_FREE(M);
        return NULL;
    }

    // initialize vectors
    for (int32_t i = 0; i < ncand; i++)
    {
        // set to -1 to indicate no candidate selected
        M->col_index[i] = -1;
        // set each of the column pointers to be NULL in case that the memory
        // allocation for first few columns fails and slip_delete_column can
        // simply free NULL pointer(s) for the rest of columns
        M->columns[i] = NULL;
    }

    // allocate memeory for and initialize each column
    for (int32_t i = 0; i < ncand; i++)
    {
        // initialize the column size as nrow
        M->columns[i] = slip_initialize_column(nrow);

        // if there is no enough memory, free the whole matrix
        if (M->columns[i] == NULL)
        {
            slip_delete_cand_columns(&M);
            return NULL;
        }
    }

    return M;
}

