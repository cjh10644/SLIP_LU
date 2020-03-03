#include "SLIP_LU_internal.h"


// This function initialize and create one column for slip_cand_columns
slip_column* slip_initialize_column
(   
    int32_t nrow
)   
{
    // nrow has been checked in slip_initialize_cand_columns

    // allocate memory for slip_column
    slip_column* col = (slip_column*) SLIP_malloc(sizeof(slip_column));
    if (col == NULL)
    {
        return col;
    }

    col->last_trial_x  = 0;
    col->last_trial_bs = 0;
    col->nz            = 0;
    col->nz_mpz        = 0;
    col->max_mpz       = 0;
    col->real_lnz      = 0;
    col->unz           =-1;

    col->i  = (int32_t*) SLIP_malloc(nrow * sizeof(int32_t));
    if (!col->i)  {return NULL;}
    col->bs = (size_t*) SLIP_malloc(nrow * sizeof(size_t));
    if (!col->bs)
    {
        SLIP_FREE(col->i);
        return NULL;
    }

    col->h  = (int32_t*) SLIP_malloc(nrow * sizeof(int32_t));;
    if (!col->h)
    {
        SLIP_FREE(col->i);
        SLIP_FREE(col->bs);
        return NULL;
    }

    // only allocate memory for the mpz_t array without initializing each entry.
    // Entries are initialized for index < max_mpz
    col->x  = (mpz_t*) SLIP_malloc(nrow * sizeof(mpz_t));
    if (!col->x)
    {
        SLIP_FREE(col->i);
        SLIP_FREE(col->bs);
        SLIP_FREE(col->h);
        return NULL;
    }

    return col;
}

