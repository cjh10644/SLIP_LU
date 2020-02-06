#include "test.h"


// This function initialize and create one column for slip_column_of_M
slip_column* slip_initialize_column
(   
    int32_t nzmax
)   
{   
    slip_column* col = (slip_column*) SLIP_malloc(sizeof(slip_column));;
    col->last_trial_x  = 0;
    col->last_trial_bs = 0;
    col->nz            = 0;
    col->nzmax         = nzmax;

    col->i  = (int32_t*) SLIP_calloc(nzmax, sizeof(int32_t));
    if (!col->i)  {return NULL;}
    col->bs = (int32_t*) SLIP_calloc(nzmax, sizeof(int32_t));
    if (!col->bs)
    {
        SLIP_FREE(col->i);
        return NULL;
    }
    col->x  = SLIP_create_mpz_array(nzmax);
    if (!col->x)
    {
        SLIP_FREE(col->i);
        SLIP_FREE(col->bs);
        return NULL;
    }
    return col;
}

