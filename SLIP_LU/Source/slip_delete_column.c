#include "SLIP_LU_internal.h"

// This function delete the slip_column struct and set the pointer to NULL

void slip_delete_column
(
    slip_column **col
)
{
    if (col  == NULL || (*col) == NULL) {return ;}
            

    SLIP_FREE((*col)->i);
    SLIP_FREE((*col)->bs);
    SLIP_delete_mpz_array(&( (*col)->x ), (*col)->max_mpz );
    SLIP_FREE((*col)->h);
    SLIP_FREE(*col);
}
