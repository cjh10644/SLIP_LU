#include "SLIP_LU_internal.h"

// This function delete the slip_cand_columns struct

void slip_delete_cand_columns
(   
    slip_cand_columns **M
)   
{
    if (M == NULL || (*M) == NULL) {return ;}
    
    for (int32_t i = 0; i < (*M)->n; i++)
    {
        slip_delete_column( &( (*M)->columns[i] ) );
    }
    SLIP_FREE((*M)->col_index);
    SLIP_FREE((*M));
}
