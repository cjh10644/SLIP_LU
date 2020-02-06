#include "test.h"

// This function delete the slip_columns_of_M struct

void slip_delete_columnsofM
(   
    slip_columns_of_M **M
)   
{
    if (M == NULL || (*M) == NULL) {return ;}
    
    for (int32_t i = 0; i < (*M)->n; i++)
    {
        slip_delete_column( &( (*M)->columns[i] ) );
    }
    SLIP_FREE((*M));
}
