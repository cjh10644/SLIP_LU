#include "test.h"

// This function updates the array of the bit size and their indices in
// candidate columns that need to be ordered and will be used to determine the
// best pivot column. We only need to compare the entries below k-th row since
// the pivots for the first (k-1) rows have been found
void slip_get_new_bs
(
    size_t *bs,                // the array of estimated bit size to be ordered
    int32_t *index,             // indices of the entries in col to be ordered
    int32_t *index_size,        // current number of the found indices
    slip_column *col,           // the candidate column
    int32_t *pinv,              // inverse of row permutation
    int32_t cand                // the index of candidate column in original A
)
{
    int32_t k = col->last_trial_bs, nz = col->nz, base = cand*n;
    int32_t *Ai = col->i, *Axb = col->bs;
    for (int32_t i = 0; i < nz; i++)
    {
        if (pinv[Ai[i]] >= k)
        {
            bs[(*index_size)] = Axb[i];
            index[(*index_size)++] = base+i;
        }
    }
}
