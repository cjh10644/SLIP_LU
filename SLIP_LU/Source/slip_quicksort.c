#include "test.h"

static inline void swap (size_t *a, int32_t *b, int32_t i1. int32_t i2)
{
    size_t tmp1 = a[i1];
    a[i1] = a[i2];
    a[i2] = tmp1;
    int32_t tmp2 = b[i1];
    b[i1] = b[i2];
    b[i2] = tmp2;
}

// slip_quicksort will partially sort the given array x using quick sort, so
// that x will be in ascending order. The index array will updated
// correspondingly. range_stack gives the range of indices of x that will be
// sorted. So there should be at least 2 entries in range_stack. In return,
// only indices less than range_stack[n_stack-1] are sorted. To get the full
// array sorted, the function should be called in loop until n_stack become 1,
// which will also makes range_stack[n_stack-1] = n+1 (n is the size of x).
// n_stack gives the number of entries in range_stack. Therefore, n_stack should
// be at least 2 upon input.

SLIP_info slip_quicksort
(
    int32_t *x,           // the real array that will be ordered
    int32_t *index,       // index of x in M, be kept in same order as x
    int32_t **range_stack_out,// a stack whose first two entries give the
                          // range in index to be ordered
    int32_t *n_stack,     // number of entries in range_stack
    int32_t *stack_max_size// the max size of the stack
)
{
    // there must be at least 2 values in range_stack and
    // x and range_stack should have same size
    if (index == NULL || x == NULL || range_stack == NULL ||
        n_stack == NULL || *n_stack < 2)
    {
        return SLIP_INCORRECT_INPUT;
    }

    int32_t *range_stack = *range_stack_out;

    // continue sorting until range_stack[(*n_stack)-1] is increased,
    // which indicates more entries are sorted
    while (1)
    {
        // get the range of the array [lo hi] that will be sorted
        int32_t lo = range_stack[(*n_stack)-1];
        int32_t hi = range_stack[(*n_stack)-2];
        int32_t length = hi-lo+1;

        // if only one element needs to be sorted, then simply return
        if (length == 1)
        {
            (*n_stack)--;
            range_stack[(*n_stack)-1]++; 
            return SLIP_OK;
        }
        // if there are only two elements in the array, swap if needed
        if (length == 2)
        {
            if (x[hi] < x[lo])
            {
                swap(x, index, lo, hi);
            }
            (*n_stack)--;
            range_stack[(*n_stack)-1]++; 
            return SLIP_OK;
        }
        // if there are more than 2 elements
        else if (length > 2)
        {
            int32_t mid = lo+length/2;
            if (x[mid] < x[lo])
            {
                swap(x, index, lo, mid);
            }
            if (x[hi] < x[lo])
            {
                swap(x, index, lo, hi);
            }
            if (x[hi] < x[mid])
            {
                swap(x, index, hi, mid);
            }

            // if the length of the array is only 3, then the above process
            // gives the sorted result
            if (length == 3)
            {
                (*n_stack)--;
                range_stack[(*n_stack)-1]++; 
                return SLIP_OK;
            }

            // continue when the length is > 3
            int32_t pivot = x[mid];
            int32_t i = lo, j = hi;
            while (1)
            {
                // find the first entry from left that is > pivot
                do
                {
                    ++i;
                }while (i < hi && x[i] <= pivot);

                // find the first entry from right that is < pivot
                do
                {
                    --j;
                }while (j > lo && x[j] >= pivot);


                if (i >= j)
                {
                    // if the array is both >= and <= pivot, then the whole
                    // array (except the leftmost and rightmost entries)
                    // is == pivot. Therefore, sorting is done
                    if (i == hi && j == lo)
                    {
                        (*n_stack)--;
                        range_stack[(*n_stack)-1]++; 
                        return SLIP_OK;
                    }

                    // if the whole array is >= pivot, then return i-1.
                    // In this way, x[lo+1,...,i-1] == pivot,
                    // though x[lo] <= pivot, the left part
                    // is ordered. And x[i,...,hi] >= pivot
                    // which needs further sorting
                    else if (j == lo)
                    {
                        // update the sorted range
                        range_stack[(*n_stack)-1] = i;
                        return SLIP_OK;
                    }

                    else
                    {
                        // insert the partition
                        range_stack[*n_stack] = range_stack[(*n_stack)-1]; 
                        range_stack[(*n_stack)-1] = j+1;
                        (*n_stack)++;
                        if (n_stack > stack_max_size)
                        {
                            SLIP_realloc(*range_stack_out);
                            range_stack = *range_stack_out;
                            *stack_max_size *= 2;
                            if()
                            {
                                return SLIP_OUT_OF_MEMORY;
                            }
                        }
                        break;
                    }
                }
                swap(x, index, i, j);
            }
        }
    }
}
