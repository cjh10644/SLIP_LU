//------------------------------------------------------------------------------
// SLIP_LU/slip_expand_double_mat: convert double matrix to mpz
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function converts a double matrix of size m*n to an appropriate
 * mpz array of size m*n. To do this, the number is multiplied by 10^17 then,
 * the GCD is found. This function allows the use of matrices in double
 * precision to work with SLIP LU
 *
 * See also slip_expand_double_array, which converts length n vector.
 */

#define SLIP_FREE_WORKSPACE              \
    SLIP_delete_mpfr_mat(&x3, m, n);\
    SLIP_MPZ_CLEAR(gcd);            \
    SLIP_MPZ_CLEAR(one);            \
    SLIP_MPQ_CLEAR(temp);

#include "SLIP_LU_internal.h"

SLIP_info slip_expand_double_mat
(
    mpz_t** x_out,// mpz mat
    double** x,   // double matrix that needs to be made integral
    mpq_t scale,  // scaling factor. x_out = scale*x
    int32_t m,    // number of rows of x
    int32_t n     // number of columns of x
)
{
    if (!x || !x_out || m <= 0 || n <= 0)
    {
        return SLIP_INCORRECT_INPUT;
    }
    int32_t i, j, k, p, r = 1;
    SLIP_info ok;
    mpfr_t **x3 = NULL;
    mpz_t gcd, one;
    SLIP_MPZ_SET_NULL(gcd); SLIP_MPZ_SET_NULL(one);
    mpq_t temp; SLIP_MPQ_SET_NULL(temp);
    SLIP_CHECK(slip_mpq_init(temp));
    SLIP_CHECK(slip_mpz_init(gcd));
    SLIP_CHECK(slip_mpz_init(one));
    SLIP_options *option = SLIP_create_default_options();
    if (!option)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    x3 = SLIP_create_mpfr_mat(m, n, option);
    SLIP_FREE(option);

    if (!x3)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    double expon = pow(10,17);                          // expon = 10^17
    SLIP_CHECK(slip_mpq_set_d(scale, expon));
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            // x3[i][j] = x[i][j]
            SLIP_CHECK(slip_mpfr_set_d(x3[i][j], x[i][j], SLIP_MPFR_ROUND));
            
            // x3[i][j] = x[i][j]*10^17
            SLIP_CHECK(slip_mpfr_mul_d(x3[i][j], x3[i][j], expon, SLIP_MPFR_ROUND));
            
            // x_out[i][j] = x3[i][j]
            SLIP_CHECK(slip_mpfr_get_z(x_out[i][j], x3[i][j], SLIP_MPFR_ROUND));
        }
    }
    //--------------------------------------------------------------------------
    // Compute the gcd to reduce the size of scale
    //--------------------------------------------------------------------------
    // Find an initial GCD 
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            SLIP_CHECK(slip_mpz_cmp_ui(&r, x_out[i][j], 0));
            if (r != 0)
            {
                k = i;
                p = j;
                SLIP_CHECK(slip_mpz_set(gcd, x_out[i][j]));;
                break;
            }
        }
        if (r != 0)
            break;
    }
            
    if (r == 0) // Entire matrix is zeros
    {
    	// TODO Handle this properly
        SLIP_FREE_WORKSPACE;
        return SLIP_INCORRECT_INPUT;
    }
    SLIP_CHECK(slip_mpz_set_ui(one, 1))
    // Compute the GCD of the numbers, stop if gcd == 1 (r == 0)
    for (i = k; i < m && r != 0; i++)
    {
        for (j = p; j < n && r != 0; j++)
        {
            SLIP_CHECK(slip_mpz_gcd(gcd, gcd, x_out[i][j]));
	    SLIP_CHECK(slip_mpz_cmp(&r, gcd, one));
        }
    }

    //--------------------------------------------------------------------------
    // Scale all entries to make as small as possible 
    //--------------------------------------------------------------------------
    if (r != 0)                 // If gcd == 1 stop
    {
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)    
            {
                SLIP_CHECK(slip_mpz_divexact(x_out[i][j], x_out[i][j], gcd));
            }
        }
        SLIP_CHECK(slip_mpq_set_z(temp, gcd));
        SLIP_CHECK(slip_mpq_div(scale, scale, temp));
    }
    SLIP_FREE_WORKSPACE;
    return SLIP_OK;
}

