//------------------------------------------------------------------------------
// SLIP_LU/SLIP_solve_double: solve Ax=b, returning solution as double matrix
//------------------------------------------------------------------------------

// SLIP_LU: (c) 2019, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
// Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
// SLIP_LU/License for the license.

//------------------------------------------------------------------------------

/* This code utilizes the SLIP LU factorization. 
 * Soln is output as double matrix.
 *
 * Input/Output arguments:
 *
 * x_doub:  A double** array which is unitialized on input. On output, contains 
 *          the solution to the linear system Ax=b in double precision.
 *
 * A:       User's input matrix. It must be populated prior to calling this function
 *
 * S:       Symbolic analysis struct. It is important that this data structure
 *          is both allocated and populated prior to calling this function. 
 *          That is, this function MUST be preceded by a call to the function
 *          SLIP_LU_analyze.
 *
 * b:       Collection of right hand side vectors. Must be populated prior to factorization
 *
 * option:  Struct containing various command parameters for the factorization
 *
 */

#define SLIP_FREE_WORKSPACE                 \
    SLIP_delete_mpq_mat(&x_mpq, n, numRHS); \
    SLIP_delete_sparse(&L);                 \
    SLIP_delete_sparse(&U);                 \
    SLIP_FREE(pinv);                        \
    SLIP_delete_mpz_array(&rhos, n);

# include "SLIP_LU_internal.h"

SLIP_info SLIP_solve_double 
(   
    double **x_doub,        // Solution vector stored as an double
    SLIP_sparse *A,         // Compressed column form full precision matrix A
    SLIP_LU_analysis *S,    // Column ordering
    SLIP_dense *b,          // Right hand side vectrors
    SLIP_options *option    // Control parameters
)
{
    //-------------------------------------------------------------------------
    // Check input
    //-------------------------------------------------------------------------
    if (!x_doub || !A || !A->p || !A->i || !A->x ||
        !S || !S->q || !b || !b->x || !option)
    {
        return SLIP_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    /* Declare memory */
    //--------------------------------------------------------------------------
    int32_t *pinv, n = A->n, numRHS = b->n;
    SLIP_info ok;
    mpq_t **x_mpq = SLIP_create_mpq_mat(n, numRHS);
    SLIP_sparse* L = SLIP_create_sparse();
    SLIP_sparse* U = SLIP_create_sparse();
    pinv = (int32_t*) SLIP_malloc(n* sizeof(int32_t));    
    mpz_t* rhos = SLIP_create_mpz_array(n);
    if (!x_mpq || !L || !U || !pinv || !rhos)
    {
        SLIP_FREE_WORKSPACE;
        return SLIP_OUT_OF_MEMORY;
    }
    //--------------------------------------------------------------------------
    // LU factorization
    //--------------------------------------------------------------------------
    SLIP_CHECK(SLIP_LU_factorize(L, U, A, S, rhos, pinv, option));

    //--------------------------------------------------------------------------
    // FB Substitution
    //--------------------------------------------------------------------------
    SLIP_CHECK(SLIP_LU_solve(x_mpq, b, rhos, L, U, pinv));

    SLIP_CHECK(SLIP_permute_x(x_mpq, n, numRHS, S));

#if 0
    // check solution (debugging only)
    SLIP_CHECK(SLIP_check_solution(A, x_mpq, b));
#endif

    //--------------------------------------------------------------------------
    // Scale solution
    //--------------------------------------------------------------------------

    SLIP_CHECK(SLIP_scale_x(x_mpq, A, b));

    //--------------------------------------------------------------------------
    // Output, free memory
    //--------------------------------------------------------------------------
    SLIP_CHECK(SLIP_get_double_soln(x_doub, x_mpq, n, numRHS));

    SLIP_FREE_WORKSPACE;
    return ok;
}

