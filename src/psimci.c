/*********************************************************************
 * Program: psimci.c
 * Author:  Mauricio Caceres Bravo <caceres@nber.org>
 * Created: Sun Feb 12 19:28:43 EST 2017
 * Updated: Tue Feb 14 16:59:12 EST 2017
 * Purpose: Stata plugin to simulate a CI under H0: b = 0 for a
 *          treatment effect given a regression specification.
 * Note:    See stata.com/plugins for more on Stata plugins
 *********************************************************************/

/**
 * @file psimci.c
 * @author Mauricio Caceres bravo
 * @date 12 Mar 2017
 * @brief Stata plugin to simulate a CI for a placebo treatment.
 *
 * See the documentation for simci.ado (e.g. help simci from Stata)
 *
 * @see http://www.stata.com/plugins
 */

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_vector.h>
#include "psimci.h"
#include "stplugin.h"
#include "stutils.c"

/**
 * @brief Main function call will execute as Stata plugin
 *
 * The function takes the first variable in @argc as the dependent
 * variable and the next k - 1 variables are covariates. @argc contains
 * the comma options passed by Stata. Currently just the proportion
 * randomzed and the number of simulations to run. See documentation for
 * simci.ado or sim_ci below for moew.
 *
 * @param argc List of variables to use
 * @param argv Comma options from Stata
 * @return Modified variables in Stata
 * @see Documentation for simci.ado
 * @warning This is meant to be run from simci.ado and not by itself
 */
STDLL stata_call(int argc, char *argv[])
{

    // Initialize the variables to use
    ST_int      i, j ;
    ST_double   z ;
    ST_retcode  rc ;

    // Get P and number of reps. Note the 0-based indexing! So the
    // functiona ssumes P and reps were the 1st and 3nd argument.
    double P    = strtod (argv[0], NULL);
    int    reps = strtod (argv[1], NULL);

    const size_t n = SF_in2();
    const int    k = SF_nvars();

    // If too few variables (at least 2 for regressio), exit
    if (k < 2) {
        return (102) ;
    }

    // Initialize GSL elements where to store data
    gsl_matrix *X  = gsl_matrix_alloc (n, k + 1);
    gsl_vector *y  = gsl_vector_alloc (n);

    // Not sure if there is another way to read data vs the double loop.
    // Note: Careful with the 0-based indexing!
    for (i = SF_in1(); i <= SF_in2(); i++) {
        if (SF_ifobs(i)) {

            // Variables 2 through k are covariates
            for (j = 2; j <= k; j++) {
                // Note we leave the first column empty
                if ( (rc = SF_vdata(j, i, &z)) ) return(rc);
                gsl_matrix_set (X, i - 1, j - 1, z);
            }

            // Note we add the constant
            gsl_matrix_set (X,  i - 1, k, 1.0);

            // Variable 1 is the dependent variable
            if ( (rc = SF_vdata(1, i, &z)) ) return(rc);
            gsl_vector_set (y,  i - 1, z);
        }
    }

    // Now we call the simulation function and output the results into b, mu
    gsl_vector *b  = gsl_vector_alloc (reps);
    gsl_vector *mu = gsl_vector_alloc (reps);
    sim_ci (X, y, P, reps, b, mu);

    // Put the results into a local Stata variable
    char obuf[15], bchar[16 * reps], muchar[16 * reps];
    strcpy (bchar,  "");
    strcpy (muchar, "");
    for (int r = 0; r < reps; r++) {
        sprintf (obuf, " %15.9f", gsl_vector_get (b,  r));
        strcat  (bchar, obuf);
        sprintf (obuf, " %15.9f", gsl_vector_get (mu, r));
        strcat  (muchar, obuf);
    }

    // Note stata locals are an illusion---they are globals prepended
    // by _ that get destroyed when entering/exiting a given space.
    SF_macro_save ("_b",  bchar);
    SF_macro_save ("_mu", muchar);

    // Cleanup
    gsl_matrix_free (X);
    gsl_vector_free (y);
    gsl_vector_free (b);
    gsl_vector_free (mu);

    return (0);
}

/**
 * @brief Simulate a confidence interval given X, y
 *
 * The idea is to simulate a non-parametric CI based on placebo
 * assignments of a treatment variable. The program assigns
 * treatment at random, hence a null effect, to individuals or
 * clusters, optionally stratifying by any number of variables (or
 * the means thereof, in the case of clusters). Consider
 *
 *     Y_ij = a + b T_j + g X_ij + e_ij
 *
 * There are C = J choose PJ ways to treat the clusters (or C =
 * P choose PN in the case of individuals). If we computed b_ols
 * for c = 1, ..., C we would know the exact distribution of our
 * estimator, conditional on the data being representative of the
 * study data. C is typically intractably large, hence we simulate
 * K draws with sum(T_jk = PJ) and run
 *
 *     Y_ij = a + b_k T_jk + g X_ij + e_ij
 *
 * Let Fhat be the empirical cdf of b_k; a valid 1 - alpha CI for
 * b is given by
 *
 *     CI(1 - a) = [Fhat^-1(a / 2), Fhat^-1(1 - a / 2)]
 *
 * The function takes the @X as the covariate matrix, which
 * must have k + 1 columns with the first column free, @y as
 * the dependent variable, and outputs the results to @b, @mu
 *
 * @param X Covariate matrix with first column blank
 * @param y Dependent variable
 * @param P Proportion in treatment
 * @param reps Number of reprtitions
 * @param b Vector of length @reps; will output coefficients here
 * @param mu Vector of length @reps; will output control means here
 * @return Modified @b, @mu with coefficients and means
 * @see Documentation for simci.ado
 */
int sim_ci (const gsl_matrix * X,
            const gsl_vector * y,
            const double P,
            const int reps,
            gsl_vector * b,
            gsl_vector * mu)
{

    const size_t n  = X->size1;
    const int k     = X->size2;
    const int np    = ceil(n * P);
    const int nc    = n - np;
    double *sy      = malloc (sizeof(double));

    gsl_vector *ones = gsl_vector_alloc (n);
    gsl_vector_set_all (ones, 1.0);
    gsl_blas_ddot (ones, y, sy);

    // Set the random seed based on the time of day (seconds)
    srand (time(NULL));
    gsl_rng *rng = gsl_rng_alloc (gsl_rng_default);
    gsl_rng_set (rng, rand());

    // Get vector of 1s and 0s
    gsl_vector *T = gsl_vector_alloc (n);
    gsl_vector_set_zero (T);
    for (int i = 0; i < np; i++) {
        gsl_vector_set (T, i, 1.0);
    }

    // Initialize elements for parallel loop
    gsl_vector *Tp ;
    gsl_matrix *Xp ;
    int nloops ;
    double *sty ;

    // Get the number of threads available to OMP
    st_printf("Parallelizing simulation; %d threads found:\n",
              get_omp_num_threads());

    // Parallelize execution: Note We need a copy of Xp and Tp for each
    // thread since they will be modified at each iteration y does not
    // change, so it's shared.
    #pragma omp parallel private(Xp, Tp, nloops, sty) shared(y, b, sy)
    {
        nloops = 0;

        // Allocate to each therad their own copy
        Tp  = gsl_vector_alloc (n);
        Xp  = gsl_matrix_alloc (n, k);
        sty = malloc (sizeof(double));

        gsl_vector_memcpy (Tp, T);
        gsl_matrix_memcpy (Xp, X);

        // Parallel for loop through simulation
        #pragma omp for
        for (int r = 0; r < reps; r++) {
            // 1. Shuffle treatment
            // 2. Set as first column of covariate matrix
            // 3. Get mean of y over controls
            // 4. Store coefficient/mean
            // 5. Repeat 1-4
            // 6. ...
            // 7. Profit?
            gsl_ran_shuffle (rng, Tp->data, n, sizeof(size_t));
            gsl_matrix_set_col (Xp, 0, Tp);
            gsl_vector_set (b, r, sim_ols(Xp, y));
            gsl_blas_ddot (Tp, y, sty);
            gsl_vector_set (mu, r, (sy - sty) / nc);
            ++nloops;
        }

        // I want to print a pretty message saying how many iterations
        // each thread completed. Since threads finish on their own,
        // messages would be print at disparate times. However, one can
        // specify "critical" code which is executed only after all
        // threads are done running.
        #pragma omp critical
        {
            st_printf("\tThread %d performed %d simulations.\n",
                      omp_get_thread_num(), nloops);
        }

        // Cleanup
        gsl_matrix_free (Xp);
        gsl_vector_free (Tp);
    }

    // Cleanup
    gsl_vector_free (T);
    gsl_rng_free (rng);

    return (0);
}

/**
 * @brief Number of threads available to OMP
 *
 * Short wrapper to get number of threads available to OMP
 *
 * @return Number of threads available to OMP
 */
int get_omp_num_threads()
{
    int thread_id;
    int nthreads = 0;

    #pragma omp parallel private(thread_id) shared(nthreads)
    {
        thread_id = omp_get_thread_num();
        #pragma omp critical
        {
            nthreads = thread_id > nthreads? thread_id: nthreads;
        }
    }
    nthreads++;

    return (nthreads);
}

/**
 * @brief Wrapper to run a linear regression
 *
 * All I want is the first coefficient of a linear regression. For
 *
 *     Y = X beta
 *
 * I want (X' X)^-1 X' Y. GSL has solvers for a system of the form
 *
 *     Ax = b
 *
 * Where A is a symmetric matrix. Take A = X' X and b = X' y, then
 * we can use any number of routines to find x (especially since A
 * is now symmetric).
 *
 * @param X A n by k gsl matrix containing covariates.
 * @param y A n by 1 gsl vector containing the dependent variable
 * @return The first coefficient of a linear regression.
 * @warning This is meant to be run within the main loop of stata_call
 */
double sim_ols(const gsl_matrix * X, const gsl_vector * y)
{

    // Allocate memory to express the system as Ax = b
    gsl_matrix *A = gsl_matrix_alloc (X->size2, X->size2);
    gsl_vector *b = gsl_vector_alloc (X->size2);
    gsl_vector *x = gsl_vector_alloc (X->size2);

    // Set A = X' X and b = X' y
    gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, A);
    gsl_blas_dgemv (CblasTrans, 1.0, X, y, 0.0, b);

    // Cholesky decomposition
    gsl_linalg_cholesky_decomp1 (A);
    gsl_linalg_cholesky_solve (A, b, x);

    // You don't have to use Cholesky; a number of methods are available
    //
    // int s;
    // gsl_permutation * P = gsl_permutation_alloc (X->size2);
    // gsl_vector * tau    = gsl_vector_alloc (X->size2);
    //
    // Householder
    // gsl_linalg_HH_solve (A, b, x);
    //
    // LU decomposition
    // gsl_linalg_LU_decomp (A, P, &s);
    // gsl_linalg_LU_solve (A, P, b, x);
    // gsl_permutation_free (P);
    //
    // QR decomposition
    // gsl_linalg_QR_decomp (A, tau);
    // gsl_linalg_QR_solve (A, tau, b, x);
    // gsl_vector_free (tau);

    // Free up space
    gsl_matrix_free (A);
    gsl_vector_free (b);

    return (gsl_vector_get(x, 0));
}

/**
 * @brief Get pctile of a function
 *
 * Basic wrapper to get the @pctile percentile of a function.
 *
 * @param x n by 1 gsl vector whose percentile we want.
 * @param pctile Percentile
 * @return @pctile percentile of x
 */
double sim_pctile(gsl_vector * x, double pctile)
{
    gsl_sort_vector (x);
    int n = x->size;
    int i = floor(n * pctile);
    double qq = gsl_vector_get (x, i);
    if (i / n == pctile) {
        qq = (qq + gsl_vector_get (x, i + 1)) / 2;
    }
    return (qq);
}
