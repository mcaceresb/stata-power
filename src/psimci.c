/*********************************************************************
 * Program: psimci.c
 * Author:  Mauricio Caceres Bravo <caceres@nber.org>
 * Created: Sun Feb 12 19:28:43 EST 2017
 * Updated: Sun Feb 12 19:28:43 EST 2017
 * Purpose: Stata plugin to simulate a CI under H0: b = 0 for a
 *          treatment effect given a regression specification.
 * Note:    See stata.com/plugins for more on Stata plugins
 *********************************************************************/

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_vector.h>
#include "stplugin.h"

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

double sim_ols(const gsl_matrix * X, const gsl_vector * y);
double pctile(gsl_vector * x, double pctile);

/**
 * @brief Main function call will execute as Stata plugin
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
 * The function takes the first variable in @argc as thep dependent
 * variable and the next k - 2 variables are covariates. k - 1 is where
 * the program will store the coefficients and k is where the program
 * will store the simulated control means.
 *
 * @argc contains the comma options passed by Stata. Currently just the
 * proportion randomzed and the number of simulations to run.
 *
 * @param argc List of variables to use
 * @param argv Comma options from Stata
 * @return Modified variables in Stata
 * @see Documentation for simci.ado
 * @warning This is meant to be run from simci.ado and not by itself
 */
STDLL stata_call(int argc, char *argv[])
{

    // Set up
    // -------

    // Initialize ALL THE THINGS!
    ST_int      i, j ;
    ST_double   z ;
    ST_retcode  rc ;

    double *mu;
    int    r, thread_id, nloops;
    char   buf[72], pbuf[72] ;

    double P    = strtod (argv[0], NULL); // Proportion randomized
    double reps = strtod (argv[1], NULL); // Number of repetitions

    // Since communication with Stata is primitive, the last 2 variables
    // specified is where we store the output. Since the number of
    // repetitions can be > than the number of observations, we get
    // the number of relevant observations, SF_in2, and the number of
    // variables sans the last 2.
    const size_t n = SF_in2();
    const int    k = SF_nvars() - 2;

    int np = ceil(n * P);
    int R  = floor(reps);

    // Set the random seed based on the time of day (seconds)
    srand(time(NULL));
    gsl_rng *rng = gsl_rng_alloc (gsl_rng_default);
    gsl_rng_set (rng, rand());

    // If too few variables (at least 2 for regressio), exit
    if (k < 2) {
        return (102) ;
    }

    // Initialize GSL elements
    gsl_matrix *X  = gsl_matrix_alloc (n, k + 1);
    gsl_vector *y  = gsl_vector_alloc (n);
    gsl_vector *T  = gsl_vector_alloc (n);

    // Get the data
    // ------------

    // Not sure if there is another way than the double loop
    for(i = SF_in1(); i <= SF_in2(); i++) {
        if (SF_ifobs(i)) {
            for(j = 2; j <= k; j++) {
                if ( (rc = SF_vdata(j, i, &z)) ) return(rc);
                gsl_matrix_set (X,  i - 1, j - 1, z);
            }
            if ( (rc = SF_vdata(1, i, &z)) ) return(rc);
            gsl_vector_set (y,  i - 1, z);
            gsl_matrix_set (X,  i - 1, k, 1.0);
            gsl_vector_set (T, i - 1, i <= np? 1 : 0);
        }
    }

    // Parallelize ALL THE THINGS!
    // ---------------------------

    // Get the number of threads
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
    sprintf(buf, "Parallelizing simulation; %d threads found:\n",
            nthreads);
    SF_display (buf);

    // Initialize variables for parallel loop execution
    gsl_vector *Tp ;
    gsl_matrix *Xp ;

    // Initialize parallel environment: Note We need a copy of Xp and
    // Tp for each thread since they will be modified at each iteration
    // y does not change, so it's shared.

    // Note: To store in a vector, you would first initialize:
    //
    //     gsl_vector *b = gsl_vector_alloc (R);
    //
    // This would be shared among all threads. Now inside the loop
    //
    //     gsl_vector_set(b, r, sim_ols(Xp, y));
    //
    // The CI is given by: pctile(b, 0.025), pctile(b, 0.975)

    #pragma omp parallel private(Xp, Tp, thread_id, nloops, pbuf, mu) shared(y)
    {
        thread_id = omp_get_thread_num();
        nloops    = 0;

        // Allocate memory for each object in each thread
        mu = malloc (sizeof(double));
        Tp = gsl_vector_alloc (n);
        Xp = gsl_matrix_alloc (n, k + 1);

        // Copy treatment vector and main matrix into each thread
        gsl_vector_memcpy (Tp, T);
        gsl_matrix_memcpy (Xp, X);

        // Parallel for loop!
        #pragma omp for
        for(r = 0; r < R; r++) {
            // 1. Shuffle treatment
            // 2. Set as first column of covariate matrix
            // 3. Get mean of y over controls
            // 4. Store coefficient in Stata
            // 5. Score mean in stata
            // 6. Repeat 1-5
            // 7. ...
            // 8. Profit?
            gsl_ran_shuffle (rng, Tp->data, n, sizeof(size_t));
            gsl_matrix_set_col (Xp, 0, Tp);
            gsl_blas_ddot (Tp, y, mu);
            // gsl_blas_ddot(1.0 - Tp, y, mu);
            rc  = SF_vstore (k + 1, r + 1, sim_ols(Xp, y));
            rc  = SF_vstore (k + 2, r + 1, *mu);
            ++nloops;
        }

        gsl_matrix_free (Xp);
        gsl_vector_free (Tp);

        // I want to print a pretty message saying how many iterations
        // each thread completed. Since threads finish on their own,
        // messages would be print at disparate times. However, one can
        // specify "critical" code which is executed only after all
        // threads are done running.
        #pragma omp critical
        {
            sprintf(pbuf, "\tThread %d performed %d simulations.\n",
                    thread_id, nloops);
            SF_display (pbuf);
        }
    }

    // Free memory
    gsl_matrix_free (X);
    gsl_vector_free (T);
    gsl_vector_free (y);
    gsl_rng_free (rng);

    return (0);
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

    // You don't have to use Cholesky; a number of methods are available
    //
    // int s;
    // gsl_permutation * P = gsl_permutation_alloc (X->size2);
    // gsl_vector * tau    = gsl_vector_alloc (X->size2);

    // Householder
    // gsl_linalg_HH_solve (A, b, x);

    // LU decomposition
    // gsl_linalg_LU_decomp (A, P, &s);
    // gsl_linalg_LU_solve (A, P, b, x);
    // gsl_permutation_free (P);

    // QR decomposition
    // gsl_linalg_QR_decomp (A, tau);
    // gsl_linalg_QR_solve (A, tau, b, x);
    // gsl_vector_free (tau);

    // Cholesky decomposition
    gsl_linalg_cholesky_decomp1 (A);
    gsl_linalg_cholesky_solve (A, b, x);

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
double pctile(gsl_vector * x, double pctile)
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
