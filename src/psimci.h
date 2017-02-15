#ifndef PSIMCI
#define PSIMCI

int get_omp_num_threads();

double sim_ols (
    const gsl_matrix * X,
    const gsl_vector * y
);

double sim_pctile (
    gsl_vector * x,
    double pctile
);

int sim_ci (
    const gsl_matrix * X,
    const gsl_vector * y,
    const double P,
    const int reps,
    gsl_vector * b,
    gsl_vector * mu
);

#endif
