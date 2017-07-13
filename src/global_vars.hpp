#include <Rcpp.h>
using namespace Rcpp;

/* Vars for simulator only */
// extern NumericVector timevector;
// extern double timestep;
// extern double vol;

/* Vars for simulator and model props files */
extern int nspecies;
extern int nreactions;

// Vars that change each simulation iteration
extern NumericVector calcium;
extern unsigned int ntimepoint;
extern double *amu;
extern unsigned long long int *x;
