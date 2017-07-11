#include <Rcpp.h>
using namespace Rcpp;

// Global variables that are used by the general simulator and the model properties file
extern NumericVector timevector;
extern NumericVector calcium;
extern double timestep;
extern double vol;

extern unsigned int ntimepoint;
extern double *amu;
extern unsigned long long int *x;

extern int nspecies;
extern int nreactions;