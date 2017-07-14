#include <Rcpp.h>
using namespace Rcpp;

// Declare variables that are necessary for _model and _sim cpp files
extern NumericVector calcium;
extern unsigned int ntimepoint;
extern double *amu;
extern unsigned long long int *x;

// Declare variables and functions that are necessary for _sim cpp files
extern int nspecies;
extern int nreactions;
void calculate_amu();
void update_system(unsigned int rIndex);