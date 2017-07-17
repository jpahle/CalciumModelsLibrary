#include <Rcpp.h>
using namespace Rcpp;

// Declare variables that are necessary for model files and simulator
extern NumericVector calcium;
extern unsigned int ntimepoint;
extern double *amu;
extern unsigned long long int *x;

// Declare variables and functions that are necessary only for simulator
extern int nspecies;
extern int nreactions;

void calculate_amu();
void update_system(unsigned int rIndex);

// void calmodulin_calculate_amu();
// void calmodulin_update_system(unsigned int rIndex);
// void camkii_calculate_amu();
// void camkii_update_system(unsigned int rIndex);
// void pkc_calculate_amu();
// void pkc_update_system(unsigned int rIndex);