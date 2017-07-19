#include "CalciumModelsLibrary_types.h"
#include <Rcpp.h>
using namespace Rcpp;


// Declare variables for model and simulator functions
extern NumericVector calcium;
extern unsigned int ntimepoint;
extern double *amu;
extern unsigned long long int *x;

extern int nspecies;
extern int nreactions;

// Declare functions for simulator function
void pkc_init();
void pkc_calculate_amu();
void pkc_update_system(unsigned int rIndex);

// Declare functions for R script calling the simulation function
R_init_ptr make_pkc_init();
R_amu_ptr make_pkc_calculate_amu();
R_stM_ptr make_pkc_update_system();