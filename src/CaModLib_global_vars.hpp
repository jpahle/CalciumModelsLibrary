#include "CalciumModelsLibrary_types.h"
#include <Rcpp.h>
using namespace Rcpp;


// Declare variables that are necessary for model files and simulator
extern NumericVector calcium;
extern unsigned int ntimepoint;
extern double *amu;
extern unsigned long long int *x;

extern int nspecies;
extern int nreactions;

// Declare variables and functions that are necessary only for simulator

// void calculate_amu();
// void update_system(unsigned int rIndex);

void pkc_init();
void pkc_calculate_amu();
void pkc_update_system(unsigned int rIndex);
XPtr<void (*)()> make_pkc_calculate_amu();