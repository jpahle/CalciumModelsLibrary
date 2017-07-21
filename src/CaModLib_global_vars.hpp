#include <Rcpp.h>
using namespace Rcpp;


// Declare shared variables
extern NumericVector calcium;
extern unsigned int ntimepoint;
extern double *amu;
extern unsigned long long int *x;
extern int nspecies;
extern int nreactions;


// Declare simulator function (prototyp)
extern NumericMatrix simulator(NumericVector param_time,
                   NumericVector param_calcium,
                   double param_timestep,
                   double param_vol,
                   NumericVector param_init_conc);

                   
// Declare subfunctions used by simulatior
extern void pkc_calculate_amu();
extern void pkc_update_system(unsigned int rIndex);