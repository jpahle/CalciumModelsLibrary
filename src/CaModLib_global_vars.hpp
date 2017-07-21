#include <Rcpp.h>
using namespace Rcpp;


// Declare simulator function (prototyp)
extern NumericMatrix simulator(NumericVector param_time,
                   NumericVector param_calcium,
                   double param_timestep,
                   double param_vol,
                   NumericVector param_init_conc);

                   
// Declare subfunctions used by simulatior
// extern void calculate_amu();
// extern void update_system(unsigned int rIndex);