#include <Rcpp.h>
using namespace Rcpp;


// Declare simulator function (prototyp)
extern NumericMatrix simulator(DataFrame param_input_df,
                   NumericVector param_sim_params,
                   double param_vol,
                   NumericVector param_init_conc);