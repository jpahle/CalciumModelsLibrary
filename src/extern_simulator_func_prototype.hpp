#include <Rcpp.h>
using namespace Rcpp;


// Declare simulator function (prototyp)
extern NumericMatrix simulator(DataFrame param_input_df,
                   NumericVector param_sim_params,
                   List param_model_params);