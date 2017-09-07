#include <Rcpp.h>
using namespace Rcpp;


// Declare simulator function (prototyp)
extern NumericMatrix simulator(DataFrame user_input_df,
                   NumericVector user_sim_params,
                   List user_model_params,
                   std::map <std::string, double> model_params);