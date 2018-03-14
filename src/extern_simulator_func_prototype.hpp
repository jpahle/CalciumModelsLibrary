#include <Rcpp.h>
using namespace Rcpp;


// Declare simulator function (prototype)
extern DataFrame simulator(DataFrame user_input_df,
                           List user_sim_params,
                           NumericVector default_vols,
                           NumericVector default_init_conc);