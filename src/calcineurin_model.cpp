#include <string>
#include <Rcpp.h>
using namespace Rcpp;



//********************************/* R EXPORT OPTIONS */********************************

// 1. USER INPUT for new models: Change value of the macro variable MODEL_NAME to the name of the new model.
#define MODEL_NAME calcineurin
// include the simulation function with macros (#define statements) that make it model specific (based on MODEL_NAME)
#include "simulator.cpp"
// Global variables
static std::map <std::string, double> model_params;
// 2. USER INPUT for new models: Change the name of the wrapper function to sim_<MODEL_NAME> and the names of the internally called functions to init_<MODEL_NAME> and simulator_<MODEL_NAME>.
//' Calcineurin Model R Wrapper Function (exported to R)
//'
//' This function calls the internal C++ simulator function to simulate the Calcineurin model (Fisher 2006). 
//' @param
//' @return
//' @examples
//' sim_calcineurin()
//' @export
// [[Rcpp::export]]
NumericMatrix sim_calcineurin(DataFrame user_input_df,
                   NumericVector user_sim_params,
                   List user_model_params) {

  // Provide default model parameters
  model_params = init_calcineurin();
  // Model Parameter Check: replace model parameters with user-supplied values if necessary
  List user_model_params_names = user_model_params.names();
  for (int i=0; i < user_model_params_names.size(); i++) {
    std::string current_name = user_model_params_names[i];
    // the map.find("key") returns map.end() if it doesn't find the key 
    if (model_params.find(current_name) != model_params.end()) {
      model_params[current_name] = user_model_params[current_name];     
    }
  }
  // Return result of the simulation
  return simulator_calcineurin(user_input_df,
                   user_sim_params,
                   user_model_params);
   
}



//********************************/* MODEL DEFINITION */********************************
// 3. USER INPUT for new models: define model parameters, number of species, 
// number of reactions, propensity equations and update_system function 

// Default model parameters
std::map <std::string, double> init() {
  // Model dimensions
  nspecies = 2;
  nreactions = 2;
  
  // Propensity equation parameters
  std::map <std::string, double> default_params;
  
  default_params["k_on"] = 1;
  default_params["k_off"] = 1;
  default_params["p"] = 3.0;
  
  return default_params;
}

// Propensity calculation:
// Calculates the propensities of all Calcineurin model reactions and stores them in the vector amu.
void calculate_amu() {
  
  // Look up model parameters in array 'model_params' initially
  double k_on = model_params["k_on"];
  double k_off = model_params["k_off"];
  double p = model_params["p"];  
  
  amu[0] = k_on * pow((double)calcium[ntimepoint],(double)p) * x[0];
  amu[1] = amu[0] + k_off * x[1];
}

// System update:
// Changes the system state (updates the particle numbers) by instantiating a chosen reaction.
void update_system(unsigned int rIndex) {
  switch (rIndex) {
  case 0:   // Forward
    x[0]--;
    x[1]++;
    break;
  case 1:   // Backward
    x[0]++;
    x[1]--;
    break;
    printf("\nError in updateSystem(): rIndex (%u) out of range!\n", rIndex);
    exit(-1);
  }
}