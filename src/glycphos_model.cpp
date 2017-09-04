#include <string>
#include <Rcpp.h>
using namespace Rcpp;



//********************************/* R EXPORT OPTIONS */********************************

// 1. USER INPUT for new models: Change value of the macro variable MODEL_NAME to the name of the new model.
#define MODEL_NAME glycphos
// include the simulation function with macros (#define statements) that make it model specific (based on MODEL_NAME)
#include "simulator.cpp"
// Global variables
static std::map <std::string, double> model_params;
// 2. USER INPUT for new models: Change the name of the wrapper function to sim_<MODEL_NAME> and the names of the internally called functions to init_<MODEL_NAME> and simulator_<MODEL_NAME>.
//' Glycogen phosphorylase Model R Wrapper Function (exported to R)
//'
//' This function calls the internal C++ simulator function to simulate the glycogen phosphorylase model (Gall 2000). 
//' @param
//' @return
//' @examples
//' sim_glycphos()
//' @export
// [[Rcpp::export]]
NumericMatrix sim_glycphos(DataFrame user_input_df,
                   NumericVector user_sim_params,
                   List user_model_params) {

  // Provide default model parameters
  model_params = init_glycphos();
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
  return simulator_glycphos(user_input_df,
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
  
  default_params["VpM1"] = 1.5;
  default_params["VpM2"] = 0.6;
  default_params["alpha"] = 9;
  default_params["gamma"] = 9;
  default_params["K11"] = 0.1;
  default_params["Kp2"] = 0.2;
  // it is not necessary to convert glucose, Ka1, Ka2, Ka5 and Ka6 into particle numbers because the units cancel
  default_params["Ka1_conc"] = 1e-7;
  default_params["Ka2_conc"] = 1e-7;
  default_params["Ka5_conc"] = 500;
  default_params["Ka6_conc"] = 600;
  default_params["gluc_conc"] = 1e-7; // in Gall 2000 model fixed at 10mM
  
  return default_params;
}

// Propensity calculation:
// Calculates the propensities of all glycogen phosphorylase model reactions and stores them in the vector amu.
void calculate_amu() {
  
  // Look up model parameters in array 'double = model_params' initially
  double VpM1 = model_params["VpM1"];
  double VpM2 = model_params["VpM2"];
  double alpha = model_params["alpha"];
  double gamma = model_params["gamma"];
  double K11 = model_params["K11"];
  double Kp2 = model_params["Kp2"];
  double Ka1_conc = model_params["Ka1_conc"];
  double Ka2_conc = model_params["Ka2_conc"];
  double Ka5_conc = model_params["Ka5_conc"];
  double Ka6_conc = model_params["Ka6_conc"];
  double gluc_conc = model_params["gluc_conc"];
  
  double activeFraction = x[1]/(x[0]+x[1]);
  
  // divide VpM1 and VpM2 by 60 to convert the units from min^-1 to s^-1
  amu[0] = (VpM1 / 60.0 * (1.0 + gamma * pow((double)calcium[ntimepoint],4) / (pow(Ka5_conc,4) + pow((double)calcium[ntimepoint],4))) * ( 1.0 - activeFraction)) / ((K11 / (1.0 + pow((double)calcium[ntimepoint],4) / pow(Ka6_conc,4))) + 1.0 - activeFraction) * (x[0]+x[1]);
  amu[1] = amu[0] + ((VpM2 / 60.0 * (1.0 + alpha * (gluc_conc) / (Ka1_conc + gluc_conc)) * activeFraction) / (Kp2 / (1 + gluc_conc / Ka2_conc) + activeFraction) * (x[0]+x[1]));
}

// System update:
// Changes the system state (updates the particle numbers) by instantiating a chosen reaction.
void update_system(unsigned int rIndex) {
  switch (rIndex) {
  case 0:   // Forward: glycogen phosphorylase kinase, regulated by calcium
    x[0]--;
    x[1]++;
    break;
  case 1:   // Backward: phosphatase, regulated by glucose
    x[0]++;
    x[1]--;
    break;
    printf("\nError in updateSystem(): rIndex (%u) out of range!\n", rIndex);
    exit(-1);
  }
}