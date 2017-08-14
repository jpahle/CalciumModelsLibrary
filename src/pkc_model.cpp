#include <Rcpp.h>
using namespace Rcpp;



//********************************/* R EXPORT OPTIONS */********************************

// 1. USER INPUT for new models: Change value of the macro variable MODEL_NAME to the name of the new model.
#define MODEL_NAME pkc
// include the simulation function with macros (#define statements) that make it model specific (based on MODEL_NAME)
#include "simulator.cpp"
// Global variables
static std::map <std::string, double> model_params;
// 2. USER INPUT for new models: Change the name of the wrapper function to sim_<MODEL_NAME> and the names of the internally called functions to init_<MODEL_NAME> and simulator_<MODEL_NAME>.
//' PKC Model R Wrapper Function (exported to R)
//'
//' This function calls the internal C++ simulator function to simulate the PKC model. 
//' @param
//' @return
//' @examples
//' sim_pkc()
//' @export
// [[Rcpp::export]]
NumericMatrix sim_pkc(DataFrame user_input_df,
                   NumericVector user_sim_params,
                   List user_model_params) {

  // Provide default model parameters
  model_params = init_pkc();
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
  return simulator_pkc(user_input_df,
                   user_sim_params,
                   user_model_params);
   
}



//********************************/* MODEL DEFINITION */********************************
// 3. USER INPUT for new models: define model parameters, number of species, 
// number of reactions, propensity equations and update_system function 

// Default model parameters
std::map <std::string, double> init() {
  // Model dimensions
  nspecies = 11;
  nreactions = 20;
  
  // Propensity equation parameters
  std::map <std::string, double> default_params;
  
  default_params["k1"] = 1;
  default_params["k2"] = 50;
  default_params["k3"] = 1.2e-7;
  default_params["k4"] = 0.1;
  default_params["k5"] = 1.2705;
  default_params["k6"] = 3.5026;
  default_params["k7"] = 1.2e-7;
  default_params["k8"] = 0.1;
  default_params["k9"] = 1;
  default_params["k10"] = 0.1;
  default_params["k11"] = 2;
  default_params["k12"] = 0.2;
  default_params["k13"] = 0.0006;
  default_params["k14"] = 0.5;
  default_params["k15"] = 7.998e-6;
  default_params["k16"] = 8.6348;
  default_params["k17"] = 6e-7;
  default_params["k18"] = 0.1;
  default_params["k19"] = 1.8e-5;
  default_params["k20"] = 2;
  default_params["AA"] = 11000;  // given as conc. remains fixed throughout the simulation
  default_params["DAG"] = 5000;  // given as conc. remains fixed throughout the simulation
  
  return default_params;
}

// Propensity calculation:
// Calculates the propensities of all PKC model reactions and stores them in the vector amu.
void calculate_amu() {
  
  // Look up model parameters in array 'double = model_params' initially
  double k1 = model_params["k1"];
  double k2 = model_params["k2"];
  double k3 = model_params["k3"];
  double k4 = model_params["k4"];
  double k5 = model_params["k5"];
  double k6 = model_params["k6"];
  double k7 = model_params["k7"];
  double k8 = model_params["k8"];
  double k9 = model_params["k9"];
  double k10 = model_params["k10"];
  double k11 = model_params["k11"];
  double k12 = model_params["k12"];
  double k13 = model_params["k13"];
  double k14 = model_params["k14"];
  double k15 = model_params["k15"];
  double k16 = model_params["k16"];
  double k17 = model_params["k17"];
  double k18 = model_params["k18"];
  double k19 = model_params["k19"];
  double k20 = model_params["k20"];
  double AA = model_params["AA"];
  double DAG = model_params["DAG"];
  
  amu[0] = k1 * x[0];
  amu[1] = amu[0] + k2 * x[5];
  amu[2] = amu[1] + k3 * AA * (double)x[0]; /* AA given as conc., hence, no scaling */
  amu[3] = amu[2] + k4 * x[6];
  amu[4] = amu[3] + k5 * x[1];
  amu[5] = amu[4] + k6 * x[7];
  amu[6] = amu[5] + k7 * AA * (double)x[1];  /* AA given as conc., hence, no scaling */
  amu[7] = amu[6] + k8 * x[8];
  amu[8] = amu[7] + k9 * x[2];
  amu[9] = amu[8] + k10 * x[9];
  amu[10] = amu[9] + k11 * x[3];
  amu[11] = amu[10] + k12 * x[4];
  amu[12] = amu[11] + calcium[ntimepoint] * k13 * (double)x[0]; /* Ca given as conc., hence, no scaling */
  amu[13] = amu[12] + k14 * x[1];
  amu[14] = amu[13] + k15 * DAG * (double)x[1]; /* DAG given as conc., hence, no scaling */
  amu[15] = amu[14] + k16 * x[2];
  amu[16] = amu[15] + k17 * DAG * (double)x[0]; /* DAG given as conc., hence, no scaling */
  amu[17] = amu[16] + k18 * x[10];
  amu[18] = amu[17] + k19 * AA * (double)x[10];  /* AA given as conc., hence, no scaling */
  amu[19] = amu[18] + k20 * x[3];
}

// System update:
// Changes the system state (updates the particle numbers) by instantiating a chosen reaction.
void update_system(unsigned int rIndex) {
  switch (rIndex) {
  case 0:   /* R1 */
    x[0]--;
    x[5]++;
    break;
  case 1:
    x[5]--;
    x[0]++;
    break;
  case 2:   /* R2 */
    x[0]--;
    x[6]++;
    break;
  case 3:
    x[6]--;
    x[0]++;
    break;
  case 4:  /* R3 */
    x[1]--;
    x[7]++;
    break;
  case 5:
    x[7]--;
    x[1]++;
    break;
  case 6:  /* R4 */
    x[1]--;
    x[8]++;
    break;
  case 7:
    x[8]--;
    x[1]++;
    break;
  case 8: /* R5 */
    x[2]--;
    x[9]++;
    break;
  case 9:
    x[9]--;
    x[2]++;
    break;
  case 10:/* R6 */
    x[3]--;
    x[4]++;
    break;
  case 11:
    x[4]--;
    x[3]++;
    break;
  case 12:/* R7 */
    x[0]--;
    x[1]++;
    break;
  case 13:
    x[1]--;
    x[0]++;
    break;
  case 14:/* R8 */
    x[1]--;
    x[2]++;
    break;
  case 15:
    x[2]--;
    x[1]++;
    break;
  case 16:/* R9 */
    x[0]--;
    x[10]++;
    break;
  case 17:
    x[10]--;
    x[0]++;
    break;
  case 18:/* R10 */
    x[10]--;
    x[3]++;
    break;
  case 19:
    x[3]--;
    x[10]++;
    break;
    printf("\nError in updateSystem(): rIndex (%u) out of range!\n", rIndex);
    exit(-1);
  }
}