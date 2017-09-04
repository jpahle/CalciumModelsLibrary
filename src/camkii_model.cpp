#include <Rcpp.h>
using namespace Rcpp;



//********************************/* R EXPORT OPTIONS */********************************

// 1. USER INPUT for new models: Change value of the macro variable MODEL_NAME to the name of the new model.
#define MODEL_NAME camkii
// include the simulation function with macros (#define statements) that make it model specific (based on MODEL_NAME)
#include "simulator.cpp"
// Global variables
static std::map <std::string, double> model_params;
// 2. USER INPUT for new models: Change the name of the wrapper function to sim_<MODEL_NAME> and the names of the internally called functions to init_<MODEL_NAME> and simulator_<MODEL_NAME>.
//' CamKII Model R Wrapper Function (exported to R)
//'
//' This function calls the internal C++ simulator function to simulate the CamKII model. 
//' @param
//' @return
//' @examples
//' sim_camkii()
//' @export
// [[Rcpp::export]]
NumericMatrix sim_camkii(DataFrame user_input_df,
                   NumericVector user_sim_params,
                   List user_model_params) {

  // Provide default model parameters
  model_params = init_camkii();
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
  return simulator_camkii(user_input_df,
                   user_sim_params,
                   user_model_params);
   
}



//********************************/* MODEL DEFINITION */********************************
// 3. USER INPUT for new models: define model parameters, number of species, 
// number of reactions, propensity equations and update_system function 

// Default model parameters
std::map <std::string, double> init() {
  // Model dimensions
  nspecies = 5;
  nreactions = 10;
  
  // Propensity equation parameters
  std::map <std::string, double> default_params;
  
  default_params["a"] = -0.22;
  default_params["b"] = 1.826;
  default_params["c"] = 0.1;
  default_params["k_IB"] = 0.01;
  default_params["k_BI"] = 0.8;
  default_params["k_PT"] = 1;
  default_params["k_TP"] = 1e-12;
  default_params["k_TA"] = 0.0008;
  default_params["k_AT"] = 0.01;
  default_params["k_AA"] = 0.29;
  default_params["c_B"] = 0.75;
  default_params["c_P"] = 1;
  default_params["c_T"] = 0.8;
  default_params["c_A"] = 0.8;
  default_params["camT"] = 1000;
  default_params["Kd"] = 1000;
  default_params["Vm_phos"] = 0.005;
  default_params["Kd_phos"] = 0.3;
  default_params["totalC"] = 40;
  default_params["h"] = 4.0;
  
  return default_params;
}

// Propensity calculation:
// Calculates the propensities of all CamKII model reactions and stores them in the vector amu.
void calculate_amu() {
  
  // Look up model parameters in array 'model_params' initially
  double a = model_params["a"];
  double b = model_params["b"];
  double c = model_params["c"];
  double k_IB = model_params["k_IB"];
  double k_BI = model_params["k_BI"];
  double k_PT = model_params["k_PT"];
  double k_TP = model_params["k_TP"];
  double k_TA = model_params["k_TA"];
  double k_AT = model_params["k_AT"];
  double k_AA = model_params["k_AA"];
  double c_B = model_params["c_B"];
  double c_P = model_params["c_P"];
  double c_T = model_params["c_T"];
  double c_A = model_params["c_A"];
  double camT = model_params["camT"];
  double Kd = model_params["Kd"];
  double Vm_phos = model_params["Vm_phos"];
  double Kd_phos = model_params["Kd_phos"];
  double totalC = model_params["totalC"];
  double h = model_params["h"];
  
  amu[0] = x[0] * ((k_IB * camT * pow((double)calcium[ntimepoint],(double)h)) / (pow((double)calcium[ntimepoint],(double)h) + pow((double)Kd,(double)h)));
  amu[1] = amu[0] + k_BI * x[1];
  
  double activeSubunits = (x[1] + x[2] + x[3] + x[4]) / totalC;
  double prob =  a * activeSubunits + b*(pow((double)activeSubunits,(double)2)) + c*(pow((double)activeSubunits,(double)3));
  amu[2] = amu[1] +  totalC * k_AA * prob * ((c_B * x[1]) / pow((double)totalC,(double)2)) * (2*c_B*x[1] + c_P*x[2] + c_T*x[3]+ c_A*x[4]);
  
  amu[3] = amu[2] + k_PT * x[2];
  amu[4] = amu[3] + k_TP * x[3] * pow((double)calcium[ntimepoint],(double)h);
  amu[5] = amu[4] + k_TA * x[3];
  amu[6] = amu[5] + k_AT * x[4] * (camT - ((camT * pow((double)calcium[ntimepoint],(double)h)) / (pow((double)calcium[ntimepoint],(double)h) + pow((double)Kd,(double)h))));
  amu[7] = amu[6] + ((Vm_phos * x[2]) / (Kd_phos + (x[2] / totalC)));
  amu[8] = amu[7] + ((Vm_phos * x[3]) / (Kd_phos + (x[3] / totalC)));
  amu[9] = amu[8] + ((Vm_phos * x[4]) / (Kd_phos + (x[4] / totalC)));
}

// System update:
// Changes the system state (updates the particle numbers) by instantiating a chosen reaction.
void update_system(unsigned int rIndex) {
  switch (rIndex) {
  case 0:   // Ca-Cam binding
    x[0]--;
    x[1]++;
    break;
  case 1:   // binding reverse
    x[0]++;
    x[1]--;
    break;
  case 2:   // phosphorylation
    x[1]--;
    x[2]++;
    break;
  case 3:   // trapping
    x[2]--;
    x[3]++;
    break;
  case 4:   // trapping reverse
    x[2]++;
    x[3]--;
    break;
  case 5:   // autonomous
    x[3]--;
    x[4]++;
    break;
  case 6:   // autonomous reverse
    x[3]++;
    x[4]--;
    break;
  case 7:   // phosphatase on W_P
    x[1]++;
    x[2]--;
    break;
  case 8:   // phosphatase on W_T
    x[1]++;
    x[3]--;
    break;
  case 9:   // phosphatase on W_A
    x[0]++;
    x[4]--;
    break;
    printf("\nError in updateSystem(): rIndex (%u) out of range!\n", rIndex);
    exit(-1);
  }
}