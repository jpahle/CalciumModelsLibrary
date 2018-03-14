#include <string>
#include <Rcpp.h>
using namespace Rcpp;


//********************************/* R EXPORT OPTIONS */********************************

// 1. USER INPUT for new models: Change value of the macro variable MODEL_NAME to the name of the new model.
#define MODEL_NAME camkii
// include the simulation function with macros (#define statements) that make it model specific (based on MODEL_NAME)
#include "simulator.cpp"
// Global variables
static std::map <std::string, double> prop_params_map;
// 2. USER INPUT for new models: Change the name of the wrapper function to sim_<MODEL_NAME> and the names of the internally called functions to init_<MODEL_NAME> and simulator_<MODEL_NAME>.
//' CamKII Model R Wrapper Function (exported to R)
//'
//' This function compares user-supplied parameters to defaults parameter values, overwrites the defaults if neccessary, and calls the internal C++ simulation function for the camkii model.
//' @param user_input_df A Dataframe: the input Calcium time series (with at least two columns: "time" in s and "Ca" in nmol/l).
//' @param user_sim_params A List: contains values for the simulation end ("endTime") and its timesteps ("timestep").
//' @param user_model_params A List: the model specific parameters. Can contain up to three different vectors named "vols" (model volumes), "init_conc" (initial conditions) and "params" (propensity equation parameters). 
//' @section Default Parameters of the CamKII Model:
//' Default Volumes: 
//' * vol = 5e-15
//' 
//' Default Initial Conditions:
//' * W_I = 40
//' * W_B = 0
//' * W_P = 0
//' * W_T = 0
//' * W_A = 0
//' 
//' Default Reaction Parameters:
//' * a = -0.22
//' * b = 1.826
//' * c = 0.1
//' * k_IB = 0.01
//' * k_BI = 0.8
//' * k_PT = 1
//' * k_TP = 1e-12
//' * k_TA = 0.0008
//' * k_AT = 0.01
//' * k_AA = 0.29
//' * c_B = 0.75
//' * c_P = 1
//' * c_T = 0.8
//' * c_A = 0.8
//' * camT = 1000
//' * Kd = 1000
//' * Vm_phos = 0.005
//' * Kd_phos = 0.3
//' * totalC = 40
//' * h = 4.0
//' @md
//' @return the result of calling the model specific version of the function "simulator" 
//' @examples
//' sim_camkii()
//' @export
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
DataFrame sim_camkii(DataFrame user_input_df,
                     List user_sim_params,
                     List user_model_params) {

  // READ INPUT
  // Provide default model parameters list
  List default_model_params = init_camkii();
  // Extract default vectors from list
  NumericVector default_vols = default_model_params["vols"];
  NumericVector default_init_conc = default_model_params["init_conc"];
  NumericVector default_params = default_model_params["params"];
  // Extract vectors from user supplied list 
  // If they exist: definition with the default vector gets overwritten
  NumericVector user_vols = default_vols;
  if (user_model_params.containsElementNamed("vols")) {
      user_vols = user_model_params["vols"];
  } else {
    Rcout << "Default volume(s) have been used." << std::endl;
  } 
  NumericVector user_init_conc = default_init_conc;
  if (user_model_params.containsElementNamed("init_conc")) {
    user_init_conc = user_model_params["init_conc"];
  } else {
    Rcout << "Default initial condition(s) have been used." << std::endl;
  }
  NumericVector user_params = default_params;
  if (user_model_params.containsElementNamed("params")) {
    user_params = user_model_params["params"];
  } else {
    Rcout << "Default reaction parameter(s) have been used." << std::endl;
  }
  // UPDATE DEFAULTS 
  // Replace entries in default_model_params with user-supplied values if necessary
  // 1.) Volumes update:
  CharacterVector user_vols_names = user_vols.names();
  for (int i = 0; i < user_vols_names.length(); i++) {
    std::string current_vol_name = as<std::string>(user_vols_names[i]);
    if (default_vols.containsElementNamed((current_vol_name).c_str())) {
      // update default values
      default_vols[current_vol_name] = user_vols[current_vol_name];    
    } else {
      Rcout << "No such index! Default values have been used. Check input parameter vectors." << std::endl;
    }
  }
  // 2.) Initial conditions update:
  CharacterVector user_init_conc_names = user_init_conc.names();
  for (int i = 0; i < user_init_conc_names.length(); i++) {
    std::string current_init_conc_name = as<std::string>(user_init_conc_names[i]);
    if (default_init_conc.containsElementNamed((current_init_conc_name).c_str())) {
      // update default values
      default_init_conc[current_init_conc_name] = user_init_conc[current_init_conc_name];    
    } else {
      Rcout << "No such index! Default values have been used. Check input parameter vectors." << std::endl;
    }
  }
  // 3.) Propensity equation parameters update:
  CharacterVector user_params_names = user_params.names();
  for (int i = 0; i < user_params_names.length(); i++) {
    std::string current_param_name = as<std::string>(user_params_names[i]);
    if (default_params.containsElementNamed((current_param_name).c_str())) {
      // update default values
      default_params[current_param_name] = user_params[current_param_name];    
    } else {
      Rcout << "No such index! Default values have been used. Check input parameter vectors." << std::endl;
    }
  }
  // Put propensity reaction parameters in a map (for function calculate_amu)
  // (take parameters from vector "default_params" which contains the updated values)
  CharacterVector default_params_names = default_params.names();
  for (int n = 0; n < default_params.length(); n++) {
    std::string current_param_name = as<std::string>(default_params_names[n]);
    prop_params_map[current_param_name] = default_params[current_param_name];  
  }
  // RUN SIMULATION
  // Return result of the included, model-specific copy of the function "simulator" 
  return simulator_camkii(user_input_df,
                          user_sim_params,
                          default_vols,
                          default_init_conc);
   
}



//********************************/* MODEL DEFINITION */********************************
// 3. USER INPUT for new models: define model parameters, number of species, 
// number of reactions, propensity equations and update_system function 

// Default model parameters
List init() {
  // Model dimensions
  nspecies = 5;
  nreactions = 10;
  
  // Default volume(s)
  NumericVector vols = NumericVector::create(
    _["vol"] = 5e-15
  );
  // Default initial conditions
  NumericVector init_conc = NumericVector::create(
    _["W_I"] = 40,
    _["W_B"] = 0,
    _["W_P"] = 0,
    _["W_T"] = 0,
    _["W_A"]= 0
  );
  // Default propensity equation parameters
  NumericVector params = NumericVector::create(
    _["a"] = -0.22,
    _["b"] = 1.826,
    _["c"] = 0.1,
    _["k_IB"] = 0.01,
    _["k_BI"] = 0.8,
    _["k_PT"] = 1,
    _["k_TP"] = 1e-12,
    _["k_TA"] = 0.0008,
    _["k_AT"] = 0.01,
    _["k_AA"] = 0.29,
    _["c_B"] = 0.75,
    _["c_P"] = 1,
    _["c_T"] = 0.8,
    _["c_A"] = 0.8,
    _["camT"] = 1000,
    _["Kd"] = 1000,
    _["Vm_phos"] = 0.005,
    _["Kd_phos"] = 0.3,
    _["totalC"] = 40,
    _["h"] = 4.0
  );
    
  // Combine and return all vectors in a default_params list
  return List::create(
    _["vols"] = vols,
    _["init_conc"] = init_conc,
    _["params"] = params
  );
  
}

// Propensity calculation:
// Calculates the propensities of all CamKII model reactions and stores them in the vector amu.
void calculate_amu() {
  
  // Look up model parameters in array 'model_params' initially
  double a = prop_params_map["a"];
  double b = prop_params_map["b"];
  double c = prop_params_map["c"];
  double k_IB = prop_params_map["k_IB"];
  double k_BI = prop_params_map["k_BI"];
  double k_PT = prop_params_map["k_PT"];
  double k_TP = prop_params_map["k_TP"];
  double k_TA = prop_params_map["k_TA"];
  double k_AT = prop_params_map["k_AT"];
  double k_AA = prop_params_map["k_AA"];
  double c_B = prop_params_map["c_B"];
  double c_P = prop_params_map["c_P"];
  double c_T = prop_params_map["c_T"];
  double c_A = prop_params_map["c_A"];
  double camT = prop_params_map["camT"];
  double Kd = prop_params_map["Kd"];
  double Vm_phos = prop_params_map["Vm_phos"];
  double Kd_phos = prop_params_map["Kd_phos"];
  double totalC = prop_params_map["totalC"];
  double h = prop_params_map["h"];
  
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

// Stoichiometric matrix
NumericMatrix get_stM() {
  
  // initialize stoich matrix (with zeroes)
  NumericMatrix stM(nspecies, nreactions);
  // create stoich matrix row vectors
  NumericVector stM_row1 = {-1, 1, 0, 0, 0, 0, 0, 0, 0, 1};
  NumericVector stM_row2 = {1, -1, -1, 0, 0, 0, 0, 1, 1, 0};
  NumericVector stM_row3 = {0, 0, 1, -1, 1, 0, 0, -1, 0, 0};
  NumericVector stM_row4 = {0, 0, 0, 1, -1, -1, 1, 0, -1, 0};
  NumericVector stM_row5 = {0, 0, 0, 0, 0, 1, -1, 0, 0, -1};
  // fill rows of stoich matrix
  stM(0, _) = stM_row1;
  stM(1, _) = stM_row2;
  stM(2, _) = stM_row3;
  stM(3, _) = stM_row4;
  stM(4, _) = stM_row5;
  
  return stM;  
}