#include <string>
#include <Rcpp.h>
using namespace Rcpp;


//********************************/* R EXPORT OPTIONS */********************************

// 1. USER INPUT for new models: Change value of the macro variable MODEL_NAME to the name of the new model.
#define MODEL_NAME calmodulin
// include the simulation function with macros (#define statements) that make it model specific (based on MODEL_NAME)
#include "simulator.cpp"
// Global variables
static std::map <std::string, double> prop_params_map;
// 2. USER INPUT for new models: Change the name of the wrapper function to sim_<MODEL_NAME> and the names of the internally called functions to init_<MODEL_NAME> and simulator_<MODEL_NAME>.
//' Calmodulin Model R Wrapper Function (exported to R)
//'
//' This function compares user-supplied parameters to defaults parameter values, overwrites the defaults if neccessary, and calls the internal C++ simulation function for the Calmodulin model.
//' @param user_input_df A Dataframe: the input Calcium time series (with at least two columns: "time" in s and "Ca" in nmol/l).
//' @param user_sim_params A List: contains values for the simulation end ("endTime") and its timesteps ("timestep").
//' @param user_model_params A List: the model specific parameters. Can contain up to three different vectors named "vols" (model volumes), "init_conc" (initial conditions) and "params" (propensity equation parameters). 
//' @section Default Parameters of the Calmodulin Model:
//' Default Volumes: 
//' * vol = 5e-14
//' 
//' Default Initial Conditions:
//' * Prot_inact = 5
//' * Prot_act = 0
//' 
//' Default Reaction Parameters:
//' * k_on = 0.025
//' * k_off = 0.005
//' * Km = 1.0
//' * h = 4.0
//' @md
//' @return the result of calling the model specific version of the function "simulator" 
//' @examples
//' sim_calmodulin()
// [[Rcpp::plugins("cpp11")]]
//' @export
// [[Rcpp::export]]
DataFrame sim_calmodulin(DataFrame user_input_df,
                   List user_sim_params,
                   List user_model_params) {

  // READ INPUT
  // Provide default model parameters list
  List default_model_params = init_calmodulin();
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
  return simulator_calmodulin(user_input_df,
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
  nspecies = 2;
  nreactions = 2;
    
  // Default volume(s)
  NumericVector vols = NumericVector::create(
    _["vol"] = 5e-14
  );
  // Default initial conditions
  NumericVector init_conc = NumericVector::create(
    _["Prot_inact"] = 5,
    _["Prot_act"] = 0
  );
  // Default propensity equation parameters
  NumericVector params = NumericVector::create(
    _["k_on"] = 0.025,
    _["k_off"] = 0.005,
    _["Km"] = 1.0,
    _["h"] = 4.0
  );
    
  // Combine and return all vectors in a default_params list
  return List::create(
    _["vols"] = vols,
    _["init_conc"] = init_conc,
    _["params"] = params
  );
}

// Propensity calculation
// Calculates the propensities of all Calmodulin model reactions and stores them in the vector amu.
void calculate_amu() {
  
  // Look up model parameters in array 'prop_params_map' initially
  // (contains updated default parameters from vector default_params)
  double k_on = prop_params_map["k_on"];
  double k_off = prop_params_map["k_off"];
  double Km = prop_params_map["Km"];
  double h = prop_params_map["h"];
  
  amu[0] = ((k_on * pow((double)calcium[ntimepoint],(double)h)) / (pow((double)Km,(double)h) + pow((double)calcium[ntimepoint],(double)h))) * x[0];
  amu[1] = amu[0] + k_off * x[1];
    
}


// Stoichiometric matrix
//              R1   R2
// Prot_inact   -1    1
// Prot_act      1   -1
NumericMatrix get_stM() {
  
  // initialize stoich matrix (with zeroes)
  NumericMatrix stM(nspecies, nreactions);
  // create stoich matrix row vectors
  NumericVector stM_row1 = {-1,  1};
  NumericVector stM_row2 = { 1, -1};
  // fill rows of stoich matrix
  stM(0, _) = stM_row1;
  stM(1, _) = stM_row2;
  
  return stM;  
}