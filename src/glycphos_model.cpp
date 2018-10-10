#include <string>
#include <Rcpp.h>
using namespace Rcpp;


//********************************/* R EXPORT OPTIONS */********************************

// 1. USER INPUT for new models: Change value of the macro variable MODEL_NAME to the name of the new model.
#define MODEL_NAME glycphos
// include the simulation function with macros (#define statements) that make it model specific (based on MODEL_NAME)
#include "simulator.cpp"
// Global variables
static std::map <std::string, double> prop_params_map;
// 2. USER INPUT for new models: Change the name of the wrapper function to sim_<MODEL_NAME> and the names of the internally called functions to init_<MODEL_NAME> and simulator_<MODEL_NAME>.
//' Glycphos Model R Wrapper Function (exported to R)
//'
//' This function compares user-supplied parameters to defaults parameter values, overwrites the defaults if neccessary, and calls the internal C++ simulation function for the glycphos model.
//' @param user_input_df A Dataframe: the input Calcium time series (with at least two columns: "time" in s and "Ca" in nmol/l).
//' @param user_sim_params A List: contains values for the simulation end ("endTime") and its timesteps ("timestep").
//' @param user_model_params A List: the model specific parameters. Can contain up to three different vectors named "vols" (model volumes), "init_conc" (initial conditions) and "params" (propensity equation parameters). 
//' @section Default Parameters of the Glycogen Phosphorylase Model:
//' Default Volumes: 
//' * vol = 5e-14
//' 
//' Default Initial Conditions:
//' * Prot_inact = 5
//' * Prot_act = 0
//' 
//' Default Reaction Parameters:
//' * VpM1 = 1.5 (in min^-1)
//' * VpM2 = 0.6 (in min^-1)
//' * alpha = 9
//' * gamma = 9
//' * K11 = 0.1
//' * Kp2 = 0.2
//' * Ka1_conc = 1e7
//' * Ka2_conc = 1e7
//' * Ka5_conc = 500
//' * Ka6_conc = 500
//' * gluc_conc = 1e7 (in Gall 2000 model fixed at 10mM)
//' @md
//' @return the result of calling the model specific version of the function "simulator" 
//' @examples
//' sim_glycphos()
// [[Rcpp::plugins("cpp11")]]
//' @export
// [[Rcpp::export]]
DataFrame sim_glycphos(DataFrame user_input_df,
                       List user_sim_params,
                       List user_model_params) {

  // READ INPUT
  // Provide default model parameters list
  List default_model_params = init_glycphos();
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
  return simulator_glycphos(user_input_df,
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
    _["Prot_inact"] = 5.0,
    _["Prot_act"] = 0
  );
  // Default propensity equation parameters
  NumericVector params = NumericVector::create(
    
    _["VpM1"] = 1.5, // in min^-1
    _["VpM2"] = 0.6, // in min^-1
    _["alpha"] = 9,
    _["gamma"] = 9,
    _["K11"] = 0.1,
    _["Kp2"] = 0.2,
    // it is not necessary to convert glucose, Ka1, Ka2, Ka5 and Ka6 into particle numbers because the units cancel
    _["Ka1_conc"] = 1e7,
    _["Ka2_conc"] = 1e7,
    _["Ka5_conc"] = 500,
    _["Ka6_conc"] = 500,
    _["gluc_conc"] = 1e7 // in Gall 2000 model fixed at 10mM
  );
    
  // Combine and return all vectors in a default_params list
  return List::create(
    _["vols"] = vols,
    _["init_conc"] = init_conc,
    _["params"] = params
  );
  
}

// Propensity calculation:
// Calculates the propensities of all glycogen phosphorylase model reactions and stores them in the vector amu.
void calculate_amu() {
  
  // Look up model parameters in array 'double = model_params' initially
  double VpM1 = prop_params_map["VpM1"];
  double VpM2 = prop_params_map["VpM2"];
  double alpha = prop_params_map["alpha"];
  double gamma = prop_params_map["gamma"];
  double K11 = prop_params_map["K11"];
  double Kp2 = prop_params_map["Kp2"];
  double Ka1_conc = prop_params_map["Ka1_conc"];
  double Ka2_conc = prop_params_map["Ka2_conc"];
  double Ka5_conc = prop_params_map["Ka5_conc"];
  double Ka6_conc = prop_params_map["Ka6_conc"];
  double gluc_conc = prop_params_map["gluc_conc"];
  
  
  double total = x[0] + x[1];
  double activeFraction = x[1]/total;
  
  double Ca_conc_pow4 = calcium[ntimepoint] * calcium[ntimepoint] * calcium[ntimepoint] * calcium[ntimepoint];
  
  double Ka5_conc_pow4 = Ka5_conc * Ka5_conc * Ka5_conc * Ka5_conc;
  double Ka6_conc_pow4 = Ka6_conc * Ka6_conc * Ka6_conc * Ka6_conc;   

  // divide VpM1 and VpM2 by 60 to convert the units from min^-1 to s^-1
  amu[0] = (VpM1 / 60.0 * (1.0 + gamma * Ca_conc_pow4 / (Ka5_conc_pow4 + Ca_conc_pow4)) * ( 1.0 - activeFraction)) / ((K11 / (1.0 + Ca_conc_pow4 / Ka6_conc_pow4)) + 1.0 - activeFraction) * total;
  amu[1] = amu[0] + ((VpM2 / 60.0 * (1.0 + alpha * gluc_conc / (Ka1_conc + gluc_conc)) * activeFraction) / (Kp2 / (1 + gluc_conc / Ka2_conc) + activeFraction) * total);
}

// Stoichiometric matrix
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