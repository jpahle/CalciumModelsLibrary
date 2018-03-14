#include <string>
#include <Rcpp.h>
using namespace Rcpp;


//********************************/* R EXPORT OPTIONS */********************************

// 1. USER INPUT for new models: Change value of the macro variable MODEL_NAME to the name of the new model.
#define MODEL_NAME pkc
// include the simulation function with macros (#define statements) that make it model specific (based on MODEL_NAME)
#include "simulator.cpp"
// Global variables
static std::map <std::string, double> prop_params_map;
// 2. USER INPUT for new models: Change the name of the wrapper function to sim_<MODEL_NAME> and the names of the internally called functions to init_<MODEL_NAME> and simulator_<MODEL_NAME>.
//' PKC Model R Wrapper Function (exported to R)
//'
//' This function compares user-supplied parameters to defaults parameter values, overwrites the defaults if neccessary, and calls the internal C++ simulation function for the pkc model.
//' @param user_input_df A Dataframe: the input Calcium time series (with at least two columns: "time" in s and "Ca" in nMol/l).
//' @param user_sim_params A List: contains values for the simulation end ("endTime") and its timesteps ("timestep").
//' @param user_model_params A List: the model specific parameters. Can contain up to three different vectors named "vols" (model volumes), "init_conc" (initial conditions) and "params" (propensity equation parameters). 
//' @section Default Parameters of the Protein Kinase C Model:
//' Default Volumes: 
//' * vol = 1e-15
//' 
//' Default Initial Conditions:
//' * PKC_inact = 1000
//' * CaPKC = 0
//' * DAGCaPKC = 0
//' * AADAGPKC_inact = 0
//' * AADAGPKC_act = 0
//' * PKCbasal = 20
//' * AAPKC = 0
//' * CaPKCmemb = 0
//' * AACaPKC = 0
//' * DAGPKCmemb = 0
//' * DAGPKC = 0
//' 
//' Default Reaction Parameters:
//' * k1 = 1
//' * k2 = 50
//' * k3 = 1.2e-7
//' * k4 = 0.1
//' * k5 = 1.2705
//' * k6 = 3.5026
//' * k7 = 1.2e-7
//' * k8 = 0.1
//' * k9 = 1
//' * k10 = 0.1
//' * k11 = 2
//' * k12 = 0.2
//' * k13 = 0.0006
//' * k14 = 0.5
//' * k15 = 7.998e-6
//' * k16 = 8.6348
//' * k17 = 6e-7
//' * k18 = 0.1
//' * k19 = 1.8e-5
//' * k20 = 2
//' * AA = 11000
//' * DAG = 5000
//' @md
//' @return the result of calling the model specific version of the function "simulator" 
//' @examples
//' sim_pkc()
//' @export
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
DataFrame sim_pkc(DataFrame user_input_df,
                  List user_sim_params,
                  List user_model_params) {

  // READ INPUT
  // Provide default model parameters list
  List default_model_params = init_pkc();
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
  return simulator_pkc(user_input_df,
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
  nspecies = 11;
  nreactions = 20;
  
  // Default volume(s)
  NumericVector vols = NumericVector::create(
    _["vol"] = 1e-15
  );
  // Default initial conditions
  NumericVector init_conc = NumericVector::create(
    _["PKC_inact"] = 1000,
    _["CaPKC"] = 0,
    _["DAGCaPKC"] = 0,
    _["AADAGPKC_inact"] = 0,
    _["AADAGPKC_act"] = 0,
    _["PKCbasal"] = 20,
    _["AAPKC"] = 0,
    _["CaPKCmemb"] = 0,
    _["AACaPKC"] = 0,
    _["DAGPKCmemb"] = 0,
    _["DAGPKC"] = 0
  );
  // Default propensity equation parameters
  NumericVector params = NumericVector::create(
    _["k1"] = 1,
    _["k2"] = 50,
    _["k3"] = 1.2e-7,
    _["k4"] = 0.1,
    _["k5"] = 1.2705,
    _["k6"] = 3.5026,
    _["k7"] = 1.2e-7,
    _["k8"] = 0.1,
    _["k9"] = 1,
    _["k10"] = 0.1,
    _["k11"] = 2,
    _["k12"] = 0.2,
    _["k13"] = 0.0006,
    _["k14"] = 0.5,
    _["k15"] = 7.998e-6,
    _["k16"] = 8.6348,
    _["k17"] = 6e-7,
    _["k18"] = 0.1,
    _["k19"] = 1.8e-5,
    _["k20"] = 2
  );
  // The ::create function can only take 20 arguments at once (Rcpp problem)
  // All additional parameters need to be added via the slow push_back (x, name) method called on params vector
  params.push_back(11000, "AA"); // given as conc. remains fixed throughout the simulation
  params.push_back(5000, "DAG"); // given as conc. remains fixed throughout the simulation
    
  // Combine and return all vectors in a default_params list
  return List::create(
    _["vols"] = vols,
    _["init_conc"] = init_conc,
    _["params"] = params
  );
  
}

// Propensity calculation:
// Calculates the propensities of all PKC model reactions and stores them in the vector amu.
void calculate_amu() {
  
  // Look up model parameters in array 'double = model_params' initially
  double k1 = prop_params_map["k1"];
  double k2 = prop_params_map["k2"];
  double k3 = prop_params_map["k3"];
  double k4 = prop_params_map["k4"];
  double k5 = prop_params_map["k5"];
  double k6 = prop_params_map["k6"];
  double k7 = prop_params_map["k7"];
  double k8 = prop_params_map["k8"];
  double k9 = prop_params_map["k9"];
  double k10 = prop_params_map["k10"];
  double k11 = prop_params_map["k11"];
  double k12 = prop_params_map["k12"];
  double k13 = prop_params_map["k13"];
  double k14 = prop_params_map["k14"];
  double k15 = prop_params_map["k15"];
  double k16 = prop_params_map["k16"];
  double k17 = prop_params_map["k17"];
  double k18 = prop_params_map["k18"];
  double k19 = prop_params_map["k19"];
  double k20 = prop_params_map["k20"];
  double AA = prop_params_map["AA"];
  double DAG = prop_params_map["DAG"];
  
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


// Stoichiometric matrix
NumericMatrix get_stM() {
  
  // initialize stoich matrix (with zeroes)
  NumericMatrix stM(nspecies, nreactions);
  // create stoich matrix row vectors
  NumericVector stM_row1 = {-1, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, -1, 1, 0, 0};
  NumericVector stM_row2 = {0, 0, 0, 0, -1, 1, -1, 1, 0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0};
  NumericVector stM_row3 = {0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0};
  NumericVector stM_row4 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 1, -1};
  NumericVector stM_row5 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0};
  NumericVector stM_row6 = {1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  NumericVector stM_row7 = {0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  NumericVector stM_row8 = {0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  NumericVector stM_row9 = {0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  NumericVector stM_row10 = {0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  NumericVector stM_row11 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1};
  // fill rows of stoich matrix
  stM(0, _) = stM_row1;
  stM(1, _) = stM_row2;
  stM(2, _) = stM_row3;
  stM(3, _) = stM_row4;
  stM(4, _) = stM_row5;
  stM(5, _) = stM_row6;
  stM(6, _) = stM_row7;
  stM(7, _) = stM_row8;
  stM(8, _) = stM_row9;
  stM(9, _) = stM_row10;
  stM(10, _) = stM_row11;
  
  return stM;  
}