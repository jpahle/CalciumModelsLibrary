#include <string>
#include <Rcpp.h>
using namespace Rcpp;


//********************************/* R EXPORT OPTIONS */********************************

// 1. USER INPUT for new models: Change value of the macro variable MODEL_NAME to the name of the new model.
#define MODEL_NAME ano
// include the simulation function with macros (#define statements) that make it model specific (based on MODEL_NAME)
#include "simulator.cpp"
// Global variables
static std::map <std::string, double> prop_params_map;
// 2. USER INPUT for new models: Change the name of the wrapper function to sim_<MODEL_NAME> and the names of the internally called functions to init_<MODEL_NAME> and simulator_<MODEL_NAME>.
//' Ano1 Model R Wrapper Function (exported to R)
//'
//' This function compares user-supplied parameters to defaults parameter values, overwrites the defaults if neccessary, and calls the internal C++ simulation function for the ano model.
//' @param user_input_df A Dataframe: the input Calcium time series (with at least two columns: "time" in s and "Ca" in nmol/l).
//' @param user_sim_params A NumericVector: contains values for the simulation end ("endTime") and its timesteps ("timestep").
//' @param user_model_params A List: the model specific parameters. Can contain up to three different vectors named "vols" (model volumes), "init_conc" (initial conditions) and "params" (propensity equation parameters). 
//' @section Default Parameters of the Ano1 Model:
//' Default Volumes: 
//' * vol = 1e-11
//' 
//' Default Initial Conditions:
//' * Cl_ext = 300
//' * C = 100
//' * C_c = 0
//' * C_1 = 0
//' * C_1c = 0
//' * C_2 = 0
//' * C_2c = 0
//' * O = 0
//' * O_c = 0
//' * O_1 = 0
//' * O_1c = 0
//' * O_2 = 0
//' * O_2c = 0
//' 
//' Default Reaction Parameters:
//' * Vm = -0.06
//' * T = 293.15
//' * a1 = 0.0077
//' * b1 = 917.1288
//' * k01 = 0.5979439
//' * k02 = 2.853
//' * acl1 = 1.8872
//' * bcl1 = 5955.783
//' * kccl1 = 1.143e-12
//' * kccl2 = 0.0009
//' * kocl1 = 1.1947e-06
//' * kocl2 = 3.4987
//' * za1 = 0
//' * zb1 = 0.0064
//' * zk01 = 0
//' * zk02 = 0.1684
//' * zacl1 = 0.1111
//' * zbcl1 = 0.3291
//' * zkccl1 = 0.1986
//' * zkccl2 = 0.0427
//' * zkocl1 = 0.6485 
//' * zkocl2 = 0.03 
//' * l = 41.6411 
//' * L = 0.6485 
//' * m = 0.0102 
//' * M = 0.0632 
//' * h = 0.3367 
//' * H = 14.2956 
//' @md
//' @return the result of calling the model specific version of the function "simulator" 
//' @examples
//' sim_ano()
//' @export
// [[Rcpp::export]]
NumericMatrix sim_ano(DataFrame user_input_df,
                   NumericVector user_sim_params,
                   List user_model_params) {

  // READ INPUT
  // Provide default model parameters list
  List default_model_params = init_ano();
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
  return simulator_ano(user_input_df,
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
  nspecies = 13;
  nreactions = 40;
  
  // Default volume(s)
  NumericVector vols = NumericVector::create(
    _["vol"] = 1e-11
  );
  // Default initial conditions
  NumericVector init_conc = NumericVector::create(
     _["Cl_ext"] = 300,
     _["C"] = 100,
     _["C_c"] = 0,
     _["C_1"] = 0,
     _["C_1c"] = 0,
     _["C_2"] = 0,
     _["C_2c"] = 0,
     _["O"] = 0,
     _["O_c"] = 0,
     _["O_1"] = 0,
     _["O_1c"] = 0,
     _["O_2"] = 0,
     _["O_2c"] = 0
  );
  // Default propensity equation parameters
  NumericVector params = NumericVector::create(
    _["Vm"] = -0.06,
    _["T"] = 293.15,
    _["a1"] = 0.0077,
    _["b1"] = 917.1288,
    _["k01"] = 0.5979439,
    _["k02"] = 2.853,
    _["acl1"] = 1.8872,
    _["bcl1"] = 5955.783,
    _["kccl1"] = 1.143e-12,
    _["kccl2"] = 0.0009,
    _["kocl1"] = 1.1947e-06,
    _["kocl2"] = 3.4987,
    _["za1"] = 0,
    _["zb1"] = 0.0064,
    _["zk01"] = 0,
    _["zk02"] = 0.1684,
    _["zacl1"] = 0.1111,
    _["zbcl1"] = 0.3291,
    _["zkccl1"] = 0.1986,
    _["zkccl2"] = 0.0427
  );
  // The ::create function can only take 20 arguments at once (Rcpp problem)
  // All additional parameters need to be added via the slow push_back (x, name) method called on params vector
  params.push_back(0.6485, "zkocl1");
  params.push_back(0.03, "zkocl2");
  params.push_back(41.6411, "l");
  params.push_back(0.6485, "L");
  params.push_back(0.0102, "m");
  params.push_back(0.0632, "M");
  params.push_back(0.3367, "h");
  params.push_back(14.2956, "H");
    
  // Combine and return all vectors in a default_params list
  return List::create(
    _["vols"] = vols,
    _["init_conc"] = init_conc,
    _["params"] = params
  );
}

// Propensity calculation:
// Calculates the propensities of all Ano1 model reactions and stores them in the vector amu.
void calculate_amu() {
  
  // Look up model parameters in array 'model_params' initially
  double Vm = prop_params_map["Vm"];
  double T = prop_params_map["T"];
  double a1 = prop_params_map["a1"];
  double b1 = prop_params_map["b1"];
  double k01 = prop_params_map["k01"];
  double k02 = prop_params_map["k02"];
  double acl1 = prop_params_map["acl1"];
  double bcl1 = prop_params_map["bcl1"];
  double kccl1 = prop_params_map["kccl1"];
  double kccl2 = prop_params_map["kccl2"];
  double kocl1 = prop_params_map["kocl1"];
  double kocl2 = prop_params_map["kocl2"];
  double za1 = prop_params_map["za1"];
  double zb1 = prop_params_map["zb1"];
  double zk01 = prop_params_map["zk01"];
  double zk02 = prop_params_map["zk02"];
  double zacl1 = prop_params_map["zacl1"];
  double zbcl1 = prop_params_map["zbcl1"];
  double zkccl1 = prop_params_map["zkccl1"];
  double zkccl2 = prop_params_map["zkccl2"];
  double zkocl1 = prop_params_map["zkocl1"];
  double zkocl2 = prop_params_map["zkocl2"];
  double l = prop_params_map["l"];
  double L = prop_params_map["L"];
  double m = prop_params_map["m"];
  double M = prop_params_map["M"];
  double h = prop_params_map["h"];
  double H = prop_params_map["H"];
  
  // Required constants
  // faradayConst = 96485.3329;
  // gasConst = 8.3144598; 
  double vterm;
  // vterm = 96485.3329 * model_params["Vm"] / (8.3144598 * model_params["T"]);
  vterm = 96485.3329 * Vm / (8.3144598 * T);
  
  //forward:      k1 * exp( z * vterm ) * x[0];
  //backward:     k1 * exp( -z * vterm ) * x[0];
  //forwardmod:   k1 * exp( z * vterm ) * x[mod] * x[0];
  //backwardmod:  k1 * exp( -z * vterm ) * x[0];
  //forward2mod:  k1 * exp( z * vterm ) * 2 * x[mod] * x[0];
  //backward2mod: k1 * exp( -z * vterm ) * x[0];
  //forward2rev:  k1 * exp( z * vterm ) * x[0];
  //backward2rev: k1 * exp( -z * vterm ) * 2 * x[0];
  
  // Propensity Equations (results are stored cumulative)
  amu[0] =             a1 * exp(za1 * vterm) * x[1]; //f: C - O
  amu[1] = amu[0] +    b1 * exp(-zb1 * vterm) * x[7]; //b: C - O
  amu[2] = amu[1] +    k01 * exp(zk01 * vterm) * 2 * calcium[ntimepoint] * x[1]; //f: C - Ca
  amu[3] = amu[2] +    l/L * k02 * exp(-zk02 * vterm) * x[3]; //b: C - Ca
  amu[4] = amu[3] +    kccl1 * exp(zkccl1 * vterm) * x[0] * x[1]; //f: C - Cl
  amu[5] = amu[4] +    kccl2 * exp(-zkccl2 * vterm) * x[2]; //b: C - Cl
  
  amu[6] = amu[5] +    acl1 * exp(zacl1 * vterm) * x[2]; //f: C_c - O
  amu[7] = amu[6] +    bcl1 * exp(-zbcl1 * vterm) * x[8]; //b: C_c - O
  amu[8] = amu[7] +    h/H * k01 * exp(zk01 * vterm) * 2 * calcium[ntimepoint] * x[2]; //f: C_c - Ca
  amu[9] = amu[8] +    l/L * k02 * exp(-zk02 * vterm) * x[4]; //b: C_c - Ca
  
  amu[10] = amu[9] +   l * a1 * exp(za1 * vterm) * x[3]; //f: C_1 - O
  amu[11] = amu[10] +  L * b1 * exp(-zb1 * vterm) * x[9]; //b: C_1 - O
  amu[12] = amu[11] +  k01 * exp(zk01 * vterm) * calcium[ntimepoint] * x[3]; //f: C_1 - Ca
  amu[13] = amu[12] +  l/L * 2 * k02 * exp(-zk02 * vterm) * x[5]; //b: C_1 - Ca
  amu[14] = amu[13] +  h * kccl1 * exp(zkccl1 * vterm) * x[0] * x[3]; //f: C_1 - Cl
  amu[15] = amu[14] +  H * kccl2 * exp(-zkccl2 * vterm) * x[4]; //b: C_1 - Cl
  
  amu[16] = amu[15] +  H*m*l/M * acl1 * exp(zacl1 * vterm) * x[4]; //f: C_1c - O
  amu[17] = amu[16] +  h*L * bcl1 * exp(-zbcl1 * vterm) * x[10]; //b: C_1c - O
  amu[18] = amu[17] +  h/H * k01 * exp(zk01 * vterm) * calcium[ntimepoint] * x[4]; //f: C_1c - Ca
  amu[19] = amu[18] +  l/L * 2 * k02 * exp(-zk02 * vterm) * x[6]; //b: C_1c - Ca
  
  amu[20] = amu[19] +  pow(l,2) * a1 * exp(za1 * vterm) * x[5]; //f: C_2 - O
  amu[21] = amu[20] +  pow(L,2) * b1 * exp(-zb1 * vterm) * x[11]; //b: C_2 - O
  amu[22] = amu[21] +  pow(h,2) * kccl1 * exp(zkccl1 * vterm) * x[0] * x[5]; //f: C_2 - Cl
  amu[23] = amu[22] +  pow(H,2) * kccl2 * exp(-zkccl2 * vterm) * x[6]; //b: C_2 - Cl
  
  amu[24] = amu[23] +  H*m*pow(l,2)/pow(m,2) * acl1 * exp(zacl1 * vterm) * x[6]; //f: C_2c - O
  amu[25] = amu[24] +  pow(h,2)*pow(L,2) * bcl1 * exp(-zbcl1 * vterm) * x[12]; //b: C_2c - O
  
  amu[26] = amu[25] +  k01 * exp(zk01 * vterm) * 2 * calcium[ntimepoint] * x[7]; //f: O - Ca
  amu[27] = amu[26] +  k02 * exp(-zk02 * vterm) * x[9]; //b: O - Ca
  amu[28] = amu[27] +  kocl1 * exp(zkocl1 * vterm) * x[0] * x[7]; //f: O - Cl
  amu[29] = amu[28] +  kocl2 * exp(-zkocl2 * vterm) * x[8]; //b: O - Cl
  
  amu[30] = amu[29] +  m/M * k01 * exp(zk01 * vterm) * 2 * calcium[ntimepoint] * x[8]; //f: O_c - Ca
  amu[31] = amu[30] +  k02 * exp(-zk02 * vterm) * x[10]; //b: O_c - Ca
  
  amu[32] = amu[31] +  k01 * exp(zk01 * vterm) * calcium[ntimepoint] * x[9]; //f: O_1 - Ca
  amu[33] = amu[32] +  2 * k02 * exp(-zk02 * vterm) * x[11]; //b: O_1 - Ca
  amu[34] = amu[33] +  m * kocl1 * exp(zkocl1 * vterm) * x[0] * x[9]; //f: O_1 - Cl
  amu[35] = amu[34] +  M * kocl2 * exp(-zkocl1 * vterm) * x[10]; //b: O_1 - Cl
  
  amu[36] = amu[35] +  m/M * k01 * exp(zk01 * vterm) * calcium[ntimepoint] * x[10]; //f: O_1c - Ca
  amu[37] = amu[36] +  2 * k02 * exp(-zk02 * vterm) * x[12]; //b: O_1c - Ca
  
  amu[38] = amu[37] +  pow(m,2) * kocl1 * exp(zkocl1 * vterm) * x[0] * x[11]; //f: O_2 - Cl
  amu[39] = amu[38] +  pow(M,2) * kocl2 * exp(-zkocl1 * vterm) * x[12]; //b: O_2 - Cl
  
}

// How it would look like it Rcpp would support the creation of vectors with more than 20 entries...
// (they don't: http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2011-November/003073.html
/* NumericMatrix get_stM() {
  
  // initialize stoich matrix (with zeroes)
  NumericMatrix stM(nspecies, nreactions);
  // initialize row vectors
  NumericVector stM_row1(nreactions);
  NumericVector stM_row2(nreactions);
  NumericVector stM_row3(nreactions);
  NumericVector stM_row4(nreactions);
  NumericVector stM_row5(nreactions);
  NumericVector stM_row6(nreactions);
  NumericVector stM_row7(nreactions);
  NumericVector stM_row8(nreactions);
  NumericVector stM_row9(nreactions);
  NumericVector stM_row10(nreactions);
  NumericVector stM_row11(nreactions);
  NumericVector stM_row12(nreactions);
  NumericVector stM_row13(nreactions);
  // fill stoich matrix rows
  NumericVector stM_row1  = NumericVector::create( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  NumericVector stM_row2  = NumericVector::create( -1, 1, -1, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  NumericVector stM_row3  = NumericVector::create( 0, 0, 0, 0, 1, -1, -1, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  NumericVector stM_row4  = NumericVector::create( 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, -1, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  NumericVector stM_row5  = NumericVector::create( 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  NumericVector stM_row6  = NumericVector::create( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  NumericVector stM_row7  = NumericVector::create( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  NumericVector stM_row8  = NumericVector::create( 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  NumericVector stM_row9  = NumericVector::create( 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0);
  NumericVector stM_row10 = NumericVector::create( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, -1, 1, -1, 1, 0, 0, 0, 0);
  NumericVector stM_row11 = NumericVector::create( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0);
  NumericVector stM_row12 = NumericVector::create( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, -1, 1);
  NumericVector stM_row13 = NumericVector::create( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1);
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
  stM(11, _) = stM_row12;
  stM(12, _) = stM_row13;
  
  return stM;  
} */