#include <Rcpp.h>
using namespace Rcpp;



//********************************/* R EXPORT OPTIONS */********************************

// 1. USER INPUT for new models: Change value of the macro variable MODEL_NAME to the name of the new model.
#define MODEL_NAME ano
// include the simulation function with macros (#define statements) that make it model specific (based on MODEL_NAME)
#include "simulator.cpp"
// Global variables
static std::map <std::string, double> model_params;
// 2. USER INPUT for new models: Change the name of the wrapper function to sim_<MODEL_NAME> and the names of the internally called functions to init_<MODEL_NAME> and simulator_<MODEL_NAME>.
//' Ano1 Model R Wrapper Function (exported to R)
//'
//' This function calls the internal C++ simulator function to simulate the Ano1 model. 
//' @param
//' @return
//' @examples
//' sim_ano()
//' @export
// [[Rcpp::export]]
NumericMatrix sim_ano(DataFrame user_input_df,
                   NumericVector user_sim_params,
                   List user_model_params) {

  // Provide default model parameters
  model_params = init_ano();
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
  return simulator_ano(user_input_df,
                   user_sim_params,
                   user_model_params);
   
}



//********************************/* MODEL DEFINITION */********************************
// 3. USER INPUT for new models: define model parameters, number of species, 
// number of reactions, propensity equations and update_system function 

// Default model parameters
std::map <std::string, double> init() {
  // Model dimensions
  nspecies = 14;
  nreactions = 40;
  
  // Propensity equation parameters
  std::map <std::string, double> default_params;
  
  default_params["Vm"] = -0.06;
  default_params["T"] = 293.15;
  default_params["a1"] = 0.0077;
  default_params["b1"] = 917.1288;
  default_params["k01"] = 0.5979439;
  default_params["k02"] = 2.853;
  default_params["acl1"] = 1.8872;
  default_params["bcl1"] = 5955.783;
  default_params["kccl1"] = 1.143e-12;
  default_params["kccl2"] = 0.0009;
  default_params["kocl1"] = 1.1947e-06;
  default_params["kocl2"] = 3.4987;
  default_params["za1"] = 0;
  default_params["zb1"] = 0.0064;
  default_params["zk01"] = 0;
  default_params["zk02"] = 0.1684;
  default_params["zacl1"] = 0.1111;
  default_params["zbcl1"] = 0.3291;
  default_params["zkccl1"] = 0.1986;
  default_params["zkccl2"] = 0.0427;
  default_params["zkocl1"] = 0.6485;
  default_params["zkocl2"] = 0.03;
  default_params["l"] = 41.6411;
  default_params["L"] = 0.1284;
  default_params["m"] = 0.0102;
  default_params["M"] = 0.0632;
  default_params["h"] = 0.3367;
  default_params["H"] = 14.2956;
  
  return default_params;
}

// Propensity calculation:
// Calculates the propensities of all Ano1 model reactions and stores them in the vector amu.
void calculate_amu() {
  
  // Required constants
  // faradayConst = 96485.3329;
  // gasConst = 8.3144598; 
  double vterm;
  vterm = 96485.3329 * model_params["Vm"] / (8.3144598 * model_params["T"]);
  
  //forward:      k1 * exp( z * vterm ) * x[0];
  //backward:     k1 * exp( -z * vterm ) * x[0];
  //forwardmod:   k1 * exp( z * vterm ) * x[mod] * x[0];
  //backwardmod:  k1 * exp( -z * vterm ) * x[0];
  //forward2mod:  k1 * exp( z * vterm ) * 2 * x[mod] * x[0];
  //backward2mod: k1 * exp( -z * vterm ) * x[0];
  //forward2rev:  k1 * exp( z * vterm ) * x[0];
  //backward2rev: k1 * exp( -z * vterm ) * 2 * x[0];
  
  amu[0] =             model_params["a1"] * exp(model_params["za1"] * vterm) * x[2]; //f: C - O
  amu[1] = amu[0] +    model_params["b1"] * exp(-model_params["zb1"] * vterm) * x[8]; //b: C - O
  amu[2] = amu[1] +    model_params["k01"] * exp(model_params["zk01"] * vterm) * 2 * x[0] * x[2]; //f: C - Ca
  amu[3] = amu[2] +    model_params["l"]/model_params["L"] * model_params["k02"] * exp(-model_params["zk02"] * vterm) * x[4]; //b: C - Ca
  amu[4] = amu[3] +    model_params["kccl1"] * exp(model_params["zkccl1"] * vterm) * x[1] * x[2]; //f: C - Cl
  amu[5] = amu[4] +    model_params["kccl2"] * exp(-model_params["zkccl2"] * vterm) * x[3]; //b: C - Cl
  
  amu[6] = amu[5] +    model_params["acl1"] * exp(model_params["zacl1"] * vterm) * x[3]; //f: C_c - O
  amu[7] = amu[6] +    model_params["bcl1"] * exp(-model_params["zbcl1"] * vterm) * x[9]; //b: C_c - O
  amu[8] = amu[7] +    model_params["h"]/model_params["H"] * model_params["k01"] * exp(model_params["zk01"] * vterm) * 2 * x[0] * x[3]; //f: C_c - Ca
  amu[9] = amu[8] +    model_params["l"]/model_params["L"] * model_params["k02"] * exp(-model_params["zk02"] * vterm) * x[5]; //b: C_c - Ca
  
  amu[10] = amu[9] +   model_params["l"] * model_params["a1"] * exp(model_params["za1"] * vterm) * x[4]; //f: C_1 - O
  amu[11] = amu[10] +  model_params["L"] * model_params["b1"] * exp(-model_params["zb1"] * vterm) * x[10]; //b: C_1 - O
  amu[12] = amu[11] +  model_params["k01"] * exp(model_params["zk01"] * vterm) * x[0] * x[4]; //f: C_1 - Ca
  amu[13] = amu[12] +  model_params["l"]/model_params["L"] * 2 * model_params["k02"] * exp(-model_params["zk02"] * vterm) * x[6]; //b: C_1 - Ca
  amu[14] = amu[13] +  model_params["h"] * model_params["kccl1"] * exp(model_params["zkccl1"] * vterm) * x[1] * x[4]; //f: C_1 - Cl
  amu[15] = amu[14] +  model_params["H"] * model_params["kccl2"] * exp(-model_params["zkccl2"] * vterm) * x[5]; //b: C_1 - Cl
  
  amu[16] = amu[15] +  model_params["H"]*model_params["m"]*model_params["l"]/model_params["M"] * model_params["acl1"] * exp(model_params["zacl1"] * vterm) * x[5]; //f: C_1c - O
  amu[17] = amu[16] +  model_params["h"]*model_params["L"] * model_params["bcl1"] * exp(-model_params["zbcl1"] * vterm) * x[11]; //b: C_1c - O
  amu[18] = amu[17] +  model_params["h"]/model_params["H"] * model_params["k01"] * exp(model_params["zk01"] * vterm) * x[0] * x[5]; //f: C_1c - Ca
  amu[19] = amu[18] +  model_params["l"]/model_params["L"] * 2 * model_params["k02"] * exp(-model_params["zk02"] * vterm) * x[7]; //b: C_1c - Ca
  
  amu[20] = amu[19] +  pow(model_params["l"],2) * model_params["a1"] * exp(model_params["za1"] * vterm) * x[6]; //f: C_2 - O
  amu[21] = amu[20] +  pow(model_params["L"],2) * model_params["b1"] * exp(-model_params["zb1"] * vterm) * x[12]; //b: C_2 - O
  amu[22] = amu[21] +  pow(model_params["h"],2) * model_params["kccl1"] * exp(model_params["zkccl1"] * vterm) * x[1] * x[6]; //f: C_2 - Cl
  amu[23] = amu[22] +  pow(model_params["H"],2) * model_params["kccl2"] * exp(-model_params["zkccl2"] * vterm) * x[7]; //b: C_2 - Cl
  
  amu[24] = amu[23] +  model_params["H"]*model_params["m"]*pow(model_params["l"],2)/pow(model_params["m"],2) * model_params["acl1"] * exp(model_params["zacl1"] * vterm) * x[7]; //f: C_2c - O
  amu[25] = amu[24] +  pow(model_params["h"],2)*pow(model_params["L"],2) * model_params["bcl1"] * exp(-model_params["zbcl1"] * vterm) * x[13]; //b: C_2c - O
  
  amu[26] = amu[25] +  model_params["k01"] * exp(model_params["zk01"] * vterm) * 2 * x[0] * x[8]; //f: O - Ca
  amu[27] = amu[26] +  model_params["k02"] * exp(-model_params["zk02"] * vterm) * x[10]; //b: O - Ca
  amu[28] = amu[27] +  model_params["kocl1"] * exp(model_params["zkocl1"] * vterm) * x[1] * x[8]; //f: O - Cl
  amu[29] = amu[28] +  model_params["kocl2"] * exp(-model_params["zkocl2"] * vterm) * x[9]; //b: O - Cl
  
  amu[30] = amu[29] +  model_params["m"]/model_params["M"] * model_params["k01"] * exp(model_params["zk01"] * vterm) * 2 * x[0] * x[9]; //f: O_c - Ca
  amu[31] = amu[30] +  model_params["k02"] * exp(-model_params["zk02"] * vterm) * x[11]; //b: O_c - Ca
  
  amu[32] = amu[31] +  model_params["k01"] * exp(model_params["zk01"] * vterm) * x[0] * x[10]; //f: O_1 - Ca
  amu[33] = amu[32] +  2 * model_params["k02"] * exp(-model_params["zk02"] * vterm) * x[12]; //b: O_1 - Ca
  amu[34] = amu[33] +  model_params["m"] * model_params["kocl1"] * exp(model_params["zkocl1"] * vterm) * x[1] * x[10]; //f: O_1 - Cl
  amu[35] = amu[34] +  model_params["M"] * model_params["kocl2"] * exp(-model_params["zkocl1"] * vterm) * x[11]; //b: O_1 - Cl
  
  amu[36] = amu[35] +  model_params["m"]/model_params["M"] * model_params["k01"] * exp(model_params["zk01"] * vterm) * x[0] * x[11]; //f: O_1c - Ca
  amu[37] = amu[36] +  2 * model_params["k02"] * exp(-model_params["zk02"] * vterm) * x[13]; //b: O_1c - Ca
  
  amu[38] = amu[37] +  pow(model_params["m"],2) * model_params["kocl1"] * exp(model_params["zkocl1"] * vterm) * x[1] * x[12]; //f: O_2 - Cl
  amu[39] = amu[38] +  pow(model_params["M"],2) * model_params["kocl2"] * exp(-model_params["zkocl1"] * vterm) * x[13]; //b: O_2 - Cl
}

// System update:
// Changes the system state (updates the particle numbers) by instantiating a chosen reaction.
void update_system(unsigned int rIndex) {
  switch (rIndex) {
    case 0: // C - O
        x[2]--;
        x[8]++;
        break;
    case 1:
        x[8]--;
        x[2]++;
        break;
    case 2: // C - Ca
        x[2]--;
        x[4]++;
        break;
    case 3:
        x[4]--;
        x[2]++;
        break;
    case 4: // C - Cl
        x[2]--;
        x[3]++;
        break;
    case 5:
        x[3]--;
        x[2]++;
        break;
    case 6: // C_c - O
        x[3]--;
        x[9]++;
        break;
    case 7:
        x[9]--;
        x[3]++;
        break;
    case 8: // C_c - Ca
        x[3]--;
        x[5]++;
        break;
    case 9:
        x[5]--;
        x[3]++;
        break;
    case 10: // C_1 - O
        x[4]--;
        x[10]++;
        break;
    case 11:
        x[10]--;
        x[4]++;
        break;
    case 12: // C_1 - Ca
        x[4]--;
        x[6]++;
        break;
    case 13:
        x[6]--;
        x[4]++;
        break;
    case 14: // C_1 - Cl
        x[4]--;
        x[5]++;
        break;
    case 15:
        x[5]--;
        x[4]++;
        break;
    case 16: // C_1c - O
        x[5]--;
        x[11]++;
        break;
    case 17:
        x[11]--;
        x[5]++;
        break;
    case 18: // C_1c - Ca
        x[5]--;
        x[7]++;
        break;
    case 19:
        x[7]--;
        x[5]++;
        break;
    case 20: // C_2 - O
        x[6]--;
        x[12]++;
        break;
    case 21:
        x[12]--;
        x[6]++;
        break;
    case 22: // C_2 - Cl
        x[6]--;
        x[7]++;
        break;
    case 23:
        x[7]--;
        x[6]++;
        break;
    case 24: // C_2c - O
        x[7]--;
        x[13]++;
        break;
    case 25:
        x[13]--;
        x[7]++;
        break;
    case 26: // O - Ca
        x[8]--;
        x[10]++;
        break;
    case 27:
        x[10]--;
        x[8]++;
        break;
    case 28: // O - Cl
        x[8]--;
        x[9]++;
        break;
    case 29:
        x[9]--;
        x[8]++;
        break;
    case 30: // O_c - Ca
        x[9]--;
        x[11]++;
        break;
    case 31:
        x[11]--;
        x[9]++;
        break;
    case 32: // O_1 - Ca
        x[10]--;
        x[12]++;
        break;
    case 33:
        x[12]--;
        x[10]++;
        break;
    case 34: // O_1 - Cl
        x[10]--;
        x[11]++;
        break;
    case 35:
        x[11]--;
        x[10]++;
        break;
    case 36: // O_1c - Ca
        x[11]--;
        x[13]++;
        break;
    case 37:
        x[13]--;
        x[11]++;
        break;
    case 38: // O_2 - Cl
        x[12]--;
        x[13]++;
        break;
    case 39:
        x[13]--;
        x[12]++;
        break;
        printf("\nError in updateSystem(): rIndex (%u) out of range!\n", rIndex);
        exit(-1);
    }
}