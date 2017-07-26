#include <Rcpp.h>
using namespace Rcpp;



//********************************/* R EXPORT OPTIONS */********************************

// USER INPUT for new models: Change value of the macro variable MODEL_NAME to the name of the new model.
// Model description
#define MODEL_NAME calmodulin

// DON'T CHANGE THIS BLOCK!
//#######################################################################################
#define Map_helper(x,y) x##y
#define Map(x,y) Map_helper(x,y)
#define simulator Map(simulator_, MODEL_NAME)
#define init Map(init_, MODEL_NAME)
#define calculate_amu Map(calculate_amu_, MODEL_NAME)
#define update_system Map(update_system_, MODEL_NAME)
#include "simulator.cpp"
// Placeholder init function since it is defined 
// after the Wrapper Function tries to call it
void init();
//#######################################################################################

// USER INPUT for new models: Change the name of the wrapper function to sim_<MODEL_NAME> and the names of the internally called functions to init_<MODEL_NAME> and simulator_<MODEL_NAME>.
//' Calmodulin Model Wrapper Function (exported to R)
//'
//' This function calls the internal C++ simulator function to simulate the Calmodulin model. 
//' @param
//' @return
//' @examples
//' sim_calmodulin()
//' @export
// [[Rcpp::export]]
NumericMatrix sim_calmodulin(DataFrame param_input_df,
                   NumericVector param_sim_params,
                   List param_model_params) {
  
  init_calmodulin();
  
  // Use user-supplied model parameters over default values
  
  
  return simulator_calmodulin(param_input_df,
                   param_sim_params,
                   param_model_params);
   
}



//********************************/* MODEL DEFINITION */********************************
// USER INPUT for new models: define model parameters, number of species, 
// number of reactions, propensity equations and update_system function 

// Default Model Parameters
static double k_on = 0.025;
static double k_off = 0.005;
static double Km = 1.0;
static int h = 4;

// Model dimensions
void init() {
  nspecies = 2;
  nreactions = 2;
}

// Propensity calculation:
// Calculates the propensities of all Calmodulin model reactions and stores them in the vector amu.
void calculate_amu() {
  amu[0] = ((k_on * pow((double)calcium[ntimepoint],(double)h)) / (pow((double)Km,(double)h) + pow((double)calcium[ntimepoint],(double)h))) * x[0];
  amu[1] = amu[0] + k_off * x[1];
}

// System update:
// Changes the system state (updates the particle numbers) by instantiating a chosen reaction.
void update_system(unsigned int rIndex) {
  switch (rIndex) {
  case 0:   // Activation
    x[0]--;
    x[1]++;
    break;
  case 1:   // Deactivation
    x[0]++;
    x[1]--;
    break;
    printf("\nError in updateSystem(): rIndex (%u) out of range!\n", rIndex);
    exit(-1);
  }
}