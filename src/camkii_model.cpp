#include <Rcpp.h>
using namespace Rcpp;



//********************************/* R EXPORT OPTIONS */********************************

// USER INPUT for new models: Change value of the macro variable MODEL_NAME to the name of the new model.
// Model description
#define MODEL_NAME camkii

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
//' CamKII Model Wrapper Function (exported to R)
//'
//' This function calls the internal C++ simulator function to simulate the CamKII model. 
//' @param
//' @return
//' @examples
//' sim_camkii()
//' @export
// [[Rcpp::export]]
NumericMatrix sim_camkii(DataFrame param_input_df,
                   NumericVector param_sim_params,
                   List param_model_params) {
  
  init_camkii();
  return simulator_camkii(param_input_df,
                   param_sim_params,
                   param_model_params);
   
}



//********************************/* MODEL DEFINITION */********************************
// USER INPUT for new models: define model parameters, number of species, 
// number of reactions, propensity equations and update_system function 

// Default Model Parameters
static double a = -0.22;
static double b = 1.826;
static double c = 0.1;
static double k_IB = 0.01;
static double k_BI = 0.8;
static double k_PT = 1;
static double k_TP = 1e-12;
static double k_TA = 0.0008;
static double k_AT = 0.01;
static double k_AA = 0.29;
static double c_B = 0.75;
static double c_P = 1;
static double c_T = 0.8;
static double c_A = 0.8;
static double camT = 1000;
static double Kd = 1000;
static double Vm_phos = 0.005;
static double Kd_phos = 0.3;
static double totalC = 40;
static int h = 4;

// Model dimensions
void init() {
  nspecies = 5;
  nreactions = 10;
}

// Propensity calculation:
// Calculates the propensities of all CamKII model reactions and stores them in the vector amu.
void calculate_amu() {
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