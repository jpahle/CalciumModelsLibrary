#include "CaModLib_global_vars.hpp"
#include <Rcpp.h>
using namespace Rcpp;

/* Global variables */
static int nspecies = 5;
static int nreactions = 10;
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
static double c_I = 0;
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

// Define (initialize) variables declared in header
// Pointers *amu and *x are initialized with the adress of filler variables
// (These values are overwritten by the simulation)
static NumericVector calcium(1);
static unsigned int ntimepoint = 0;
static double amu_init_value = 0.1;
static double *amu = &amu_init_value;
static unsigned long long int x_init_value = 1000000;
static unsigned long long int *x = &x_init_value; 


//' Propensity Calculation
//'
//' Calculates the propensities of all CamKII model reactions and stores them in the vector amu.
//'
//' @param None
//' @return None
//' @examples
//' calculate_amu() 
//' @export
// [[Rcpp::export]]
void camkii_calculate_amu() {

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

//' System Update
//'
//' Changes the system state (updates the particle numbers) by instantiating a chosen reaction.
//'
//' @param rIndex An unsigned integer: the id of the chosen reaction.
//' @return void
//' @examples
//' update_system() 
//' @export
// [[Rcpp::export]]
void camkii_update_system(unsigned int rIndex) {
  
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