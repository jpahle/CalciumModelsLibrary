#include <Rcpp.h>
using namespace Rcpp;


/* Global variables */
static int nspecies = 2;
static int nreactions = 2;
static double k_on = 0.025;
static double k_off = 0.005;
static double Km = 1.0;
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
//' Calculates the propensities of all Calmodulin model reactions and stores them in the vector amu.
//'
//' @param None
//' @return None
//' @examples
//' calculate_amu() 
//' @export
// [[Rcpp::export]]
void calmodulin_calculate_amu() {

  amu[0] = ((k_on * pow((double)calcium[ntimepoint],(double)h)) / (pow((double)Km,(double)h) + pow((double)calcium[ntimepoint],(double)h))) * x[0];
  amu[1] = amu[0] + k_off * x[1];
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
void calmodulin_update_system(unsigned int rIndex) {
  
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