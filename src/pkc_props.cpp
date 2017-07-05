#include <Rcpp.h>
using namespace Rcpp;

//' Calculate propensity for a PKC Model reaction.
//'
//' Return the propensity of a PKC Model reaction for a given vector of particle numbers and a reaction Id. 
//' 
//' @param part_num A numeric vector: the particle numbers of the model species.
//' @param calcium A numeric vector: the calcium particle number.
//' @param rId An integer value: the id of the specified reaction for which the propensity should be calculated.
//' @return A double value (the propensity of the specified reaction).
//' @examples
//' pkc_props()
//' @export
// [[Rcpp::export]]
double pkc_props(NumericVector part_num, double calcium, int rId) {
  
  // model parameters
  double k1 = 1;
  double k2 = 50;
  double k3 = 1.2e-7;
  double k4 = 0.1;
  double k5 = 1.2705;
  double k6 = 3.5026;
  double k7 = 1.2e-7;
  double k8 = 0.1;
  double k9 = 1;
  double k10 = 0.1;
  double k11 = 2;
  double k12 = 0.2;
  double k13 = 0.0006;
  double k14 = 0.5;
  double k15 = 7.998e-6;
  double k16 = 8.6348;
  double k17 = 6e-7;
  double k18 = 0.1;
  double k19 = 1.8e-5;
  double k20 = 2;
  double AA = 11000;
  double DAG = 5000;
  
  // calculate propensity of selected reaction
  NumericVector x = part_num;
  double prop;
  
  switch (rId) {
  case 0:
  prop = k1 * x[0];
  case 1:
    prop = k2 * x[5];
    break;
  case 2:
    prop = k3 * AA * (double)x[0]; /* AA given as conc., hence, no scaling */
    break;
  case 3:
    prop = k4 * x[6];
    break;
  case 4:
    prop = k5 * x[1];
    break;
  case 5:
    prop = k6 * x[7];
    break;
  case 6:
    prop = k7 * AA * (double)x[1];  /* AA given as conc., hence, no scaling */
    break;
  case 7:
    prop = k8 * x[8];
    break;
  case 8:
    prop = k9 * x[2];
    break;
  case 9:
    prop = k10 * x[9];
    break;
  case 10:
    prop = k11 * x[3];
    break;
  case 11:
    prop = k12 * x[4];
    break;
  case 12:
    prop = calcium * k13 * (double)x[0]; /* Ca given as conc., hence, no scaling */
    break;
  case 13:
    prop = k14 * x[1];
    break;
  case 14:
    prop = k15 * DAG * (double)x[1]; /* DAG given as conc., hence, no scaling */
    break;
  case 15:
    prop = k16 * x[2];
    break;
  case 16:
    prop = k17 * DAG * (double)x[0]; /* DAG given as conc., hence, no scaling */
    break;
  case 17:
    prop = k18 * x[10];
    break;
  case 18:
    prop = k19 * AA * (double)x[10];  /* AA given as conc., hence, no scaling */
    break;
  case 19:
    prop = k20 * x[3];
    break;
  default:
    printf("\nError in propensity calculation: reaction Index (%u) out of range!\n", rId);
    prop = 0;
  }
    
  return prop;
}

//' Define stoichiometric matrix of the PKC model
//' 
//' Create and return the stoichiometric matrix of the PKC model as a numeric matrix.
//' 
//' @return A numeric matrix: the stoichiometric matrix
//' @examples
//' pkc_stM()
//' @export
// [[Rcpp::export]] 
NumericMatrix pkc_stM() {

  // initialize Stoich. Matrix filled with zeroes
  NumericMatrix stM(11,20);
  // R1 forward
  stM(0,0) = -1;
  stM(5,0) = 1;
  // R1 reverse
  stM(0,1) = 1;
  stM(5,1) = -1;
  // R2 forward
  stM(0,2) = -1;
  stM(6,2) = 1;
  // R2 reverse
  stM(0,3) = 1;
  stM(6,3) = -1;
  // R3 forward
  stM(1,4) = -1;
  stM(7,4) = 1;
  // R3 reverse
  stM(1,5) = 1;
  stM(7,5) = -1;
  // R4 forward
  stM(1,6) = -1;
  stM(8,6) = 1;
  // R4 reverse
  stM(1,7) = 1;
  stM(8,7) = -1;
  // R5 forward
  stM(2,8) = -1;
  stM(9,8) = 1;
  // R5 reverse
  stM(2,9) = 1;
  stM(9,9) = -1;
  // R6 forward
  stM(3,10) = -1;
  stM(4,10) = 1;
  // R6 reverse
  stM(3,11) = 1;
  stM(4,11) = -1;
  // R7 forward
  stM(0,12) = -1;
  stM(1,12) = 1;
  // R7 reverse
  stM(0,13) = 1;
  stM(1,13) = -1;
  // R8 forward
  stM(1,14) = -1;
  stM(2,14) = 1;
  // R8 reverse
  stM(1,15) = 1;
  stM(2,15) = -1;
  // R9 forward
  stM(0,16) = -1;
  stM(10,16) = 1;
  // R9 reverse
  stM(0,17) = 1;
  stM(10,17) = -1;
  // R10 forward
  stM(3,18) = 1;
  stM(10,18) = -1;
  // R10 reverse  
  stM(3,19) = -1;
  stM(10,19) = 1;
  
  return stM;
}
