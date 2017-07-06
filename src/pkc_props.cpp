#include <Rcpp.h>
using namespace Rcpp;

//' Calculate propensities for the PKC Model.
//'
//' Return the propensity vector of the PKC Model for a given vector of particle numbers. 
//' 
//' @param part_num A numeric vector: the particle numbers of the model species.
//' @param calcium A numeric vector: the calcium particle number.
//' @param rId An integer value: the id of the specified reaction for which the propensity should be calculated.
//' @return A numeric vector containing a cumulative sum of all reaction propensities.
//' @examples
//' pkc_props()
//' @export
// [[Rcpp::export]]
NumericVector pkc_props(NumericVector part_num, double calcium) {
  
  // model parameters
  NumericVector x = part_num;
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
  NumericVector amu(20);
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
  amu[12] = amu[11] + calcium * k13 * (double)x[0]; /* Ca given as conc., hence, no scaling */
  amu[13] = amu[12] + k14 * x[1];
  amu[14] = amu[13] + k15 * DAG * (double)x[1]; /* DAG given as conc., hence, no scaling */
  amu[15] = amu[14] + k16 * x[2];
  amu[16] = amu[15] + k17 * DAG * (double)x[0]; /* DAG given as conc., hence, no scaling */
  amu[17] = amu[16] + k18 * x[10];
  amu[18] = amu[17] + k19 * AA * (double)x[10];  /* AA given as conc., hence, no scaling */
  amu[19] = amu[18] + k20 * x[3];
      
  return amu;
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
