#include <Rcpp.h>
using namespace Rcpp;

//' Calculate propensities for the CamKII Model.
//'
//' Return the propensity vector of the CamKII Model for a given vector of particle numbers. 
//' 
//' @param part_num A numeric vector: the particle numbers of the model species.
//' @param calcium A numeric vector: the calcium particle number.
//' @param rId An integer value: the id of the specified reaction for which the propensity should be calculated.
//' @return A numeric vector containing a cumulative sum of all reaction propensities.
//' @examples
//' camkii_props()
//' @export
// [[Rcpp::export]]
NumericVector camkii_props(NumericVector part_num, double calcium) {
  
  // model parameters
  NumericVector x = part_num;
  double a = -0.22;
  double b = 1.826;
  double c = 0.1;
  double k_IB = 0.01;
  double k_BI = 0.8;
  double k_PT = 1;
  double k_TP = 1e-12;
  double k_TA = 0.0008;
  double k_AT = 0.01;
  double k_AA = 0.29;
  double c_I = 0;
  double c_B = 0.75;
  double c_P = 1;
  double c_T = 0.8;
  double c_A = 0.8;
  double camT = 1000;
  double Kd = 1000;
  double Vm_phos = 0.005;
  double Kd_phos = 0.3;
  double totalC = 40;
  int h = 4;
  // calculate propensity of selected reaction
  NumericVector amu(10);
  amu[0] = x[0] * ((k_IB * camT * pow((double)calcium,(double)h)) / (pow((double)calcium,(double)h) + pow((double)Kd,(double)h)));
  amu[1] = amu[0] + k_BI * x[1];
  
  double activeSubunits = (x[1] + x[2] + x[3] + x[4]) / totalC;
  double prob =  a * activeSubunits + b*(pow((double)activeSubunits,(double)2)) + c*(pow((double)activeSubunits,(double)3));
  amu[2] = amu[1] +  totalC * k_AA * prob * ((c_B * x[1]) / pow((double)totalC,(double)2)) * (2*c_B*x[1] + c_P*x[2] + c_T*x[3]+ c_A*x[4]);
  
  amu[3] = amu[2] + k_PT * x[2];
  amu[4] = amu[3] + k_TP * x[3] * pow((double)calcium,(double)h);
  amu[5] = amu[4] + k_TA * x[3];
  amu[6] = amu[5] + k_AT * x[4] * (camT - ((camT * pow((double)calcium,(double)h)) / (pow((double)calcium,(double)h) + pow((double)Kd,(double)h))));
  amu[7] = amu[6] + ((Vm_phos * x[2]) / (Kd_phos + (x[2] / totalC)));
  amu[8] = amu[7] + ((Vm_phos * x[3]) / (Kd_phos + (x[3] / totalC)));
  amu[9] = amu[8] + ((Vm_phos * x[4]) / (Kd_phos + (x[4] / totalC)));

  return amu;
}

//' Define stoichiometric matrix of the CamKII model
//' 
//' Create and return the stoichiometric matrix of the CamKII model as a numeric matrix.
//' 
//' @return A numeric matrix: the stoichiometric matrix
//' @examples
//' camkii_stM()
//' @export
// [[Rcpp::export]] 
NumericMatrix camkii_stM() {

  // initialize Stoich. Matrix filled with zeroes
  NumericMatrix stM(5,10);
  // Ca-Cam binding
  stM(0,0) = -1;
  stM(1,0) = 1;
  // binding reverse
  stM(0,1) = 1;
  stM(1,1) = -1;
  // phosphorylation
  stM(1,2) = -1;
  stM(2,2) = 1;
  // trapping
  stM(2,3) = -1;
  stM(3,3) = 1;
  // trapping reverse
  stM(2,4) = 1;
  stM(3,4) = -1;
  // autonomous
  stM(3,5) = -1;
  stM(4,5) = 1;
  // autonomous reverse
  stM(3,6) = 1;
  stM(4,6) = -1;
  // phosphatase on W_P
  stM(1,7) = 1;
  stM(2,7) = -1;
  // phosphatase on W_T
  stM(1,8) = 1;
  stM(3,8) = -1;
  // phosphatase on W_A
  stM(0,9) = 1;
  stM(4,9) = -1;
    
  return stM;
}
