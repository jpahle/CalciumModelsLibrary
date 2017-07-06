#include <Rcpp.h>
using namespace Rcpp;

//' Calculate propensities for the Calmodulin Model.
//'
//' Return the propensity vector of the Calmodulin Model for a given vector of particle numbers. 
//' 
//' @param part_num A numeric vector: the particle numbers of the model species.
//' @param calcium A numeric vector: the calcium particle number.
//' @param rId An integer value: the id of the specified reaction for which the propensity should be calculated.
//' @return A numeric vector containing a cumulative sum of all reaction propensities.
//' @examples
//' calmodulin_props()
//' @export
// [[Rcpp::export]]
NumericVector calmodulin_props(NumericVector part_num, double calcium) {
  
  // model parameters
  NumericVector x = part_num;
  double k_on = 0.025;
  double k_off = 0.005;
  double Km = 1.0;
  int h = 4;
  // calculate cumulative propensity vector (amu[0], amu[0] + amu[1], etc.)
  NumericVector amu(2);
  amu[0] = ((k_on * pow((double)calcium,(double)h)) / (pow((double)Km,(double)h) + pow((double)calcium,(double)h))) * x[0];
  amu[1] = amu[0] + k_off * x[1];
  
  return amu;  
}

//' Define stoichiometric matrix of the model
//' 
//' Create and return the stoichiometric matrix of the model as a numeric matrix.
//' 
//' @return A numeric matrix: the stoichiometric matrix
//' @examples
//' calmodulin_stM()
//' @export
// [[Rcpp::export]] 
NumericMatrix calmodulin_stM() {

  NumericMatrix stM(2,2);
  stM(0,0) = -1;
  stM(0,1) = 1;
  stM(1,0) = 1;
  stM(1,1) = -1;
  
  return stM;
}
