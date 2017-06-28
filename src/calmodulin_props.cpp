#include <Rcpp.h>
using namespace Rcpp;

//' Calculate propensity for a Calmodulin Model reaction.
//'
//' Return the propensity of a Calmodulin Model reaction for a given vector of particle numbers and a reaction Id. 
//' 
//' @param x A numeric vector: the particle numbers of the model species.
//' @param calcium A numeric vector: the calcium particle number.
//' @param rId An integer value: the id of the specified reaction for which the propensity should be calculated.
//' @return A double value (the propensity of the specified reaction).
//' @examples
//' calmodulin_props()
//' @export
// [[Rcpp::export]]
double calmodulin_props(NumericVector x, double calcium, int rId) {
  
  // model parameters
  double k_on = 0.025;
  double k_off = 0.005;
  double Km = 1.0;
  double E0_conc = 5;
  int h = 4;
  
  // calculate propensity of selected reaction
  double amu;
  switch (rId) {
  case 1:
    amu = ((k_on * pow((double)calcium,(double)h)) / (pow((double)Km,(double)h) + pow((double)calcium,(double)h))) * x[0];
  case 2:
    amu =  k_off * x[1] + amu[0];
  default:
    amu = 0;
  }
  return amu;
}

