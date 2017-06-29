#include <Rcpp.h>
using namespace Rcpp;

//' Calculate propensity for a Calmodulin Model reaction.
//'
//' Return the propensity of a Calmodulin Model reaction for a given vector of particle numbers and a reaction Id. 
//' 
//' @param part_num A numeric vector: the particle numbers of the model species.
//' @param calcium A numeric vector: the calcium particle number.
//' @param rId An integer value: the id of the specified reaction for which the propensity should be calculated.
//' @return A double value (the propensity of the specified reaction).
//' @examples
//' calmodulin_props()
//' @export
// [[Rcpp::export]]
double calmodulin_props(NumericVector part_num, double calcium, int rId) {
  
  // model parameters
  NumericVector x = part_num;
  double k_on = 0.025;
  double k_off = 0.005;
  double Km = 1.0;
  int h = 4;
  
  // calculate propensity of selected reaction
  double a;
  switch (rId) {
  case 0:
    a = ((k_on * pow((double)calcium,(double)h)) / (pow((double)Km,(double)h) + pow((double)calcium,(double)h))) * x[0];
    break;
  case 1:
    a =  k_off * x[1];
    break;
  default:
    printf("\nError in propensity calculation: reaction Index (%u) out of range!\n", rId);
    a = 0;
  }
  return a;
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
