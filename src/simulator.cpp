#include <Rcpp.h>
using namespace Rcpp;

//' Stochastic Simulatior (Gillespie).
//'
//' Simulate a model using an implementation of Gillespie's Direct Method Stochastic Simulation Algorithm.
//'
//' @param time A numeric vector: the times of the observations.
//' @param calcium_conc A numeric vector: the concentration of cytosolic calcium [nmol].
//' @param init_conc A numeric vector: the initial concentrations of the model species.
//' @param calc_props A function: returns a vector of model propensities for a given vector of concentrations.
//' @param provide_stM A function: returns a matrix of the stoichiometric coefficients of the reaction system.
//' @param dt A numeric, the time interval between two output samples.
//' @param vol A numeric, the volume of the system [l].
//' @return A dataframe with time and the active protein time series as columns.
//' @examples
//' simulator()
//' @export
// [[Rcpp::export]]
NumericMatrix simulator(NumericVector time,
                         NumericVector calcium,
                         NumericVector init_conc,
                         Function calc_props,
                         Function provide_stM,
                         double dt,
                         double vol
                         ) {

  /* VARIABLES */
  // fetch stoichiometric matrix
  NumericMatrix stM = provide_stM();
  // create empty vector for propensities (length = no. of cols of the stoich matrix)
  NumericVector amu(stM.ncol());
  // particle number from concentration (nmol/l) factor
  double f;
  f = 6.0221415e14*vol;
  // create initial particle number vector x from initial conc vector ic
  NumericVector ic = init_conc;
  NumericVector x(ic.length());
  int i;
  for (i=0; i < ic.length(); i++) {
    x[i] = (int)floor(ic[i]*f);  
  }
  // get calcium particle numbers
  //NumericVector calcium;
  //int j;
  //for (j=0; j < calcium_conc.length(); j++) {
  //  calcium[j] = (double)floor(calcium_conc[j]*f);
  //}
  //Rcout << "x: " << x << std::endl;
  
  // create algorithm parameters 
  int ntimepoint;
  double r2;
  int noutput;
  double currentTime;
  double outputTime;
  double startTime;
  double endTime;
  double tau;
  unsigned int rIndex;
  int xID;
  // initial algorithm time settings
  startTime = time[0];
  endTime = time[time.length()-1];
  currentTime = startTime;
  noutput = 0;
  ntimepoint = 0;
  outputTime = currentTime;
  int nintervals = (int)floor((endTime-startTime)/dt+0.5)+1;
  // create return value matrix
  NumericMatrix retval(nintervals, 2+x.length());
  
  /* SIMULATION ALGORITHM */
  while (currentTime < endTime) {
    R_CheckUserInterrupt();
    // calcium particle number at current timepoint (ntimepoint)
    double curr_ca;
    curr_ca = calcium[ntimepoint];
    // calculate propensity amu for every reaction
    int n;
    for (n=0; n < stM.ncol(); n++) {
      double a = as<double>(calc_props(x, curr_ca, n));
      if (n == 0) {
        amu[n] = a;
      } else {
        amu[n] = a + amu[n-1];
      }
    }
    
    Rcout << "amu: " << amu << std::endl;
    
    // terminate simulation if sum of propensities is <= 0
    //if (a0 <= 0) {
    //  exit(-1);
    //}
    
    //Rcout << "a0: " << a0 << std::endl;
    
    // calculate time step tau
    tau = - log(runif(1)[0])/amu[amu.length()-1];
    
    //Rcout << "tau: " << tau << std::endl;
    
    // check if reaction time exceeds time until the next observation
    // if true -> update output directly
    if ((currentTime+tau)>=time[ntimepoint+1]) {
      currentTime = time[ntimepoint+1];
      while ((currentTime > outputTime)&&(outputTime < endTime)) {
        retval(noutput, 0) = outputTime;
        retval(noutput, 1) = calcium[ntimepoint];
        for (xID=2; xID < 2+x.length(); xID++) {
          retval(noutput, xID) = x[xID-2]/f;
        }
        noutput++;
        outputTime += dt;
      }
      ntimepoint++;
    } else {
      // if false -> select reaction to fire
      r2 = amu[amu.length()-1] * runif(1)[0];
      
      //Rcout << "amu: " << amu << std::endl;
      //Rcout << "r2: " << r2 << std::endl;
      
      rIndex = 0;
      for (rIndex=0; amu[rIndex] < r2; rIndex++);
      // propagate time
      currentTime += tau;
      // update output
      while ((currentTime > outputTime)&&(outputTime < endTime)) {
        retval(noutput, 0) = outputTime;
        retval(noutput, 1) = calcium[ntimepoint];
        for (xID=2; xID < 2+x.length(); xID++) {
          retval(noutput, xID) = x[xID-2]/f;
        }
        noutput++;
        outputTime += dt;
      }
      
      Rcout << "rIndex: " << rIndex << std::endl;
      
      // update system state
      int species;
      for (species=0; species < stM.nrow();species++) {
        x[species] += stM(species,rIndex);
      }
      
      //Rcout << "rIndex: " << rIndex << std::endl;
      Rcout << "x: " << x << " at " << ntimepoint << std::endl;
      
    }
  }
  // update output
  while (outputTime <= endTime) {
    retval(noutput, 0) = outputTime;
    retval(noutput, 1) = calcium[ntimepoint];
    for (xID=2; xID < 2+x.length(); xID++) {
      retval(noutput, xID) = x[xID-2]/f;
    }
    noutput++;
    outputTime += dt;
  }

  return retval;
}
