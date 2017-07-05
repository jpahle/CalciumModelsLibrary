#include <Rcpp.h>
using namespace Rcpp;

//' Stochastic Simulatior (Gillespie).
//'
//' Simulate a model using an implementation of Gillespie's Direct Method Stochastic Simulation Algorithm.
//'
//' @param time A numeric vector: the times of the observations.
//' @param calcium_conc A numeric vector: the concentration of cytosolic calcium [nmol].
//' @param init_conc A numeric vector: the initial concentrations of the model species.
//' @param calc_props A function: calculates and returns the propensity of a selected reaction given a vector of current particle numbers.
//' @param provide_stM A function: returns a matrix of the stoichiometric coefficients of the reaction system.
//' @param dt A numeric, the time interval between two output samples.
//' @param vol A numeric, the volume of the system [l].
//' @return A dataframe with time, calcium and the active protein time series as columns.
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
  // check user supplied timestep: if too low -> exit simulation
  if (dt < 0.00005) {
    Rcout << "Fatal error: Value of dt too low! Timesteps below the threshold of 0.00005 cause large rounding errors. Simulation aborted!" << std::endl;
    return retval;
  }
    
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
    // calculate time step tau
    tau = - log(runif(1)[0])/amu[amu.length()-1];
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
      // update system state
      int species;
      for (species=0; species < stM.nrow();species++) {
        x[species] += stM(species,rIndex);
      }
    }
  }
  // update output
  while (floor(outputTime*10000) <= floor(endTime*10000)) {
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
