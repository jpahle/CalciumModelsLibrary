#include "global_vars.hpp"
#include <Rcpp.h>
using namespace Rcpp;

/* Global Variables */
NumericVector calcium;
unsigned int ntimepoint;
double *amu;
unsigned long long int *x;


//' Couple a simulated PKC protein to a given calcium time series.
//'
//' Takes a calcium time series and simulates the Ca-dependent protein PKC.
//'
//' Implementation based on the PKC Model by Manninnen (2006)
//'
//' @param param_time A numeric vector: the times of the observations.
//' @param param_calcium A numeric vector: the concentrations of cytosolic calcium [nmol/l].
//' @param param_timestep A double, the time interval between two output samples.
//' @param param_vol A double, the volume of the system [l].
//' @param param_init_conc A numeric vector: the initial concentrations of the model species [nmol/l].
//' @param calculate_amu A function: provides parameter values and result of propensity equations.
//' @param update_system A function: contains stoichiometric effects of all reactions.
//' @return A dataframe with time and the active protein time series as columns.
//' @examples
//' simulator2()
//' @export
// [[Rcpp::export]]
NumericMatrix simulator2(NumericVector param_time,
                   NumericVector param_calcium,
                   double param_timestep,
                   double param_vol,
                   NumericVector param_init_conc,
                   Function calc_propensities,
                   Function update_system_state
                   ) {
  
  /* get R random generator state */
  GetRNGstate();
  
  /* Memory allocation */
  amu = (double *)Calloc(nreactions, double);
  x = (unsigned long long int *)Calloc(nspecies, unsigned long long int);

  
  /* Variables */
  NumericVector timevector = param_time;
  calcium = param_calcium;
  double timestep = param_timestep;
  // check user supplied timestep: if too low -> exit simulation
  if (timestep < 0.00005) {
    Rcout << "Fatal error: Value of dt too low! Timesteps below the threshold of 0.00005 cause large rounding errors. Simulation aborted!" << std::endl;
    return 0;
  }
  double vol = param_vol;
  // Particle number to concentration (nmol/l) factor
  double f;
  f = 6.0221415e14*vol;
  // Initial particle number vector x from initial conc vector ic
  NumericVector ic = param_init_conc;
  int i;
  for (i=0; i < ic.length(); i++) {
    x[i] = (unsigned long long int)floor(ic[i]*f);  
  }
  int xID;
  // Control variables
  ntimepoint = 0;
  int noutput;
  noutput = 0;
  // Variables for random steps
  double tau;
  double r2;
  unsigned int rIndex;
  // Time variables
  double startTime;
  double endTime;
  double currentTime;
  double outputTime;
  startTime = timevector[0];
  endTime = timevector[timevector.length()-1];
  currentTime = startTime;
  outputTime = currentTime;
  // Return value
  int nintervals = (int)floor((endTime-startTime)/timestep+0.5)+1;
  NumericMatrix retval(nintervals, 2+nspecies);

  
  /* Simulation loop */
  while (currentTime < endTime) {
    // R_CheckUserInterrupt();
    // calculate propensity amu for every reaction
    calc_propensities();
    // calculate time step tau
    tau = - log(runif(1)[0])/amu[nreactions-1];
    // check if reaction time exceeds time until the next observation
    // and update output
    if ((currentTime+tau)>=timevector[ntimepoint+1]) {
      currentTime = timevector[ntimepoint+1];
      while ((currentTime > outputTime)&&(outputTime < endTime)) {
        retval(noutput, 0) = outputTime;
        retval(noutput, 1) = calcium[ntimepoint];
        for (xID=2; xID < 2+nspecies; xID++) {
          retval(noutput, xID) = x[xID-2]/f;
        }
        noutput++;
        outputTime += timestep;
      }
      ntimepoint++;
    } else {
      // select reaction to fire
      r2 = amu[nreactions-1] * runif(1)[0];
      rIndex = 0;
      for (rIndex=0; amu[rIndex] < r2; rIndex++);
      // propagate time
      currentTime += tau;
      // update output
      while ((currentTime > outputTime)&&(outputTime < endTime)) {
        retval(noutput, 0) = outputTime;
        retval(noutput, 1) = calcium[ntimepoint];
        for (xID=2; xID < 2+nspecies; xID++) {
          retval(noutput, xID) = x[xID-2]/f;
        }
        noutput++;
        outputTime += timestep;
      }
      // update system state
      update_system_state(rIndex);
    }
  }
  // update output
  while (outputTime <= endTime) {
    retval(noutput, 0) = outputTime;
    retval(noutput, 1) = calcium[ntimepoint];
    for (xID=2; xID < 2+nspecies; xID++) {
      retval(noutput, xID) = x[xID-2]/f;
    }
    noutput++;
    outputTime += timestep;
  }
 
  // Free dyn. allocated pointers
  Free(amu);
  Free(x);
    
  /* send random generator state back to R*/
  PutRNGstate();
  
  return retval;
}
