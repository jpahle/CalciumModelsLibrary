#include "CaModLib_global_vars.hpp"
#include <Rcpp.h>
using namespace Rcpp;


/* Global variables */
NumericVector timevector;
double timestep;
double vol;
double PKCinact0_conc;
double PKCbasal0_conc;


//' Couple a simulated Ca-dependent protein to a given calcium time series.
//'
//' Takes a calcium time series and simulates the coupled Ca-dependent protein.
//'
//' @param param_time A numeric vector: the times of the observations.
//' @param param_calcium A numeric vector: the concentrations of cytosolic calcium [nmol/l].
//' @param param_timestep A numeric, the time interval between two output samples.
//' @param param_vol A numeric, the volume of the system [l].
//' @param param_init_conc A numeric vector: the initial concentrations of model species [nmol/l].
//' @return A dataframe with time and the active protein time series as columns.
//' @examples
//' simulator()
//' @export
// [[Rcpp::export]]
NumericMatrix simulator(NumericVector param_time,
                   NumericVector param_calcium,
                   double param_timestep,
                   double param_vol,
                   NumericVector param_init_conc) {

  // get R random generator state
  GetRNGstate();
  
  /* VARIABLES */
    
  // get parameter values from arguments
  timevector = param_time;
  calcium = param_calcium;
  timestep = param_timestep;
  vol = param_vol;
  NumericVector ic = param_init_conc;
  // particle number to concentration (nmol/l) factor
  double f;
  f = 6.0221415e14*vol;
  // Control variables
  int noutput;
  ntimepoint = 0;
  noutput = 0;
  int xID;
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
  NumericMatrix retval(nintervals, nspecies+2); // nspecies+2 because time and calcium
  // Memory allocation for propensity and particle number pointers
  amu = (double *)Calloc(nreactions, double);
  x = (unsigned long long int *)Calloc(nspecies, unsigned long long int);
  // Initial particle numbers
  int i;
  for (i=0; i < ic.length(); i++) {
    x[i] = (unsigned long long int)floor(ic[i]*f);  
  }
  
  /* SIMULATION LOOP */
  while (currentTime < endTime) {
    R_CheckUserInterrupt();
    // calculate propensity amu for every reaction
    calculate_amu();
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
      update_system(rIndex);
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
