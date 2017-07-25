#include "extern_simulator_func_prototype.hpp"
#include <Rcpp.h>
using namespace Rcpp;


// Global shared variables
extern NumericVector timevector;
extern double timestep;
extern double vol;
extern NumericVector calcium;
extern unsigned int ntimepoint;
extern double *amu;
extern unsigned long long int *x;
extern int nspecies;
extern int nreactions;


extern void calculate_amu();
extern void update_system(unsigned int rIndex);


//' Couple a simulated Ca-dependent protein to a given calcium time series.
//'
//' Takes a calcium time series and simulates the coupled Ca-dependent protein.
//'
//' @param param_input_df A data frame: contains the times of the observations (column "time") and the cytosolic calcium concentration [nmol/l] (column "Ca").
//' @param param_sim_params A numeric vector: contains all simulation parameters (timestep = the time interval between two output samples, endTime = the time at which to end the simulation).
//' @param param_vol A numeric: the volume of the system [l].
//' @param param_init_conc A numeric vector: the initial concentrations of model species [nmol/l].
//' @return A dataframe with time and the active protein time series as columns.
//' @examples
//' simulator()
NumericMatrix simulator(DataFrame param_input_df,
                   NumericVector param_sim_params,
                   double param_vol,
                   NumericVector param_init_conc) {

  // get R random generator state
  GetRNGstate();
  
  /* VARIABLES */
  // get parameter values from arguments
  // timevector = param_time;
  // calcium = param_calcium;
  timevector = param_input_df["time"];
  calcium = param_input_df["Ca"];
  timestep = param_sim_params["timestep"];
  vol = param_vol;
  NumericVector ic = param_init_conc;
  // particle number <-> concentration (nmol/l) factor (n/f = c <=> c*f = n)
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
  //endTime = timevector[timevector.length()-1];
  endTime = param_sim_params["endTime"];
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
    // pkc_calculate_amu();
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
      // pkc_update_system(rIndex);
      update_system(rIndex);
    }
  }
  // update output
  while (floor(outputTime*10000) <= floor(endTime*10000)) {
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
