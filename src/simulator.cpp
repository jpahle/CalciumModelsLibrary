#include "extern_simulator_func_prototype.hpp"
#include <Rcpp.h>
using namespace Rcpp;


// The general simulator file is included in every C++ model file -> if they provide a MODEL_NAME, all functions will be renamed in a model specific fashion
#ifdef MODEL_NAME
  #define Map_helper(x,y) x##y
  #define Map(x,y) Map_helper(x,y)
  #define simulator Map(simulator_, MODEL_NAME)
  #define init Map(init_, MODEL_NAME)
  #define calculate_amu Map(calculate_amu_, MODEL_NAME)
  #define get_stM Map(get_stM_, MODEL_NAME)
  
  // Placeholder init function since the R Wrapper Function tries to call it before its 'real' definition in the C++ model file
  List init();
#endif


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
// Global shared functions
extern void calculate_amu();
extern NumericMatrix get_stM();
extern void update_system(unsigned int rIndex);


//' Stochastic Simulator (Gillespie's Direct Method).
//'
//' Simulate a calcium dependent protein coupled to an input calcium time series using an implementation of Gillespie's Direct Method SSA.
//'
//' @param user_input_df A data frame: contains the times of the observations (column "time") and the cytosolic calcium concentration [nmol/l] (column "Ca").
//' @param user_sim_params A numeric vector: contains all simulation parameters ("timestep": the time interval between two output samples, "endTime": the time at which to end the simulation).
//' @param default_vols A numeric vector: contains updated default values of all volumes [l].
//' @param default_init_conc A numeric vector: contains updated default values of all initial concentrations [nmol/l].
//' @return A dataframe with time and the active protein time series as columns.
//' @examples
//' simulator()
DataFrame simulator(DataFrame user_input_df,
                   NumericVector user_sim_params,
                   NumericVector default_vols,
                   NumericVector default_init_conc) {

  // get R random generator state
  GetRNGstate();
  
  /* VARIABLES */
  // Get parameter values from arguments
  timevector = user_input_df["time"];
  calcium = user_input_df["Ca"];
  timestep = user_sim_params["timestep"];
  // currently only works with ONE volume!
  // (code takes only first volume entry)
  vol = default_vols[0];
  NumericVector ic = default_init_conc; 
  // Particle number <-> concentration (nmol/l) factor (n/f = c <=> c*f = n)
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
  endTime = user_sim_params["endTime"];
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
    // Calculate propensity amu for every reaction
    calculate_amu();
    // Calculate time step tau
    tau = - log(runif(1)[0])/amu[nreactions-1];
    // Check if reaction time exceeds time until the next observation 
    if ((currentTime+tau)>=timevector[ntimepoint+1]) {
      // Set current simulation time to next timepoint in input calcium time series
      currentTime = timevector[ntimepoint+1];
      // Update output
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
      // Select reaction to fire
      r2 = amu[nreactions-1] * runif(1)[0];
      rIndex = 0;
      for (rIndex=0; amu[rIndex] < r2; rIndex++);
      // Propagate time
      currentTime += tau;
     // Update output
      while ((currentTime > outputTime)&&(outputTime < endTime)) {
        retval(noutput, 0) = outputTime;
        retval(noutput, 1) = calcium[ntimepoint];
        for (xID=2; xID < 2+nspecies; xID++) {
          retval(noutput, xID) = x[xID-2]/f;
        }
        noutput++;
        outputTime += timestep;
      }
      // Update system state
      // get stoich matrix
      NumericMatrix stM = get_stM();
      // get selected reaction (a stM column vector) 
      NumericVector selected_re = stM(_, rIndex);
      // add stoich coefficients (either -1, 0 or 1) of column vector to x
      for (int k = 0; k < nspecies; k++) {
        x[k] = x[k] + selected_re[k];
      }
    }
  }
  // Update output
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
  
  // Send random generator state back to R
  PutRNGstate();
  
  // Convert NumericMatrix retval to DataFrame
  DataFrame df_retval(retval);
  
  return df_retval;
}
