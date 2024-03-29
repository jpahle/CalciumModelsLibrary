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
extern double f;
// Global shared functions
extern void calculate_amu();
extern NumericMatrix get_stM();
extern void update_system(unsigned int rIndex);
 

//' Stochastic Simulator (Gillespie's Direct Method).
//'
//' Simulate a calcium dependent protein coupled to an input calcium time series using an implementation of Gillespie's Direct Method SSA.
//'
//' @param user_input_df A data frame: contains the times of the observations (column "time") and the cytosolic calcium concentration [nmol/l] (column "Ca").
//' @param user_sim_params A List: contains parameters defining the simulation output times 
//'                        (can either be a) a user supplied vector with sim output time points or b) parameters to generate an evenly spaced sim output times vector: 
//'                        "timestep": the time interval between two output samples, "endTime": the time at which to end the simulation and its output).
//' @param default_vols A numeric vector: contains updated default values of all volumes [l].
//' @param default_init_conc A numeric vector: contains updated default values of all initial concentrations [nmol/l].
//' @return A dataframe with time and the active protein time series as columns.
//' @examples
//' simulator()
DataFrame simulator(DataFrame user_input_df,
                    List user_sim_params,
                    NumericVector default_vols,
                    NumericVector default_init_conc) {

  // get R random generator state
  GetRNGstate();
  
  
  
  /* VARIABLES */
  // ------------ Read input calcium signal data frame ------------
  calcium = user_input_df["Ca"];
  timevector = user_input_df["time"];
  //  ------------ Define sim output times: ------------
  // 1.) sim output times can be generated from timestep and endTime (evenly spaced)
  // (use default sim output params if none are supplied by user)
  int timestep_set = 0;
  if (user_sim_params.containsElementNamed("timestep")) {
    timestep = user_sim_params["timestep"];
    timestep_set = 1;
  } else {
    timestep = 0.01;
    timestep_set = 1;
  }
  double endTime;
  int endTime_set = 0;
  if (user_sim_params.containsElementNamed("endTime")) {
    endTime = user_sim_params["endTime"];
    endTime_set = 1;
  } else {
    endTime = 100;
    endTime_set = 1;
  }
  // 2.) sim output times can be supplied as vector by user (even or unevenly spaced) 
  // initialize the sim output times vector (so that it and its length exist, even if its not used)
  NumericVector user_output_times_vector(1);
  int user_output_times_set = 0;
  if (user_sim_params.containsElementNamed("outputTimes")) {
    user_output_times_vector = user_sim_params["outputTimes"];
    // set flag to use custom user supplied sim output time vector (if available)
    user_output_times_set = 1;
  }
  // ------------ Memory allocation for propensity and particle number pointers ------------
  amu = (double *)calloc(nreactions, sizeof(double));
  x = (unsigned long long int *)calloc(nspecies, sizeof(unsigned long long int));
  // ------------ Conversion from concentration (nmol/l) to particle numbers (factor: n/f = c <=> c*f = n) ------------
  vol = default_vols[0];
  // initial concentration vector
  NumericVector ic = default_init_conc; 
  // conversion factor
  f = 6.0221415e14*vol;
  int i;
  for (i=0; i < ic.length(); i++) {
    x[i] = (unsigned long long int)floor(ic[i]*f);  
  }
  // ------------ Control variables ------------
  int noutput;
  noutput = 0;
  ntimepoint = 0;
  int xID;
  // ------------ Variables for random steps ------------
  double tau;
  double r2;
  unsigned int rIndex;
  // ------------ Time variables ------------
  double startTime;
  double currentTime;
  double outputTime;
  startTime = timevector[0];
  currentTime = startTime;
  outputTime = currentTime;
  // ------------ Define return value (numeric matrix; no. of rows = no. of output time point; no. of cols. = time + ca + no. of species) ------------
  // 1.) timestep and endTime are used to generate a number (nintervals) of evenly spaced intervals 
  int nintervals = 0;
  if (timestep_set == 1 && endTime_set == 1) {
    nintervals = (int)floor((endTime-startTime)/timestep+0.5)+1;
  }
  // 2.) take number of intervals from user supplied sim output times vector (can be unevenly spaced -> different timestep lengths)
  // define the timestep_vector outside of if statement so that it is available to the simulation loop
  // since user_output_times_vector is initialized with length 1 earlier, both the vector and its length exist
  // therefore: even if timestep_vector is not used, the definition does not fail (without init. of user_output_times_vector, length would be negative/undefined -> error!)
  NumericVector timestep_vector(user_output_times_vector.length());
  if (user_output_times_set == 1) {
    // no. of intervals is equal to the length of the sim output vector (intervals can be of different sizes)
    nintervals = user_output_times_vector.length();
    // endTime taken from sim output times vector (indexing starts with 0 -> last entry: vector length - 1)
    endTime = user_output_times_vector[user_output_times_vector.length()-1];
    // vector with timesteps derived from sim output time vector
    // For a vector a = [1,2,3,10,87,...], the intervals between its items are given by a[2:end] - a[1:(end-1)]
    int id;
    for (id=0; id < timestep_vector.length(); id++) {
      timestep_vector[id] = fabs(user_output_times_vector[id+1] - user_output_times_vector[id]);
    }
  }
  NumericMatrix retval(nintervals, nspecies+2); // nspecies+2 because time and calcium
  
  
  
  
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
        if (user_output_times_set == 1) {
          outputTime += timestep_vector[noutput];
        } else {
          outputTime += timestep;
        }
        noutput++;
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
        if (user_output_times_set == 1) {
          outputTime += timestep_vector[noutput];
        } else {
          outputTime += timestep;
        }
        noutput++;
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
    if (user_output_times_set == 1) {
      outputTime += timestep_vector[noutput];
    } else {
      outputTime += timestep;
    }
    noutput++;
  }
     
  // Free dyn. allocated pointers
  free(amu);
  free(x);
  
  // Send random generator state back to R
  PutRNGstate();
  
  // Convert NumericMatrix retval to DataFrame
  DataFrame df_retval(retval);
  
  return df_retval;
}