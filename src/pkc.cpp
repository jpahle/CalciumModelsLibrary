#include <Rcpp.h>
using namespace Rcpp;

//' Global variables
NumericVector timevector;
NumericVector calcium;
double timestep;
double vol;
double k1;
double k2;
double k3;
double k4;
double k5;
double k6;
double k7;
double k8;
double k9;
double k10;
double k11;
double k12;
double k13;
double k14;
double k15;
double k16;
double k17;
double k18;
double k19;
double k20;
double AA;   // given as conc. remains fixed throughout the simulation
double DAG;  // given as conc. remains fixed throughout the simulation
double PKCinact0_conc;
double PKCbasal0_conc;

unsigned int ntimepoint;
double *amu;
unsigned long long int *x;

//' Propensity Calculation
//'
//' Calculates the propensities of all PKC model reactions and stores them in the vector amu.
//'
//' @param None
//' @return None
//' @examples
//' calculate_amu() 
//' @export
// [[Rcpp::export]]
void calculate_amu() {

  amu[0] = k1 * x[0];
  amu[1] = amu[0] + k2 * x[5];
  amu[2] = amu[1] + k3 * AA * (double)x[0]; /* AA given as conc., hence, no scaling */
  amu[3] = amu[2] + k4 * x[6];
  amu[4] = amu[3] + k5 * x[1];
  amu[5] = amu[4] + k6 * x[7];
  amu[6] = amu[5] + k7 * AA * (double)x[1];  /* AA given as conc., hence, no scaling */
  amu[7] = amu[6] + k8 * x[8];
  amu[8] = amu[7] + k9 * x[2];
  amu[9] = amu[8] + k10 * x[9];
  amu[10] = amu[9] + k11 * x[3];
  amu[11] = amu[10] + k12 * x[4];
  amu[12] = amu[11] + calcium[ntimepoint] * k13 * (double)x[0]; /* Ca given as conc., hence, no scaling */
  amu[13] = amu[12] + k14 * x[1];
  amu[14] = amu[13] + k15 * DAG * (double)x[1]; /* DAG given as conc., hence, no scaling */
  amu[15] = amu[14] + k16 * x[2];
  amu[16] = amu[15] + k17 * DAG * (double)x[0]; /* DAG given as conc., hence, no scaling */
  amu[17] = amu[16] + k18 * x[10];
  amu[18] = amu[17] + k19 * AA * (double)x[10];  /* AA given as conc., hence, no scaling */
  amu[19] = amu[18] + k20 * x[3];
}

//' System Update
//'
//' Changes the system state (updates the particle numbers) by instantiating a chosen reaction.
//'
//' @param rIndex An unsigned integer: the id of the chosen reaction.
//' @return void
//' @examples
//' update_system() 
//' @export
// [[Rcpp::export]]
void update_system(unsigned int rIndex) {
  
  switch (rIndex) {
  case 0:   /* R1 */
    x[0]--;
    x[5]++;
    break;
  case 1:
    x[5]--;
    x[0]++;
    break;
  case 2:   /* R2 */
    x[0]--;
    x[6]++;
    break;
  case 3:
    x[6]--;
    x[0]++;
    break;
  case 4:  /* R3 */
    x[1]--;
    x[7]++;
    break;
  case 5:
    x[7]--;
    x[1]++;
    break;
  case 6:  /* R4 */
    x[1]--;
    x[8]++;
    break;
  case 7:
    x[8]--;
    x[1]++;
    break;
  case 8: /* R5 */
    x[2]--;
    x[9]++;
    break;
  case 9:
    x[9]--;
    x[2]++;
    break;
  case 10:/* R6 */
    x[3]--;
    x[4]++;
    break;
  case 11:
    x[4]--;
    x[3]++;
    break;
  case 12:/* R7 */
    x[0]--;
    x[1]++;
    break;
  case 13:
    x[1]--;
    x[0]++;
    break;
  case 14:/* R8 */
    x[1]--;
    x[2]++;
    break;
  case 15:
    x[2]--;
    x[1]++;
    break;
  case 16:/* R9 */
    x[0]--;
    x[10]++;
    break;
  case 17:
    x[10]--;
    x[0]++;
    break;
  case 18:/* R10 */
    x[10]--;
    x[3]++;
    break;
  case 19:
    x[3]--;
    x[10]++;
    break;
    printf("\nError in updateSystem(): rIndex (%u) out of range!\n", rIndex);
    exit(-1);
  }
}

//' Couple a simulated PKC protein to a given calcium time series.
//'
//' Takes a calcium time series and simulates the Ca-dependent protein PKC.
//'
//' Implementation based on the PKC Model by Manninnen (2006)
//'
//' @param timeseries A numeric vector: the times of the observations.
//' @param _calcium A numeric vector: the concentrations of cytosolic calcium [nmol].
//' @param _dt A numeric, the time interval between two output samples.
//' @param _vol A numeric, the volume of the system [l].
//' @param k1-20 Double values, reaction rate parameters [1/s].
//' @param _AA A double, concentration of AA [nmol].
//' @param _DAG A double, concentration of DAG [nmol].
//' @param _PKCinact0_conc A double, initial concentration of inactive PKC [nmol].
//' @param _PKCbasal0_conc A double, initial concentration of basal PKC [nmol].
//' @return A dataframe with time and the active protein time series as columns.
//' @examples
//' pkc()
//' @export
// [[Rcpp::export]]
NumericMatrix pkc(NumericVector param_time,
                   NumericVector param_calcium,
                   double param_timestep,
                   double param_vol,
                   double param_k1,
                   double param_k2,
                   double param_k3,
                   double param_k4,
                   double param_k5,
                   double param_k6,
                   double param_k7,
                   double param_k8,
                   double param_k9,
                   double param_k10,
                   double param_k11,
                   double param_k12,
                   double param_k13,
                   double param_k14,
                   double param_k15,
                   double param_k16,
                   double param_k17,
                   double param_k18,
                   double param_k19,
                   double param_k20,
                   double param_AA,
                   double param_DAG,
                   double param_PKCinact0_conc,
                   double param_PKCbasal0_conc) {

  /* VARIABLES */
  
  /* get parameter values from arguments */
  timevector = param_time;
  calcium = param_calcium;
  vol = param_vol;
  timestep = param_timestep;
  k1 = param_k1;
  k2 = param_k2;
  k3 = param_k3;
  k4 = param_k4;
  k5 = param_k5;
  k6 = param_k6;
  k7 = param_k7;
  k8 = param_k8;
  k9 = param_k9;
  k10 = param_k10;
  k11 = param_k11;
  k12 = param_k12;
  k13 = param_k13;
  k14 = param_k14;
  k15 = param_k15;
  k16 = param_k16;
  k17 = param_k17;
  k18 = param_k18;
  k19 = param_k19;
  k20 = param_k20;
  AA = param_AA;
  DAG = param_DAG;
  PKCinact0_conc = param_PKCinact0_conc;
  PKCbasal0_conc = param_PKCbasal0_conc;
  
  // Model specifications
  unsigned int nreactions = 20;
  unsigned int nspecies = 11;
  
  // particle number to concentration (nmol/l) factor
  double f;
  f = 6.0221415e14*vol;
  
  // Control variables
  int noutput;
  ntimepoint = 0;
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
  NumericMatrix retval(nintervals, 13);
  
  /* Memory allocation */
  amu = (double *)Calloc(nreactions, double);
  x = (unsigned long long int *)Calloc(nspecies, unsigned long long int);
  
  // Initial conditions
  x[0] = (unsigned long long int)floor(PKCinact0_conc*f);
  x[1] = 0;
  x[2] = 0;
  x[3] = 0;
  x[4] = 0;
  x[5] = (unsigned long long int)floor(PKCbasal0_conc*f);
  x[6] = 0;
  x[7] = 0;
  x[8] = 0;
  x[9] = 0;
  x[10] = 0;
    
  // Simulation loop
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
        retval(noutput, 1) = x[0]/f;
        retval(noutput, 2) = x[1]/f;
        retval(noutput, 3) = x[2]/f;
        retval(noutput, 4) = x[3]/f;
        retval(noutput, 5) = x[4]/f;
        retval(noutput, 6) = x[5]/f;
        retval(noutput, 7) = x[6]/f;
        retval(noutput, 8) = x[7]/f;
        retval(noutput, 9) = x[8]/f;
        retval(noutput, 10) = x[9]/f;
        retval(noutput, 11) = x[10]/f;
        retval(noutput, 12) = calcium[ntimepoint];
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
        retval(noutput, 1) = x[0]/f;
        retval(noutput, 2) = x[1]/f;
        retval(noutput, 3) = x[2]/f;
        retval(noutput, 4) = x[3]/f;
        retval(noutput, 5) = x[4]/f;
        retval(noutput, 6) = x[5]/f;
        retval(noutput, 7) = x[6]/f;
        retval(noutput, 8) = x[7]/f;
        retval(noutput, 9) = x[8]/f;
        retval(noutput, 10) = x[9]/f;
        retval(noutput, 11) = x[10]/f;
        retval(noutput, 12) = calcium[ntimepoint];
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
    retval(noutput, 0) = outputTime;
    retval(noutput, 1) = x[0]/f;
    retval(noutput, 2) = x[1]/f;
    retval(noutput, 3) = x[2]/f;
    retval(noutput, 4) = x[3]/f;
    retval(noutput, 5) = x[4]/f;
    retval(noutput, 6) = x[5]/f;
    retval(noutput, 7) = x[6]/f;
    retval(noutput, 8) = x[7]/f;
    retval(noutput, 9) = x[8]/f;
    retval(noutput, 10) = x[9]/f;
    retval(noutput, 11) = x[10]/f;
    retval(noutput, 12) = calcium[ntimepoint];
    noutput++;
    outputTime += timestep;
  }
  // Free dyn. allocated pointers
  Free(amu);
  Free(x);
  
  return retval;
  }
