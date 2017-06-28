#include <Rcpp.h>
using namespace Rcpp;

//' Stochastic Simulatior (Gillespie).
//'
//' Simulate a model using an implementation of Gillespie's SSA.
//'
//' @param time A numeric vector: the times of the observations.
//' @param calcium A numeric vector: the particle numbers of cytosolic calcium [nmol].
//' @param init_conc A numeric vector: the initial concentrations of the model species.
//' @param model_props A function: returns a vector of model propensities for a given vector of concentrations.
//' @param a numeric matrix: the stoichiometric coefficients of the reaction system.
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
                         NumericMatrix stM,
                         double dt,
                         double vol,
                         ) {

  /* VARIABLES */
    // create empty vector for propensities (length = no. of cols of the stoich matrix)
  NumericVector amu(stM.ncol());
  // particle number from concentration (nmol/l) factor
  double f;
  f = 6.0221415e14*vol;
  // create initial particle number vector x
  NumericVector x(init_conc.length())
  for (i=0; i < init_con.length(); i++) {
    x[i] = (int)floor(init_conc*f)  
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
    for (n=1; n < stM.ncol(); n++) {
      amu[n] = calc_props(x, curr_ca, n)
    }
    //amu[0] = ((k_on * pow((double)calcium[ntimepoint],(double)h)) / (pow((double)Km,(double)h) + pow((double)calcium[ntimepoint],(double)h))) * x[0];
    //amu[1] =  k_off * x[1] + amu[0];
    // calculate time step tau
    tau = - log(runif(1)[0])/amu[1];
    // check if reaction time exceeds time until the next observation
    // and update output
    if ((currentTime+tau)>=time[ntimepoint+1]) {
      currentTime = time[ntimepoint+1];
      while ((currentTime > outputTime)&&(outputTime < endTime)) {
        retval(noutput, 0) = outputTime;
        retval(noutput, 1) = x[0]/f;
        retval(noutput, 2) = x[1]/f;
        retval(noutput, 3) = calcium[ntimepoint];
        noutput++;
        outputTime += dt;
      }
      ntimepoint++;
    } else {
      // select reaction to fire
      r2 = amu[1] * runif(1)[0];
      rIndex = 0;
      for (rIndex=0; amu[rIndex] < r2; rIndex++);
      // propagate time
      currentTime += tau;
      // update output
      while ((currentTime > outputTime)&&(outputTime < endTime)) {
        retval(noutput, 0) = outputTime;
        retval(noutput, 1) = x[0]/f;
        retval(noutput, 2) = x[1]/f;
        retval(noutput, 3) = calcium[ntimepoint];
        noutput++;
        outputTime += dt;
      }
      // update system state
      switch (rIndex) {
      case 0: // protein activation
        x[0]--;
        x[1]++;
        break;
      case 1: // protein deactivation
        x[0]++;
        x[1]--;
        break;
      default:
        printf("\nError in updateSystem(): rIndex (%u) out of range!\n", rIndex);
        exit(-1);
      }
    }
  }
  // update output
  while (outputTime <= endTime) {
    retval(noutput,0) = outputTime;
    retval(noutput,1) = x[0]/f;
    retval(noutput,2) = x[1]/f;
    retval(noutput,3) = calcium[ntimepoint];
    noutput++;
    outputTime += dt;
  }

  return retval;
}
