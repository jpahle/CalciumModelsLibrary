#include <Rcpp.h>
using namespace Rcpp;

//' Couple a simulated calmodulin protein to a given calcium time series.
//'
//' Takes a calcium time series and simulates a Ca-dependent protein, e.g.
//' Calmodulin.
//'
//' This is the deterministic equivalent of the protein activation process.
//' dE_act/dt = E_inact * (k_act * Ca_cyt^p) / (K_m^p + Ca_cyt^p) - k_inact * E_act
//'
//' @param time A numeric vector: the times of the observations.
//' @param calcium A numeric vector: the concentrations of cytosolic calcium [nmol].
//' @param dt A numeric, the time interval between two output samples.
//' @param vol A numeric, the volume of the system [l].
//' @param k_on A numeric, the on rate of the protein.
//' @param Km A numeric, the half-maximal constant of the Hill activation.
//' @param k_off A numeric, the off rate of the protein.
//' @param E0_conc A numeric, the total concentration of protein [nmol].
//' @param h An integer, the cooperativity of the Hill activation.
//' @return A dataframe with time and the active protein time series as columns.
//' @examples
//' calmodulin()
//' @export
// [[Rcpp::export]]
NumericMatrix calmodulin(NumericVector time,
                         NumericVector calcium,
                         double dt,
                         double vol,
                         double k_on,
                         double Km,
                         double k_off,
                         double E0_conc,
                         int h) {

  /* VARIABLES */
  int ntimepoint;
  NumericVector amu(2);
  IntegerVector x(2);

  /* particle number to concentration (nmol/l) factor */
  /* const double f = 6.0221415e14*CellVol; */
  double f;

  double r2;
  int noutput;
  double currentTime;
  double outputTime;
  double startTime;
  double endTime;
  double tau;
  unsigned int rIndex;

  startTime = time[0];
  endTime = time[time.length()-1];

  f = 6.0221415e14*vol;

  x[0] = (int)floor(E0_conc*f);
  x[1] = 0;

  currentTime = startTime;
  noutput = 0;
  ntimepoint = 0;
  outputTime = currentTime;
  int nintervals = (int)floor((endTime-startTime)/dt+0.5)+1;
  NumericMatrix retval(nintervals, 4);

  while (currentTime < endTime) {
    R_CheckUserInterrupt();
    // calculate propensity amu for every reaction (modified Larsen 2003)
    amu[0] = ((k_on * pow((double)calcium[ntimepoint],(double)h)) / (pow((double)Km,(double)h) + pow((double)calcium[ntimepoint],(double)h))) * x[0];
    amu[1] =  k_off * x[1] + amu[0];
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
//  return DataFrame::create(Named("time")=v,
//                           Named("protein")=s);
}
