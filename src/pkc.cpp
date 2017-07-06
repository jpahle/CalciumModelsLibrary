#include <Rcpp.h>
using namespace Rcpp;

//' Couple a simulated PKC protein to a given calcium time series.
//'
//' Takes a calcium time series and simulates the Ca-dependent protein PKC.
//'
//' Implementation based on the PKC Model by Manninnen (2006)
//'
//' @param time A numeric vector: the times of the observations.
//' @param calcium A numeric vector: the concentrations of cytosolic calcium [nmol].
//' @param dt A numeric, the time interval between two output samples.
//' @param vol A numeric, the volume of the system [l].
//' @param k1-20 Double values, reaction rate parameters [1/s].
//' @param AA A double, concentration of AA [nmol].
//' @param DAG A double, concentration of DAG [nmol].
//' @param PKCinact0_conc A double, initial concentration of inactive PKC [nmol].
//' @param PKCbasal0_conc A double, initial concentration of basal PKC [nmol].
//' @return A dataframe with time and the active protein time series as columns.
//' @examples
//' pkc()
//' @export
// [[Rcpp::export]]
NumericMatrix pkc(NumericVector time,
                   NumericVector calcium,
                   double dt,
                   double vol,
                   double k1,
                   double k2,
                   double k3,
                   double k4,
                   double k5,
                   double k6,
                   double k7,
                   double k8,
                   double k9,
                   double k10,
                   double k11,
                   double k12,
                   double k13,
                   double k14,
                   double k15,
                   double k16,
                   double k17,
                   double k18,
                   double k19,
                   double k20,
                   double AA,
                   double DAG,
                   double PKCinact0_conc,
                   double PKCbasal0_conc) {

  /* VARIABLES */
  int ntimepoint;
  NumericVector amu(10);
  IntegerVector x(5);

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

  x[0] = (int)floor(PKCinact0_conc*f);
  x[1] = 0;
  x[2] = 0;
  x[3] = 0;
  x[4] = 0;
  x[5] = (int)floor(PKCbasal0_conc*f);
  x[6] = 0;
  x[7] = 0;
  x[8] = 0;
  x[9] = 0;
  x[10] = 0;

  currentTime = startTime;
  noutput = 0;
  ntimepoint = 0;
  outputTime = currentTime;
  int nintervals = (int)floor((endTime-startTime)/dt+0.5)+1;
  NumericMatrix retval(nintervals, 13);

  while (currentTime < endTime) {
    R_CheckUserInterrupt();
    // calculate propensity amu for every reaction
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

    // calculate time step tau
    tau = - log(runif(1)[0])/amu[9];
    // check if reaction time exceeds time until the next observation
    // and update output
    if ((currentTime+tau)>=time[ntimepoint+1]) {
      currentTime = time[ntimepoint+1];
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
        outputTime += dt;
      }
      ntimepoint++;
    } else {
      // select reaction to fire
      r2 = amu[9] * runif(1)[0];
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
        outputTime += dt;
      }
      // update system state
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
    
    // Debugging
    Rcout << "outputTime: \n" << outputTime << std::endl;
    
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
    outputTime += dt;
  }

  return retval;
//  return DataFrame::create(Named("time")=v,
//                           Named("protein")=s);
}
