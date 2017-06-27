#include <Rcpp.h>
using namespace Rcpp;

//' Couple a simulated CamKII protein to a given calcium time series.
//'
//' Takes a calcium time series and simulates the Ca-dependent protein CamKII.
//'
//' Implementation based on the CamKII Model by Dupont (2003, 2010)
//'
//' @param time A numeric vector: the times of the observations.
//' @param calcium A numeric vector: the concentrations of cytosolic calcium [nmol].
//' @param dt A numeric, the time interval between two output samples.
//' @param vol A numeric, the volume of the system [l].
//' @param a A numeric, the first parameter of the cubic function calculating k_AA.
//' @param b A numeric, the second parameter of the cubic function calculating k_AA.
//' @param c A numeric, the third parameter of the cubic function calculating k_AA.
//' @param k_IB A numeric, the rate of association of CaM to a non-phosphorylated subunit.
//' @param k_BI A numeric, the rate of dissociation of CaM from a non-phosphorylated subunit.
//' @param k_PT A numeric, the ate of dissociation of Ca from CaM bound to a phosphorylated subunit.
//' @param k_TP A numeric, the rate of association of Ca to CaM bound to a phosphorylated subunit.
//' @param k_TA A numeric, the rate of dissociation of CaM from a phosphorylated subunit.
//' @param k_AT A numeric, the rate of association of CaM to a phosphorylated subunit.
//' @param k_AA A numeric, the phenomenological rate constant for autophosphorylation.
//' @param c_I A numeric, the coefficient of kinase activity of inactive subunits.
//' @param c_B A numeric, the coefficient of kinase activity of bound subunits.
//' @param c_P A numeric, the coefficient of kinase activity of phosphorylated subunits.
//' @param c_T A numeric, the coefficient of kinase activity of trapped subunits.
//' @param c_A A numeric, the coefficient of kinase activity of autonomous subunits.
//' @param camT A numeric, the total concentration of calmodulin [nmol].
//' @param Kd A numeric, the half maximal concentration characterizing Ca binding to CaM.
//' @param Vm_phos A numeric, the maximal velocity of the phosphatase reactions.
//' @param Kd_phos A numeric, the half maximal concentration characterizing phosphorylated subunits being dephosphorylated.
//' @param totalC A numeric, the total fraction of subunits that can phosphorylate or be phosphorylated.
//' @param Wi_conc A numeric, the initial concentration of the inactive fraction [nmol].
//' @param Wb_conc A numeric, the initial concentration of the bound fraction [nmol].
//' @param Wp_conc A numeric, the initial concentration of the phosphorylated fraction [nmol].
//' @param Wt_conc A numeric, the initial concentration of the trapped fraction [nmol].
//' @param Wa_conc A numeric, the initial concentration of the autonomous fraction [nmol].
//' @param h An integer, the Hill coefficient of Ca binding.
//' @return A dataframe with time and the active protein time series as columns.
//' @examples
//' camkii()
//' @export
// [[Rcpp::export]]
NumericMatrix camkii(NumericVector time,
                         NumericVector calcium,
                         double dt,
                         double vol,
                         double a,
                         double b,
                         double c,
                         double k_IB,
                         double k_BI,
                         double k_PT,
                         double k_TP,
                         double k_TA,
                         double k_AT,
                         double k_AA,
                         double c_I,
                         double c_B,
                         double c_P,
                         double c_T,
                         double c_A,
                         double camT,
                         double Kd,
                         double Vm_phos,
                         double Kd_phos,
                         double totalC,
                         double Wi_conc,
                         double Wb_conc,
                         double Wp_conc,
                         double Wt_conc,
                         double Wa_conc,
                         int h) {

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

  x[0] = (int)floor(Wi_conc*f);
  x[1] = (int)floor(Wb_conc*f);
  x[2] = (int)floor(Wp_conc*f);
  x[3] = (int)floor(Wt_conc*f);
  x[4] = (int)floor(Wa_conc*f);

  currentTime = startTime;
  noutput = 0;
  ntimepoint = 0;
  outputTime = currentTime;
  int nintervals = (int)floor((endTime-startTime)/dt+0.5)+1;
  NumericMatrix retval(nintervals, 7);

  while (currentTime < endTime) {
    R_CheckUserInterrupt();
    // calculate propensity amu for every reaction (modified Dupont 2003, 2010)
    amu[0] = x[0] * ((k_IB * camT * pow((double)calcium[ntimepoint],(double)h)) / (pow((double)calcium[ntimepoint],(double)h) + pow((double)Kd,(double)h)));
    amu[1] = amu[0] + k_BI * x[1];
    double activeSubunits = (x[1] + x[2] + x[3] + x[4]) / totalC;
    double prob =  a * activeSubunits + b*(pow((double)activeSubunits,(double)2)) + c*(pow((double)activeSubunits,(double)3));
    amu[2] = amu[1] + totalC * k_AA * prob * ((c_B * x[1]) / pow((double)totalC,(double)2)) * (2*c_B*x[1] + c_P*x[2] + c_T*x[3]+ c_A*x[4]);
    amu[3] = amu[2] + k_PT * x[2];
    amu[4] = amu[3] + k_TP * x[3] * pow((double)calcium[ntimepoint],(double)h);
    amu[5] = amu[4] + k_TA * x[3];
    amu[6] = amu[5] + k_AT * x[4] * (camT - ((camT * pow((double)calcium[ntimepoint],(double)h)) / (pow((double)calcium[ntimepoint],(double)h) + pow((double)Kd,(double)h))));
    amu[7] = amu[6] + ((Vm_phos * x[2]) / (Kd_phos + (x[2] / totalC)));
    amu[8] = amu[7] + ((Vm_phos * x[3]) / (Kd_phos + (x[3] / totalC)));
    amu[9] = amu[8] + ((Vm_phos * x[4]) / (Kd_phos + (x[4] / totalC)));

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
        retval(noutput, 6) = calcium[ntimepoint];
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
        retval(noutput, 6) = calcium[ntimepoint];
		noutput++;
        outputTime += dt;
      }
      // update system state
      switch (rIndex) {
      case 0: // Ca-Cam binding
        x[0]--;
        x[1]++;
        break;
      case 1: // binding reverse
        x[0]++;
        x[1]--;
        break;
      case 2: // phosphorylation
        x[1]--;
        x[2]++;
        break;
      case 3: // trapping
        x[2]--;
        x[3]++;
        break;
      case 4: // trapping reverse
        x[2]++;
        x[3]--;
        break;
      case 5: // autonomous
        x[3]--;
        x[4]++;
        break;
      case 6: // autonomous reverse
        x[3]++;
        x[4]--;
        break;
      case 7: // phosphatase on W_P
        x[2]--;
        x[1]++;
        break;
      case 8: // phosphatase on W_T
        x[3]--;
        x[1]++;
        break;
      case 9: // phosphatase on W_A
        x[4]--;
        x[1]++;
        break;
      default:
        printf("\nError in updateSystem(): rIndex (%u) out of range!\n", rIndex);
        exit(-1);
      }
    }
  }
  // update output
  while (outputTime <= endTime) {
    retval(noutput, 0) = outputTime;
	retval(noutput, 1) = x[0]/f;
	retval(noutput, 2) = x[1]/f;
	retval(noutput, 3) = x[2]/f;
	retval(noutput, 4) = x[3]/f;
	retval(noutput, 5) = x[4]/f;
	retval(noutput, 6) = calcium[ntimepoint];
	noutput++;
    outputTime += dt;
  }

  return retval;
//  return DataFrame::create(Named("time")=v,
//                           Named("protein")=s);
}
