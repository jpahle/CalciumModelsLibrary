#include <Rcpp.h>
using namespace Rcpp;


// Definitions of simulator variables
// Linker will look here for the definitions of 
// variables that have been externally declared in simulator.cpp 
NumericVector timevector;
double timestep;
double vol;
NumericVector calcium;
unsigned int ntimepoint;
double *amu;
unsigned long long int *x;
int nspecies;
int nreactions;


// Empty placeholder functions
// Since the simulator function 'blueprint' in simulator.cpp is also compiled (Rcpp Issue, it doesn't need to be compiled) we create these placeholders to satisfy the compiler.
// Necessary because excluding simulator.cpp from the compilation process is not possible with the general g++ compiler provided by Rtools.
// These functions are never used since define statements in the model file rename the functions provided by the model file and expected in the included simulator by adding the "_MODEL_NAME" suffix.  
void calculate_amu() {
}
void update_system(unsigned int rIndex) {
}