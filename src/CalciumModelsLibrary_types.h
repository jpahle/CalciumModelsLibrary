#include <Rcpp.h>
using namespace Rcpp;


/* Types for pointers to functions that return the default model parameter set */
typedef void (*init_ptr)();
// wrap the C++ func pointer in an R external pointer
typedef XPtr<init_ptr> R_init_ptr;


/* Types for pointers to functions that calculate the propensities */
typedef void (*amu_ptr)();
// wrap the C++ func pointer in an R external pointer
typedef XPtr<amu_ptr> R_amu_ptr;


/* Types for pointers to functions that update the system state */
typedef void (*stM_ptr)(unsigned int rIndex);
// wrap the C++ func pointer in an R external pointer
typedef XPtr<stM_ptr> R_stM_ptr;