#include <Rcpp.h>
using namespace Rcpp;

typedef void (*amu_ptr)();
// wrap the C++ func pointer in an R external pointer
typedef XPtr<amu_ptr> R_amu_ptr;