#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This test makes sure that every Lobatto
// shape function starting with the third 
// is zero at x = -1 and zero at x = 1. 

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int main(int argc, char* argv[])
{
  // Maximum index of Lobatto function tested.
  int max_test_poly_degree = MAX_P;
  int ok = 1;

  // Maximum allowed error.
  double max_allowed_error = 1e-10;

  // Loop over polynomial degrees, starting with 2.
  for (int n=2; n<max_test_poly_degree; n++) {
    double val_left = lobatto_val_ref(-1.0, n);
    double val_right = lobatto_val_ref(1.0, n);
      
    if (fabs(val_left) > max_allowed_error) {
      info("n = %d, val_left = %g", n, val_left); 
      ok = 0;
    }
    if (fabs(val_right) > max_allowed_error) {
      info("n = %d, val_right = %g", n, val_right); 
      ok = 0;
    }
  }
  if (ok) {
      info("Success!");
      return ERROR_SUCCESS;
  } else {
      info("Failure!");
      return ERROR_FAILURE;
  }
}
