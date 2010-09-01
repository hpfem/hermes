#include "hermes1d.h"

#include "legendre.h"
#include "lobatto.h"
#include "quad_std.h"

// This test makes sure that every Lobatto
// shape function starting with the third 
// is zero at x = -1 and zero at x = 1. 

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int main(int argc, char* argv[])
{
  // maximum index of Lobatto function tested
  int max_test_poly_degree = MAX_P;
  int ok = 1;

  // maximum allowed error
  double max_allowed_error = 1e-10;

  // loop over polynomial degrees, starting with 1
  for (int n=2; n<max_test_poly_degree; n++) {
    double val_left = lobatto_val_ref(-1.0, n);
    double val_right = lobatto_val_ref(1.0, n);
      
    if (fabs(val_left) > max_allowed_error) {
      printf("n = %d, val_left = %g\n", n, val_left); 
      ok = 0;
    }
    if (fabs(val_right) > max_allowed_error) {
      printf("n = %d, val_right = %g\n", n, val_right); 
      ok = 0;
    }
  }
  if (ok) {
      printf("Success!\n");
      return ERROR_SUCCESS;
  } else {
      printf("Failure!\n");
      return ERROR_FAILURE;
  }
}
