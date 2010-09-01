#include <hermes2d.h>
#include <shapeset_common.h>

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

// This test makes sure that the coefficient of 
// every kernel function in shape funtion are correct.
// And make sure the parameters will not be changed by 
// any accident. 

int main(int argc, char* argv[])
{
  // Maximum index of Kernel function tested.
  int max_test_poly_degree = 10;
  int ok = 1;

  // Maximum allowed error.
  double max_allowed_error = 1e-15;
  double x = -1.0;

  double val_kernel[max_test_poly_degree];
  memset (val_kernel, 0, sizeof(double));
  double val_test[max_test_poly_degree];
  memset (val_test, 0, sizeof(double));

  // For future use it is convenient to docompose 
  // the high-order Lobatto shape function l_{2}(x), 
  // l_{3}(x)... into products of the form
  //
  //          l_{k}(x) = l_{0}(x)*l_{1}(x)*phi_{k-2}(x), (2 <= k)    
  // 
  // Since all Lobatto shape functions l_{k}(x), k <= 2, vanish 
  // at -1.0 and 1.0, the kernel function phi_{k-2}, k = 2, 3, ... 
  // are polynomials of the order k-2. 
  //
  // More information please refer to Pavel's 2003 book on page 55.
  val_test[0] = (-2.0 * sqrt(3.0 / 2.0)); 
  val_test[1] = (-2.0 * sqrt(5.0 / 2.0) * (x));
  val_test[2] = (-1.0 / 2.0 * sqrt(7.0 / 2.0) * (5 * (x) * (x) - 1));
  val_test[3] = (-1.0 / 2.0 * sqrt(9.0 / 2.0) * (7 * (x) * (x) - 3) * (x));
  val_test[4] = (-1.0 / 4.0 * sqrt(11.0 / 2.0) * (21 * (x) * (x) * (x) * (x) - 14 * (x) * (x) + 1));
  val_test[5] = (-1.0 / 4.0 * sqrt(13.0 / 2.0) * ((33 * (x) * (x) - 30) * (x) * (x) + 5) * (x));
  val_test[6] = (-1.0 / 32.0 * sqrt(15.0 / 2.0) * (((429 * (x) * (x) - 495) * (x) * (x) + 135) * (x) * (x) - 5));
  val_test[7] = (-1.0 / 32.0 * sqrt(17.0 / 2.0) * (((715 * (x) * (x) - 1001) * (x) * (x) + 385) * (x) * (x) - 35) * (x));
  val_test[8] = (-1.0 / 64.0 * sqrt(19.0 / 2.0) * ((((2431 * (x) * (x) - 4004) * (x) * (x) + 2002) * (x) * (x) - 308) * (x) * (x) + 7));
  val_test[9] = (-1.0 / 64.0 * sqrt(21.0 / 2.0) * ((((4199 * (x) * (x) - 7956) * (x) * (x) + 4914) * (x) * (x) - 1092) * (x) * (x) + 63) * (x));

  val_kernel[0] = phi0(x); 
  val_kernel[1] = phi1(x); 
  val_kernel[2] = phi2(x); 
  val_kernel[3] = phi3(x); 
  val_kernel[4] = phi4(x); 
  val_kernel[5] = phi5(x); 
  val_kernel[6] = phi6(x); 
  val_kernel[7] = phi7(x); 
  val_kernel[8] = phi8(x); 
  val_kernel[9] = phi9(x); 

  // Loop over polynomial degrees, starting with 0.
  for (int n = 0; n < max_test_poly_degree; n++) {

      if (fabs(val_test[n] - val_kernel[n]) > max_allowed_error) {
        printf("n = %d, val_test = %g\n", n, val_test[n]);
        return ERROR_FAILURE;
      }
  }
      printf("Success!\n");
      return ERROR_SUCCESS;
}

