#include "hermes1d.h"

#include "legendre.h"
#include "quad_std.h"

// This test makes sure that all Legendre 
// polynomials (starting with the linear 
// one) integrated from -1 to 1 numerically
// with all quadrature rules of same or 
// higher order, give zero.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int main(int argc, char* argv[])
{
  // maximum poly degree of Legendre polynomials tested
  int max_test_poly_degree = MAX_P;
  
  // precalculating the values of Legendre polynomials
  // and their derivatives at all possible quadrature
  // points in (-1,1)
  precalculate_legendre_1d();

  // maximum allowed error
  double max_allowed_error = 1e-12;

  // loop over polynomial degrees, starting with 1
  double max_actual_error = 0; 
  for (int poly_deg=1; poly_deg < max_test_poly_degree + 1; poly_deg++) {
    // integrating the Legendre polynomial of degree 'poly_deg'
    // from -1 to 1 using Gauss quadratures of orders 1, 2, ...
    // MAX_P
    for (int quad_order=poly_deg; quad_order < max_test_poly_degree + 1; quad_order++) {
      int num_pts = g_quad_1d_std.get_num_points(quad_order);
      double2 *quad_tab = g_quad_1d_std.get_points(quad_order);
      double val = 0;
      for (int i=0; i<num_pts; i++) {
        //double point_i = quad_tab[i][0];
        double weight_i = quad_tab[i][1];
        //val += legendre_val_ref(point_i, poly_deg) * weight_i;
        val += legendre_val_ref_tab[quad_order][i][poly_deg] * weight_i;
      }
      printf("poly_deg = %d, quad_order = %d, integral = %g\n",
             poly_deg, quad_order, val);      
      if (max_actual_error > max_allowed_error) {
        printf("Failure!\n");
        return ERROR_FAILURE;
      }
      if (fabs(val) > max_actual_error) max_actual_error = fabs(val);
    }
  }

  printf("Success!\n");
  return ERROR_SUCCESS;
}
