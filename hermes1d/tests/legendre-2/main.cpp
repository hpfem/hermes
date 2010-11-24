#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This test makes sure that Legendre polynomials
// are orthonormal. It may take a lot of time.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int main(int argc, char* argv[])
{
  // Maximum poly degree of Legendre polynomials tested.
  int max_test_poly_degree = MAX_P;

  // Precalculating the values of Legendre polynomials
  // and their derivatives at all possible quadrature
  // points in (-1,1).
  precalculate_legendre_1d();

  // Maximum allowed error.
  double max_allowed_error = 1e-12;

  // Loop over polynomial degrees, starting with 1.
  double max_actual_error = 0; 
  for (int poly_deg_1=0; poly_deg_1 < max_test_poly_degree+1; poly_deg_1++) {
    for (int poly_deg_2=0; poly_deg_2 < max_test_poly_degree+1; poly_deg_2++) {
      // Integrating the product of Legendre polynomials of degrees 
      // 'poly_deg_1' and 'poly_deg_2'from -1 to 1 using Gauss quadratures 
      // of order poly_deg_1 + poly_deg_2.
      for (int quad_order = poly_deg_1 + poly_deg_2; 
           quad_order < 2*max_test_poly_degree + 1; quad_order++) {
        int num_pts = g_quad_1d_std.get_num_points(quad_order);
        double2 *quad_tab = g_quad_1d_std.get_points(quad_order);
        double val = 0;
        for (int i=0; i < num_pts; i++) {
          double weight_i = quad_tab[i][1];
          val += legendre_val_ref_tab[quad_order][i][poly_deg_1] * 
                 legendre_val_ref_tab[quad_order][i][poly_deg_2] * weight_i;
        }
        double val_final;
        if (poly_deg_1 == poly_deg_2) val_final = val - 1.0;
        else val_final = val;  
        info("poly_deg_1 = %d, poly_deg_2 = %d, quad_order = %d, val_final = %g", 
                poly_deg_1, poly_deg_2, quad_order, val_final);      
        if (max_actual_error > max_allowed_error) {
          info("Failure!");
          return ERROR_FAILURE;
        }
        if (fabs(val_final) > max_actual_error) {
          max_actual_error = fabs(val_final);
        }
      }
    }
  }

  info("Success!");
  return ERROR_SUCCESS;
}
