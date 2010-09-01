#include "hermes1d.h"

// This test makes sure that the Legendre 
// polynomials are compatible with their 
// derivatives (these are defined  
// independently of each other.

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
  double max_allowed_error = 1e-10;

  // loop over polynomial degrees, starting with 1
  double max_actual_error = 0; 
  for (int poly_deg=1; poly_deg < max_test_poly_degree; poly_deg++) {
    // using the formula (2n+1)P_n(x) = P'_{n+1} - P'_{n-1}
    for (int quad_order=poly_deg; quad_order < MAX_QUAD_ORDER; quad_order++) {
      int num_pts = g_quad_1d_std.get_num_points(quad_order);
      double2 *quad_tab = g_quad_1d_std.get_points(quad_order);
      for (int i=0; i < num_pts; i++) {
        //double point_i = quad_tab[i][0];
        //double val_1 = (2.*poly_deg + 1.) 
        //               * legendre_val_ref(point_i, poly_deg);
        //double der_1 = legendre_der_ref(point_i, poly_deg + 1);
        //double der_2 = legendre_der_ref(point_i, poly_deg - 1);
        double val_1 = (2.*poly_deg + 1.) 
                       * legendre_val_ref_tab[quad_order][i][poly_deg];
        double der_1 = legendre_der_ref_tab[quad_order][i][poly_deg+1];
        double der_2 = legendre_der_ref_tab[quad_order][i][poly_deg-1];
        val_1 *= leg_norm_const_ref(poly_deg);
        der_1 *= leg_norm_const_ref(poly_deg + 1);
        der_2 *= leg_norm_const_ref(poly_deg - 1);
        double check_val = fabs(val_1 - (der_1 - der_2));
          printf("poly_deg = %d, i = %d, quad_order = %d, val_1 = %g, der_1 = %g, der_2 = %g, check_val = %g\n", poly_deg, i, quad_order, val_1, der_1, der_2, check_val);      
        if (check_val > max_allowed_error) {
          printf("Failure!\n");
          return ERROR_FAILURE;
        }
      }
    }
  }

  printf("Success!\n");
  return ERROR_SUCCESS;
}
