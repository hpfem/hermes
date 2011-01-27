#include <hermes2d.h>

// Class that represents the semi-analytic solution on the rectangle [0,1]x[0,1].
// Contains coordinates of the Gauss-Kronrod quadrature nodes and the corresponding
// function values and quadrature weights, as read from the file supplied to the
// constructor.
class SemiAnalyticSolution
{  
  std::vector<long double> x; //  x coordinate.
  std::vector<long double> y; //  y coordinate.
  std::vector<long double> u; //  u(x,y).
  std::vector<long double> w; //  Gauss-Kronrod quadrature weights.
  
  unsigned long int n;  // Number of quadrature points (size of vectors x,y,u,w).
  
  bool NA;  // Indicates that the exact solution could not be read from the input file.
  
  public:
    SemiAnalyticSolution(std::string file);
    
    // The following two functions use the quadrature specified by the input file to
    // calculate integrals over the rectangle [0,1]x[0,1].
    //
    // Calculates L2-norm of the exact solution.
    double get_l2_norm();     
    //
    // Calculates L2-norm of the relative difference between
    // the exact solution and the one computed by Hermes.  
    double get_l2_rel_err(Solution *approx_sln);  
};
    