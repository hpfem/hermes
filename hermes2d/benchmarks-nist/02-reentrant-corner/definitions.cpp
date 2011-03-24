#include "weakform_library/laplace.h"

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, int param) 
        : ExactSolutionScalar(mesh) {
    if (param == 0) {
      OMEGA = ((5.0 * M_PI)/ 4.0);
      ALPHA = (M_PI/ OMEGA);
     }
    else if (param == 1) {
      OMEGA = ((3.0 * M_PI)/ 2.0);
      ALPHA = (M_PI/ OMEGA);
    }
    else if (param == 2) {
      OMEGA = ((7.0 * M_PI)/ 4.0);
      ALPHA = (M_PI/ OMEGA);
    }
    else{
      OMEGA = (2.0 * M_PI);
      ALPHA = (M_PI/ OMEGA);
    }
  };

  double get_angle(double y, double x) const {
    double theta = atan2(y, x);
    if (theta < 0)
      theta += 2 * M_PI;
    return theta;
  };

  virtual double value(double x, double y) const {
    return (pow(sqrt(x*x + y*y), ALPHA) * sin(ALPHA * get_angle(y, x)));
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    double a = sqrt(x*x + y*y);
    double b = pow(a, (ALPHA - 1.0));
    double c = pow(a, ALPHA);
    double d = ((y*y)/(x*x) + 1.0 );

    dx = (((ALPHA* x* sin(ALPHA * get_angle(y,x)) *b)/a) 
         - ((ALPHA *y *cos(ALPHA * get_angle(y, x)) * c)/(pow(x, 2.0) *d)));
    dy = (((ALPHA* cos(ALPHA* get_angle(y, x)) *c)/(x * d)) 
         + ((ALPHA* y* sin(ALPHA* get_angle(y, x)) *b)/a));
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(10);
  }

  double OMEGA;
  double ALPHA;
};

