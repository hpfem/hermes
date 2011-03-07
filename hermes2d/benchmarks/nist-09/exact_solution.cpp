class ExactSolutionNIST09 : public ExactSolution1D
{
public:
  ExactSolutionNIST09(Mesh* mesh, int user_parameter) : ExactSolution1D(mesh) {
    switch(user_parameter) {
    case 0:
      alpha = 20;
      x_loc = -0.05;
      y_loc = -0.05;
      r_zero = 0.7;
      break;
    case 1:
      alpha = 1000;
      x_loc = -0.05;
      y_loc = -0.05;
      r_zero = 0.7;
      break;
    case 2:
      alpha = 1000;
      x_loc = 1.5;
      y_loc = 0.25;
      r_zero = 0.92;
      break;
    case 3:
      alpha = 50;
      x_loc = 0.5;
      y_loc = 0.5;
      r_zero = 0.25;
      break;
    default:   // The same as 0.
      alpha = 20;
      x_loc = -0.05;
      y_loc = -0.05;
      r_zero = 0.7;
      break;
    }
  };

  template<typename Real>
  Real rhs(Real x, Real y) {  
    Real a = pow(x - x_loc, 2);
    Real b = pow(y - y_loc, 2);
    Real c = sqrt(a + b);
    Real d = ((alpha*x - (alpha * x_loc)) * (2*x - (2 * x_loc)));
    Real e = ((alpha*y - (alpha * y_loc)) * (2*y - (2 * y_loc)));
    Real f = (pow(alpha*c - (alpha * r_zero), 2) + 1.0);
    Real g = (alpha * c - (alpha * r_zero));

    return ((alpha/(c * f)) - (d/(2 * pow(a + b, 1.5) * f)) - ((alpha * d * g)/((a + b) * pow(f, 2))) + 
           (alpha/(c * f)) - (e/(2 * pow(a + b, 1.5) * f)) - ((alpha * e * g)/((a + b) * pow(f, 2))));
  }

  double fn(double x, double y) {
    return atan(alpha * (sqrt(pow(x - x_loc, 2) + pow(y - y_loc, 2)) - r_zero));
  };

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    double a = pow(x - x_loc, 2);
    double b = pow(y - y_loc, 2);
    double c = sqrt(a + b);
    double d = (alpha*x - (alpha * x_loc));
    double e = (alpha*y - (alpha * y_loc));
    double f = (pow(alpha*c - (alpha * r_zero), 2) + 1.0);

    dx = (d/(c * f));
    dy = (e/(c * f));

    return fn(x, y);
  };

  // Members.
  double alpha;
  double x_loc;
  double y_loc;
  double r_zero;
};