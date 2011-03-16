class ExactSolutionNIST06 : public ExactSolutionScalar
{
public:
  ExactSolutionNIST06(Mesh* mesh, double epsilon) : ExactSolutionScalar(mesh), epsilon(epsilon) {};

  double fn(double x, double y) {
    return (1 - exp(-(1-x)/epsilon)) * (1 - exp(-(1-y)/epsilon))
      * cos(M_PI * (x + y)); 
  };

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    dx = -M_PI*(1 - exp(-(1 - x)/epsilon))*(1 - exp(-(1 - y)/epsilon))*sin(M_PI*(x + y)) - (1 - exp(-(1 - y)/epsilon))*cos(M_PI*(x + y))*exp(-(1 - x)/epsilon)/epsilon;
    dy = -M_PI*(1 - exp(-(1 - x)/epsilon))*(1 - exp(-(1 - y)/epsilon))*sin(M_PI*(x + y)) - (1 - exp(-(1 - x)/epsilon))*cos(M_PI*(x + y))*exp(-(1 - y)/epsilon)/epsilon;
    return fn(x, y);
  };

  // Right-hand side.
  template<typename Real>
  Real rhs(Real x, Real y) {
    return -epsilon*(-2*pow(M_PI,2)*(1 - exp(-(1 - x)/epsilon))*(1 - exp(-(1 - y)/epsilon))*cos(M_PI*(x + y)) + 2*M_PI*(1 - exp(-(1 - x)/epsilon))*exp(-(1 - y)/epsilon)*sin(M_PI*(x + y))/epsilon + 2*M_PI*(1 - exp(-(1 - y)/epsilon))*exp(-(1 - x)/epsilon)*sin(M_PI*(x + y))/epsilon - (1 - exp(-(1 - y)/epsilon))*cos(M_PI*(x + y))*exp(-(1 - x)/epsilon)/pow(epsilon,2) - (1 - exp(-(1 - x)/epsilon))*cos(M_PI*(x + y))*exp(-(1 - y)/epsilon)/pow(epsilon,2)) - 3*M_PI*(1 - exp(-(1 - x)/epsilon))*(1 - exp(-(1 - y)/epsilon))*sin(M_PI*(x + y)) - 2*(1 - exp(-(1 - y)/epsilon))*cos(M_PI*(x + y))*exp(-(1 - x)/epsilon)/epsilon - (1 - exp(-(1 - x)/epsilon))*cos(M_PI*(x + y))*exp(-(1 - y)/epsilon)/epsilon;
  }

  // Members.
  double epsilon;
};