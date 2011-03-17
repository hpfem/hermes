class ExactSolutionPoisson : public ExactSolutionScalar
{
public:
  ExactSolutionPoisson(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  // Exact solution.
  double fn(double x, double y) {
    return x*x +y*y;
  }

  // Exact solution with derivatives.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    dx = 2*x;
    dy = 2*y;
    return x*x + y*y;
  };
};
