class ExactSolutionNIST08 : public ExactSolution1D
{
public:
  ExactSolutionNIST08(Mesh* mesh, double alpha) : ExactSolution1D(mesh), alpha(alpha) {};

  double fn(double x, double y) {
    double r = sqrt(x*x + y*y);
    return sin(1/(alpha + r));
  };

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    double r = sqrt(x*x + y*y);
    double h = 1/(alpha + r);
    dx = -cos(h) * h * h * x / r;
    dy = -cos(h) * h * h * y / r;
    return fn(x, y);
  };

  // Members.
  double alpha;
};