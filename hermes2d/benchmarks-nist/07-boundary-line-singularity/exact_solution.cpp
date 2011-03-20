class ExactSolutionNIST07 : public ExactSolutionScalar
{
public:
  ExactSolutionNIST07(Mesh* mesh, double alpha) : ExactSolutionScalar(mesh), alpha(alpha) {};

  double fn(double x, double y) {
    return pow(x, alpha);
  };

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    dx = (alpha/(pow(x, 0.4)));
    dy = 0;
    return fn(x, y);
  };

  // Members.
  double alpha;
};