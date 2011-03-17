class ExactSolutionPoisson : public ExactSolutionScalar
{
public:
  ExactSolutionPoisson(Mesh* mesh, double slope) : ExactSolutionScalar(mesh), slope(slope) {};

  // Exact solution.
  double fn(double x, double y) {
    return atan(slope * (sqrt(sqr(x-1.25) + sqr(y+0.25)) - M_PI/3));
  }

  // Exact solution with derivatives.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
      double t = sqrt(sqr(x-1.25) + sqr(y+0.25));
      double u = t * (sqr(slope) * sqr(t - M_PI/3) + 1);
      dx = slope * (x-1.25) / u;
      dy = slope * (y+0.25) / u;
      return fn(x, y);
  };

  // Member.
  double slope;
};
