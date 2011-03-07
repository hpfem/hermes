class ExactSolutionPoisson : public ExactSolution1D
{
public:
  ExactSolutionPoisson(Mesh* mesh, double slope) : ExactSolution1D(mesh), slope(slope) {};

  // Right-hand side.
  double fn(double x, double y) {
    return atan(slope * (sqrt(sqr(x-1.25) + sqr(y+0.25)) - M_PI/3));
  }

  // Function representing an exact one-dimension valued solution.
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