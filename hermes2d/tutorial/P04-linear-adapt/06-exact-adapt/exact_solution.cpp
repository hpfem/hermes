class ExactSolutionCustom : public ExactSolution1D
{
public:
  ExactSolutionCustom(Mesh* mesh) : ExactSolution1D(mesh) { };

  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    dx = 0.25 * pow(x*x + y*y, -0.75) * 2 * x;
    dy = 0.25 * pow(x*x + y*y, -0.75) * 2 * y;
    return pow(x*x + y*y, 0.25);
  };
};