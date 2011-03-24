class ExactSolutionCustom : public ExactSolutionScalar
{
public:
  ExactSolutionCustom(Mesh* mesh) : ExactSolutionScalar(mesh) { };

  virtual scalar value (double x, double y) {
    return pow(x*x + y*y, 0.25);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) {
    dx = 0.25 * pow(x*x + y*y, -0.75) * 2 * x;
    dy = 0.25 * pow(x*x + y*y, -0.75) * 2 * y;
  }

  virtual Ord ord(Ord x, Ord y) {
    return pow(x*x + y*y, 0.25);
  }
};