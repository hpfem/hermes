class ExactSolutionCustom : public ExactSolutionScalar
{
public:
  ExactSolutionCustom(Mesh* mesh) : ExactSolutionScalar(mesh) { };

  virtual scalar value (double x, double y) const {
    return pow(x*x + y*y, 0.25);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 0.25 * pow(x*x + y*y, -0.75) * 2 * x;
    dy = 0.25 * pow(x*x + y*y, -0.75) * 2 * y;
  }

  virtual Ord ord(Ord x, Ord y) const {
    return pow(x*x + y*y, 0.25);
  }
};