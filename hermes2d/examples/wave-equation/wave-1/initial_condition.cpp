class InitialConditionWave : public ExactSolutionScalar
{
public:
  InitialConditionWave(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) {
    return exp(-x*x - y*y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) {
  dx = exp(-x*x - y*y) * (-2*x);
  dy = exp(-x*x - y*y) * (-2*y);
  };

  virtual Ord ord(Ord x, Ord y) {
    return Ord(20);
  }
};