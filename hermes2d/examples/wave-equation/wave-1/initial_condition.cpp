class InitialConditionWave : public ExactSolutionScalar
{
public:
  InitialConditionWave(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
  dx = exp(-x*x - y*y) * (-2*x);
  dy = exp(-x*x - y*y) * (-2*y);
  return exp(-x*x - y*y);
  };
};