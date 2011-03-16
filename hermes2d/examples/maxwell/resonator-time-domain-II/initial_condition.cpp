class InitialConditionWave : public ExactSolutionVector
{
public:
  InitialConditionWave(Mesh* mesh) : ExactSolutionVector(mesh) {};

  // Function representing an exact one-dimension valued solution.
  virtual scalar2 exact_function (double x, double y, scalar2& dx, scalar2& dy) {
  dx[0] = sin(x) * sin(y);
  dx[1] = cos(x) * cos(y);
  dy[0] = -cos(x) * cos(y);
  dy[1] = -sin(x) * sin(y);
  return scalar2(-cos(x) * sin(y), sin(x) * cos(y));
  };
};
