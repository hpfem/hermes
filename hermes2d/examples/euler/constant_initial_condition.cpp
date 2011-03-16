class InitialSolutionEulerDensity : public ExactSolutionScalar
{
public:
  InitialSolutionEulerDensity(Mesh* mesh, double value) : ExactSolutionScalar(mesh), value(value) {};

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return value;
  };

  // Value.
  double value;
};
class InitialSolutionEulerDensityVelX : public ExactSolutionScalar
{
public:
  InitialSolutionEulerDensityVelX(Mesh* mesh, double value) : ExactSolutionScalar(mesh), value(value) {};

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return value;
  };

  // Value.
  double value;
};
class InitialSolutionEulerDensityVelY : public ExactSolutionScalar
{
public:
  InitialSolutionEulerDensityVelY(Mesh* mesh, double value) : ExactSolutionScalar(mesh), value(value) {};

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return value;
  };

  // Value.
  double value;
};
class InitialSolutionEulerDensityEnergy : public ExactSolutionScalar
{
public:
  InitialSolutionEulerDensityEnergy(Mesh* mesh, double value) : ExactSolutionScalar(mesh), value(value) {};

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return value;
  };

  // Value.
  double value;
};
