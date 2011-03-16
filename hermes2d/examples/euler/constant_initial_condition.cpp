class InitialSolutionEulerDensity : public ExactSolution1D
{
public:
  InitialSolutionEulerDensity(Mesh* mesh, double value) : ExactSolution1D(mesh), value(value) {};

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return value;
  };

  // Value.
  double value;
};
class InitialSolutionEulerDensityVelX : public ExactSolution1D
{
public:
  InitialSolutionEulerDensityVelX(Mesh* mesh, double value) : ExactSolution1D(mesh), value(value) {};

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return value;
  };

  // Value.
  double value;
};
class InitialSolutionEulerDensityVelY : public ExactSolution1D
{
public:
  InitialSolutionEulerDensityVelY(Mesh* mesh, double value) : ExactSolution1D(mesh), value(value) {};

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return value;
  };

  // Value.
  double value;
};
class InitialSolutionEulerDensityEnergy : public ExactSolution1D
{
public:
  InitialSolutionEulerDensityEnergy(Mesh* mesh, double value) : ExactSolution1D(mesh), value(value) {};

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return value;
  };

  // Value.
  double value;
};
