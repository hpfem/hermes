class InitialSolutionVoltage : public ExactSolutionScalar
{
public:
  InitialSolutionVoltage(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  // Function representing an exact scalar valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    // y^2 function for the domain.
    //return (y+100e-6) * (y+100e-6) / (40000e-12);
    return 0.0;
  };
};

class InitialSolutionConcentration : public ExactSolutionScalar
{
public:
  InitialSolutionConcentration(Mesh* mesh, double C0) : ExactSolutionScalar(mesh), C0(C0) {};

  // Function representing an exact scalar valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return C0;
  };

  double C0;
};

class InitialSolutionU1 : public ExactSolutionScalar
{
public:
  InitialSolutionU1(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  // Function representing an exact scalar valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return 0.0;
  };
};

class InitialSolutionU2 : public ExactSolutionScalar
{
public:
  InitialSolutionU2(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  // Function representing an exact scalar valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return 0.0;
  };
};