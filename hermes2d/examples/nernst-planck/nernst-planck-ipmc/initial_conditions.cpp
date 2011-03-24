class InitialSolutionVoltage : public ExactSolutionScalar
{
public:
  InitialSolutionVoltage(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) {
    return 0.0;
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) {
    return Ord(0);
  }
};

class InitialSolutionConcentration : public ExactSolutionScalar
{
public:
  InitialSolutionConcentration(Mesh* mesh, double C0) : ExactSolutionScalar(mesh), C0(C0) {};

  virtual scalar value (double x, double y) {
    return C0;
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) {
    return Ord(0);
  }

  double C0;
};

class InitialSolutionU1 : public ExactSolutionScalar
{
public:
  InitialSolutionU1(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) {
    return 0.0;
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) {
    return Ord(0);
  }
};

class InitialSolutionU2 : public ExactSolutionScalar
{
public:
  InitialSolutionU2(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) {
    return 0.0;
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) {
    return Ord(0);
  }
};