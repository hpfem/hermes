class InitialSolutionVoltage : public ExactSolutionScalar
{
public:
  InitialSolutionVoltage(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) const {
    return 0.0;
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(0);
  }
};

class InitialSolutionConcentration : public ExactSolutionScalar
{
public:
  InitialSolutionConcentration(Mesh* mesh, double C0) : ExactSolutionScalar(mesh), C0(C0) {};

  virtual scalar value (double x, double y) const {
    return C0;
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(0);
  }

  double C0;
};