class InitialSolutionEulerDensity : public ExactSolutionScalar
{
public:
  InitialSolutionEulerDensity(Mesh* mesh, double constant) : ExactSolutionScalar(mesh), constant(constant) {};

  virtual scalar value (double x, double y) const {
    return constant; 
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return constant;
  }

  // Value.
  double constant;
};

class InitialSolutionEulerDensityVelX : public ExactSolutionScalar
{
public:
  InitialSolutionEulerDensityVelX(Mesh* mesh, double constant) : ExactSolutionScalar(mesh), constant(constant) {};

  virtual scalar value (double x, double y) const {
    return constant; 
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return constant;
  }

  // Value.
  double constant;
};

class InitialSolutionEulerDensityVelY : public ExactSolutionScalar
{
public:
  InitialSolutionEulerDensityVelY(Mesh* mesh, double constant) : ExactSolutionScalar(mesh), constant(constant) {};

  virtual scalar value (double x, double y) const {
    return constant; 
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return constant;
  }

  // Value.
  double constant;
};

class InitialSolutionEulerDensityEnergy : public ExactSolutionScalar
{
public:
  InitialSolutionEulerDensityEnergy(Mesh* mesh, double constant) : ExactSolutionScalar(mesh), constant(constant) {};

  virtual scalar value (double x, double y) const {
    return constant; 
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return constant;
  }

  // Value.
  double constant;
};

class InitialSolutionConcentration : public ExactSolutionScalar
{
public:
  InitialSolutionConcentration(Mesh* mesh, double constant) : ExactSolutionScalar(mesh), constant(constant) {};

  virtual scalar value (double x, double y) const {
    return constant; 
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return constant;
  }

  // Value.
  double constant;
};