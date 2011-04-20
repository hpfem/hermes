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

class InitialSolutionEulerDensityVelY_LShape : public ExactSolutionScalar
{
public:
  InitialSolutionEulerDensityVelY_LShape(Mesh* mesh, double constant) : ExactSolutionScalar(mesh), constant(constant) {};

  virtual scalar value (double x, double y) const {


    if(y < 0)
      return constant;
    else
      if(x < -0.1)
        return 0;
      else {
        double x_tan = x + 0.1;
        double y_tan = y;

        double angle;
        if(x_tan < 1E-6)
          angle = M_PI / 2;
        else
          angle = std::atan(y_tan / x_tan);
        return constant * (1 - (angle) / (M_PI / 2));
      }
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

class InitialSolutionEulerDensityEnergy_LShape : public ExactSolutionScalar
{
public:
  InitialSolutionEulerDensityEnergy_LShape(Mesh* mesh, InitialSolutionEulerDensity* rho, InitialSolutionEulerDensityVelX* rho_v_x, InitialSolutionEulerDensityVelY_LShape* rho_v_y, double p, double kappa) 
        : ExactSolutionScalar(mesh), rho(rho), rho_v_x(rho_v_x), rho_v_y(rho_v_y), p(p), kappa(kappa) {};

  virtual scalar value (double x, double y) const {
    return QuantityCalculator::calc_energy(rho->value(x, y), rho_v_x->value(x, y), rho_v_y->value(x, y), P_EXT, KAPPA); 
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 0;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(1);
  }

  // Value.
  InitialSolutionEulerDensity* rho;
  InitialSolutionEulerDensityVelX* rho_v_x;
  InitialSolutionEulerDensityVelY_LShape* rho_v_y;
  double p;
  double kappa;
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

/* Essential boundary condition for the coupled problem. */

class ConcentrationTimedepEssentialBC : public EssentialBoundaryCondition {
public:
  ConcentrationTimedepEssentialBC(std::string marker, double constant, double startup_time) 
           : EssentialBoundaryCondition(Hermes::vector<std::string>()), startup_time(startup_time), constant(constant)
  {
    markers.push_back(marker);
  }

  ~ConcentrationTimedepEssentialBC() {};

  inline EssentialBCValueType get_value_type() const { 
    return EssentialBoundaryCondition::BC_FUNCTION; 
  }

  virtual scalar value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
  {
    if(this->get_current_time() < startup_time)
      return 0.0;
    else
      if(this->get_current_time() < 2 * startup_time)
        return ((this->get_current_time() - startup_time) / startup_time) * constant;
      else
        return constant;
  }

  double startup_time;
  double constant;
};