/// Essential boundary condition for the coupled problem.
class ConcentrationTimedepEssentialBC : public EssentialBoundaryCondition<double> {
public:
  ConcentrationTimedepEssentialBC(std::string marker, double constant, double startup_time) 
           : EssentialBoundaryCondition<double>(Hermes::vector<std::string>()), startup_time(startup_time), constant(constant)
  {
    markers.push_back(marker);
  }

  ~ConcentrationTimedepEssentialBC() {};

  inline EssentialBCValueType get_value_type() const { 
    return EssentialBoundaryCondition<double>::BC_FUNCTION; 
  }

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const
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

/// Class for Heating induced vortex, linear initial condition.
class InitialSolutionLinearProgress : public ExactSolutionScalar<double>
{
public:
  InitialSolutionLinearProgress(MeshSharedPtr mesh, double max, double min, double size) : ExactSolutionScalar<double>(mesh), max(max), min(min), size(size) {};

  virtual double value (double x, double y) const {
    return min + ((max - min) * (size - y) / size); 
  };

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 0;
    dy = - (max - min) / size;
  };

  virtual Ord ord(double x, double y) const {
    return Ord(1);
  }

  MeshFunction<double>* clone() const { if(this->get_type() == HERMES_SLN) return Solution<double>::clone(); else return new InitialSolutionLinearProgress(mesh, max, min, size); }

  // Value.
  double max, min, size;
};