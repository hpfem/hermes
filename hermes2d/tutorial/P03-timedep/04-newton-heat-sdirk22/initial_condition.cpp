class InitialSolutionHeatTransfer : public ExactSolutionScalar
{
public:
  InitialSolutionHeatTransfer(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  // Function representing an exact one-dimension valued solution.
  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = (y+10)/10.;
    dy = (x+10)/10.;
  };

  virtual scalar value (double x, double y) const {
    return (x+10)*(y+10)/100.;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return (x+10)*(y+10)/100.;
  }
};