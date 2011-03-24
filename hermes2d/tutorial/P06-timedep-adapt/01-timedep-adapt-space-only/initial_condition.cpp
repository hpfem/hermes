class InitialSolutionHeatTransfer : public ExactSolutionScalar
{
public:
  InitialSolutionHeatTransfer(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) {
    return (x+10)*(y+10)/100.;
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) {
    dx = (y+10)/10.;
    dy = (x+10)/10.;
  }

  virtual Ord ord(Ord x, Ord y) {
    return (x+10)*(y+10)/100.;
  }
};