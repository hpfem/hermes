class InitialSolutionHeatTransfer : public ExactSolution1D
{
public:
  InitialSolutionHeatTransfer(Mesh* mesh) : ExactSolution1D(mesh) {};

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    dx = (y+10)/10.;
    dy = (x+10)/10.;
    return (x+10)*(y+10)/100. + 2;
  };
};