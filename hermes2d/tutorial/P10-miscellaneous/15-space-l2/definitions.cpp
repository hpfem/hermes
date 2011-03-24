// Exact solution.
class CustomExactSolution : public ExactSolutionScalar
{
  public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar value (double x, double y) const {
    return - pow(x, 4) * pow(y, 5); 
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 0; // Not needed for L2 projections.
    dy = 0; // Not needed for L2 projections.
  };

  virtual Ord ord(Ord x, Ord y) const {
    return - pow(x, 4) * pow(y, 5);
  }
};