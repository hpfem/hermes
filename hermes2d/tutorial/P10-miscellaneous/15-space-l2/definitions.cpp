// Exact solution.
class CustomExactSolution : public ExactSolutionScalar
{
  public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionScalar(mesh) {};

  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return - pow(x, 4) * pow(y, 5); 
    dx = 0; // not needed for L2-projection
    dy = 0; // not needed for L2-projection
  };
};