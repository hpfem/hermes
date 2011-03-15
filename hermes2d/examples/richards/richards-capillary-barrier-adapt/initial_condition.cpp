class InitialSolutionRichards : public ExactSolution1D
{
public:
  InitialSolutionRichards(Mesh* mesh, double value) : ExactSolution1D(mesh), value(value) {};

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    return value;
  };

  // Value.
  double value;
};
