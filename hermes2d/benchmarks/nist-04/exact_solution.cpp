class ExactSolutionNIST04 : public ExactSolutionScalar
{
public:
  ExactSolutionNIST04(Mesh* mesh, double ALPHA_P, double X_LOC, double Y_LOC) : ExactSolutionScalar(mesh), ALPHA_P(ALPHA_P), X_LOC(X_LOC), Y_LOC(Y_LOC) {};

  double fn(double x, double y) {
    return exp(-ALPHA_P * (pow((x - X_LOC), 2) + pow((y - Y_LOC), 2)));
  };

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    double a = -ALPHA_P * ( (x - X_LOC) * (x - X_LOC) + (y - Y_LOC) * (y - Y_LOC));

    dx = -exp(a) * (2 * ALPHA_P * (x - X_LOC));
    dy = -exp(a) * (2 * ALPHA_P * (y - Y_LOC));

    return fn(x, y);
  };

  // Members.
  double ALPHA_P;
  double X_LOC;
  double Y_LOC;
};