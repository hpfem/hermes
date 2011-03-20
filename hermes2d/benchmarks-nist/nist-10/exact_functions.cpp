class ExactFunctionNIST10
{
public:
  ExactFunctionNIST10(double k, double alpha) : k(k), alpha(alpha) {};

  double fn(double x, double y) {
    if (x <= 0) 
      return cos(k * y);
    else 
      return cos(k * y) + pow(x, alpha);
  }

  // Members.
  double k, alpha;
};

// Right-hand side.
class RightHandSideNIST10 : public ExactFunctionNIST10
{
public:
  RightHandSideNIST10(double k, double alpha) : ExactFunctionNIST10(k, alpha) {};

  double value(double x, double y) {
    if (x < 0)
      return fn(x, y) * k * k;
    else 
      return fn(x, y) * k * k - alpha *(alpha - 1) * pow(x, alpha - 2.) - k * k * pow(x, alpha);
  }
};

// Exact solution.
class ExactSolutionNIST10 : public ExactSolutionScalar, public ExactFunctionNIST10
{
public:
  ExactSolutionNIST10(Mesh* mesh, double k, double alpha) : ExactSolutionScalar(mesh), ExactFunctionNIST10(k, alpha) {};

  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    if (x <= 0) 
      dx = 0;
    else 
      dx = alpha * pow(x, alpha - 1);
    dy = -sin(k * y) * k;

    return fn(x, y);
  };
};