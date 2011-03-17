// Exact solution to the 1D problem -u'' + K*K*u = K*K in (-1,1) with zero Dirichlet BC.
class MyExactFunction : public ExactSolutionScalar
{
public:
  MyExactFunction(Mesh* mesh, double K) : ExactSolutionScalar(mesh), K(K) {};

  double uhat(double x) {
    return 1. - (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
  }
  double duhat_dx(double x) {
    return -K * (exp(K*x) - exp(-K*x)) / (exp(K) + exp(-K));
  }
  double dduhat_dxx(double x) {
    return -K*K * (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
  }

  // Member.
  double K;
};

// Right-hand side.
class MyRightHandSide : public MyExactFunction
{
public:
  MyRightHandSide(Mesh* mesh, double K) : MyExactFunction(mesh, K) {};

  double rhs(double x, double y) {
    return -(dduhat_dxx(x)*uhat(y) + uhat(x)*dduhat_dxx(y)) + K*K*uhat(x)*uhat(y);
  }
};

// Exact solution.
class MyExactSolution : public MyExactFunction
{
public:
  MyExactSolution(Mesh* mesh, double K) : MyExactFunction(mesh, K) {};

  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    dx = duhat_dx(x) * uhat(y);
    dy = uhat(x) * duhat_dx(y);
    return uhat(x) * uhat(y);
  };
};





