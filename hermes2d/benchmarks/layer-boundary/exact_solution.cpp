class ExactSolutionPerturbedPoisson : public ExactSolution1D
{
public:
  ExactSolutionPerturbedPoisson(Mesh* mesh, double K) : ExactSolution1D(mesh), K(K) {};

  // Exact solution to the 1D problem -u'' + K*K*u = K*K in (-1,1) with zero Dirichlet BC.
  double uhat(double x) {
    return 1. - (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
  }
  double duhat_dx(double x) {
    return -K * (exp(K*x) - exp(-K*x)) / (exp(K) + exp(-K));
  }
  double dduhat_dxx(double x) {
    return -K*K * (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
  }

  // Right-hand side.
  double rhs(double x, double y) {
    return -(dduhat_dxx(x)*uhat(y) + uhat(x)*dduhat_dxx(y)) + K*K*uhat(x)*uhat(y);
  }

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    dx = duhat_dx(x) * uhat(y);
    dy = uhat(x) * duhat_dx(y);
    return uhat(x) * uhat(y);
  };

  // Member.
  double K;
};