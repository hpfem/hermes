class ExactSolutionFitzHughNagumo1 : public ExactSolutionScalar
{
public:
  ExactSolutionFitzHughNagumo1(Mesh* mesh, double sigma, double d_u) : ExactSolutionScalar(mesh), sigma(sigma),
                                                                    d_u(d_u) {
  }

  virtual scalar value (double x, double y) const {
    return U(x)*U(y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = dUdt(x)*U(y);
    dy = U(x)*dUdt(y);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }

  double U(double t) {
    return cos(M_PI*t/2);
  }
  double dUdt(double t) {
    return -sin(M_PI*t/2)*(M_PI/2.);
  }
  double ddUdtt(double t) {
    return -cos(M_PI*t/2)*(M_PI/2.)*(M_PI/2.);
  }

  // Members.
  double sigma;
  double d_u;
};


class ExactSolutionFitzHughNagumo2 : public ExactSolutionScalar
{
public:
  ExactSolutionFitzHughNagumo2(Mesh* mesh, double K, double d_v) : ExactSolutionScalar(mesh), K(K), d_v(d_v) {
  }
  virtual scalar value (double x, double y) const {
    return V(x)*V(y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = dVdt(x)*V(y);
    dy = V(x)*dVdt(y);
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }

  double V(double t) {
    return 1. - (exp(K*t) + exp(-K*t))/(exp(K) + exp(-K));
  }
  double dVdt(double t) {
    return -K*(exp(K*t) - exp(-K*t))/(exp(K) + exp(-K));
  }
  double ddVdtt(double t) {
    return -K*K*(exp(K*t) + exp(-K*t))/(exp(K) + exp(-K));
  }
  
  // Members.
  double K;
  double d_v;
};