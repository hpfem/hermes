// Exact solution u(x,y) = U(x)*U(y) and its derivatives.
double U(double t) {
  return cos(M_PI*t/2);
}
double dUdt(double t) {
  return -sin(M_PI*t/2)*(M_PI/2.);
}
double ddUdtt(double t) {
  return -cos(M_PI*t/2)*(M_PI/2.)*(M_PI/2.);
}
static double uexact(double x, double y, double& dx, double& dy)
{
  dx = dUdt(x)*U(y);
  dy = U(x)*dUdt(y);
  return U(x)*U(y);
}

// Exact solution v(x,y) = V(x)*V(y).
double V(double t) {
  return 1. - (exp(K*t) + exp(-K*t))/(exp(K) + exp(-K));
}
double dVdt(double t) {
  return -K*(exp(K*t) - exp(-K*t))/(exp(K) + exp(-K));
}
double ddVdtt(double t) {
  return -K*K*(exp(K*t) + exp(-K*t))/(exp(K) + exp(-K));
}
static double vexact(double x, double y, double& dx, double& dy)
{
  dx = dVdt(x)*V(y);
  dy = V(x)*dVdt(y);
  return V(x)*V(y);
}

// Right-hand side functions g_1 and g_2.
double g_1(double x, double y)
{
  double Laplace_u = ddUdtt(x)*U(y) + U(x)*ddUdtt(y);
  double u = U(x)*U(y);
  double v = V(x)*V(y);
  return -D_u*D_u * Laplace_u - u + SIGMA*v;
}

double g_2(double x, double y)
{
  double Laplace_v = ddVdtt(x)*V(y) + V(x)*ddVdtt(y);
  double u = U(x)*U(y);
  double v = V(x)*V(y);
  return -D_v*D_v * Laplace_v - u + v;
}
