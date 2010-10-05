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

// Exact solution u(x,y) to the 2D problem is defined as the
// Cartesian product of the 1D solutions.
static double sol_exact(double x, double y, double z, 
                        double& dx, double& dy, double& dz)
{
  dx = duhat_dx(x) * uhat(y) * uhat(z);
  dy = uhat(x) * duhat_dx(y) * uhat(z);
  dz = uhat(x) * uhat(y) * duhat_dx(z);
  return uhat(x) * uhat(y) * uhat(z);
}

// Right-hand side.
double rhs(double x, double y, double z) {
  return -(  dduhat_dxx(x)*uhat(y)*uhat(z) 
           + uhat(x)*dduhat_dxx(y)*uhat(z) 
           + uhat(x)*uhat(y)*dduhat_dxx(z)
	     ) + K*K*uhat(x)*uhat(y)*uhat(z);
}
