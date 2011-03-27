// bilinear form for the Jacobi matrix 
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...basis function
// v...test function
// u_prev...previous solution
double jacobian_vol(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i]*dvdx[i]*weights[i];
  }
  return val;
};

double residual_vol(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  double val = 0;
  int comp = 0;    // solution component
  int si = 0;      // solution index (only 0 is relevant for this example)
  for(int i = 0; i<num; i++) {
    val += (du_prevdx[si][comp][i]*dvdx[i] + f(x[i])*v[i])*weights[i];
  }
  return val;
};

double residual_surf_right(double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
                           double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM],
                           double v, double dvdx, void *user_data)
{
    // FIXME: Later, the value 'Val_neum_right' will enter through user_data,
    // not as a global variable
    // NOTE: the minus sign here is due to the fact that the surface
    // integral -\int_{\partial \Omega} \partial u/\partial nu times v
    // has a negative sign in front of it. But the convention for 
    // defining weak forms in Hermes1D is that all of them are with 
    // positive signs. 
    // NOTE: the Neumann boundary condition deals with the outer
    // normal derivative. At the right end of the interval (a, b), this 
    // is equal to the x-derivative, but at the left end point of (a, b)
    // this is minus one times the x-derivative. In other words, if you 
    // want your solution to be increasing with slope 1 at 'b', the normal
    // derivative at 'b' needs to be 1. If you want the solution to be 
    // increasing with slope 1 at 'a', then the normal derivative at 'a'
    // needs to be -1.   
    return -Val_neum_right * v; 
}
