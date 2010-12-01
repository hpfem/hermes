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
    val += (du_prevdx[si][comp][i]*dvdx[i] - f(x[i])*v[i])*weights[i];
  }
  return val;
};

double jacobian_surf_left(double x, double u, double dudx,
        double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{
  return (1/Val_newton_alpha_left)*u*v;
}

double jacobian_surf_right(double x, double u, double dudx,
        double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{
  return (1/Val_newton_alpha_right)*u*v;
}

double residual_surf_left(double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
  int comp = 0;    // solution component
  int si = 0;      // solution index (only 0 is relevant for this example)
  return ((u_prev[si][comp] - Val_newton_beta_left)/Val_newton_alpha_left) * v; 
}

double residual_surf_right(double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
  int comp = 0;    // solution component
  int si = 0;      // solution index (only 0 is relevant for this example)
  return ((u_prev[si][comp] - Val_newton_beta_right)/Val_newton_alpha_right) * v; 
}
