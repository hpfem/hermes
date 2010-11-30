// bilinear form for the Jacobi matrix 
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...basis function
// v...test function
// u_prev...previous solution
double jacobian_vol_inner(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  int m = 0;
  double val = 0;
  int comp = 0;    // solution component
  for(int i = 0; i<num; i++) {
    val += (D[comp][0] * dudx[i] * dvdx[i] + Sa[comp][m] * u[i] * v[i]) * weights[i]; // inner core
  }
  return val;
}
double jacobian_vol_outer(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  int m = 1;
  double val = 0;
  int comp = 0;    // solution component
  for(int i = 0; i<num; i++) {
    val += (D[comp][m] * dudx[i] * dvdx[i] + Sa[comp][m] * u[i] * v[i]) * weights[i]; // outer core
  }
  return val;
}
double jacobian_vol_reflector(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                void *user_data)
{
  int m = 2;
  double val = 0;
  int comp = 0;    // solution component
  for(int i = 0; i<num; i++) {
    val += (D[comp][m] * dudx[i] * dvdx[i] + Sa[comp][m] * u[i] * v[i]) * weights[i]; // reflector
  }
  return val;
}

double residual_vol_inner(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  int m = 0;
  double val = 0;
  int comp = 0;    // solution component
  int last_newton = 0, last_global = 1;   // solution indices
  for(int i = 0; i<num; i++) {
    val += (D[comp][m] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sa[comp][m] * u_prev[last_newton][comp][i] * v[i] // inner core
    				- chi[comp] / K_EFF * nSf[comp][m] * u_prev[last_global][comp][i] * v[i]) * weights[i];
  }
  return val;
}
double residual_vol_outer(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  int m = 1;
  double val = 0;
  int comp = 0;    // solution component
  int last_newton = 0, last_global = 1;   // solution indices
  for(int i = 0; i<num; i++) {
    val += (D[comp][m] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sa[comp][m] * u_prev[last_newton][comp][i] * v[i] // outer core
    				- chi[comp] / K_EFF * nSf[comp][m] * u_prev[last_global][comp][i] * v[i]) * weights[i];
  }
  return val;
}
double residual_vol_reflector(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double *v, double *dvdx, void *user_data)
{
  int m = 2;
  double val = 0;
  int comp = 0;    // solution component
  int last_newton = 0, last_global = 1;   // solution indices
  for(int i = 0; i<num; i++) {
    val += (D[comp][m] * du_prevdx[last_newton][comp][i] * dvdx[i] + Sa[comp][m] * u_prev[last_newton][comp][i] * v[i] // reflector
    				- chi[comp] / K_EFF * nSf[comp][m] * u_prev[last_global][comp][i] * v[i]) * weights[i];
  }
  return val;
}

double residual_surf_left(double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
  int m = 0;
  int comp = 0;    // solution component
  return D[comp][m] * Val_neumann_left * v; 
}

double jacobian_surf_right(double x, double u, double dudx,
        double v, double dvdx, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], void *user_data)
{
  return Val_albedo_right * u * v;
}

double residual_surf_right(double x, double u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
        double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM], double v,
        double dvdx, void *user_data)
{
  int comp = 0;    // solution component
  int last_newton = 0, last_global = 1;   // solution indices
  return Val_albedo_right * u_prev[last_newton][comp] * v; 
}
