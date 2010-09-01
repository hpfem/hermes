double get_ctrl_alpha(double t) 
{
  // FIXME: this implementation is highly inefficient
  for (int i = 0; i<N_ctrl-1; i++) {
    if (time_ctrl[i] <= t && t <= time_ctrl[i+1]) {
      // return linear interpolant
      return alpha_ctrl[i] + 
	(t - time_ctrl[i])*(alpha_ctrl[i+1] - alpha_ctrl[i])/(time_ctrl[i+1] - time_ctrl[i]);
    }
  }
  error("Internal: time interval not found in get_ctrl_alpha().");
}

double get_ctrl_zeta(double t) 
{
  // FIXME: this implementation is highly inefficient
  for (int i = 0; i<N_ctrl-1; i++) {
    if (time_ctrl[i] <= t && t <= time_ctrl[i+1]) {
      // return linear interpolant
      return zeta_ctrl[i] + 
	(t - time_ctrl[i])*(zeta_ctrl[i+1] - zeta_ctrl[i])/(time_ctrl[i+1] - time_ctrl[i]);
    }
  }
  error("Internal: time interval not found in get_ctrl_zeta().");
}

double jacobian_0_0(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i] * v[i] * weights[i];
  }
  return val;
};

double jacobian_0_2(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* phi = u_prev[si][3];
  double* theta = u_prev[si][4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -u[i] * cos(phi[i]) * cos(theta[i]) * v[i] * weights[i];
  }
  return val;
};

double jacobian_0_3(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* vel = u_prev[si][2];
  double* phi = u_prev[si][3];
  double* theta = u_prev[si][4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += vel[i] * sin(phi[i]) * u[i] * cos(theta[i]) * v[i] * weights[i];
  }
  return val;
};

double jacobian_0_4(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* vel = u_prev[si][2];
  double* phi = u_prev[si][3];
  double* theta = u_prev[si][4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += vel[i] * cos(phi[i]) * sin(theta[i]) * u[i] * v[i] * weights[i];
  }
  return val;
};

double jacobian_1_1(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i] * v[i] * weights[i];
  }
  return val;
};

double jacobian_1_2(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* phi = u_prev[si][3];
  double* theta = u_prev[si][4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -u[i] * cos(phi[i]) * sin(theta[i]) * v[i] * weights[i];
  }
  return val;
};

double jacobian_1_3(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* vel = u_prev[si][2];
  double* phi = u_prev[si][3];
  double* theta = u_prev[si][4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += vel[i] * sin(phi[i]) * u[i] * sin(theta[i]) * v[i] * weights[i];
  }
  return val;
};

double jacobian_1_4(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* vel = u_prev[si][2];
  double* phi = u_prev[si][3];
  double* theta = u_prev[si][4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -vel[i] * cos(phi[i]) * cos(theta[i]) * u[i] * v[i] * weights[i];
  }
  return val;
};

double jacobian_2_2(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i] * v[i] * weights[i];
  }
  return val;
};

double jacobian_3_3(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i] * v[i] * weights[i];
  }
  return val;
};

double jacobian_4_2(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* phi = u_prev[si][3];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -u[i] * sin(phi[i]) * v[i] * weights[i];
  }
  return val;
};

double jacobian_4_3(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* vel = u_prev[si][2];
  double* phi = u_prev[si][3];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += -vel[i] * cos(phi[i]) * u[i] * v[i] * weights[i];
  }
  return val;
};

double jacobian_4_4(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i] * v[i] * weights[i];
  }
  return val;
};

double residual_0(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* dxposdt = du_prevdx[si][0];
  double* vel = u_prev[si][2];
  double* phi = u_prev[si][3];
  double* theta = u_prev[si][4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (dxposdt[i] - vel[i] * cos(phi[i]) * cos(theta[i])) * v[i] * weights[i];
  }
  return val;
};

double residual_1(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* dyposdt = du_prevdx[si][1];
  double* vel = u_prev[si][2];
  double* phi = u_prev[si][3];
  double* theta = u_prev[si][4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (dyposdt[i] - vel[i] * cos(phi[i]) * sin(theta[i])) * v[i] * weights[i];
  }
  return val;
};

double residual_2(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* dveldt = du_prevdx[si][2];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (dveldt[i] - get_ctrl_alpha(x[i])) * v[i] * weights[i];
  }
  return val;
};


double residual_3(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* dphidt = du_prevdx[si][3];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (dphidt[i] - get_ctrl_zeta(x[i])) * v[i] * weights[i];
  }
  return val;
};

double residual_4(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data)
{
  int si = 0;      // solution index (only 0 is relevant for this example)
  // renaming for better readability
  double* vel = u_prev[si][2];
  double* phi = u_prev[si][3];
  double* dthetadt = du_prevdx[si][4];

  double val = 0;
  for(int i = 0; i<num; i++) {
    val += (dthetadt[i] - vel[i] * sin(phi[i])) * v[i] * weights[i];
  }
  return val;
};
