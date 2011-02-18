double stac_jacobian_vol(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext)
{
  Func<double>* K_sln = u_ext[0];             // Stage increment from the Runge-Kutta method.
  Func<double>* sln_prev_time = ext->fn[0];   // Previous time step solution.
                                              // To obtain current solution, add these two together.
  double result = 0;
  for (int i = 0; i < n; i++) {
    double sln_val_i = sln_prev_time->val[i] + K_sln->val[i];
    result += -wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) * lambda(e->x[i], e->y[i], sln_val_i);
  }

  return result / HEATCAP / RHO;
}

Ord stac_jacobian_vol_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                          Geom<Ord> *e, ExtData<Ord> *ext)
{
  // This assumes that lambda is at most a cubic polynomial.
  return u->dx[0] * v->dx[0] * pow(e->x[0], 3);
}

double stac_residual_vol(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext)
{
  Func<double>* K_sln = u_ext[0];             
  Func<double>* sln_prev_time = ext->fn[0];   

  // This is a temporary workaround. The stage time t_n + h * c_i
  // can be accessed via u_stage_time->val[0];
  // In this particular case the stage time is not needed as 
  // the form does not depend explicitly on time.
  //Func<double>* u_stage_time = ext->fn[0]; 
  //double current_time = u_stage_time->val[0];

  double result = 0;
  for (int i = 0; i < n; i++) {
    double sln_val_i = sln_prev_time->val[i] + K_sln->val[i];
    double sln_dx_i = sln_prev_time->dx[i] + K_sln->dx[i];
    double sln_dy_i = sln_prev_time->dy[i] + K_sln->dy[i];
    result += -wt[i] * ( sln_dx_i * v->dx[i] +  sln_dy_i * v->dy[i]) * lambda(e->x[i], e->y[i], sln_val_i);	       
  }

  return result / HEATCAP / RHO;
}

Ord stac_residual_vol_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                          Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

double stac_jacobian_top(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext)
{
  Func<double>* K_sln = u_ext[0];             
  Func<double>* sln_prev_time = ext->fn[0];   

  double result = 0;
  for (int i = 0; i < n; i++) {
    double sln_val_i = sln_prev_time->val[i] + K_sln->val[i];
    result += wt[i] * u->val[i] * v->val[i] * lambda(e->x[i], e->y[i], sln_val_i);
  }

  return - result / HEATCAP / RHO * ALPHA_TOP;
}

Ord stac_jacobian_top_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                          Geom<Ord> *e, ExtData<Ord> *ext)
{
  return u->val[0] * v->val[0] * pow(e->x[0], 3);
}

double stac_residual_top(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext)
{
  Func<double>* K_sln = u_ext[0];
  Func<double>* sln_prev_time = ext->fn[0];

  double result = 0;
  for (int i = 0; i < n; i++) {
    double sln_val_i = sln_prev_time->val[i] + K_sln->val[i];
    result += wt[i] * (TEMP_EXT_TOP - sln_val_i) * v->val[i] * lambda(e->x[i], e->y[i], sln_val_i);		       
  }

  return result / HEATCAP / RHO * ALPHA_TOP;
}

Ord stac_residual_top_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                          Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

double stac_jacobian_bottom(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext)
{
  Func<double>* K_sln = u_ext[0];
  Func<double>* sln_prev_time = ext->fn[0];

  double result = 0;
  for (int i = 0; i < n; i++) {
    double sln_val_i = sln_prev_time->val[i] + K_sln->val[i];
    result += wt[i] * u->val[i] * v->val[i] * lambda(e->x[i], e->y[i], sln_val_i);
  }

  return - result / HEATCAP / RHO * ALPHA_BOTTOM;
}

Ord stac_jacobian_bottom_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                             Geom<Ord> *e, ExtData<Ord> *ext)
{
  return u->val[0] * v->val[0] * pow(e->x[0], 3);
}

double stac_residual_bottom(int n, double *wt, Func<double> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<double> *ext)
{
  Func<double>* K_sln = u_ext[0];
  Func<double>* sln_prev_time = ext->fn[0];

  // This is a temporary workaround. The stage time t_n + h * c_i
  // can be accessed via u_stage_time->val[0];
  Func<double>* u_stage_time = ext->fn[1]; 
  
  double stage_time = u_stage_time->val[0];

  double result = 0;
  for (int i = 0; i < n; i++) {
    double sln_val_i = sln_prev_time->val[i] + K_sln->val[i];
    result += wt[i] * (T_fire(e->x[i], stage_time) - sln_val_i) * v->val[i] * lambda(e->x[i], e->y[i], sln_val_i);		       
  }

  return result / HEATCAP / RHO * ALPHA_BOTTOM;
}

Ord stac_residual_bottom_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                             Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}


