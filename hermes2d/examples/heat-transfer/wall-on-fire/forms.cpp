template<typename Real, typename Scalar>
Scalar stac_jacobian_vol(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += -wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  }

  return result * LAMBDA / HEATCAP / RHO;
}

template<typename Real, typename Scalar>
Scalar stac_residual_vol(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* K_sln = u_ext[0];
  Func<Scalar>* sln_prev_time = ext->fn[0];

  // This is a temporary workaround. The stage time t_n + h * c_i
  // can be accessed via u_stage_time->val[0];
  // In this particular case the stage time is not needed as 
  // the form does not depend explicitly on time.
  //Func<Scalar>* u_stage_time = ext->fn[0]; 
  //Scalar current_time = u_stage_time->val[0];

  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    Scalar sln_dx_i = sln_prev_time->dx[i] + K_sln->dx[i];
    Scalar sln_dy_i = sln_prev_time->dy[i] + K_sln->dy[i];
    result += -wt[i] * ( sln_dx_i * v->dx[i] +  sln_dy_i * v->dy[i]);	       
  }

  return result * LAMBDA / HEATCAP / RHO;
}

template<typename Real, typename Scalar>
Scalar stac_jacobian_top(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  return - LAMBDA / HEATCAP / RHO * ALPHA_TOP * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar stac_residual_top(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* K_sln = u_ext[0];
  Func<Scalar>* sln_prev_time = ext->fn[0];

  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    Scalar sln_val_i = sln_prev_time->val[i] + K_sln->val[i];
    result += wt[i] * (TEMP_EXT_TOP - sln_val_i) * v->val[i];		       
  }

  return LAMBDA / HEATCAP / RHO * ALPHA_TOP * result;
}

template<typename Real, typename Scalar>
Scalar stac_jacobian_bottom(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  return - LAMBDA / HEATCAP / RHO * ALPHA_BOTTOM * int_u_v<Real, Scalar>(n, wt, u, v);
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
    result += wt[i] * (T_fire(e->x[i], stage_time) - sln_val_i) * v->val[i];		       
  }

  return LAMBDA / HEATCAP / RHO * ALPHA_BOTTOM * result;
}

Ord stac_residual_bottom_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                             Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}


