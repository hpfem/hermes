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
Scalar stac_jacobian_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                          Geom<Real> *e, ExtData<Scalar> *ext)
{
  return - LAMBDA / HEATCAP / RHO * ALPHA * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar stac_residual_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                          Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* K_sln = u_ext[0];
  Func<Scalar>* sln_prev_time = ext->fn[0];

  // This is a temporary workaround. The stage time t_n + h * c_i
  // can be accessed via u_stage_time->val[0];
  Func<Scalar>* u_stage_time = ext->fn[1]; 
  
  Scalar stage_time = u_stage_time->val[0];
  Real stage_ext_temp = temp_ext<Real>(stage_time);

  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    Scalar sln_val_i = sln_prev_time->val[i] + K_sln->val[i];
    result += wt[i] * (stage_ext_temp - sln_val_i) * v->val[i];		       
  }

  return LAMBDA / HEATCAP / RHO * ALPHA * result;
}
