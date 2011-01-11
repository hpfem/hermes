template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                     Geom<Real> *e, ExtData<Scalar> *ext)
{
  return - LAMBDA / HEATCAP / RHO * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                          Geom<Real> *e, ExtData<Scalar> *ext)
{
  return - LAMBDA / HEATCAP / RHO * ALPHA * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                   Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* u_prev = u_ext[0];

  // This is a temporary workaround. The stage time t_n + h * c_i
  // can be accessed via u_stage_time->val[0];
  // In this particular case the stage time is not needed as 
  // the form does not depend explicitly on time.
  Func<Scalar>* u_stage_time = ext->fn[0]; 
  
  Scalar current_time = u_stage_time->val[0];

  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += -wt[i] * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i]);		       
  }

  return result * LAMBDA / HEATCAP / RHO;
}


template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                        Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* u_prev = u_ext[0];

  // This is a temporary workaround. The stage time t_n + h * c_i
  // can be accessed via u_stage_time->val[0];
  // In this particular case the stage time is not needed as 
  // the form does not depend explicitly on time.
  Func<Scalar>* u_stage_time = ext->fn[0]; 
  
  Scalar current_time = u_stage_time->val[0];
  Real current_exterior_temperature = temp_ext<Real>(current_time);

  Scalar result = 0;
  for (int i = 0; i < n; i++) {
    result += -wt[i] * u_prev->val[i] * v->val[i];		       
  }

  return LAMBDA / HEATCAP / RHO * ALPHA * 
         (current_exterior_temperature * int_v<Real, Scalar>(n, wt, v) + result);
}

