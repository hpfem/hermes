// Jacobian matrix.
template<typename Real, typename Scalar>
Scalar stac_jacobian(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                     Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* u_prev = u_ext[0];

  // This is a temporary workaround. The stage time t_n + h * c_i
  // can be accessed via u_stage_time->val[0];
  // In this particular case the stage time is not needed as 
  // the form does not depend explicitly on time.
  Func<Scalar>* u_stage_time = ext->fn[0]; 

  // Stationary part of the Jacobian matrix (time derivative term left out).
  Scalar result1 = 0, result2 = 0;
  for (int i = 0; i < n; i++) {
    result1 += -wt[i] * dlam_du(u_prev->val[i]) * u->val[i] * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i]);
    result2 += -wt[i] * lam(u_prev->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
  }

  return result1 + result2;
}

// Residual vector
template<typename Real, typename Scalar>
Scalar stac_residual(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                     Geom<Real> *e, ExtData<Scalar> *ext)
{
  Func<Scalar>* u_prev = u_ext[0];

  // This is a temporary workaround. The stage time t_n + h * c_i
  // can be accessed via u_stage_time->val[0];
  // In this particular case the stage time is not needed as 
  // the form does not depend explicitly on time.
  Func<Scalar>* u_stage_time = ext->fn[0]; 

  // Stationary part of the residual (time derivative term left out).
  Scalar result1 = 0, result2 = 0;
  for (int i = 0; i < n; i++) {
    result1 = result1 - wt[i] * lam(u_prev->val[i]) * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i]);
    result2 = result2 + wt[i] * heat_src(e->x[i], e->y[i]) * v->val[i];
  }

  return result1 + result2;
}
