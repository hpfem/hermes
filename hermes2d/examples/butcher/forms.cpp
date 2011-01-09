// Jacobian matrix.
template<typename Real, typename Scalar>
Scalar jac(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
           Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev = u_ext[0];

  // This is a temporary workaround. The stage time t_n + h * c_i
  // can be accessed via u_stage_time->val[0];
  // In this particular case the stage time is not needed as 
  // the form does not depend explicitly on time.
  Func<Scalar>* u_stage_time = ext->fn[0]; 

  // Stationary part of the Jacobian matrix (time derivative term left out).
  for (int i = 0; i < n; i++) {
    result += -wt[i] * (dlam_du(u_prev->val[i]) * u->val[i] *
                       (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
                       + lam(u_prev->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  }

  return result;
}

// Fesidual vector
template<typename Real, typename Scalar>
Scalar res(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
           Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev = u_ext[0];

  // This is a temporary workaround. The stage time t_n + h * c_i
  // can be accessed via u_stage_time->val[0];
  // In this particular case the stage time is not needed as 
  // the form does not depend explicitly on time.
  Func<Scalar>* u_stage_time = ext->fn[0]; 

  // Stationary part of the residual (time derivative term left out).
  for (int i = 0; i < n; i++) {
    result += -wt[i] * (lam(u_prev->val[i]) * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
		       - heat_src(e->x[i], e->y[i]) * v->val[i]);
  }

  return result;
}
