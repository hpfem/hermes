// Jacobian for the implicit Euler time discretization
template<typename Real, typename Scalar>
Scalar J_euler(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] / TAU + dlam_du(u_prev_newton->val[i]) * u->val[i] *
                       (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                       + lam(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  return result;
}

// Residuum for the implicit Euler time discretization
template<typename Real, typename Scalar>
Scalar F_euler(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = u_ext[0];
  Func<Scalar>* u_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((u_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / TAU +
                      lam(u_prev_newton->val[i]) * (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
		      - heat_src(e->x[i], e->y[i]) * v->val[i]);
  return result;
}

// Jacobian for the Crank-Nicolson time discretization
template<typename Real, typename Scalar>
Scalar J_cranic(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] / TAU +
                       0.5 * dlam_du(u_prev_newton->val[i]) * u->val[i] * (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                       + 0.5 * lam(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  return result;
}

// Residuum for the Crank-Nicolson time discretization
template<typename Real, typename Scalar>
Scalar F_cranic(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = u_ext[0];
  Func<Scalar>* u_prev_time = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * ((u_prev_newton->val[i] - u_prev_time->val[i]) * v->val[i] / TAU
                       + 0.5 * lam(u_prev_newton->val[i]) * (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                       + 0.5 * lam(u_prev_time->val[i]) * (u_prev_time->dx[i] * v->dx[i] + u_prev_time->dy[i] * v->dy[i])
                       - heat_src(e->x[i], e->y[i]) * v->val[i]);
  return result;
}
