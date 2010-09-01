// Jacobian matrix.
template<typename Real, typename Scalar>
Scalar jac(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev_newton = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (u->val[i] * v->val[i] / TAU + dlam_du(u_prev_newton->val[i]) * u->val[i] *
                       (u_prev_newton->dx[i] * v->dx[i] + u_prev_newton->dy[i] * v->dy[i])
                       + lam(u_prev_newton->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
  return result;
}

// Fesidual vector
template<typename Real, typename Scalar>
Scalar res(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
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
