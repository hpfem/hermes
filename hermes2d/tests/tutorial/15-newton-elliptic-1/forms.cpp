// Heat sources (can be a general function of 'x' and 'y').
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Jacobian matrix
template<typename Real, typename Scalar>
Scalar jac(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  //Func<Scalar>* u_prev = ext->fn[0];
  Func<Scalar>* u_prev = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (dlam_du(u_prev->val[i]) * u->val[i] * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
                       + lam(u_prev->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));

  return result;
}

// Fesidual vector
template<typename Real, typename Scalar>
Scalar res(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  //Func<Scalar>* u_prev = ext->fn[0];
  Func<Scalar>* u_prev = u_ext[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (lam(u_prev->val[i]) * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
		       - heat_src(e->x[i], e->y[i]) * v->val[i]);
  return result;
}
