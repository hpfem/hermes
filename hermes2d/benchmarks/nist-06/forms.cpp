template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar val = 0;
  for (int i=0; i < n; i++) {
    val = val + wt[i] * EPSILON * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
    val = val + wt[i] * (2*u->dx[i] + u->dy[i]) * v->val[i];
  }

  return val;
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}

// integration order for the linear_form
Ord linear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

