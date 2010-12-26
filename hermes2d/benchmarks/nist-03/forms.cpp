// bilinear forms
template<typename Real, typename Scalar>
Scalar bilinear_form_0_0(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (A * u->dx[i] * v->dx[i] + B * u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_0_1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (C * u->dx[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_0(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (B * u->dx[i] * v->dx[i] + A * u->dy[i] * v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (C * u->dx[i] * v->dy[i]);
  return result;
}

