// first equation:
template<typename Real, typename Scalar>
Scalar bilinear_form_sym_0_0(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * ( (c_TT/TAU) * u->val[i] * v->val[i] + d_TT * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]) );
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_sym_0_1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * ( d_Tw * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]) );
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf_0_0_ext(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * ( k_TT * u->val[i] * v->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_0(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * ( (c_TT/TAU) * ext->fn[0]->val[i] * v->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_0_ext(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * ( k_TT * TEMP_EXTERIOR * v->val[i] );
  return result;
}

// second equation:
template<typename Real, typename Scalar>
Scalar bilinear_form_sym_1_0(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * ( d_wT * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]) );
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_sym_1_1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * ( (c_ww/TAU) * u->val[i] * v->val[i] + d_ww * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]) );
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf_1_1_ext(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * ( k_ww * u->val[i] * v->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * ( (c_ww/TAU) * ext->fn[0]->val[i] * v->val[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_1_ext(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * e->x[i] * ( k_ww * MOIST_EXTERIOR * v->val[i] );
  return result;
}
