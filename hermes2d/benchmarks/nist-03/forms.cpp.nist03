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
Scalar bilinear_form_1_1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * (B * u->dx[i] * v->dx[i] + A * u->dy[i] * v->dy[i]);
  return result;
}



Ord bilinear_form_0_0_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

Ord bilinear_form_0_1_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

Ord bilinear_form_1_1_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// Right hand side.
template<typename Real>
Real rhs_u(Real x, Real y)
{
  return A * dudxdudx(x, y) + B * dudydudy(x, y) + C * dvdxdvdy(x, y);
}

template<typename Real, typename Scalar>
Scalar linear_form_u(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, rhs_u, v, e);
}

template<typename Real>
Real rhs_v(Real x, Real y)
{
  return B * dvdxdvdx(x, y) + A * dvdydvdy(x, y) + C * dudxdudy(x, y);
}

template<typename Real, typename Scalar>
Scalar linear_form_v(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, rhs_v, v, e);
}

// integration order for linear_form_0
Ord linear_form_0_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}

// integration order for linear_form_1
Ord linear_form_1_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(30);
}