// Integrals
template<typename Real, typename Scalar>
Scalar int_a_dx_b_dy_c_dz(double a, double b, double c, int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e) 
{
  Scalar Integral = 0.0;
  for (int i = 0; i < n; i++)
    Integral += wt[i] * (a * u->dx[i] * v->dx[i] + b * u->dy[i] * v->dy[i] + c * u->dz[i] * v->dz[i]);
  return Integral;
}

template<typename Real, typename Scalar>
Scalar int_a_dudx_dvdy_b_dudy_dvdx(double a, double b, int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e) 
{
  Scalar Integral = 0.0;
  for (int i = 0; i < n; i++)
    Integral += wt[i] * (a * u->dx[i] * v->dy[i] + b * u->dy[i] * v->dx[i]);
  return Integral;
}

template<typename Real, typename Scalar>
Scalar int_a_dudx_dvdz_b_dudz_dvdx(double a, double b, int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e) 
{
  Scalar Integral = 0.0;
  for (int i = 0; i < n; i++)
    Integral += wt[i] * (a * u->dx[i] * v->dz[i] + b * u->dz[i] * v->dx[i]);
  return Integral;
}

template<typename Real, typename Scalar>
Scalar int_a_dudy_dvdz_b_dudz_dvdy(double a, double b, int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e) 
{
  Scalar Integral = 0.0;
  for (int i = 0; i < n; i++)
    Integral += wt[i] * (a * u->dy[i] * v->dz[i] + b * u->dz[i] * v->dy[i]);
  return Integral;
}

// 1. equation
template<typename Real, typename Scalar>
Scalar bilinear_form_0_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) 
{
  return int_a_dx_b_dy_c_dz<Real, Scalar>(lambda + 2*mu, mu, mu, n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_0_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) 
{
  return int_a_dudx_dvdy_b_dudy_dvdx<Real, Scalar>(lambda, mu, n, wt, v, u, e);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_0_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) 
{
  return int_a_dudx_dvdz_b_dudz_dvdx<Real, Scalar>(lambda, mu, n, wt, v, u, e);
}

template<typename Real, typename Scalar>
Scalar surf_linear_form_x(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) 
{
  Scalar res = 0.0;
  for (int i = 0; i < n; i++)
    res += wt[i] * (f_x * v->val[i]);
  return res;
}

// 2. equation
template<typename Real, typename Scalar>
Scalar bilinear_form_1_1(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) 
{
  return int_a_dx_b_dy_c_dz<Real, Scalar>(mu, lambda + 2*mu, mu, n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_1_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) 
{
  return int_a_dudy_dvdz_b_dudz_dvdy<Real, Scalar>(lambda, mu, n, wt, v, u, e);
}

template<typename Real, typename Scalar>
Scalar surf_linear_form_y(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) 
{
  Scalar res = 0.0;
  for (int i = 0; i < n; i++)
    res += wt[i] * (f_y * v->val[i]);
  return res;
}

// 3. equation
template<typename Real, typename Scalar>
Scalar bilinear_form_2_2(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) 
{
  return int_a_dx_b_dy_c_dz<Real, Scalar>(mu, mu, lambda + 2*mu, n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar surf_linear_form_z(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) 
{
  Scalar res = 0.0;
  for (int i = 0; i < n; i++)
    res += wt[i] * (f_z * v->val[i]);
  return res;
}
