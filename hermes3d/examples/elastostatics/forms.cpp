// Integrals
template<typename real, typename scalar>
scalar int_a_dx_b_dy_c_dz(double a, double b, double c, int n, double *wt, fn_t<real> *u, fn_t<real> *v, geom_t<real> *e) 
{
  scalar Integral = 0.0;
  for (int i = 0; i < n; i++)
    Integral += wt[i] * (a * u->dx[i] * v->dx[i] + b * u->dy[i] * v->dy[i] + c * u->dz[i] * v->dz[i]);
  return Integral;
}

template<typename real, typename scalar>
scalar int_a_dudx_dvdy_b_dudy_dvdx(double a, double b, int n, double *wt, fn_t<real> *u, fn_t<real> *v, geom_t<real> *e) 
{
  scalar Integral = 0.0;
  for (int i = 0; i < n; i++)
    Integral += wt[i] * (a * u->dx[i] * v->dy[i] + b * u->dy[i] * v->dx[i]);
  return Integral;
}

template<typename real, typename scalar>
scalar int_a_dudx_dvdz_b_dudz_dvdx(double a, double b, int n, double *wt, fn_t<real> *u, fn_t<real> *v, geom_t<real> *e) 
{
  scalar Integral = 0.0;
  for (int i = 0; i < n; i++)
    Integral += wt[i] * (a * u->dx[i] * v->dz[i] + b * u->dz[i] * v->dx[i]);
  return Integral;
}

template<typename real, typename scalar>
scalar int_a_dudy_dvdz_b_dudz_dvdy(double a, double b, int n, double *wt, fn_t<real> *u, fn_t<real> *v, geom_t<real> *e) 
{
  scalar Integral = 0.0;
  for (int i = 0; i < n; i++)
    Integral += wt[i] * (a * u->dy[i] * v->dz[i] + b * u->dz[i] * v->dy[i]);
  return Integral;
}

// 1. equation
template<typename real, typename scalar>
scalar bilinear_form_0_0(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data) 
{
  return int_a_dx_b_dy_c_dz<real, scalar>(lambda + 2*mu, mu, mu, n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar bilinear_form_0_1(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data) 
{
  return int_a_dudx_dvdy_b_dudy_dvdx<real, scalar>(lambda, mu, n, wt, v, u, e);
}

template<typename real, typename scalar>
scalar bilinear_form_0_2(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data) 
{
  return int_a_dudx_dvdz_b_dudz_dvdx<real, scalar>(lambda, mu, n, wt, v, u, e);
}

template<typename real, typename scalar>
scalar surf_linear_form_0(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data) 
{
  return 0.0;
}

// 2. equation
template<typename real, typename scalar>
scalar bilinear_form_1_1(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data) 
{
  return int_a_dx_b_dy_c_dz<real, scalar>(mu, lambda + 2*mu, mu, n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar bilinear_form_1_2(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data) 
{
  return int_a_dudy_dvdz_b_dudz_dvdy<real, scalar>(lambda, mu, n, wt, v, u, e);
}

template<typename real, typename scalar>
scalar surf_linear_form_1(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data) 
{
  return 0.0;
}

// 3. equation
template<typename real, typename scalar>
scalar bilinear_form_2_2(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data) 
{
  return int_a_dx_b_dy_c_dz<real, scalar>(mu, mu, lambda + 2*mu, n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar surf_linear_form_2(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data) 
{
  scalar res = 0.0;
  for (int i = 0; i < n; i++)
    res += wt[i] * (f * v->fn[i]);
  return res;
}
