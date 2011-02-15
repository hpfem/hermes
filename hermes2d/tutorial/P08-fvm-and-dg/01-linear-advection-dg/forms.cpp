template<typename Real>
Real F(Real x, Real y)
{
  return 0;
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += -wt[i] * u->val[i] * calculate_a_dot_v<Real>(e->x[i], e->y[i], v->dx[i], v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real,Scalar>(n, wt, F<Real>, v, e);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_interface(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++) {
    Real a_dot_n = calculate_a_dot_v<Real>(e->x[i], e->y[i], e->nx[i], e->ny[i]);
    Real jump_v = v->get_val_central(i) - v->get_val_neighbor(i);
    result += wt[i] * upwind_flux<Real,Scalar>(u->get_val_central(i), u->get_val_neighbor(i), a_dot_n) * jump_v;
  }
  
  return result;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_boundary(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++) {
    Real x = e->x[i], y = e->y[i];
    Real a_dot_n = calculate_a_dot_v<Real>(x, y, e->nx[i], e->ny[i]);
    result += wt[i] * upwind_flux<Real,Scalar>(u->val[i], 0, a_dot_n) * v->val[i];
  }
  
  return result;
}

template<typename Real, typename Scalar>
Scalar linear_form_boundary(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++) {
    Real x = e->x[i], y = e->y[i];
    Real a_dot_n = calculate_a_dot_v<Real>(x, y, e->nx[i], e->ny[i]);
    result += -wt[i] * upwind_flux<Real,Scalar>(0, g<Real,Scalar>(e->edge_marker,x,y), a_dot_n) * v->val[i];
  }
  
  return result;
}
