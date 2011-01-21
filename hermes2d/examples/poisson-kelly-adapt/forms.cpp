template<typename Real, typename Scalar>
Scalar bilinear_form_1(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
               Geom<Real> *e, ExtData<Scalar> *ext)
{
  return EPS_1 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_2(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
               Geom<Real> *e, ExtData<Scalar> *ext)
{
  return EPS_2 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_outer(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return CONST_GAMMA_OUTER * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf_bottom(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return H * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf_bottom(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return T0 * H * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar kelly_interface_estimator(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, 
                                 Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.;
  
  double EPS_C = e->elem_marker == MATERIAL_1 ? EPS_1 : EPS_2;
  double EPS_N = e->get_neighbor_marker() == MATERIAL_1 ? EPS_1 : EPS_2;
  
  for (int i = 0; i < n; i++)
    result += wt[i] * sqr( e->nx[i] * (EPS_C*u->get_dx_central(i) - EPS_N*u->get_dx_neighbor(i)) +
                           e->ny[i] * (EPS_C*u->get_dy_central(i) - EPS_N*u->get_dy_neighbor(i))  );
 
  return result;  // Multiplication by element diameter will be done automatically by the KellyTypeAdapt class.
                  // This allows to call this function only once for each interface and add the correctly scaled
                  // result to the total error estimate for both elements sharing the interface.
}

template<typename Real, typename Scalar>
Scalar kelly_newton_boundary_estimator(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, 
                                       Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.;
  
  double EPS_C = e->elem_marker == MATERIAL_1 ? EPS_1 : EPS_2;
  
  for (int i = 0; i < n; i++)
    result += wt[i] * sqr( H*u->val[i] - H*T0 - e->nx[i] * EPS_C * u->dx[i] - e->ny[i] * EPS_C * u->dy[i] );
  
  return e->diam * result;
}

template<typename Real, typename Scalar>
Scalar kelly_neumann_boundary_estimator(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, 
                                        Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.;
  
  double EPS_C = e->elem_marker == MATERIAL_1 ? EPS_1 : EPS_2;
  
  for (int i = 0; i < n; i++)
    result += wt[i] * sqr( CONST_GAMMA_OUTER - e->nx[i] * EPS_C * u->dx[i] - e->ny[i] * EPS_C * u->dy[i] );
  
  return e->diam * result;
}

template<typename Real, typename Scalar>
Scalar kelly_zero_neumann_boundary_estimator(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, 
                                             Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0.;
  
  double EPS_C = e->elem_marker == MATERIAL_1 ? EPS_1 : EPS_2;
  
  for (int i = 0; i < n; i++)
    result += wt[i] * sqr( e->nx[i] * EPS_C * u->dx[i] + e->ny[i] * EPS_C * u->dy[i] );
  
  return e->diam * result;
}