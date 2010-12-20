/* This form is not used now.
template<typename Real, typename Scalar>
Scalar bilinear_form_vol(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  double a;
  switch (e->elem_marker) {
    case SOUTH_EAST: a = A_SE; break;
    case NORTH_EAST: a = A_NE; break;
    case NORTH_WEST: a = A_NW; break;
    case SOUTH_WEST: a = A_SW; break;
    default: if (e->elem_marker >= 0) error("Unknown element marker %d detected.", e->elem_marker);
  }

  // For automatic quadrature order calculation (e->elem_marker should be -9999).
  if (e->elem_marker < 0) a = 1.0; 

  return a * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}
*/

scalar bilinear_form_vol_SE(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
{
  return A_SE * int_grad_u_grad_v<double, scalar>(n, wt, u, v);
}

scalar bilinear_form_vol_NE(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
{
  return A_NE * int_grad_u_grad_v<double, scalar>(n, wt, u, v);
}

scalar bilinear_form_vol_NW(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
{
  return A_NW * int_grad_u_grad_v<double, scalar>(n, wt, u, v);
}

scalar bilinear_form_vol_SW(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
{
  return A_SW * int_grad_u_grad_v<double, scalar>(n, wt, u, v);
}

Ord bilinear_form_vol_Ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return int_grad_u_grad_v<Ord, Ord>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_vol(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return RHS * int_v<Real, Scalar>(n, wt, v);
}

/* This form is not used now.

template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  double normal_der;
  double a;

  switch (e->elem_marker) {
    case SOUTH_EAST: a = A_SE; break;
    case NORTH_EAST: a = A_NE; break;
    case NORTH_WEST: a = A_NW; break;
    case SOUTH_WEST: a = A_SW; break;
    default: if (e->elem_marker >= 0) error("Unknown element marker %d detected.", e->elem_marker);
  }

  // For automatic qua drature order calculation (e->elem_marker should be -9999).
  if (e->elem_marker < 0) a = 1.0; 

  if (e->edge_marker == BDY_VERTICAL) normal_der = -1;
  else {
    if (e->edge_marker == BDY_TOP) normal_der = 1;
    else if (e->edge_marker >= 0) error("Unknown edge marker %d detected.", e->edge_marker);
  }

  // For automatic quadrature order calculation (e->edge_marker should be -8888).
  if (e->edge_marker < 0) normal_der = 1.0; 

  return a * normal_der * int_v<Real, Scalar>(n, wt, v);
}
*/

double linear_form_surf_VERTICAL_SE(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  return - A_SE * int_v<double, double>(n, wt, v);
}
double linear_form_surf_VERTICAL_NE(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  return - A_NE * int_v<double, double>(n, wt, v);
}
double linear_form_surf_VERTICAL_NW(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  return - A_NW * int_v<double, double>(n, wt, v);
}
double linear_form_surf_VERTICAL_SW(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  return - A_SW * int_v<double, double>(n, wt, v);
}

double linear_form_surf_TOP_NE(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  return A_NE * int_v<double, double>(n, wt, v);
}
double linear_form_surf_TOP_NW(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  return A_NW * int_v<double, double>(n, wt, v);
}


Ord linear_form_surf_Ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return A_SW * int_v<Ord, Ord>(n, wt, v);
}
