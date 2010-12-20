// Weak form (volumetric, left-hand side).
template<typename Real, typename Scalar>
Scalar bilinear_form_vol(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  int elem_marker = e->elem_marker;
  double c1, c2, c3, c4;
  if (elem_marker < 0) {
    // This is for Order calculation only:
    c1 = c2 = c3 = c4 = -5555.0;
  } else {
    // FIXME: these global arrays need to be removed.
    int index = _global_mat_markers.find_index(elem_marker);
    c1 = _global_c1_array[index];
    c2 = _global_c2_array[index];
    c3 = _global_c3_array[index];
    c4 = _global_c4_array[index];
  }
  return  
    c1 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)
    + c2 * int_dudx_v<Real, Scalar>(n, wt, u, v)
    + c3 * int_dudy_v<Real, Scalar>(n, wt, u, v)
    + c4 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Weak form (volumetric, right-hand side).
template<typename Real, typename Scalar>
Scalar linear_form_vol(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext)
{
  int elem_marker = e->elem_marker;
  double c5;
  if (elem_marker < 0) {
    // This is for Order calculation only:
    c5 = -5555.0;
  } else {
    // FIXME: these global arrays need to be removed.
    int index = _global_mat_markers.find_index(elem_marker);
    c5 = _global_c5_array[index];
  }
  return c5 * int_v<Real, Scalar>(n, wt, v);
}

// Weak form (surface, left-hand side).
template<typename Real, typename Scalar>
Scalar bilinear_form_surf_newton(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                          Geom<Real> *e, ExtData<Scalar> *ext)
{
  int edge_marker = e->edge_marker;
  int elem_marker = e->elem_marker;
  double const_newton_1;
  double c1;
  if (edge_marker < 0) {
    // This is for Order calculation only:
    const_newton_1 = -5555.0;
    c1 = -5555.0;
  } else {
    // FIXME: these global arrays need to be removed.
    if (_global_bdy_values_newton.size() > 0) {
      double_pair newton_pair = _global_bdy_values_newton[_global_bc_types->find_index_newton(edge_marker)];
      const_newton_1 = newton_pair.first;
    }
    else error("Internal in ModuleBasic: bilinear_form_surf_newton() should not have been called.");
    // FIXME: these global arrays need to be removed.
    c1 = _global_c1_array[_global_mat_markers.find_index(elem_marker)];
  }
  return c1 * const_newton_1 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Weak form (surface, neumann, right-hand side).
template<typename Real, typename Scalar>
Scalar linear_form_surf_neumann(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                        Geom<Real> *e, ExtData<Scalar> *ext)
{
  int edge_marker = e->edge_marker;
  int elem_marker = e->elem_marker;
  double const_neumann;
  double c1;
  double result = 0;
  if (edge_marker < 0) {
    // This is for Order calculation only:
    const_neumann = -5555.0;
    c1 = -5555.0;
  } else {
    // FIXME: these global arrays need to be removed.
    if (_global_bdy_values_neumann.size() > 0) {
      int index = _global_bc_types->find_index_neumann(edge_marker);
      const_neumann = _global_bdy_values_neumann[index];
    }
    else error("Internal in ModuleBasic: linear_form_surf_neumann() should not have been called.");
    // FIXME: these global arrays need to be removed.
    c1 = _global_c1_array[_global_mat_markers.find_index(elem_marker)];  
  }
  return c1 * const_neumann * int_v<Real, Scalar>(n, wt, v);
}

// Weak form (surface, newton, right-hand side).
template<typename Real, typename Scalar>
Scalar linear_form_surf_newton(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                        Geom<Real> *e, ExtData<Scalar> *ext)
{
  int edge_marker = e->edge_marker;
  int elem_marker = e->elem_marker;
  double const_newton_2;
  double c1;
  double result = 0;
  if (edge_marker < 0) {
    // This is for Order calculation only:
    const_newton_2 = -5555.0;
    c1 = -5555.0;
  } else {
    // FIXME: these global arrays need to be removed.
    if (_global_bdy_values_newton.size() > 0) {
      int index = _global_bc_types->find_index_newton(edge_marker);
      double_pair newton_pair = _global_bdy_values_newton[index];
      const_newton_2 = newton_pair.second;
    }
    else error("Internal in ModuleBasic: linear_form_surf_newton() should not have been called.");
    // FIXME: these global arrays need to be removed.
    c1 = _global_c1_array[_global_mat_markers.find_index(elem_marker)];  
  }
  return c1 * const_newton_2 * int_v<Real, Scalar>(n, wt, v);
}
