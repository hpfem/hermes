Using Material and Boundary Markers: Typical Way (25-markers-typical)
---------------------------------------------------------------------

**Git reference:** Tutorial example `25-markers-typical <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P10-miscellaneous/25-markers-typical>`_. 

This example shows how to solve a problem with complex area / boundary
properties using the typical approach -- to define an individual weak form 
for each element and boundary marker and register them separately.
The equation solved is

.. math::
         - \nabla \cdot (a(x,y) \nabla u) = f \ \ \ \mbox{in}\ \Omega

where $\Omega$ is a square domain subdivided into four 
identical quadrants. The parameter $a(x,y)$ is constant 
in each of them::

    // Material markers.
    const int SOUTH_EAST = 10;
    const int NORTH_EAST = 20;
    const int NORTH_WEST = 30;
    const int SOUTH_WEST = 40;

    // Corresponding material constants.
    const double A_SE = 1.0;
    const double A_NE = 1.5;
    const double A_NW = 0.5;
    const double A_SW = 2.0;

Boundary conditions are zero Dirichlet on the bottom edge,
Neumann $\partial u / \partial n = -1$ along the vertical edges,
and Neumann $\partial u / \partial n = 1$ along the top edge::

    // Enter boundary markers.
    BCTypes bc_types;
    bc_types.add_bc_dirichlet(BDY_BOTTOM);
    bc_types.add_bc_neumann(Hermes::Tuple<int>(BDY_TOP_NE, BDY_TOP_NW));
    bc_types.add_bc_neumann(Hermes::Tuple<std::string>(BDY_VERTICAL_SE, BDY_VERTICAL_NE, BDY_VERTICAL_NW, BDY_VERTICAL_SW));

    // Enter Dirichlet boundary values.
    BCValues bc_values;
    bc_values.add_zero(BDY_BOTTOM);

The forms are registered as::

    // Initialize the weak formulation.
    WeakForm wf;
    wf.add_matrix_form(callback(bilinear_form_vol_SE), HERMES_NONSYM, SOUTH_EAST);
    wf.add_matrix_form(callback(bilinear_form_vol_NE), HERMES_NONSYM, NORTH_EAST);
    wf.add_matrix_form(callback(bilinear_form_vol_NW), HERMES_NONSYM, NORTH_WEST);
    wf.add_matrix_form(callback(bilinear_form_vol_SW), HERMES_NONSYM, SOUTH_WEST);

    wf.add_vector_form(callback(linear_form_vol));

    wf.add_vector_form_surf(callback(linear_form_surf_VERTICAL_SE), BDY_VERTICAL_SE);
    wf.add_vector_form_surf(callback(linear_form_surf_VERTICAL_NE), BDY_VERTICAL_NE);
    wf.add_vector_form_surf(callback(linear_form_surf_VERTICAL_NW), BDY_VERTICAL_NW);
    wf.add_vector_form_surf(callback(linear_form_surf_VERTICAL_SW), BDY_VERTICAL_SW);
    wf.add_vector_form_surf(callback(linear_form_surf_TOP_NE), BDY_TOP_NE);
    wf.add_vector_form_surf(callback(linear_form_surf_TOP_NW), BDY_TOP_NW);

And the forms are really simple::

    // Bilinear volume forms.
    template<typename Real, typename Scalar>
    Scalar bilinear_form_vol_SE(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    { return A_SE * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v); }
    template<typename Real, typename Scalar>
    Scalar bilinear_form_vol_NE(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    { return A_NE * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v); }
    template<typename Real, typename Scalar>
    Scalar bilinear_form_vol_NW(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    { return A_NW * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v); }
    template<typename Real, typename Scalar>
    Scalar bilinear_form_vol_SW(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    { return A_SW * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v); }

    // Linear volume forms.
    template<typename Real, typename Scalar>
    Scalar linear_form_vol(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    { return RHS * int_v<Real, Scalar>(n, wt, v); }

    // Linear surface forms.
    template<typename Real, typename Scalar>
    Scalar linear_form_surf_VERTICAL_SE(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    { return - A_SE * int_v<Real, Scalar>(n, wt, v); }
    template<typename Real, typename Scalar>
    Scalar linear_form_surf_VERTICAL_NE(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    { return - A_NE * int_v<Real, Scalar>(n, wt, v); }
    template<typename Real, typename Scalar>
    Scalar linear_form_surf_VERTICAL_NW(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    { return - A_NW * int_v<Real, Scalar>(n, wt, v); }
    template<typename Real, typename Scalar>
    Scalar linear_form_surf_VERTICAL_SW(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    { return - A_SW * int_v<Real, Scalar>(n, wt, v); }
    template<typename Real, typename Scalar>
    Scalar linear_form_surf_TOP_NE(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    { return A_NE * int_v<Real, Scalar>(n, wt, v); }
    template<typename Real, typename Scalar>
    Scalar linear_form_surf_TOP_NW(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    { return A_NW * int_v<Real, Scalar>(n, wt, v); }
