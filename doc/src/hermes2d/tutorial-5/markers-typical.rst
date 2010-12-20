Using Material and Boundary Markers in Weak Forms (35)
------------------------------------------------------

**Git reference:** Tutorial example `35-markers-typical <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/35-markers-typical>`_. 

This example shows how to solve a problem with complex area / boundary
properties using the typical approach -- to define an individual weak form 
for each element and boundary marker and register them separately.
The equation solved is::

.. math::
    \mbox{div}(a(x,y) \nabla u) = f \ \ \ \mbox{in}\ \Omega

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
    bc_types.add_bc_neumann(Hermes::Tuple<std::string>(BDY_VERTICAL_SE,BDY_VERTICAL_NE, BDY_VERTICAL_NW, BDY_VERTICAL_SW));

    // Enter Dirichlet boundary values.
    BCValues bc_values;
    bc_values.add_zero(BDY_BOTTOM);

The forms are registered as::

    // Initialize the weak formulation.
    WeakForm wf;
    wf.add_matrix_form(bilinear_form_vol_SE, bilinear_form_vol_Ord, HERMES_UNSYM, SOUTH_EAST);
    wf.add_matrix_form(bilinear_form_vol_NE, bilinear_form_vol_Ord, HERMES_UNSYM, NORTH_EAST);
    wf.add_matrix_form(bilinear_form_vol_NW, bilinear_form_vol_Ord, HERMES_UNSYM, NORTH_WEST);
    wf.add_matrix_form(bilinear_form_vol_SW, bilinear_form_vol_Ord, HERMES_UNSYM, SOUTH_WEST);

    wf.add_vector_form(callback(linear_form_vol));

    wf.add_vector_form_surf(linear_form_surf_VERTICAL_SE, linear_form_surf_Ord, BDY_VERTICAL_SE);
    wf.add_vector_form_surf(linear_form_surf_VERTICAL_NE, linear_form_surf_Ord, BDY_VERTICAL_NE);
    wf.add_vector_form_surf(linear_form_surf_VERTICAL_NW, linear_form_surf_Ord, BDY_VERTICAL_NW);
    wf.add_vector_form_surf(linear_form_surf_VERTICAL_SW, linear_form_surf_Ord, BDY_VERTICAL_SW);
    wf.add_vector_form_surf(linear_form_surf_TOP_NE, linear_form_surf_Ord, BDY_TOP_NE);
    wf.add_vector_form_surf(linear_form_surf_TOP_NW, linear_form_surf_Ord, BDY_TOP_NW);

And the forms are really simple::

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