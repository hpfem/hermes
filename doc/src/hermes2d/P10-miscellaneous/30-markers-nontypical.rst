Using Material and Boundary Markers: Nontypical Way (30-markers-nontypical)
---------------------------------------------------------------------------

**Git reference:** Tutorial example `30-markers-nontypical <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P10-miscellaneous/30-markers-nontypical>`_. 

The reader knows from before that individual weak form can be defined and 
registered for each material marker, as well as for each boundary marker.  
In some situations, however, it might be useful to access the material 
and/or boundary markers from inside the weak forms. We created this 
(a bit artificial) example to illustrate it. 

The equation solved is

.. math::
         -\nabla \cdot (a(x,y) \nabla u) = f \ \ \ \mbox{in}\ \Omega

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
    bc_types.add_bc_neumann(Hermes::Tuple<int>(BDY_VERTICAL, BDY_TOP));

    // Enter Dirichlet boundary values.
    BCValues bc_values;
    bc_values.add_zero(BDY_BOTTOM);

The material and boundary markers are accessed in the weak forms 
as follows. In particular notice that one needs to handle calls
for automatic determination of the quadrature order::

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

    template<typename Real, typename Scalar>
    Scalar linear_form_vol(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {

      return RHS * int_v<Real, Scalar>(n, wt, v);
    }

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

      // For automatic quadrature order calculation (e->elem_marker should be -9999).
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

Finally the forms are registered as::

    // Initialize the weak formulation.
    WeakForm wf;
    wf.add_matrix_form(callback(bilinear_form_vol));
    wf.add_vector_form(callback(linear_form_vol));
    wf.add_vector_form_surf(callback(linear_form_surf));





