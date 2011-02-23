Time-Dependent Problems (01-cathedral-ie)
----------------------------------

**Git reference:** Tutorial example `01-cathedral-ie <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P03-timedep/01-cathedral-ie>`_. 

Model problem
~~~~~~~~~~~~~

This section describes the implementation of a simple time-dependent
heat transfer model that describes, in a naive approximation, how the St. Vitus 
cathedral in Prague responds to changes in the surrounding air temperature
during one 24-hour cycle. The geometry is shown below:

.. image:: cathedral-ie/vitus1.png
   :align: center
   :width: 400
   :height: 500
   :alt: Model geometry and temperature distribution after 24 hours.

We will solve the standard heat transfer equation

.. math::
    :label: eqvit1

       c \varrho\frac{\partial T}{\partial t} - \lambda \Delta T = 0

equipped with a Dirichlet condition

.. math::

     T = T_{init}

on the bottom edge $\Gamma_{ground}$ and a Newton condition

.. math::

     \frac{\partial T}{\partial \nu} = \alpha(T_{ext}(t) - T)

on the rest of the boundary $\Gamma_{air}$. Here, $c$ is the heat capacity of the material,
$\varrho$ the material density, $\lambda$ the thermal conductivity,
$T_{init}$ the fixed temperature on the
ground (same as the initial temperature of the building), and $\alpha$
the heat transfer coefficient 
between the building and the surrounding air. The surrounding air temperature
$T_{ext}$ is time-dependent of the form

.. math::

     T_{ext}(t) = T_{init} + 10\sin(2\pi t/T_{final}),

where $T_{final}$ is 24 hours (translated into seconds).

Equation :eq:`eqvit1` is equipped with an initial condition of the
form

.. math::

     T(x,y,0) = T_{init}(x,y) \ \ \ \mbox{in} \ \Omega.

For simplicity we will use the implicit Euler method with a constant
time step $\tau$, which transforms equation :eq:`eqvit1` into

.. math::

     c \varrho\frac{T^{n+1} - T^n}{\tau} - \lambda \Delta T^{n+1} = 0.

Weak formulation
~~~~~~~~~~~~~~~~

The corresponding weak formulation is

.. math::

     \int_{\Omega} c \varrho\frac{T^{n+1}}{\tau}v + \int_{\Omega} \lambda \nabla T^{n+1}\cdot \nabla v + \int_{\Gamma_{air}} \alpha \lambda T^{n+1}v = \int_{\Omega} c \varrho\frac{T^{n}}{\tau}v + \int_{\Gamma_{air}} \alpha \lambda T_{ext}(t^{n+1})v.

In this example we use string boundary markers::

    // Boundary markers.
    const std::string BDY_GROUND = "Boundary ground";
    const std::string BDY_AIR = "Boundary air";

Boundary condition types are defined using the BCTypes class::

    // Enter boundary markers.
    BCTypes bc_types;
    bc_types.add_bc_dirichlet(BDY_GROUND);
    bc_types.add_bc_newton(BDY_AIR);

Values for Dirichlet boundary conditions are set via the BCValues class::

    // Enter Dirichlet boundary values.
    BCValues bc_values;
    bc_values.add_const(BDY_GROUND, TEMP_INIT);

Then the space for the temperature $T$ is set up::

    // Initialize an H1 space with default shepeset.
    H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
    int ndof = Space::get_num_dofs(&space);
    info("ndof = %d.", ndof);

Defining weak forms and accessing external functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bilinear and linear forms are defined as follows::

    template<typename Real, typename Scalar>
    Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, u, v) / TAU +
             LAMBDA * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }
  
    template<typename Real, typename Scalar>
    Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Func<Real> *temp_prev = ext->fn[0];
      return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, temp_prev, v) / TAU;
    }
  
    template<typename Real, typename Scalar>
    Scalar bilinear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return LAMBDA * ALPHA * int_u_v<Real, Scalar>(n, wt, u, v);
    }
  
    template<typename Real, typename Scalar>
    Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return LAMBDA * ALPHA * temp_ext(TIME) * int_v<Real, Scalar>(n, wt, v);
    }

Notice how the previous time level temperature is accessed:

::

      Func<Real> *temp_prev = ext->fn[0];
    
Setting initial condition
~~~~~~~~~~~~~~~~~~~~~~~~~ 

Next we need to initialize the previous time level solution tsln with the initial condition $T_{init}$.
Besides holding the finite element solution, the Solution class
can be forced to return zero, to return a constant, or to return an arbitrary function
using the methods set_zero(), set_const() and set_exact(), respectively.
Here we simply call set_const() and supply the initial temperature::

    // Set constant initial condition.
    Solution tsln(&mesh, TEMP_INIT);

Registering external functions in weak forms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The weak forms are registered as follows::

    // Initialize weak formulation.
    WeakForm wf;
    wf.add_matrix_form(callback(bilinear_form));
    wf.add_matrix_form_surf(callback(bilinear_form_surf), BDY_AIR);
    wf.add_vector_form(callback(linear_form), HERMES_ANY, &tsln);
    wf.add_vector_form_surf(callback(linear_form_surf), BDY_AIR);

Notice how the previous time level solution 'tsln' is registered. A few lines above
we saw how it is accessed from inside the weak form. 

Initializing the discrete problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next, the DiscreteProblem class and the matrix solver structures are initialized::

    // Initialize the FE problem.
    bool is_linear = true;
    DiscreteProblem dp(&wf, &space, is_linear);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

Assembling and the 'rhs_only' flag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We are now ready to start the time stepping. Since the stiffness matrix does
not depend on the solution, it only needs to be assembled once in the first time
step. For all remaining time steps it will be the same, and we just need to
re-construct the load vector. This is done via the Boolean variable rhsonly
which is set to false before the time stepping begins. For completeness, we show 
the entire time stepping loop below::

    // Time stepping:
    int ts = 1;
    bool rhs_only = false;
    do 
    {
      info("---- Time step %d, time %3.5f, ext_temp %g", ts, current_time, temp_ext(current_time));

      // First time assemble both the stiffness matrix and right-hand side vector,
      // then just the right-hand side vector.
      if (rhs_only == false) info("Assembling the stiffness matrix and right-hand side vector.");
      else info("Assembling the right-hand side vector (only).");
      dp.assemble(matrix, rhs, rhs_only);
      rhs_only = true;

      // Solve the linear system and if successful, obtain the solution.
      info("Solving the matrix problem.");
      if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &tsln);
      else error ("Matrix solver failed.\n");

      // Visualize the solution.
      char title[100];
      sprintf(title, "Time %3.2f, exterior temperature %3.5f", current_time, temp_ext(current_time));
      Tview.set_title(title);
      Tview.show(&tsln);

      // Update global time.
      current_time += time_step;

      // Increase time step counter
      ts++;
    }
    while (current_time < T_FINAL);


