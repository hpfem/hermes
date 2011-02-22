Trilinos - Coupled (05-trilinos-coupled)
-----------------------

**Git reference:** Example `05-trilinos-coupled
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P09-trilinos/05-trilinos-coupled>`_.

The purpose of this example is to show how to use Trilinos for nonlinear time-dependent coupled PDE systems.
Solved by NOX solver via Newton or JFNK, with or without preconditioning. 

Model problem
~~~~~~~~~~~~~

We solve the simplified flame
propagation problem from `P03-timedep/flame <http://hpfem.org/hermes/doc/src/hermes2d/timedep/flame.html>`_.

Weak formulation
~~~~~~~~~~~~~~~~

The code is the same as in example 19 until the definition of the weak formulation, where we
use diagonal blocks of the Jacobian for preconditioning::

    // Initialize weak formulation.
    WeakForm wf(2, JFNK ? true : false);
    if (!JFNK || (JFNK && PRECOND == 1))
    {
      wf.add_matrix_form(0, 0, callback(newton_bilinear_form_0_0), HERMES_NONSYM, HERMES_ANY, &omega_dt);
      wf.add_matrix_form_surf(0, 0, callback(newton_bilinear_form_0_0_surf), 3);
      wf.add_matrix_form(1, 1, callback(newton_bilinear_form_1_1), HERMES_NONSYM, HERMES_ANY, &omega_dc);
      wf.add_matrix_form(0, 1, callback(newton_bilinear_form_0_1), HERMES_NONSYM, HERMES_ANY, &omega_dc);
      wf.add_matrix_form(1, 0, callback(newton_bilinear_form_1_0), HERMES_NONSYM, HERMES_ANY, &omega_dt);
    }
    else if (PRECOND == 2)
    {
      wf.add_matrix_form(0, 0, callback(precond_0_0));
      wf.add_matrix_form(1, 1, callback(precond_1_1));
    }
    wf.add_vector_form(0, callback(newton_linear_form_0), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&t_prev_time_1, &t_prev_time_2, &omega));
    wf.add_vector_form_surf(0, callback(newton_linear_form_0_surf), 3);
    wf.add_vector_form(1, callback(newton_linear_form_1), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&c_prev_time_1, &c_prev_time_2, &omega));

Calculating initial condition for NOX
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next we project the initial conditions to obtain an initial coefficient vector for NOX::

  // Project the functions "t_prev_time_1" and "c_prev_time_1" on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solutions on the FE meshes.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(Hermes::Tuple<Space *>(t_space, c_space), 
                                       Hermes::Tuple<MeshFunction*>(&t_prev_time_1, &c_prev_time_1),
                                       coeff_vec);

Initializing DiscreteProblem, NOX, and preconditioner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Then we initialize the DiscreteProblem class, NOX solver, and preconditioner::

    // Initialize finite element problem.
    DiscreteProblem dp(&wf, Hermes::Tuple<Space*>(t_space, c_space));

    // Initialize NOX solver and preconditioner.
    NoxSolver solver(&dp);
    RCP<Precond> pc = rcp(new MlPrecond("sa"));
    if (PRECOND)
    {
      if (JFNK) solver.set_precond(pc);
      else solver.set_precond("Ifpack");
    }
    if (TRILINOS_OUTPUT)
      solver.set_output_flags(NOX::Utils::Error | NOX::Utils::OuterIteration |
                              NOX::Utils::OuterIterationStatusTest |
                              NOX::Utils::LinearSolverDetails);

Setting output flags
~~~~~~~~~~~~~~~~~~~~

Output flags are set as follows::

    if (TRILINOS_OUTPUT)
      solver.set_output_flags(NOX::Utils::Error | NOX::Utils::OuterIteration |
                              NOX::Utils::OuterIterationStatusTest |
                              NOX::Utils::LinearSolverDetails);

Time stepping loop
~~~~~~~~~~~~~~~~~~

The time stepping loop is as usual::

  // Time stepping loop:
  double total_time = 0.0;
  cpu_time.tick_reset();
  for (int ts = 1; total_time <= T_FINAL; ts++)
  {
    info("---- Time step %d, t = %g s", ts, total_time + TAU);

    cpu_time.tick(HERMES_SKIP);
    solver.set_init_sln(coeff_vec);
    if (solver.solve())
    {
      Solution::vector_to_solutions(solver.get_solution(), Hermes::Tuple<Space *>(t_space, c_space), 
                Hermes::Tuple<Solution *>(&t_prev_newton, &c_prev_newton));

      cpu_time.tick();
      info("Number of nonlin iterations: %d (norm of residual: %g)",
          solver.get_num_iters(), solver.get_residual());
      info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)",
          solver.get_num_lin_iters(), solver.get_achieved_tol());

      // Time measurement.
      cpu_time.tick(HERMES_SKIP);

      // Visualization.
      DXDYFilter omega_view(omega_fn, Hermes::Tuple<MeshFunction*>(&t_prev_newton, &c_prev_newton));
      rview.set_min_max_range(0.0,2.0);
      rview.show(&omega_view);
      cpu_time.tick(HERMES_SKIP);
			
      // Skip visualization time.
      cpu_time.tick(HERMES_SKIP);

      // Update global time.
      total_time += TAU;

      // Saving solutions for the next time step.
      t_prev_time_2.copy(&t_prev_time_1);
      c_prev_time_2.copy(&c_prev_time_1);
      t_prev_time_1 = t_prev_newton;
      c_prev_time_1 = c_prev_newton;
    }
    else
      error("NOX failed.");

