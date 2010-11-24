Trilinos - Coupled (44)
-----------------------

**Git reference:** Example `trilinos-coupled
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/44-trilinos-coupled>`_.

The purpose of this example is to show how to use Trilinos for nonlinear time-dependent coupled PDE systems.
Solved by NOX solver via Newton or JFNK, with or without preconditioning. We solve the simplified flame
propagation problem from `tutorial example 19 <http://hpfem.org/hermes/doc/src/hermes2d/tutorial-3.html#flame-propagation-problem-19>`_.

The code is the same as in example 19 until the definition of the weak formulation, where we
use diagonal blocks of the Jacobian for preconditioning::

    // Initialize weak formulation.
    WeakForm wf(2, JFNK ? true : false);
    if (!JFNK || (JFNK && PRECOND == 1))
    {
      wf.add_matrix_form(callback(newton_bilinear_form_0_0), HERMES_UNSYM, H2D_ANY, &omega_dt);
      wf.add_matrix_form_surf(0, 0, callback(newton_bilinear_form_0_0_surf), 3);
      wf.add_matrix_form(1, 1, callback(newton_bilinear_form_1_1), HERMES_UNSYM, H2D_ANY, &omega_dc);
      wf.add_matrix_form(0, 1, callback(newton_bilinear_form_0_1), HERMES_UNSYM, H2D_ANY, &omega_dc);
      wf.add_matrix_form(1, 0, callback(newton_bilinear_form_1_0), HERMES_UNSYM, H2D_ANY, &omega_dt);
    }
    else if (PRECOND == 2)
    {
      wf.add_matrix_form(0, 0, callback(precond_0_0));
      wf.add_matrix_form(1, 1, callback(precond_1_1));
    }
    wf.add_vector_form(0, callback(newton_linear_form_0), H2D_ANY, 
                       Tuple<MeshFunction*>(&t_prev_time_1, &t_prev_time_2, &omega));
    wf.add_vector_form_surf(0, callback(newton_linear_form_0_surf), 3);
    wf.add_vector_form(1, callback(newton_linear_form_1), H2D_ANY, 
                       Tuple<MeshFunction*>(&c_prev_time_1, &c_prev_time_2, &omega));

Next we project the initial conditions to obtain a coefficient vector::

    // Project the functions "t_iter" and "c_iter" on the FE space 
    // in order to obtain initial vector for NOX. 
    info("Projecting initial solutions on the FE meshes.");
    Vector* coeff_vec = new AVector(ndof);
    project_global(Tuple<Space *>(t_space, c_space), Tuple<int>(H2D_H1_NORM, H2D_H1_NORM), 
                   Tuple<MeshFunction*>(&t_prev_time_1, &c_prev_time_1), 
                   Tuple<Solution*>(&t_prev_time_1, &c_prev_time_1),
                   coeff_vec);

Then we initialize the FeProblem class, NOX solver, and preconditioner::

    // Initialize finite element problem.
    FeProblem fep(&wf, Tuple<Space*>(t_space, c_space));

    // Initialize NOX solver and preconditioner.
    NoxSolver solver(&fep);
    RCP<Precond> pc = rcp(new MlPrecond("sa"));
    if (PRECOND)
    {
      if (JFNK) solver.set_precond(pc);
      else solver.set_precond("Ifpack");
    }

Output flags are set as follows::

    if (TRILINOS_OUTPUT)
      solver.set_output_flags(NOX::Utils::Error | NOX::Utils::OuterIteration |
                              NOX::Utils::OuterIterationStatusTest |
                              NOX::Utils::LinearSolverDetails);

The time stepping loop is as usual. Skipping info outputs and visualization,
it has the form::

    for (int ts = 1; total_time <= 60.0; ts++)
    {
      info("---- Time step %d, t = %g s", ts, total_time + TAU);

      cpu_time.tick(HERMES_SKIP);
      solver.set_init_sln(coeff_vec->get_c_array());
      bool solved = solver.solve();
      if (solved)
      {
        double* coeffs = solver.get_solution_vector();
        t_prev_newton.set_coeff_vector(t_space, coeffs, ndof);
        c_prev_newton.set_coeff_vector(c_space, coeffs, ndof);

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


