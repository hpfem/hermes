Trilinos - Nonlinear (41)
-------------------------

**Git reference:** Example `trilinos-nonlinear 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/41-trilinos-nonlinear>`_.

The purpose of this example is to show how to use Trilinos for nonlinear PDE problems. It 
compares performance of the Newton's method in Hermes (assembling via the DiscreteProblem 
class and matrix problem solution via UMFpack) with the performance of the Trilinos/NOX 
solver (using the Hermes FeProblem class to assemble discrete problems).

This example is concerned with the nonlinear equation 

.. math ::
    - \nabla (k(u) \nabla u) = f

where

.. math ::
    k(u) = (1 + u_x^2 + u_y^2)^{-0.5}.


Boundary conditions are chosen zero Dirichlet and a manufactured exact 
solution has the form 

.. math::
    u(x, y) = (x - x^2) (y - y^2).

The Trilinos part starts by projecting the function init_cond() on the finite 
element space to generate an initial coefficient vector for the Newton's method::

    // Project the initial condition on the FE space.
    info("Projecting initial condition on the FE space.");
    // The NULL pointer means that we do not want the projection result as a Solution.
    sln_tmp = new Solution(&mesh, init_cond);
    project_global(&space, H2D_H1_NORM, sln_tmp, NULL, coeff_vec);
    delete sln_tmp;

Note that since init_cond() is zero in this case, we could have just set the initial
coefficient vector to zero as in example 40, but we want to keep the example more general.

Next we initialize the weak formulation (matrix form added only if needed), initialize
the FeProblem class, initialize the NOX solver and supply an initial coefficient vector, 
set preconditioner, and call the NOX solver to assemble and solve the discrete problem::

    // Initialize the weak formulation for Trilinos.
    WeakForm wf2(1, JFNK ? true : false);
    if (!JFNK || (JFNK && PRECOND == 1)) wf2.add_matrix_form(callback(jacobian_form_nox), HERMES_SYM);
    if (JFNK && PRECOND == 2) wf2.add_matrix_form(callback(precond_form_nox), HERMES_SYM);
    wf2.add_vector_form(callback(residual_form_nox));

    // Initialize FeProblem.
    FeProblem fep(&wf2, &space);

    // Initialize the NOX solver with the vector "coeff_vec".
    info("Initializing NOX.");
    NoxSolver nox_solver(&fep);
    nox_solver.set_init_sln(coeff_vec->get_c_array());

    // Choose preconditioning.
    RCP<Precond> pc = rcp(new MlPrecond("sa"));
    if (PRECOND)
    {
      if (JFNK) nox_solver.set_precond(pc);
      else nox_solver.set_precond("ML");
    }

    // Solve the matrix problem using NOX.
    info("Assembling by FeProblem, solving by NOX.");
    bool solved = nox_solver.solve();

The solution coefficient vector is extracted from NOX as in example 40, and 
a Solution is created and visualized as usual.
