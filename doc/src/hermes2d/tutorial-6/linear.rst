Trilinos - Linear (40)
----------------------

**Git reference:** Example `trilinos-linear 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/40-trilinos-linear>`_.

The purpose of this example is to show how to use Trilinos for linear PDE problems. 
The problem is solved in two ways for comparison purposes: First, assembling is done 
by the LinearProblem class in Hermes and the discrete problem is solved using UMFpack. 
Second, assembling is done using the FeProblem class in Hermes and the discrete problem 
is solved using the Trilinos NOX solver (using Newton's method or JFNK, with or 
without preconditioning). Note: Assembling in the LinearProblem and FeProblem
is the same, these two classes have only marginal difference and they are going 
to be merged soon. The PDE solved is 

.. math::
    -\Delta u = f

with an exact solution 

.. math::
    u(x,y) = x^2 + y^2.

The first part (LinearProblem + UMFpack) needs not be discussed here. The Trilinos 
part begins with defining a zero initial coefficient vector for the NOX solver::

    // Set initial vector for NOX to zero. Alternatively, you can obtain 
    // an initial vector by projecting init_cond() on the FE space, see below.
    coeff_vec->set_zero();

Alternatively, one can obtain an initial vector by projecting the function 
init_cond() on the FE space::

    // Project the initial condition on the FE space. 
    info("Projecting initial solution on the FE mesh.");
    // The NULL pointer means that we do not want the projection result as a Solution.
    Solution* sln_tmp = new Solution(mesh, init_cond);
    project_global(&space, H2D_H1_NORM, sln_tmp, NULL, coeff_vec);
    delete sln_tmp;

The latter approach is not relevant in this simple example (where moreover init_cond()
is a zero function) but this projection may be practical in more difficult problems 
where one wants to start from a nonzero initial condition. Next, perform the following 
steps:

(1) Initialize the weak formulation::

      // Initialize the weak formulation for Trilinos.
      WeakForm wf2(1, JFNK ? true : false);
      wf2.add_matrix_form(callback(jacobian_form), H2D_SYM);
      wf2.add_vector_form(callback(residual_form));

(2) Initialize the FeProblem class::
 
      // Initialize FeProblem.
      FeProblem* fep = new FeProblem(&wf2, &space);

(3) Initialize NOX solver::

      // Initialize NOX solver.
      NoxSolver* nox_solver = new NoxSolver(fep);

(4) Supply the initial coefficient vector coeff_vec::

      nox_solver->set_init_sln(coeff_vec->get_c_array());

(5) If the user wants preconditioning, set a preconditioner::

      // Choose preconditioning.
      RCP<Precond> pc = rcp(new MlPrecond("sa"));
      if (PRECOND)
      {
        if (JFNK) nox_solver->set_precond(pc);
        else nox_solver->set_precond("ML");
      }

See NOX documentation for more preconditioning choices.

Note: RCP stands for Reference Counted Pointer (a sophisticated smart pointer
management system). This is a very powerful tool provided by the Teuchos library, 
whose usage greatly reduces memory corruption related segfaults.

Now we are ready to call the NOX solver to assemble the discrete problem and solve it::

    // Assemble and solve using NOX.
    bool solved = nox_solver->solve();

The solution vector is extracted from NOX and turned into a Solution as follows::

    double *coeffs = nox_solver->get_solution_vector();
    sln_nox.set_coeff_vector(&space, coeffs, ndof);

That's it! 
