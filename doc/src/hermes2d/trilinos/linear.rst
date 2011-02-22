Trilinos - Linear (01-trilinos-linear)
----------------------

**Git reference:** Example `01-trilinos-linear 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P09-trilinos/01-trilinos-linear>`_.

The purpose of this example is to show how to use Trilinos for linear PDE problems. 
The problem is solved in two ways for comparison purposes: First, assembling is done 
by the DiscreteProblem class in Hermes and the discrete problem is solved using UMFpack. 
Second, assembling is done using the DiscreteProblem class in Hermes and the discrete problem 
is solved using the Trilinos NOX solver (using Newton's method or JFNK, with or 
without preconditioning).

Model problem
~~~~~~~~~~~~~

The PDE solved is 

.. math::
    -\Delta u = f

with an exact solution 

.. math::
    u(x,y) = x^2 + y^2.

The first part (DiscreteProblem + UMFpack) needs not be discussed. 

Initializing weak formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    // Initialize the weak formulation for Trilinos.
    WeakForm wf2(1, JFNK ? true : false);
    wf2.add_matrix_form(callback(jacobian_form), HERMES_SYM);
    wf2.add_vector_form(callback(residual_form));

Initializing DiscreteProblem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
::

    // Initialize DiscreteProblem.
    is_linear = false;
    DiscreteProblem dp2(&wf2, &space, is_linear);

Calculating initial condition for NOX
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    // Set initial vector for NOX.
    info("Projecting to obtain initial vector for the Newton's method.");
    scalar* coeff_vec = new scalar[ndof] ;
    Solution* init_sln = new ExactSolution(&mesh, init_cond);
    OGProjection::project_global(&space, init_sln, coeff_vec);
    delete init_sln;

Initializing NOX
~~~~~~~~~~~~~~~~

::

    // Initialize the NOX solver with the vector "coeff_vec".
    info("Initializing NOX.");
    NoxSolver nox_solver(&dp2);
    nox_solver.set_init_sln(coeff_vec);

Setting a preconditioner
~~~~~~~~~~~~~~~~~~~~~~~~

::

    // Choose preconditioning.
    RCP<Precond> pc = rcp(new MlPrecond("sa"));
    if (PRECOND)
    {
      if (JFNK) nox_solver.set_precond(pc);
      else nox_solver.set_precond("ML");
    }

See NOX documentation for more preconditioning choices.

Reference Counted Pointers (RCP)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Note: RCP stands for Reference Counted Pointer (a sophisticated smart pointer
management system). This is a very powerful tool provided by the Teuchos library, 
whose usage greatly reduces memory corruption related segfaults. For more reading 
we refer `here <http://trilinos.sandia.gov/packages/docs/r5.0/packages/teuchos/doc/html/group__RefCountPtr__stuff.html>`_.

Calling NOX
~~~~~~~~~~~

Now we are ready to call the NOX solver to assemble the discrete problem and solve it::

    // Assemble and solve using NOX.
    Solution sln2;
    if (nox_solver.solve())
    {
      Solution::vector_to_solution(nox_solver.get_solution(), &space, &sln2);

      info("Number of nonlin iterations: %d (norm of residual: %g)", 
        nox_solver.get_num_iters(), nox_solver.get_residual());
      info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
        nox_solver.get_num_lin_iters(), nox_solver.get_achieved_tol());
    }
    else error("NOX failed");

Comparing relative errors
~~~~~~~~~~~~~~~~~~~~~~~~~

::

    // Compare relative errors.
    Solution ex;
    ex.set_exact(&mesh, &exact);
    rel_err_1 = calc_abs_error(&sln1, &ex, HERMES_H1_NORM) / calc_norm(&ex, HERMES_H1_NORM) * 100;
    info("Solution 1 (%s):  exact H1 error: %g (time %g s)", MatrixSolverNames[matrix_solver].c_str(), rel_err_1, time1);
    rel_err_2 = calc_abs_error(&sln2, &ex, HERMES_H1_NORM) / calc_norm(&ex, HERMES_H1_NORM) * 100;
    info("Solution 2 (NOX): exact H1 error: %g (time %g + %g = %g [s])", rel_err_2, proj_time, time2, proj_time+time2);

That's it! 

Sample results
~~~~~~~~~~~~~~

You should see the following result:

.. image:: 40/1.png
   :align: center
   :width: 800
   :alt: Sample result
