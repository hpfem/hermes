Trilinos - Nonlinear (02-trilinos-nonlinear)
-------------------------

**Git reference:** Example `02-trilinos-nonlinear 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P09-trilinos/02-trilinos-nonlinear>`_.

The purpose of this example is to show how to use Trilinos for nonlinear PDE problems. It 
compares performance of the Newton's method in Hermes (assembling via the DiscreteProblem 
class and matrix problem solution via UMFpack) with the performance of the Trilinos/NOX 
solver (using the Hermes DiscreteProblem class to assemble discrete problems).

Model problem
~~~~~~~~~~~~~

This example is concerned with the nonlinear equation 

.. math ::
    - \nabla (k(u) \nabla u) = f

where

.. math ::
    k(u) = (1 + u_x^2 + u_y^2)^{-0.5}.


Boundary conditions are chosen zero Dirichlet.

Manufactured exact solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We have manufactured an exact solution has the form 

.. math::
    u(x, y) = (x - x^2) (y - y^2).

Again, we do not find it necessary to discuss the Hermes/UMFpack part.

Calculating initial condition for NOX
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Trilinos part starts by projecting the function init_cond() on the finite 
element space to generate an initial coefficient vector for the Newton's method::

    // Project the initial condition on the FE space.
    info("Projecting initial condition on the FE space.");
    sln_tmp = new Solution(&mesh, init_cond);
    OGProjection::project_global(&space, sln_tmp, coeff_vec, matrix_solver);
    delete sln_tmp;

Note that since init_cond() is zero in this case, we could have just set the initial
coefficient vector to zero, but we want to keep the example more general.

Initializating weak forms for Newton or JFNK
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next we initialize the weak formulation (matrix form added only if needed), initialize
the DiscreteProblem class, initialize the NOX solver and supply an initial coefficient vector, 
set preconditioner::

    // Initialize the weak formulation for Trilinos.
    WeakForm wf2(1, JFNK ? true : false);
    if (!JFNK || (JFNK && PRECOND == 1)) wf2.add_matrix_form(callback(jacobian_form_nox), HERMES_SYM);
    if (JFNK && PRECOND == 2) wf2.add_matrix_form(callback(precond_form_nox), HERMES_SYM);
    wf2.add_vector_form(callback(residual_form_nox));

Initializing DiscreteProblem, NOX, and preconditioner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Then we initialize the DiscreteProblem class, NOX solver, and set preconditioner::

    // Initialize DiscreteProblem.
    DiscreteProblem dp2(&wf2, &space);

    // Initialize the NOX solver with the vector "coeff_vec".
    info("Initializing NOX.");
    NoxSolver nox_solver(&dp2);
    nox_solver.set_init_sln(coeff_vec);

    // Choose preconditioning.
    RCP<Precond> pc = rcp(new MlPrecond("sa"));
    if (PRECOND)
    {
      if (JFNK) nox_solver.set_precond(pc);
      else nox_solver.set_precond("ML");
    }

Calling NOX
~~~~~~~~~~~

Next we call the NOX solver to assemble and solve the discrete problem::

    // Solve the nonlinear problem using NOX.
    info("Assembling by DiscreteProblem, solving by NOX.");
    Solution sln2;
    if (nox_solver.solve())
    {
      Solution::vector_to_solution(nox_solver.get_solution(), &space, &sln2);
      info("Number of nonlin iterations: %d (norm of residual: %g)", 
           nox_solver.get_num_iters(), nox_solver.get_residual());
      info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
           nox_solver.get_num_lin_iters(), nox_solver.get_achieved_tol());
    }
    else
      error("NOX failed.");


The solution coefficient vector is extracted from NOX as in example 40, and 
a Solution is created and visualized as usual.

Sample results
~~~~~~~~~~~~~~

You should see the following result:

.. image:: 41/1.png
   :align: center
   :width: 800
   :alt: Sample result
