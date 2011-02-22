Trilinos - Adapt (04-trilinos-adapt)
---------------------

**Git reference:** Example `04-trilinos-adapt
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P09-trilinos/04-trilinos-adapt>`_.

The purpose of this example is to show how to use Trilinos while adapting mesh.
Solved by NOX solver, either using Newton's method or JFNK, with or without 
preconditioning. 

Model problem
~~~~~~~~~~~~~

The underlying problem is benchmark 
`layer-internal <http://hpfem.org/hermes/doc/src/hermes3d/benchmarks/layer-interior.html>`_.

One little difference vs. benchmark "layer-internal" is that we'll be solving the 
finite element problem both on the coarse and fine meshes in each adaptivity step.

Initializations on coarse mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the beginning of each adaptivity step we initialize the DiscreteProblem class,
NOX solver, and preconditioner on the coarse mesh::

    info("---- Adaptivity step %d:", as);
   
    // Initialize finite element problem.
    DiscreteProblem dp(&wf, &space);

    // Initialize NOX solver.
    NoxSolver solver(&dp);

    // Choose preconditioner.
    RCP<Precond> pc = rcp(new MlPrecond("sa"));
    if (PRECOND)
    {
      if (JFNK) solver.set_precond(pc);
      else solver.set_precond("ML");
    }

Assembling and solving on coarse mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Then we assemble and solve on coarse mesh, and convert the resulting 
coefficient vector into a Solution. Skipping info outputs and 
visualization, this reads::

    // Assemble on coarse mesh and solve the matrix problem using NOX.
    int ndof = Space::get_num_dofs(&space);
    info("Coarse mesh problem (ndof: %d): Assembling by DiscreteProblem, solving by NOX.", ndof);
    if (solver.solve())
    {
      Solution::vector_to_solution(solver.get_solution(), &space, &sln);
    }
    else
      error("NOX failed on coarse mesh.");

Creating reference mesh and space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next we create a uniformly refined mesh and H1 space on it::

    // Create uniformly refined reference mesh.
    Mesh rmesh; rmesh.copy(&mesh); 
    rmesh.refine_all_elements();

    // Create reference FE space.
    H1Space rspace(&rmesh, &bc_types, &bc_values, P_INIT);
    int order_increase = 1;
    rspace.copy_orders(&space, order_increase); // Increase orders by one.

Initializations on reference mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Then the DiscreteProblem, NOX solver and preconditioner are initialized
on the fine mesh::

    // Initialize FE problem on reference mesh.
    DiscreteProblem ref_dp(&wf, &rspace);

    // Initialize NOX solver.
    NoxSolver ref_solver(&ref_dp);
    if (PRECOND)
    {
      if (JFNK) ref_solver.set_precond(pc);
      else ref_solver.set_precond("ML");
    }

Assembling and solving on reference mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fine mesh problem is solved and the solution coefficient vector converted
into a Solution. Again, skipping info outputs and visualization this reads::

    // Assemble on fine mesh and solve the matrix problem using NOX.
    ndof = Space::get_num_dofs(&rspace);
    info("Fine mesh problem (ndof: %d): Assembling by DiscreteProblem, solving by NOX.", ndof);
    if (ref_solver.solve())
    {
      Solution::vector_to_solution(ref_solver.get_solution(), &rspace, &ref_sln);
    }
    else
      error("NOX failed on fine mesh.");

The rest
~~~~~~~~

Hence now we have a pair of solutions to guide automatic hp-adaptivity, and 
we proceed as in benchmark "layer-internal".


