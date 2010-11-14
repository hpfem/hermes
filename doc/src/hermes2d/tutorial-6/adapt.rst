Trilinos - Adapt (43)
---------------------

**Git reference:** Example `trilinos-adapt
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/43-trilinos-adapt>`_.

The purpose of this example is to show how to use Trilinos while adapting mesh.
Solved by NOX solver, either using Newton's method or JFNK, with or without 
preconditioning. The underlying problem is benchmark 
`layer-internal <http://hpfem.org/hermes/doc/src/hermes2d/benchmarks.html#interior-layer-elliptic>`_.

One little difference vs. benchmark "layer-internal" is that we'll be solving the 
finite element problem both on the coarse and fine meshes in each adaptivity step.
So, at the beginning of each adaptivity step we initialize the FeProblem class,
NOX solver, and preconditioner on the coarse mesh::

    info("---- Adaptivity step %d:", as);
   
    // Initialize finite element problem.
    FeProblem fep(&wf, &space);

    // Initialize NOX solver.
    NoxSolver solver(&fep);

    // Choose preconditioner.
    RCP<Precond> pc = rcp(new MlPrecond("sa"));
    if (PRECOND)
    {
      if (JFNK) solver.set_precond(pc);
      else solver.set_precond("ML");
    }

Then we assemble and solve on coarse mesh, and convert the resulting 
coefficient vector into a Solution. Skipping info outputs and 
visualization, this reads::

    // Assemble on coarse mesh and solve the matrix problem using NOX.
    bool solved = solver.solve();
    if (solved)
    {
      double* coeffs = solver.get_solution_vector();
      sln.set_coeff_vector(&space, coeffs, ndof);
    }
    else
      error("NOX failed on coarse mesh.");

Next we create a uniformly refined mesh and H1 space on it::

    // Create uniformly refined reference mesh.
    Mesh rmesh; rmesh.copy(&mesh); 
    rmesh.refine_all_elements();
    // Reference FE space.
    H1Space rspace(&rmesh, bc_types, essential_bc_values, P_INIT);
    int order_increase = 1;
    rspace.copy_orders(&space, order_increase); // increase orders by one

Then the FeProblem, NOX solver and preconditioner are initialized
on the fine mesh::

    // Initialize FE problem on reference mesh.
    FeProblem ref_fep(&wf, &rspace);

    // Initialize NOX solver.
    NoxSolver ref_solver(&ref_fep);
    if (PRECOND)
    {
      if (JFNK) ref_solver.set_precond(pc);
      else ref_solver.set_precond("ML");
    }

Fine mesh problem is solved and the solution coefficient vector converted
into a Solution. Again, skipping info outputs and visualization this reads::

    // Assemble on fine mesh and solve the matrix problem using NOX.
    solved = ref_solver.solve();
    if (solved)
    {
      double* s = ref_solver.get_solution_vector();
      ref_sln.set_coeff_vector(&rspace, coeffs, ndof);
    }
    else
      error("NOX failed on fine mesh.");

Hence now we have the two solutions to guide automatic hp-adaptivity and 
to adapt the mesh, we proceed as in benchmark "layer-internal".


