Newton's Method and Adaptivity (01-newton-adapt)
-----------------------------------

**Git reference:** Tutorial example `01-newton-adapt 
<http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P05-nonlinear-adapt/01-newton-adapt>`_.

We still keep the simple model problem

.. math::

    -\nabla \cdot (\lambda(u)\nabla u) - f(x,y) = 0 \ \ \ \mbox{in } \Omega = (-10,10)^2,

equipped with nonhomogeneous Dirichlet boundary conditions 

.. math::

    u(x, y) = (x+10)(y+10)/100 \ \ \ \mbox{on } \partial \Omega,

but this time it will be solved using automatic adaptivity. 

Initial solve on coarse mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We solve the nonlinear problem on the coarse mesh once, to obtain a good starting 
point for the Newton's method on the reference mesh.

First we calculate an initial coefficient vector $\bfY_0$ by projecting 
the initial condition on the coarse mesh::

    // Project the initial condition on the FE space to obtain initial 
    // coefficient vector for the Newton's method.
    info("Projecting initial condition to obtain initial vector on the coarse mesh.");
    scalar* coeff_vec_coarse = new scalar[Space::get_num_dofs(&space)] ;
    Solution* init_sln = new Solution(&mesh, init_cond);
    OGProjection::project_global(&space, init_sln, coeff_vec_coarse, matrix_solver); 
    delete init_sln;

Then the Newton's loop is run::

    // Newton's loop on the coarse mesh. This is needed to obtain a good 
    // starting point for the Newton's method on the reference mesh.
    bool verbose = true;
    if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, 
        NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

The resulting coefficient vector is translated into a coarse mesh solution::

    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solution(coeff_vec_coarse, &space, &sln);

Adaptivity loop
~~~~~~~~~~~~~~~

At the beginning of the adaptivity loop, we construct a reference space (and mesh)::

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = construct_refined_space(&space);

Next, the discrete problem on the reference mesh is initialized::

    // Initialize discrete problem on the reference mesh.
    DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space, is_linear);

We initialize matrix solver::

    // Initialize matrix solver.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

Calculating initial coefficient vector on the reference mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This step is very important since the reference space is large, and the 
quality of the initial coefficient vector matters a lot. In the first 
adaptivity step, we use the coarse mesh solution (that's why we have 
computer it), and in all other steps we use the previous fine mesh 
solution::

    // Calculate initial coefficient vector on the reference mesh.
    scalar* coeff_vec = new scalar[Space::get_num_dofs(ref_space)];
    if (as == 1) 
    {
      // In the first step, project the coarse mesh solution.
      info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
      OGProjection::project_global(ref_space, &sln, coeff_vec, matrix_solver);
    }
    else 
    {
      // In all other steps, project the previous fine mesh solution.
      info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
      OGProjection::project_global(ref_space, &ref_sln, coeff_vec, matrix_solver);
    }

Newton's loop on the reference mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next we run the Newton's method on the reference mesh::

    // Newton's loop on the fine mesh.
    info("Solving on fine mesh:");
    if (!solve_newton(coeff_vec, dp, solver, matrix, rhs, 
		      NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

    // Translate the resulting coefficient vector into the Solution ref_sln.
    Solution::vector_to_solution(coeff_vec, ref_space, &ref_sln);

Error estimation
~~~~~~~~~~~~~~~~

Now we have the desired solution pair to guide automatic adaptivity, and we can calculate 
the error estimates::

    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    Adapt* adaptivity = new Adapt(&space, HERMES_H1_NORM);
    bool solutions_for_adapt = true;
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt, 
                         HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

Adapting the coarse mesh
~~~~~~~~~~~~~~~~~~~~~~~~

Then we adapt the coarse mesh, and project the fine mesh solution on the new
coarse mesh::

    // If err_est_rel too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting the coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (Space::get_num_dofs(&space) >= NDOF_STOP) 
      {
        done = true;
        break;
      }
      
      // Project last fine mesh solution on the new coarse mesh
      // to obtain new coars emesh solution.
      info("Projecting reference solution on new coarse mesh for error calculation.");
      OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 

      // View the coarse mesh solution.
      sview.show(&sln);
      oview.show(&space);
    }

Sample results
~~~~~~~~~~~~~~

In our experience, the Newton's loop on the new coarse mesh can be skipped since this 
does not affect convergence and one saves some CPU time. This is illustrated in the 
following convergence comparison:

Convergence in the number of DOF (with and without Newton solve on the new coarse mesh):

.. image:: newton-adapt/conv_dof_compar.png
   :align: center
   :width: 600
   :height: 400
   :alt: DOF convergence graph for tutorial example 01-newton-adapt.

Convergence in CPU time (with and without Newton solve on coarse mesh):

.. image:: newton-adapt/conv_cpu_compar.png
   :align: center
   :width: 600
   :height: 400
   :alt: CPU convergence graph for tutorial example 01-newton-adapt.

