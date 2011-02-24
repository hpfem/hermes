Adapting Mesh to an Exact Function (06-exact-adapt)
---------------------------------------

**Git reference:** Tutorial example `06-exact-adapt <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/P04-linear-adapt/06-exact-adapt>`_. 

This technique can be useful, for example, when a time-dependent proces
starts from a complicated initial condition that would not be represented
with sufficient accuracy on a coarse initial mesh. 

As usual, the adaptivity algorithm expects a pair of solutions on the 
coarse and globally refined meshes. So the adaptivity loop begins with 
refining the coarse mesh::

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = construct_refined_space(&space);

Instead of calculating a solution on the fine mesh, we set the exact 
function::

    // Assign the function f() to the fine mesh.
    info("Assigning f() to the reference mesh.");
    bool is_linear = true;
    ref_sln.set_exact(ref_space->get_mesh(), f);

The coarse mesh solution is obtained by projecting 'sln_fine'::

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 

Error estimates are calculated as usual::

    // Calculate element errors and total error estimate.
    info("Calculating exact error."); 
    Adapt* adaptivity = new Adapt(&space, HERMES_H1_NORM);
    // Note: the error estimate is now equal to the exact error.
    bool solutions_for_adapt = true;
    double err_exact_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt, 
                           HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

Mesh adaptation is standard as well::

    // If err_exact_rel too large, adapt the mesh.
    if (err_exact_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      
      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
    if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;

Sample solution and mesh are shown below:

.. image:: adapt_exact/img.png
   :align: center
   :width: 800
   :alt: Resulting solution and mesh

