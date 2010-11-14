Adapting Mesh to an Exact Function (33)
---------------------------------------

**Git reference:** Tutorial example `33-exact-adapt <http://git.hpfem.org/hermes.git/tree/HEAD:/hermes2d/tutorial/33-exact-adapt>`_. 

This technique can be useful, for example, when a time-dependent proces
starts from a complicated initial condition that would not be represented
with sufficient accuracy on a coarse initial mesh. 

As usual, the adaptivity algorithm expects a pair of solutions on the 
coarse and globally refined meshes. So the adaptivity loop begins with 
refining the coarse mesh::

    // Construct globally refined reference mesh
    // and setup reference space.
    Mesh *ref_mesh = new Mesh();
    ref_mesh->copy(space.get_mesh());
    ref_mesh->refine_all_elements();
    Space* ref_space = space.dup(ref_mesh);
    int order_increase = 1;
    ref_space->copy_orders(&space, order_increase);

Instead of calculating a solution on the fine mesh, we set the exact 
function::

    // Assign the function f() to the fine mesh.
    sln_fine.set_exact(ref_mesh, f);

The coarse mesh solution is obtained by projecting 'sln_fine'::

    // Project the function f() on the coarse mesh.
    project_global(&space, H2D_H1_NORM, &sln_file, &sln_coarse);

Error estimates are calculated as usual::

    // Calculate element errors and total error estimate.
    info("Calculating error.");
    Adapt hp(&space, H2D_H1_NORM);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est_rel = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;

Mesh adaptation is standard as well::

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else {
      info("Adapting coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (get_num_dofs(&space) >= NDOF_STOP) done = true;
    }

Sample solution and mesh are shown below:

.. image:: 33/img.png
   :align: center
   :width: 800
   :alt: Resulting solution and mesh

