Hermes::vector<MeshFunctionSharedPtr<double> > prev_slns(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e);
EulerEquationsWeakFormStabilization wf_stabilization(prev_rho);

if(SHOCK_CAPTURING && SHOCK_CAPTURING_TYPE == FEISTAUER)
  wf.set_stabilization(prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, NU_1, NU_2);

// Solver.
LinearSolver<double> solver(&wf, spaces);

#pragma region 6. Time stepping loop.
int iteration = 0;
for(double t = 0.0; t < TIME_INTERVAL_LENGTH; t += time_step_n)
{
  Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f.", iteration++, t);

#pragma region 6.1. Periodic global derefinements.
  if (iteration > 1 && iteration % UNREF_FREQ == 0 && REFINEMENT_COUNT > 0) 
  {
    Hermes::Mixins::Loggable::Static::info("Global mesh derefinement.");
    REFINEMENT_COUNT = 0;

    space_rho->unrefine_all_mesh_elements(true);

    space_rho->adjust_element_order(-1, P_INIT);
    space_rho_v_x->adjust_element_order(-1, P_INIT);
    space_rho_v_y->adjust_element_order(-1, P_INIT);
    space_e->adjust_element_order(-1, P_INIT);
    Space<double>::assign_dofs(spaces);
  }
#pragma endregion

#pragma region 7. Adaptivity loop.
  int as = 1; int ndofs_prev = 0; bool done = false;
  do
  {
    // Info.
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);
    // Set the current time step.
    wf.set_current_time_step(time_step_n);

#pragma region 7.1. Construct globally refined reference mesh and setup reference space.
    int order_increase = CAND_LIST == H2D_HP_ANISO ? 1 : 0;

    Mesh::ReferenceMeshCreator refMeshCreatorFlow(mesh);
    MeshSharedPtr ref_mesh = refMeshCreatorFlow.create_ref_mesh();

    Space<double>::ReferenceSpaceCreator refSpaceCreatorRho(space_rho, ref_mesh, order_increase);
    SpaceSharedPtr<double> ref_space_rho = refSpaceCreatorRho.create_ref_space();
    Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVx(space_rho_v_x, ref_mesh, order_increase);
    SpaceSharedPtr<double> ref_space_rho_v_x = refSpaceCreatorRhoVx.create_ref_space();
    Space<double>::ReferenceSpaceCreator refSpaceCreatorRhoVy(space_rho_v_y, ref_mesh, order_increase);
    SpaceSharedPtr<double> ref_space_rho_v_y = refSpaceCreatorRhoVy.create_ref_space();
    Space<double>::ReferenceSpaceCreator refSpaceCreatorE(space_e, ref_mesh, order_increase);
    SpaceSharedPtr<double> ref_space_e = refSpaceCreatorE.create_ref_space();
    Hermes::vector<SpaceSharedPtr<double>  > ref_spaces(ref_space_rho, ref_space_rho_v_x, ref_space_rho_v_y, ref_space_e);
    solver.set_spaces(ref_spaces);
    
    if(ndofs_prev != 0)
      if(Space<double>::get_num_dofs(ref_spaces) == ndofs_prev)
        selector.set_error_weights(2.0 * selector.get_error_weight_h(), 1.0, 1.0);
      else
        selector.set_error_weights(1.0, 1.0, 1.0);

    ndofs_prev = Space<double>::get_num_dofs(ref_spaces);

    // Project the previous time level solution onto the new fine mesh
    Hermes::Mixins::Loggable::Static::info("Projecting the previous time level solution onto the new fine mesh.");
    OGProjection<double>::project_global(ref_spaces, prev_slns, prev_slns);
#pragma endregion

    if(SHOCK_CAPTURING && SHOCK_CAPTURING_TYPE == FEISTAUER)
    {
      SpaceSharedPtr<double> ref_space_stabilization(new L2Space<double>(ref_mesh, 0));
      int mesh_size = ref_mesh->get_num_active_elements();
      DiscreteProblem<double> dp_stabilization(&wf_stabilization, ref_space_stabilization);
      dp_stabilization.set_space(ref_space_stabilization);
      dp_stabilization.assemble(rhs_stabilization);
      if(!wf.discreteIndicator)
      {
        wf.set_discreteIndicator(new bool[mesh_size], mesh_size);
        memset(wf.discreteIndicator, 0, mesh_size * sizeof(bool));
      }
      Element* e;
      for_all_active_elements(e, ref_space_stabilization->get_mesh())
      {
        AsmList<double> al;
        ref_space_stabilization->get_element_assembly_list(e, &al);
        if(rhs_stabilization->get(al.get_dof()[0]) >= 1)
          wf.discreteIndicator[e->id] = true;
      }
    }

    // Solve the problem.
    solver.solve();

#pragma region *. Get the solution with optional shock capturing.
    if(!SHOCK_CAPTURING)
      Solution<double>::vector_to_solutions(solver.get_sln_vector(), ref_spaces, rslns);
    else
    {
      if(SHOCK_CAPTURING_TYPE == KRIVODONOVA)
      {
        FluxLimiter* flux_limiter = new FluxLimiter(FluxLimiter::Krivodonova, solver.get_sln_vector(), ref_spaces);
        flux_limiter->limit_according_to_detector();
        flux_limiter->get_limited_solutions(rslns);
        delete flux_limiter;
      }

      if(SHOCK_CAPTURING_TYPE == KUZMIN)
      {
        PostProcessing::VertexBasedLimiter limiter(ref_spaces, solver.get_sln_vector(), 1);
        limiter.get_solutions(rslns);
      }
    }
#pragma endregion

    // Calculate time step according to CFL condition.
    CFL.calculate(rslns, (ref_spaces)[0]->get_mesh(), time_step_n);

#pragma region 7.2. Project to coarse mesh -> error estimation -> space adaptivity
    // Project the fine mesh solution onto the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
    OGProjection<double>::project_global(spaces, rslns, slns, Hermes::vector<NormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");
    errorCalculator.calculate_errors(slns, rslns);
    double err_est_rel_total = errorCalculator.get_total_error_squared() * 100;

    // Report results.
    Hermes::Mixins::Loggable::Static::info("err_est_rel: %g%%", err_est_rel_total);

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < adaptivityErrorStop(iteration))
      done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
      done = adaptivity.adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector, &selector, &selector));
      REFINEMENT_COUNT++;
      as++;
    }
#pragma endregion

#pragma region 7.3. Visualization and saving on disk.
    if(done && (iteration - 1) % EVERY_NTH_STEP == 0)
    {
      // Hermes visualization.
      if(HERMES_VISUALIZATION)
      {        
        Mach_number->reinit();
        pressure->reinit();
        pressure_view.show(pressure, 1);
        Mach_number_view.show(Mach_number, 1);
        order_view.show((ref_spaces)[0]);
      }
      // Output solution in VTK format.
      if(VTK_VISUALIZATION)
      {
        pressure->reinit();
        Linearizer lin(FileExport);
        char filename[40];
        sprintf(filename, "Pressure-%i.vtk", iteration - 1);
        lin.save_solution_vtk(pressure, filename, "Pressure", false);
        sprintf(filename, "VelocityX-%i.vtk", iteration - 1);
        lin.save_solution_vtk(prev_rho_v_x, filename, "VelocityX", false);
        sprintf(filename, "VelocityY-%i.vtk", iteration - 1);
        lin.save_solution_vtk(prev_rho_v_y, filename, "VelocityY", false);
        sprintf(filename, "Rho-%i.vtk", iteration - 1);
        lin.save_solution_vtk(prev_rho, filename, "Rho", false);
      }
    }
#pragma endregion
  }
  while (done == false);
#pragma endregion

  // Copy the solutions into the previous time level ones.
  prev_rho->copy(rsln_rho);
  prev_rho_v_x->copy(rsln_rho_v_x);
  prev_rho_v_y->copy(rsln_rho_v_y);
  prev_e->copy(rsln_e);
}
#pragma endregion

pressure_view.close();
Mach_number_view.close();

return 0;