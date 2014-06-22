#pragma region 1. Load mesh and initialize spaces.
    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    MeshReaderH2D mloader;
    mloader.load(MESH_FILENAME, mesh);

    // Perform initial mesh refinements.
    for (int i = 0; i < INIT_REF_NUM; i++) 
      mesh->refine_all_elements(0, true);

    // Initialize boundary condition types and spaces with default shapesets.
    SpaceSharedPtr<double> space_rho(new L2Space<double>(mesh, P_INIT));
    SpaceSharedPtr<double> space_rho_v_x(new L2Space<double>(mesh, P_INIT));
    SpaceSharedPtr<double> space_rho_v_y(new L2Space<double>(mesh, P_INIT));
    SpaceSharedPtr<double> space_e(new L2Space<double>(mesh, P_INIT));
    Hermes::vector<SpaceSharedPtr<double> > spaces(space_rho, space_rho_v_x, space_rho_v_y, space_e);
    int ndof = Space<double>::get_num_dofs(spaces);
    Hermes::Mixins::Loggable::Static::info("ndof: %d", ndof);
  #pragma endregion

  #pragma region 2. Initialize solutions.
    MeshFunctionSharedPtr<double> sln_rho(new Solution<double>(mesh));
    MeshFunctionSharedPtr<double> sln_rho_v_x(new Solution<double> (mesh));
    MeshFunctionSharedPtr<double> sln_rho_v_y(new Solution<double> (mesh));
    MeshFunctionSharedPtr<double> sln_e(new Solution<double> (mesh));
    Hermes::vector<MeshFunctionSharedPtr<double> > slns(sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e);

    MeshFunctionSharedPtr<double> rsln_rho(new Solution<double>(mesh));
    MeshFunctionSharedPtr<double> rsln_rho_v_x(new Solution<double> (mesh));
    MeshFunctionSharedPtr<double> rsln_rho_v_y(new Solution<double> (mesh));
    MeshFunctionSharedPtr<double> rsln_e(new Solution<double> (mesh));
    Hermes::vector<MeshFunctionSharedPtr<double> > rslns(rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e);
  #pragma endregion

  #pragma region 3. Filters for visualization of Mach number, pressure + visualization setup.
    MeshFunctionSharedPtr<double>  Mach_number(new MachNumberFilter(rslns, KAPPA));
    MeshFunctionSharedPtr<double>  pressure(new PressureFilter(rslns, KAPPA));

    ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
    ScalarView Mach_number_view("Mach number", new WinGeom(650, 0, 600, 300));
    ScalarView eview("Error - density", new WinGeom(0, 330, 600, 300));
    ScalarView eview1("Error - momentum", new WinGeom(0, 660, 600, 300));
    OrderView order_view("Orders", new WinGeom(650, 330, 600, 300));
  #pragma endregion

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);
  
  Vector<double>* rhs_stabilization = create_vector<double>(HermesCommonApi.get_integral_param_value(matrixSolverType));
  
  #pragma region 4. Adaptivity setup.
    // Initialize refinement selector.
    L2ProjBasedSelector<double> selector(CAND_LIST);
    selector.set_dof_score_exponent(2.0);

    //selector.set_error_weights(1.0, 1.0, 1.0);

    // Error calculation.
    DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 4);
    // Stopping criterion for an adaptivity step.
    AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
    Adapt<double> adaptivity(spaces, &errorCalculator, &stoppingCriterion);
  #pragma endregion