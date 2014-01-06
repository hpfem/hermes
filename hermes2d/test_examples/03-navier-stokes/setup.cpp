

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Initial mesh refinements.
  mesh->refine_towards_boundary(BDY_OBSTACLE, 2, false);
  mesh->refine_towards_boundary(BDY_TOP, 2, true);     // '4' is the number of levels,
  mesh->refine_towards_boundary(BDY_BOTTOM, 2, true);  // 'true' stands for anisotropic refinements.
  mesh->refine_all_elements();

  // Initialize boundary conditions.
  EssentialBCNonConst bc_left_vel_x(BDY_LEFT, VEL_INLET, H, STARTUP_TIME);
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_other_vel_x(Hermes::vector<std::string>(BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE), 0.0);
  Hermes::Hermes2D::EssentialBCs<double> bcs_vel_x(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc_left_vel_x, &bc_other_vel_x));
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_vel_y(Hermes::vector<std::string>(BDY_LEFT, BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE), 0.0);
  Hermes::Hermes2D::EssentialBCs<double> bcs_vel_y(&bc_vel_y);
  Hermes::Hermes2D::EssentialBCs<double> bcs_pressure;

  // Spaces for velocity components and pressure.
  SpaceSharedPtr<double> xvel_space(new H1Space<double>(mesh, &bcs_vel_x, P_INIT_VEL));
  SpaceSharedPtr<double> yvel_space(new H1Space<double>(mesh, &bcs_vel_y, P_INIT_VEL));
#ifdef PRESSURE_IN_L2
  SpaceSharedPtr<double> p_space(new L2Space<double> (mesh, P_INIT_PRESSURE));
#else
  SpaceSharedPtr<double> p_space(new H1Space<double> (mesh, &bcs_pressure, P_INIT_PRESSURE));
#endif
  Hermes::vector<SpaceSharedPtr<double> > spaces(xvel_space, yvel_space, p_space);

  // Calculate and report the number of degrees of freedom.
  int ndof = Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(xvel_space, yvel_space, p_space));

  // Define projection norms.
  NormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  NormType p_proj_norm = HERMES_L2_NORM;
#else
  NormType p_proj_norm = HERMES_H1_NORM;
#endif
  Hermes::vector<NormType> proj_norms(vel_proj_norm, vel_proj_norm, p_proj_norm);

  // Solutions for the Newton's iteration and time stepping.
  MeshFunctionSharedPtr<double> xvel_prev_time(new ConstantSolution<double> (mesh, 0.0));
  MeshFunctionSharedPtr<double> yvel_prev_time(new ConstantSolution<double> (mesh, 0.0));
  MeshFunctionSharedPtr<double> p_prev_time(new ConstantSolution<double> (mesh, 0.0));
  Hermes::vector<MeshFunctionSharedPtr<double> > sln_prev_time(xvel_prev_time, yvel_prev_time, p_prev_time);
  double* coeff_vec = new double[Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(xvel_space, yvel_space, p_space))];

  // Project the initial condition on the FE space to obtain initial coefficient vector for the Newton's method.
  OGProjection<double>::project_global(spaces, sln_prev_time, coeff_vec, proj_norms);

  // Initialize weak formulation.
  WeakFormNSNewton wf(STOKES, RE, TAU, xvel_prev_time, yvel_prev_time);
  UExtFunctionSharedPtr<double> fn_0(new CustomUExtFunction(0));
  UExtFunctionSharedPtr<double> fn_1(new CustomUExtFunction(1));
  wf.set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(xvel_prev_time, yvel_prev_time));
  wf.set_u_ext_fn(Hermes::vector<UExtFunctionSharedPtr<double> >(fn_0, fn_1));

  // Initialize views.
  Views::VectorView vview("velocity[m/s]", new Views::WinGeom(0, 0, 750, 240));
  Views::ScalarView pview("pressure[Pa]", new Views::WinGeom(0, 290, 750, 240));
  vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
  if (HERMES_VISUALIZATION)
    pview.show_mesh(true);