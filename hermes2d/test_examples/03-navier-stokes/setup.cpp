

// Load the mesh.
MeshSharedPtr mesh(new Mesh);
MeshReaderH2DXML mloader;
mloader.load("domain.xml", mesh);

// Initial mesh refinements.
for (int i = 0; i < INIT_REF; i++)
    mesh->refine_all_elements();

// Initialize boundary conditions.
EssentialBCNonConst bc_left_vel_x(BDY_TOP, VEL_INLET, H, STARTUP_TIME);
Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_other_vel_x(std::vector<std::string>({ BDY_LEFT, BDY_RIGHT }), 0.0);
Hermes::Hermes2D::EssentialBCs<double> bcs_vel_x({ &bc_left_vel_x, &bc_other_vel_x });
Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_vel_y(std::vector<std::string>({ BDY_LEFT, BDY_RIGHT }), 0.0);
Hermes::Hermes2D::EssentialBCs<double> bcs_vel_y(&bc_vel_y);
Hermes::Hermes2D::EssentialBCs<double> bcs_pressure;

// Spaces for velocity components and pressure.
SpaceSharedPtr<double> xvel_space(new H1Space<double>(mesh, &bcs_vel_x, P_INIT_VEL));
SpaceSharedPtr<double> yvel_space(new H1Space<double>(mesh, &bcs_vel_y, P_INIT_VEL));
#ifdef PRESSURE_IN_L2
SpaceSharedPtr<double> p_space(new L2Space<double>(mesh, P_INIT_PRESSURE));
#else
SpaceSharedPtr<double> p_space(new H1Space<double>(mesh, &bcs_pressure, P_INIT_PRESSURE));
#endif
std::vector<SpaceSharedPtr<double> > spaces({ xvel_space, yvel_space, p_space });

// Calculate and report the number of degrees of freedom.
int ndof = Space<double>::get_num_dofs(spaces);

// Define projection norms.
NormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
NormType p_proj_norm = HERMES_L2_NORM;
#else
NormType p_proj_norm = HERMES_H1_NORM;
#endif
std::vector<NormType> proj_norms({ vel_proj_norm, vel_proj_norm, p_proj_norm });

// Solutions for the Newton's iteration and time stepping.
MeshFunctionSharedPtr<double> xvel_prev_time(new ConstantSolution<double>(mesh, 0.0));
MeshFunctionSharedPtr<double> yvel_prev_time(new ConstantSolution<double>(mesh, 0.0));
MeshFunctionSharedPtr<double> p_prev_time(new ConstantSolution<double>(mesh, 0.0));
std::vector<MeshFunctionSharedPtr<double> > sln_prev_time = { xvel_prev_time, yvel_prev_time, p_prev_time };
double* coeff_vec = new double[Space<double>::get_num_dofs(spaces)];

// Project the initial condition on the FE space to obtain initial coefficient vector for the Newton's method.
OGProjection<double>::project_global(spaces, sln_prev_time, coeff_vec, proj_norms);

// Initialize weak formulation.
WeakFormSharedPtr<double> wf(new WeakFormNSNewton(STOKES, RE, TAU, xvel_prev_time, yvel_prev_time));
UExtFunctionSharedPtr<double> fn_0(new CustomUExtFunction(0));
UExtFunctionSharedPtr<double> fn_1(new CustomUExtFunction(1));
wf->set_ext({ xvel_prev_time, yvel_prev_time });
wf->set_u_ext_fn({fn_0, fn_1});

// Initialize views.
Views::VectorView vview("velocity[m/s]", new Views::WinGeom(0, 0, 350, 240));
Views::ScalarView pview("pressure[Pa]", new Views::WinGeom(0, 290, 350, 240));
vview.fix_scale_width(80);
//pview.set_min_max_range(-0.9, 1.0);
pview.fix_scale_width(80);
pview.show_mesh(true);
