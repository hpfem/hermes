#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Teuchos;
using namespace RefinementSelectors;

// This example solves the compressible Euler equations coupled with an advection-diffution equation
// using a basic piecewise-constant finite volume method for the flow and continuous FEM for the concentration
// being advected by the flow.
//
// Equations: Compressible Euler equations, perfect gas state equation, advection-diffusion equation.
//
// Domains: Various
//
// BC: Various.
//
// IC: Various.
//
// The following parameters can be changed:

// Visualization.
const bool HERMES_VISUALIZATION = false;          // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const unsigned int EVERY_NTH_STEP = 1;            // Set visual output for every nth step.

// Use of preconditioning.
const bool PRECONDITIONING = true;
const double NOX_LINEAR_TOLERANCE = 1E-1;
const double NOX_NONLINEAR_TOLERANCE = 1E-1;
unsigned NOX_MESSAGE_TYPE = NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration | NOX::Utils::InnerIteration | NOX::Utils::Parameters | NOX::Utils::Details | NOX::Utils::LinearSolverDetails;

// Shock capturing.
bool SHOCK_CAPTURING = false;

// Quantitative parameter of the discontinuity detector.
double DISCONTINUITY_DETECTOR_PARAM = 1;

const int P_INIT_FLOW = 0;                        // Polynomial degree for the Euler equations (for the flow).
const int P_INIT_CONCENTRATION = 1;               // Polynomial degree for the concentration.
double time_step = 1E-2;                          // Time step.

// Adaptivity.
const int UNREF_FREQ = 5;                         // Every UNREF_FREQth time step the mesh is unrefined.
int REFINEMENT_COUNT = 0;                         // Number of mesh refinements between two unrefinements.
                                                  // The mesh is not unrefined unless there has been a refinement since
                                                  // last unrefinement.
const double THRESHOLD = 0.1;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //  error is processed. If more elements have similar errors, refine
                                                  //  all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //  than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                 //  than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST_FLOW = H2D_H_ANISO,      // Predefined list of element refinement candidates. Possible values are
      CAND_LIST_CONCENTRATION = H2D_HP_ANISO;     // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1;                        // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;                     // Adaptivity process stops when the number of degrees of freedom grows over
                                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
                                                  // Matrix solver for orthogonal projections.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

unsigned int INIT_REF_NUM_FLOW = 3;               // Number of initial uniform mesh refinements of the mesh for the flow.
unsigned int INIT_REF_NUM_CONCENTRATION = 3;      // Number of initial uniform mesh refinements of the mesh for the concentration.
unsigned int INIT_REF_NUM_CONCENTRATION_BDY = 0;  // Number of initial mesh refinements of the mesh for the concentration towards the 
                                                  // part of the boundary where the concentration is prescribed.
// Equation parameters.
const double P_EXT = 2.5;                               // Exterior pressure (dimensionless).
const double RHO_EXT = 1.0;                             // Inlet density (dimensionless).   
const double V1_EXT = 0.0;                              // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.25;                             // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;                               // Kappa.
const double CONCENTRATION_EXT = 1.0;                   // Concentration on the boundary.

const double EPSILON = 0.01;                           // Diffusivity.

// Boundary markers.
const std::string BDY_INLET = "Inflow";
const std::string BDY_OUTLET = "Outflow";
const std::string BDY_SOLID_WALL = "Solid Wall";
const std::string BDY_DIRICHLET_CONCENTRATION = "Inflow";
Hermes::vector<std::string> BDY_NATURAL_CONCENTRATION = Hermes::vector<std::string>("Outflow", "Solid Wall");

// Weak forms.
#include "../forms_implicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh basemesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);

  // Initialize the meshes.
  Mesh mesh_flow, mesh_concentration;
  mesh_flow.copy(&basemesh);
  mesh_concentration.copy(&basemesh);

  for(unsigned int i = 0; i < INIT_REF_NUM_CONCENTRATION; i++)
    mesh_concentration.refine_all_elements();

  mesh_concentration.refine_towards_boundary(BDY_DIRICHLET_CONCENTRATION, INIT_REF_NUM_CONCENTRATION_BDY);
  
  for(unsigned int i = 0; i < INIT_REF_NUM_FLOW; i++)
    mesh_flow.refine_all_elements();

  // Initialize boundary condition types and spaces with default shapesets.
  // For the concentration.
  EssentialBCs bcs_concentration;

  bcs_concentration.add_boundary_condition(new DefaultEssentialBCConst(BDY_DIRICHLET_CONCENTRATION, CONCENTRATION_EXT));
  
  L2Space space_rho(&mesh_flow, P_INIT_FLOW);
  L2Space space_rho_v_x(&mesh_flow, P_INIT_FLOW);
  L2Space space_rho_v_y(&mesh_flow, P_INIT_FLOW);
  L2Space space_e(&mesh_flow, P_INIT_FLOW);
  // Space for concentration.
  H1Space space_c(&mesh_concentration, &bcs_concentration, P_INIT_CONCENTRATION);

  int ndof = Space::get_num_dofs(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_c));
  info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  InitialSolutionEulerDensity sln_rho(&mesh_flow, RHO_EXT);
  InitialSolutionEulerDensityVelX sln_rho_v_x(&mesh_flow, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY_LShape sln_rho_v_y(&mesh_flow, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy_LShape sln_e(&mesh_flow, &sln_rho, &sln_rho_v_x, &sln_rho_v_y, P_EXT, KAPPA);
  InitialSolutionConcentration sln_c(&mesh_concentration, 0.0);

  InitialSolutionEulerDensity prev_rho(&mesh_flow, RHO_EXT);
  InitialSolutionEulerDensityVelX prev_rho_v_x(&mesh_flow, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY_LShape prev_rho_v_y(&mesh_flow, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy_LShape prev_e(&mesh_flow, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, P_EXT, KAPPA);
  InitialSolutionConcentration prev_c(&mesh_concentration, 0.0);

  Solution rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e, rsln_c;

  // Numerical flux.
  OsherSolomonNumericalFlux num_flux(KAPPA);

  // Initialize weak formulation.
  EulerEquationsWeakFormImplicitCoupled wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL,
    BDY_INLET, BDY_OUTLET, BDY_NATURAL_CONCENTRATION, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c, PRECONDITIONING, EPSILON, true);
  
  wf.set_time_step(time_step);

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA, RHO_EXT, P_EXT);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));

  ScalarView s5("Concentration", new WinGeom(700, 400, 600, 300));
  
  /*
  ScalarView s1("1", new WinGeom(0, 0, 600, 300));
  ScalarView s2("2", new WinGeom(700, 0, 600, 300));
  ScalarView s3("3", new WinGeom(0, 400, 600, 300));
  ScalarView s4("4", new WinGeom(700, 400, 600, 300));
  ScalarView s5("Concentration", new WinGeom(350, 200, 600, 300));
  */

  // Initialize refinement selector.
  L2ProjBasedSelector l2selector(CAND_LIST_FLOW, CONV_EXP, H2DRS_DEFAULT_ORDER);
  H1ProjBasedSelector h1selector(CAND_LIST_CONCENTRATION, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Select preconditioner.
  RCP<Precond> pc = rcp(new IfpackPrecond("point-relax"));

  int iteration = 0; double t = 0;
  for(t = 0.0; t < 3.0; t += time_step) {
    info("---- Time step %d, time %3.5f.", iteration++, t);

    // Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0 && REFINEMENT_COUNT > 0) {
      REFINEMENT_COUNT = 0;
      info("Global mesh derefinement.");
      mesh_flow.unrefine_all_elements();
      mesh_concentration.unrefine_all_elements();
      /*
      space_rho.adjust_element_order(-1, P_INIT_FLOW);
      space_rho_v_x.adjust_element_order(-1, P_INIT_FLOW);
      space_rho_v_y.adjust_element_order(-1, P_INIT_FLOW);
      space_e.adjust_element_order(-1, P_INIT_FLOW);
      */
      space_c.adjust_element_order(-1, P_INIT_CONCENTRATION);
    }

    // Adaptivity loop:
    int as = 1; 
    bool done = false;
    do {
      info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      int order_increase = 0;
      Hermes::vector<Space *>* ref_spaces = Space::construct_refined_spaces(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e, &space_c), order_increase);
      (*ref_spaces)[4]->adjust_element_order(+1, P_INIT_CONCENTRATION);

      // Report NDOFs.
      info("ndof_coarse: %d, ndof_fine: %d.", 
        Space::get_num_dofs(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e, &space_c)), Space::get_num_dofs(*ref_spaces));

      // Very imporant, set the meshes for the flow as the same.
      (*ref_spaces)[1]->get_mesh()->set_seq((*ref_spaces)[0]->get_mesh()->get_seq());
      (*ref_spaces)[2]->get_mesh()->set_seq((*ref_spaces)[0]->get_mesh()->get_seq());
      (*ref_spaces)[3]->get_mesh()->set_seq((*ref_spaces)[0]->get_mesh()->get_seq());

      // Project the previous time level solution onto the new fine mesh
      // in order to obtain initial vector for NOX. 
      info("Projecting initial solution on the FE mesh.");
      scalar* coeff_vec = new scalar[Space::get_num_dofs(*ref_spaces)];
      OGProjection::project_global(*ref_spaces, Hermes::vector<MeshFunction *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c), coeff_vec);

      // Initialize the FE problem.
      DiscreteProblem dp(&wf, *ref_spaces);
      
      // Initialize NOX solver.
      NoxSolver solver(&dp, NOX_MESSAGE_TYPE, "GMRES", "Newton", NOX_LINEAR_TOLERANCE, "None", 0, 0, 1, NOX_NONLINEAR_TOLERANCE);
      solver.set_init_sln(coeff_vec);
      if(PRECONDITIONING)
        solver.set_precond(pc);

      info("Assembling by DiscreteProblem, solving by NOX.");
      if (solver.solve())
        Solution::vector_to_solutions(solver.get_solution(), *ref_spaces, 
          Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e, &rsln_c));
      else
        error("NOX failed.");
      
      info("Number of nonlin iterations: %d (norm of residual: %g)", 
        solver.get_num_iters(), solver.get_residual());
      info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
        solver.get_num_lin_iters(), solver.get_achieved_tol());
      
      if(SHOCK_CAPTURING) {
        Hermes::vector<Space*> flow_spaces((*ref_spaces)[0], (*ref_spaces)[1], (*ref_spaces)[2], (*ref_spaces)[3]);
        
        scalar* flow_solution_vector = new scalar[Space::get_num_dofs(flow_spaces)];

        OGProjection::project_global(flow_spaces, Hermes::vector<MeshFunction *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), flow_solution_vector);

        DiscontinuityDetector discontinuity_detector(flow_spaces, Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));

        std::set<int> discontinuous_elements = discontinuity_detector.get_discontinuous_element_ids(DISCONTINUITY_DETECTOR_PARAM);

        FluxLimiter flux_limiter(flow_solution_vector, flow_spaces, Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));

        flux_limiter.limit_according_to_detector(discontinuous_elements);
      }

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection::project_global(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e, &space_c), Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e, &rsln_c), 
                     Hermes::vector<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e, &sln_c), matrix_solver, 
                     Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_H1_NORM)); 

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e, &space_c), Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_H1_NORM));
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e, &sln_c),
							  Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e, &rsln_c)) * 100;

      // Report results.
      info("err_est_rel: %g%%", err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) 
        done = true;
      else {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *>(&l2selector, &l2selector, &l2selector, &l2selector, &h1selector), 
                                 THRESHOLD, STRATEGY, MESH_REGULARITY);

        REFINEMENT_COUNT++;
        if (Space::get_num_dofs(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e, &space_c)) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }

      // Save orders.
      Orderizer ord;
      char filename[40];
      sprintf(filename, "Flow-mesh-%i-%i.vtk", iteration - 1, as - 1);
      ord.save_orders_vtk((*ref_spaces)[0], filename);
      sprintf(filename, "Concentration-mesh-%i-%i.vtk", iteration - 1, as - 1);
      ord.save_orders_vtk((*ref_spaces)[4], filename);
      
      // Clean up.
      delete adaptivity;
      if(!done)
        for(unsigned int i = 0; i < ref_spaces->size(); i++)
          delete (*ref_spaces)[i]->get_mesh();


      for(unsigned int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i];
    }
    while (done == false);

    // Copy the solutions into previous time level ones.
    prev_rho.copy(&rsln_rho);
    prev_rho_v_x.copy(&rsln_rho_v_x);
    prev_rho_v_y.copy(&rsln_rho_v_y);
    prev_e.copy(&rsln_e);
    prev_c.copy(&rsln_c);

    // Visualization.
    if((iteration - 1) % EVERY_NTH_STEP == 0) {
      // Hermes visualization.
      if(HERMES_VISUALIZATION) {
        Mach_number.reinit();
        pressure.reinit();
        entropy.reinit();
        pressure_view.show(&pressure);
        entropy_production_view.show(&entropy);
        Mach_number_view.show(&Mach_number);
        s5.show(&prev_c);
        /*
        s1.show(&prev_rho);
        s2.show(&prev_rho_v_x);
        s3.show(&prev_rho_v_y);
        s4.show(&prev_e);
        s5.show(&prev_c);
        */
      }
      // Output solution in VTK format.
      if(VTK_VISUALIZATION) {
        pressure.reinit();
        Mach_number.reinit();
        Linearizer lin;
        char filename[40];
        sprintf(filename, "pressure-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&pressure, filename, "Pressure", false);
        sprintf(filename, "pressure-3D-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&pressure, filename, "Pressure", true);
        sprintf(filename, "Mach number-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&Mach_number, filename, "MachNumber", false);
        sprintf(filename, "Mach number-3D-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&Mach_number, filename, "MachNumber", true);
        sprintf(filename, "Concentration-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&prev_c, filename, "Concentration", true);
        sprintf(filename, "Concentration-3D-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&prev_c, filename, "Concentration", true);
 
      }
    }
  }
  
  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();
  s5.close();

  /*
  s1.close();
  s2.close();
  s3.close();
  s4.close();
  */

  return 0;
}
