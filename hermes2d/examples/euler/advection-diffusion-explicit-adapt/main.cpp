#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

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
const bool HERMES_VISUALIZATION = true;               // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;                  // Set to "true" to enable VTK output.
const unsigned int EVERY_NTH_STEP = 1;            // Set visual output for every nth step.

// Shock capturing.
bool SHOCK_CAPTURING = false;

// Quantitative parameter of the discontinuity detector.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;

// Stability for the concentration part.
double ADVECTION_STABILITY_CONSTANT = 1.0;
const double DIFFUSION_STABILITY_CONSTANT = 1.0;

const int P_INIT_FLOW = 0;                        // Polynomial degree for the Euler equations (for the flow).
const int P_INIT_CONCENTRATION = 1;               // Polynomial degree for the concentration.
double CFL_NUMBER = 1.0;                          // CFL value.
double time_step = 1E-5, util_time_step;          // Initial and utility time step.

// Adaptivity.
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is unrefined.
int REFINEMENT_COUNT_FLOW = 0;                         // Number of mesh refinements between two unrefinements.
                                                  // The mesh is not unrefined unless there has been a refinement since
                                                  // last unrefinement.
int REFINEMENT_COUNT_CONCENTRATION = 0;                         // Number of mesh refinements between two unrefinements.
                                                  // The mesh is not unrefined unless there has been a refinement since
                                                  // last unrefinement.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
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
const double ERR_STOP_FLOW = 1.0;                 // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const double ERR_STOP_CONCENTRATION = 10.0;        // Stopping criterion for adaptivity (rel. error tolerance between the
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
const double V1_EXT = 1.25;                             // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.0;                            // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;                               // Kappa.
const double CONCENTRATION_EXT = 0.1;                  // Concentration on the boundary.
const double CONCENTRATION_EXT_STARTUP_TIME = 0.0;     // Start time of the concentration on the boundary.

const double EPSILON = 0.01;                           // Diffusivity.

// Boundary markers.
const std::string BDY_INLET = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_BOTTOM = "3";
const std::string BDY_SOLID_WALL_TOP = "4";
const std::string BDY_DIRICHLET_CONCENTRATION = "3";
Hermes::vector<std::string> BDY_NATURAL_CONCENTRATION = Hermes::vector<std::string>("2");

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh basemesh;
  H2DReader mloader;
  mloader.load("GAMM-channel.mesh", &basemesh);

  // Initialize the meshes.
  Mesh mesh_flow, mesh_concentration;
  mesh_flow.copy(&basemesh);
  mesh_concentration.copy(&basemesh);

  for(unsigned int i = 0; i < INIT_REF_NUM_CONCENTRATION; i++)
    mesh_concentration.refine_all_elements(0, true);

  mesh_concentration.refine_towards_boundary(BDY_DIRICHLET_CONCENTRATION, INIT_REF_NUM_CONCENTRATION_BDY);
  //mesh_flow.refine_towards_boundary(BDY_DIRICHLET_CONCENTRATION, INIT_REF_NUM_CONCENTRATION_BDY);
  
  for(unsigned int i = 0; i < INIT_REF_NUM_FLOW; i++)
    mesh_flow.refine_all_elements(0, true);

  // Initialize boundary condition types and spaces with default shapesets.
  // For the concentration.
  EssentialBCs bcs_concentration;

  bcs_concentration.add_boundary_condition(new ConcentrationTimedepEssentialBC(BDY_DIRICHLET_CONCENTRATION, CONCENTRATION_EXT, CONCENTRATION_EXT_STARTUP_TIME));
  bcs_concentration.add_boundary_condition(new ConcentrationTimedepEssentialBC(BDY_SOLID_WALL_TOP, 0.0, CONCENTRATION_EXT_STARTUP_TIME));
  bcs_concentration.add_boundary_condition(new ConcentrationTimedepEssentialBC(BDY_INLET, 0.0, CONCENTRATION_EXT_STARTUP_TIME));
  
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
  InitialSolutionEulerDensityVelY sln_rho_v_y(&mesh_flow, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy sln_e(&mesh_flow, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
  InitialSolutionConcentration sln_c(&mesh_concentration, 0.0);

  InitialSolutionEulerDensity prev_rho(&mesh_flow, RHO_EXT);
  InitialSolutionEulerDensityVelX prev_rho_v_x(&mesh_flow, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY prev_rho_v_y(&mesh_flow, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy prev_e(&mesh_flow, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
  InitialSolutionConcentration prev_c(&mesh_concentration, 0.0);

  Solution rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e, rsln_c;

  // Numerical flux.
  OsherSolomonNumericalFlux num_flux(KAPPA);

  // Initialize weak formulation.
  EulerEquationsWeakFormSemiImplicitCoupled wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM,
    BDY_SOLID_WALL_TOP, BDY_INLET, BDY_OUTLET, BDY_NATURAL_CONCENTRATION, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c, EPSILON, (P_INIT_FLOW == 0 && CAND_LIST_FLOW == H2D_H_ANISO));
  
  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA, RHO_EXT, P_EXT);

  /*
  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));
  ScalarView s5("Concentration", new WinGeom(700, 400, 600, 300));
  */
  
  ScalarView s1("1-coarse", new WinGeom(0, 0, 600, 300));
  ScalarView s2("2-coarse", new WinGeom(700, 0, 600, 300));
  ScalarView s3("3-coarse", new WinGeom(0, 400, 600, 300));
  ScalarView s4("4-coarse", new WinGeom(700, 400, 600, 300));
  ScalarView s5("Concentration", new WinGeom(350, 200, 600, 300));

  ScalarView s11("1-fine", new WinGeom(0, 0, 600, 300));
  ScalarView s21("2-fine", new WinGeom(700, 0, 600, 300));
  ScalarView s31("3-fine", new WinGeom(0, 400, 600, 300));
  ScalarView s41("4-fine", new WinGeom(700, 400, 600, 300));
  ScalarView s51("Concentration-fine", new WinGeom(350, 200, 600, 300));

  ScalarView s111("1-prev", new WinGeom(0, 0, 600, 300));
  ScalarView s211("2-prev", new WinGeom(700, 0, 600, 300));
  ScalarView s311("3-prev", new WinGeom(0, 400, 600, 300));
  ScalarView s411("4-prev", new WinGeom(700, 400, 600, 300));
  ScalarView s511("Concentration-prev", new WinGeom(350, 200, 600, 300));

  // Initialize refinement selector.
  L2ProjBasedSelector l2selector_flow(CAND_LIST_FLOW, CONV_EXP, H2DRS_DEFAULT_ORDER);
  L2ProjBasedSelector l2selector_concentration(CAND_LIST_CONCENTRATION, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  // Set up Advection-Diffusion-Equation stability calculation class.
  ADEStabilityCalculation ADES(ADVECTION_STABILITY_CONSTANT, DIFFUSION_STABILITY_CONSTANT, EPSILON);

  int iteration = 0; double t = 0;
  for(t = 0.0; t < 100.0; t += time_step) {
    info("---- Time step %d, time %3.5f.", iteration++, t);

    // Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0 && (REFINEMENT_COUNT_FLOW > 0 || REFINEMENT_COUNT_CONCENTRATION > 0)) {
      info("Global mesh derefinement.");
      if(REFINEMENT_COUNT_FLOW > 0) {
        REFINEMENT_COUNT_FLOW = 0;
        space_rho.unrefine_all_mesh_elements();
        space_rho_v_x.copy_orders(&space_rho);
        space_rho_v_y.copy_orders(&space_rho);
        space_e.copy_orders(&space_rho);
      }
      if(REFINEMENT_COUNT_CONCENTRATION > 0) {
        REFINEMENT_COUNT_CONCENTRATION = 0;
        space_c.unrefine_all_mesh_elements();
        space_c.adjust_element_order(-1, P_INIT_CONCENTRATION);
      }
    }

    // Adaptivity loop:
    int as = 1; 
    bool done = false;
    do {
      info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      int order_increase = 0;
      if(CAND_LIST_FLOW == H2D_HP_ANISO)
        order_increase = 1;
      Hermes::vector<Space *>* ref_spaces = Space::construct_refined_spaces(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e, &space_c), order_increase);
      if(CAND_LIST_FLOW != H2D_HP_ANISO)
        (*ref_spaces)[4]->adjust_element_order(+1, P_INIT_CONCENTRATION);

      // Project the previous time level solution onto the new fine mesh.
      info("Projecting the previous time level solution onto the new fine mesh.");
      OGProjection::project_global(*ref_spaces, Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c), 
        Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c), matrix_solver);

      if(iteration == 1)
        prev_c.set_const((*ref_spaces)[4]->get_mesh(), 0.0);
      
      s111.show(&prev_rho);
      s211.show(&prev_rho_v_x);
      s311.show(&prev_rho_v_y);
      s411.show(&prev_e);
      s511.show(&prev_c);

      //View::wait();

      if(as > 1) {
        delete rsln_rho.get_mesh();
        delete rsln_rho_v_x.get_mesh();
        delete rsln_rho_v_y.get_mesh();
        delete rsln_e.get_mesh();
      }

      // Report NDOFs.
      info("ndof_coarse: %d, ndof_fine: %d.", 
        Space::get_num_dofs(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e, &space_c)), Space::get_num_dofs(*ref_spaces));

      // Very imporant, set the meshes for the flow as the same.
      (*ref_spaces)[1]->get_mesh()->set_seq((*ref_spaces)[0]->get_mesh()->get_seq());
      (*ref_spaces)[2]->get_mesh()->set_seq((*ref_spaces)[0]->get_mesh()->get_seq());
      (*ref_spaces)[3]->get_mesh()->set_seq((*ref_spaces)[0]->get_mesh()->get_seq());

      // Set up the solver, matrix, and rhs according to the solver selection.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Initialize the FE problem.
      DiscreteProblem dp(&wf, *ref_spaces);

      wf.set_time_step(time_step);

      // Assemble stiffness matrix and rhs.
      info("Assembling the stiffness matrix and right-hand side vector.");
      dp.assemble(matrix, rhs);

      // Solve the matrix problem.
      info("Solving the matrix problem.");
      if (solver->solve())
        Solution::vector_to_solutions(solver->get_solution(), *ref_spaces, 
          Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e, &rsln_c));
      else
      error ("Matrix solver failed.\n");
      
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

      util_time_step = time_step;
      CFL.calculate_semi_implicit(Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), &mesh_flow, util_time_step);

      time_step = util_time_step;

      ADES.calculate(Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y), &mesh_concentration, util_time_step);

      // Calculate element errors and total error estimate.
      info("Calculating error estimates.");
      Adapt adaptivity_flow(Hermes::vector<Space *>(&space_rho, &space_rho_v_x,
      &space_rho_v_y, &space_e), Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      
      double err_est_rel_total_flow = adaptivity_flow.calc_err_est(Hermes::vector<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e),
							  Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e)) * 100;

      Adapt adaptivity_concentration(&space_c, HERMES_L2_NORM);
      
      double err_est_rel_total_concentration = adaptivity_concentration.calc_err_est(&sln_c, &rsln_c) * 100;

      s1.show(&sln_rho);
      s2.show(&sln_rho_v_x);
      s3.show(&sln_rho_v_y);
      s4.show(&sln_e);
      s5.show(&sln_c);

      s11.show(&rsln_rho);
      s21.show(&rsln_rho_v_x);
      s31.show(&rsln_rho_v_y);
      s41.show(&rsln_e);
      s51.show(&rsln_c);

      //s1.wait_for_keypress();

      // Report results.
      info("Error estimate for the flow part: %g%%", err_est_rel_total_flow);

      info("Error estimate for the concentration part: %g%%", err_est_rel_total_concentration);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total_flow < ERR_STOP_FLOW && err_est_rel_total_concentration < ERR_STOP_CONCENTRATION) 
        done = true;
      else {
        info("Adapting coarse meshes.");
        if(err_est_rel_total_flow > ERR_STOP_FLOW) {
          done = adaptivity_flow.adapt(Hermes::vector<RefinementSelectors::Selector *>(&l2selector_flow, &l2selector_flow, &l2selector_flow, &l2selector_flow), 
           THRESHOLD, STRATEGY, MESH_REGULARITY);
          REFINEMENT_COUNT_FLOW++;
        }
        else
          done = true;
        if(err_est_rel_total_concentration > ERR_STOP_CONCENTRATION) {
          if(!adaptivity_concentration.adapt(&l2selector_concentration, THRESHOLD, STRATEGY, MESH_REGULARITY))
            done = false;
          REFINEMENT_COUNT_CONCENTRATION++;
        }

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
      delete solver;
      delete matrix;
      delete rhs;
      for(unsigned int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i];
    }
    while (done == false);

    // Copy the solutions into the previous time level ones.

    prev_rho.copy(&rsln_rho);
    prev_rho_v_x.copy(&rsln_rho_v_x);
    prev_rho_v_y.copy(&rsln_rho_v_y);
    prev_e.copy(&rsln_e);
    prev_c.copy(&rsln_c);
    delete rsln_rho.get_mesh();
    delete rsln_rho_v_x.get_mesh();
    delete rsln_rho_v_y.get_mesh();
    delete rsln_e.get_mesh();
    delete rsln_c.get_mesh();

    // Visualization.
    if((iteration - 1) % EVERY_NTH_STEP == 0) {
      // Hermes visualization.
      if(HERMES_VISUALIZATION) {
        /*
        Mach_number.reinit();
        pressure.reinit();
        entropy.reinit();
        pressure_view.show(&pressure);
        entropy_production_view.show(&entropy);
        Mach_number_view.show(&Mach_number);
        s5.show(&prev_c);
        */
        s1.show(&prev_rho);
        s2.show(&prev_rho_v_x);
        s3.show(&prev_rho_v_y);
        s4.show(&prev_e);
        s5.show(&prev_c);
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
  
  /*
  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();
  s5.close();
  */

  s1.close();
  s2.close();
  s3.close();
  s4.close();
  s5.close();

  return 0;
}
