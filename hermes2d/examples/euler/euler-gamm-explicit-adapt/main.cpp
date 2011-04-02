#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This example solves the compressible Euler equations using FVM and automatic h-adaptivity.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: GAMM channel, see mesh file GAMM-channel.mesh
//
// BC: Normal velocity component is zero on solid walls.
//     Subsonic state prescribed on inlet and outlet.
//
// IC: Constant subsonic state identical to inlet. 
//
// The following parameters can be changed:
// Shock capturing.
bool SHOCK_CAPTURING = true;
// Quantitative parameter of the discontinuity detector.
double DISCONTINUITY_DETECTOR_PARAM = 0.05;

const int P_INIT = 0;                                   // Initial polynomial degree.                      
const int INIT_REF_NUM = 2;                             // Number of initial uniform mesh refinements.                       
double CFL = 0.8;                                 // CFL value.
double time_step = 1E-4;                                // Time step.

// Adaptivity.
const int UNREF_FREQ = 10;                        // Every UNREF_FREQth time step the mesh is unrefined.
int REFINEMENT_COUNT = 0;                         // Number of mesh refinements between two unrefinements.
                                                  // The mesh is not unrefined unless there has been a refinement since
                                                  // last unrefinement.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;           // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
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
const double ERR_STOP = 0.5;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;                     // Adaptivity process stops when the number of degrees of freedom grows over
                                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Equation parameters.
const double P_EXT = 2.5;                               // Exterior pressure (dimensionless).
const double RHO_EXT = 1.0;                             // Inlet density (dimensionless).   
const double V1_EXT = 1.25;                             // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.0;                              // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;                               // Kappa.

// Boundary markers.
const std::string BDY_INLET = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_BOTTOM = "3";
const std::string BDY_SOLID_WALL_TOP = "4";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../constant_initial_condition.cpp"

// Refinement criterion function.
int criterion(Element * e)
{
  if(e->vn[0]->x == 0.5 && e->vn[0]->y == 0)
    return 0;
  if(e->vn[1]->x == 0.5 && e->vn[1]->y == 0)
    return 0;
  if(e->vn[0]->x == 1.5 && e->vn[0]->y == 0)
    return 0;
  if(e->vn[1]->x == 1.5 && e->vn[1]->y == 0)
    return 0;
  
  return -1;
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh basemesh;
  Mesh mesh;
  H2DReader mloader;
  mloader.load("GAMM-channel.mesh", &basemesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  basemesh.refine_towards_boundary(BDY_SOLID_WALL_BOTTOM, 2);
  basemesh.refine_towards_boundary(BDY_SOLID_WALL_TOP, 2);
  basemesh.refine_towards_boundary(BDY_INLET, 2);
  basemesh.refine_towards_boundary(BDY_OUTLET, 2);
  basemesh.refine_by_criterion(criterion, 2);
  mesh.copy(&basemesh);

  // Initialize boundary condition types and spaces with default shapesets.
  L2Space space_rho(&mesh, P_INIT);
  L2Space space_rho_v_x(&mesh, P_INIT);
  L2Space space_rho_v_y(&mesh, P_INIT);
  L2Space space_e(&mesh, P_INIT);
  int ndof = Space::get_num_dofs(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  InitialSolutionEulerDensity sln_rho(&mesh, RHO_EXT);
  InitialSolutionEulerDensityVelX sln_rho_v_x(&mesh, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY sln_rho_v_y(&mesh, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy sln_e(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
  
  InitialSolutionEulerDensity prev_rho(&mesh, RHO_EXT);
  InitialSolutionEulerDensityVelX prev_rho_v_x(&mesh, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY prev_rho_v_y(&mesh, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy prev_e(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  Solution rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e;
  
  // Numerical flux.
  OsherSolomonNumericalFlux num_flux(KAPPA);

  // Initialize weak formulation.
  EulerEquationsWeakFormExplicit wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP, 
    BDY_INLET, BDY_OUTLET, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);

  // Initialize the FE problem.
  bool is_linear = true;
  
  DiscreteProblem dp(&wf, Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), is_linear);

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA, RHO_EXT, P_EXT);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));

  /*
  ScalarView s1("1", new WinGeom(0, 0, 600, 300));
  ScalarView s2("2", new WinGeom(700, 0, 600, 300));
  ScalarView s3("3", new WinGeom(0, 400, 600, 300));
  ScalarView s4("4", new WinGeom(700, 400, 600, 300));
  */

  // Initialize refinement selector.
  L2ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  int iteration = 0; double t = 0;
  for(t = 0.0; t < 3.0; t += time_step) {
    info("---- Time step %d, time %3.5f.", iteration++, t);

    // Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0 && REFINEMENT_COUNT > 0) {
      REFINEMENT_COUNT = 0;
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      mesh.copy(&basemesh);
      space_rho.set_uniform_order(P_INIT);
      space_rho_v_x.set_uniform_order(P_INIT);
      space_rho_v_y.set_uniform_order(P_INIT);
      space_e.set_uniform_order(P_INIT);
    }
    
    // Adaptivity loop:
    int as = 1; 
    bool done = false;
    do {
      info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      // Global polynomial order increase = 0;
      int order_increase = 0;
      Hermes::vector<Space *>* ref_spaces = Space::construct_refined_spaces(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), order_increase);

      // Project the previous time level solution onto the new fine mesh.
      info("Projecting the previous time level solution onto the new fine mesh.");
      OGProjection::project_global(*ref_spaces, Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 
                     Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), matrix_solver); 

      // Assemble the reference problem.
      info("Solving on reference mesh.");
      bool is_linear = true;
      DiscreteProblem* dp = new DiscreteProblem(&wf, *ref_spaces, is_linear);
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      dp->assemble(matrix, rhs);

      // Solve the linear system of the reference problem. If successful, obtain the solutions.
      scalar* solution_vector;
      if(solver->solve()) {
        solution_vector = solver->get_solution();
        Solution::vector_to_solutions(solver->get_solution(), *ref_spaces, 
                                              Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));
      }
      else 
        error ("Matrix solver failed.\n");

      if(SHOCK_CAPTURING) {
        DiscontinuityDetector discontinuity_detector(*ref_spaces, Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));

        std::set<int> discontinuous_elements = discontinuity_detector.get_discontinuous_element_ids(DISCONTINUITY_DETECTOR_PARAM);

        FluxLimiter flux_limiter(solution_vector, *ref_spaces, Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));

        flux_limiter.limit_according_to_detector(discontinuous_elements);
      }

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection::project_global(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), 
                     Hermes::vector<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), matrix_solver, 
                     Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e),
							  Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e)) * 100;

      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
        Space::get_num_dofs(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e)), Space::get_num_dofs(*ref_spaces), err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) 
        done = true;
      else 
      {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *>(&selector, &selector, &selector, &selector), 
                                 THRESHOLD, STRATEGY, MESH_REGULARITY);

        REFINEMENT_COUNT++;
        if (Space::get_num_dofs(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e)) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
      delete adaptivity;
      if(!done)
        for(unsigned int i = 0; i < ref_spaces->size(); i++)
          delete (*ref_spaces)[i]->get_mesh();
      for(unsigned int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i];
      delete dp;
    }
    while (done == false);

    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&rsln_rho);
    prev_rho_v_x.copy(&rsln_rho_v_x);
    prev_rho_v_y.copy(&rsln_rho_v_y);
    prev_e.copy(&rsln_e);

    // Visualization.
    Mach_number.reinit();
    pressure.reinit();
    entropy.reinit();
    pressure_view.show(&pressure);
    entropy_production_view.show(&entropy);
    Mach_number_view.show(&Mach_number);

    /*
    s1.show(&prev_rho);
    s2.show(&prev_rho_v_x);
    s3.show(&prev_rho_v_y);
    s4.show(&prev_e);
    */

    delete rsln_rho.get_mesh(); 
    delete rsln_rho_v_x.get_mesh(); 
    delete rsln_rho_v_y.get_mesh();
    delete rsln_e.get_mesh();
  }
  
  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();

  /*
  s1.close();
  s2.close();
  s3.close();
  s4.close();
  */
  
  return 0;
}