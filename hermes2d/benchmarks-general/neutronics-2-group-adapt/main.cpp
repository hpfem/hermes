#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"

#include "hermes2d.h"

using namespace RefinementSelectors;

// This benchmark uses automatic adaptivity to solve a 2-group neutron diffusion equation with a known exact solution.
// The solution reflects the typical behavior observed in real cases, where one component is very smooth and the
// other more oscillating. Typical boundary conditions prescribed in real models have also been chosen.
//
// Author: Milan Hanus (University of West Bohemia, Pilsen, Czech Republic).
//
// EQUATION:
//
//    L_1 = Q_1
//    L_2 = Q_2
//
// where
//
//    L_1 = - \nabla \cdot D_1 \nabla \phi_1  
//          + \SigmaR_1 \phi_1  
//          - \SigmaS_21 \phi_2 
//          - \chi_1 ( \nu\SigmaF_1 \phi_1 + \nu\SigmaF_2 \phi_2 ) 
//    L_2 = - \nabla \cdot D_2 \nabla \phi_2  
//          + \SigmaR_2 \phi_2  
//          - \SigmaS_12 \phi_1 
//          - \chi_2 ( \nu\SigmaF_1 \phi_1 + \nu\SigmaF_2 \phi_2 ) 
//
// BC:
//
//  Homogeneous Neumann on symmetry axes (BDY_SYMMETRY),
//  homogeneous Dirichlet on zero flux boundary (BDY_FLUX),
//  homogeneous Newton on albedo boundary (BDY_GAMMA):
//    -d D_1\phi_1 / d n = GAMMA_1 \phi_1,
//    -d D_2\phi_2 / d n = GAMMA_2 \phi_2,
//    GAMMA_1 = 8*D_1, GAMMA_2 = 8*D_2.
//

// INITIALIZATION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adaptivity control:

const int P_INIT[2] =
  {1, 1};                                         // Initial polynomial orders for the individual solution components.
const int INIT_REF_NUM[2] =
  {1, 1};                                         // Initial uniform mesh refinement for the individual solution components.
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const bool MULTIMESH = true;                      // true = use multi-mesh, false = use single-mesh.
                                                  // Note: in the single mesh option, the meshes are forced to be geometrically
                                                  // the same but the polynomial degrees can still vary.
const double THRESHOLD_MULTI = 0.3;               // error threshold for element refinement (multi-mesh)
const double THRESHOLD_SINGLE = 0.7;              // error threshold for element refinement (single-mesh)                                         
const CandList CAND_LIST = H2D_HP_ANISO_H;        // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.01;                     // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference and coarse mesh solution in percent).
const int NDOF_STOP = 200000;                     // Adaptivity process stops when the number of degrees of freedom grows over
                                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPT_NUM = 60;                     // Adaptivity process stops when the number of adaptation steps grows over
                                                  // this limit.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Variables used for reporting of results
const int ERR_PLOT = 0;                           // Row in the convergence graphs for exact errors.
const int ERR_EST_PLOT = 1;                       // Row in the convergence graphs for error estimates.
const int GROUP_1 = 0;                            // Row in the DOF evolution graph for group 1.
const int GROUP_2 = 1;                            // Row in the DOF evolution graph for group 2.

// Weak forms, input data and some other utility functions.
#include "definitions.cpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;
  
  // Load the mesh.
  Mesh mesh1, mesh2;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh1);

  if (MULTIMESH) 
  {
    // Obtain meshes for the 2nd group by cloning the mesh loaded for the 1st group.
    mesh2.copy(&mesh1);
    
    // Initial uniform refinements.
    for (int i = 0; i < INIT_REF_NUM[0]; i++) mesh1.refine_all_elements();
    for (int i = 0; i < INIT_REF_NUM[1]; i++) mesh2.refine_all_elements();
  } 
  else // Use just one mesh for both groups.
    for (int i = 0; i < INIT_REF_NUM[0]; i++) mesh1.refine_all_elements();


  // Create H1 space with default shapesets.
  H1Space space1(&mesh1, P_INIT[0]);
  H1Space space2(MULTIMESH ? &mesh2 : &mesh1, P_INIT[1]);

  // Load physical data of the problem for the 4 energy groups.
  MaterialPropertyMaps matprop(2, std::set<std::string>(regions, regions+4));
  matprop.set_D(D);
  matprop.set_Sigma_r(Sr);
  matprop.set_Sigma_s(Ss);
  matprop.set_Sigma_s_nnz_structure(Ss_nnz);
  matprop.set_nuSigma_f(nSf);
  matprop.set_nu(nu);
  matprop.set_chi(chi);
  matprop.validate();
  
  std::cout << std::endl << matprop << std::endl;
  
  // Initialize the weak formulation.  
  const double a = 0, b = 1;
  CustomWeakForm wf( 
    matprop, 
    Hermes::vector<DefaultFunction*>(
      new CustomRightHandSide_g1(a, b, matprop, regions), 
      new CustomRightHandSide_g2(a, b, matprop, regions)
    ),
    "gamma", 8.
  );
  
  // Initialize coarse and reference mesh solutions and pointers to them.
  Solution sln1, sln2; 
  Solution ref_sln1, ref_sln2;
  Hermes::vector<Solution*> slns(&sln1, &sln2);
  Hermes::vector<Solution*> ref_slns(&ref_sln1, &ref_sln2);
  
  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
  //selector.set_error_weights(2.1, 0.9, sqrt(2.0));

  // Initialize views.
  ScalarView view1("Neutron flux 1", new WinGeom(0, 0, 400, 350));
  ScalarView view2("Neutron flux 2", new WinGeom(410, 0, 400, 350));
  OrderView oview1("Mesh and orders for group 1", new WinGeom(820, 0, 350, 300));
  OrderView oview2("Mesh and orders for group 2", new WinGeom(1180, 0, 350, 300));
  ScalarView view3("Error in neutron flux 1", new WinGeom(0, 405, 400, 350));
  ScalarView view4("Error in neutron flux 2", new WinGeom(410, 405, 400, 350));

  // Show meshes.
  view1.show_mesh(false); view1.set_3d_mode(true);
  view2.show_mesh(false); view2.set_3d_mode(true);
  view3.show_mesh(false); view3.set_3d_mode(true);
  view4.show_mesh(false); view4.set_3d_mode(true);

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof("Error convergence", "Degrees of freedom", "Error [%]");
  graph_dof.add_row("exact error (H1)", "b", "-", "o");
  graph_dof.add_row("est.  error (H1)", "r", "-", "s");
  graph_dof.set_log_y();
  graph_dof.show_legend();
  graph_dof.show_grid();

  SimpleGraph graph_cpu("Error convergence", "CPU time [s]", "Error [%]");
  graph_cpu.add_row("exact error (H1)", "b", "-", "o");
  graph_cpu.add_row("est.  error (H1)", "r", "-", "s");
  graph_cpu.set_log_y();
  graph_cpu.show_legend();
  graph_cpu.show_grid();
  
  PNGGraph graph_dof_evol("Evolution of NDOF", "Adaptation step", "Num. DOF");
  graph_dof_evol.add_row("group 1", "b", "-", "o");
  graph_dof_evol.add_row("group 2", "r", "-", "s");
  graph_dof_evol.set_log_y();
  graph_dof_evol.show_legend();
  graph_dof_evol.set_legend_pos("bottom right");
  graph_dof_evol.show_grid();
    
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  do
  {
    int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&space1, &space2));

    info("!---- Adaptivity step %d ---------------------------------------------", as);

    // Construct globally refined reference meshes.
    Mesh ref_mesh1, ref_mesh2;
    ref_mesh1.copy(&mesh1);
    ref_mesh1.refine_all_elements();
    if (MULTIMESH) {
      ref_mesh2.copy(&mesh2);
      ref_mesh2.refine_all_elements();
    }
    
    // Setup spaces for the reference solution.
    int order_increase = 1;
    Space *ref_space1 = space1.dup(&ref_mesh1, order_increase);
    Space *ref_space2 = space2.dup(MULTIMESH ? &ref_mesh2 : &ref_mesh1, order_increase);
    
    int ref_ndof = Space::get_num_dofs(Hermes::vector<Space *>(ref_space1, ref_space2));
    info("------------------ Reference solution; NDOF=%d -------------------", ref_ndof);
                            
    if (ref_ndof >= NDOF_STOP) break;

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Assemble the reference problem.
    info("Solving on reference mesh.");
    DiscreteProblem* dp = new DiscreteProblem(&wf, Hermes::vector<Space*>(ref_space1, ref_space2));
    
    // Time measurement.
    cpu_time.tick();
    
    // Initial coefficient vector for the Newton's method.  
    scalar* coeff_vec = new scalar[ref_ndof];
    memset(coeff_vec, 0, ref_ndof * sizeof(scalar));

    // Perform Newton's iteration.
    if (!hermes2d.solve_newton(coeff_vec, dp, solver, matrix, rhs, true, 1e-8, 10, true)) 
      error("Newton's iteration failed.");
    delete dp;
    
    Solution::vector_to_solutions(coeff_vec, Hermes::vector<Space*>(ref_space1, ref_space2), ref_slns);
    delete ref_space1;
    delete ref_space2;
    
    // Time measurement.
    cpu_time.tick();
        
    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(Hermes::vector<Space*>(&space1, &space2), ref_slns, slns, matrix_solver);
    
    // View the distribution of polynomial orders on the coarse meshes.
    info("flux1_dof=%d, flux2_dof=%d", Space::get_num_dofs(&space1), Space::get_num_dofs(&space2));
    oview1.show(&space1);
    oview2.show(&space2);                   
    
    // Calculate element errors and total error estimate.
    info("Calculating error estimate and exact error.");
    Adapt adaptivity(Hermes::vector<Space*>(&space1, &space2));

    double err_est_energ_total = adaptivity.calc_err_est(slns, ref_slns) * 100;
    
    Adapt adaptivity_proj(Hermes::vector<Space*>(&space1, &space2));

    double err_est_h1_total = adaptivity_proj.calc_err_est(slns, ref_slns, NULL, false) * 100;

    Hermes::vector<double> err_exact_h1;
    CustomExactSolution_g1 ex1(&mesh1);
    CustomExactSolution_g2 ex2(MULTIMESH ? &mesh2 : &mesh1);
    Hermes::vector<Solution*> exslns(&ex1, &ex2); 
    double error_h1 = adaptivity_proj.calc_err_exact(slns, exslns, &err_exact_h1, false) * 100;
    
    // Report results.
    cpu_time.tick(); 

    // Error w.r.t. the exact solution.
    DiffFilter err_distrib_1(Hermes::vector<MeshFunction*>(&ex1, &sln1));
    DiffFilter err_distrib_2(Hermes::vector<MeshFunction*>(&ex2, &sln2));

    info("Per-component error wrt. exact solution (H1 norm): %g%%, %g%%", 
         err_exact_h1[0] * 100, err_exact_h1[1] * 100);
    info("Total error wrt. exact solution (H1 norm): %g%%", error_h1);
    info("Total error wrt. ref. solution  (H1 norm): %g%%", err_est_h1_total);
    info("Total error wrt. ref. solution  (E norm):  %g%%", err_est_energ_total);

    view1.show(&sln1);
    view2.show(&sln2);
    //view3.show(&ex1);
    //view4.show(&ex2);
    view3.show(&err_distrib_1);
    view4.show(&err_distrib_2);

    if (ndof > 100) {
      // Add entry to DOF convergence graphs.
      graph_dof.add_values(ERR_PLOT, ndof, error_h1);
      graph_dof.add_values(ERR_EST_PLOT, ndof, err_est_h1_total);
      // Add entry to DOF evolution graphs.
      graph_dof_evol.add_values(GROUP_1, as, Space::get_num_dofs(&space1));
      graph_dof_evol.add_values(GROUP_2, as, Space::get_num_dofs(&space2));
      // Add entry to CPU convergence graphs.
      graph_cpu.add_values(ERR_PLOT, cpu_time.accumulated(), error_h1);
      graph_cpu.add_values(ERR_EST_PLOT, cpu_time.accumulated(), err_est_h1_total);
    }

    cpu_time.tick(HERMES_SKIP);
    
    // If err_est_energ_total too large, adapt the mesh.
    if (err_est_energ_total < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity.adapt( Hermes::vector<RefinementSelectors::Selector*>(&selector,&selector), 
                               MULTIMESH ? THRESHOLD_MULTI : THRESHOLD_SINGLE, STRATEGY, MESH_REGULARITY );
      
      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
      else
      {
        // Show the reference solution - the final result.
        view1.set_title("Fine mesh solution");
        view1.show_mesh(false);
        view1.show(&ref_sln1);
        view2.set_title("Fine mesh solution");
        view2.show_mesh(false);
        view2.show(&ref_sln2);
      }
    }
    
    // Clean up.
    delete [] coeff_vec;
    delete solver;
    delete matrix;
    delete rhs;
  }
  while (done == false);

  cpu_time.tick();
  verbose("Total running time: %g s", cpu_time.accumulated());

  ////////////////////////  Save plots with results.  //////////////////////////
  
  std::stringstream str;
  make_str_from_adapt_opts(str);
    
  // Save plots of final distribution of polynomial orders over each mesh.
  
  std::stringstream o1;
  o1 << "mesh_" << str.str();
  if (MULTIMESH) {
    o1 << "-1.bmp";
    std::stringstream o2;
    o2 << "mesh_" << str.str() << "-2.bmp";
    oview2.save_screenshot(o2.str().c_str(), true);
  } else {
    o1 << ".bmp";
  }
  oview1.save_screenshot(o1.str().c_str(), true);
  
  
  // Save convergence graphs.
  
  std::stringstream ccfile, cdfile;
  cdfile << "conv_dof_" << str.str() << ".dat";
  ccfile << "conv_cpu_" << str.str() << ".dat";
    
  graph_dof_evol.save("dof_evol.gp");
  graph_dof.save(cdfile.str().c_str());
  graph_cpu.save(ccfile.str().c_str());  
  
  // Wait for all views to be closed.
  View::wait();
  return 0;
};
