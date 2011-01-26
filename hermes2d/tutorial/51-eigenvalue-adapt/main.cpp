#define HERMES_REPORT_ALL
#include "hermes2d.h"
#include <stdio.h>

using namespace RefinementSelectors;
using Teuchos::RCP;
using Teuchos::rcp;
using Hermes::EigenSolver;

//  This example uses automatic adaptivity to solve the eigenproblem for the 
//  Laplace operator in a square with zero boundary conditions. Python and 
//  Pysparse must be installed. 
//
//  PDE: -Laplace u = lambda_k u,
//  where lambda_0, lambda_1, ... are the eigenvalues.
//
//  Domain: Square (0, pi)^2.
//
//  BC:  Homogeneous Dirichlet.
//
//  The following parameters can be changed:

const int NUMBER_OF_EIGENVALUES = 6;              // Desired number of eigenvalues. Maximum is 6.
int P_INIT = 2;                                   // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial mesh refinements.
double TARGET_VALUE = 2.0;                        // PySparse parameter: Eigenvalues in the vicinity of 
                                                  // this number will be computed. 
double TOL = 1e-10;                               // Pysparse parameter: Error tolerance.
int MAX_ITER = 1000;                              // PySparse parameter: Maximum number of iterations.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
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
const double CONV_EXP = 0.5;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.001;                    // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;                     // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const int BDY_BOTTOM = 1, BDY_RIGHT = 2, BDY_TOP = 3, BDY_LEFT = 4;

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  if (NUMBER_OF_EIGENVALUES > 6) error("Maximum number of eigenvalues is 6.");
  info("Desired number of eigenvalues: %d.", NUMBER_OF_EIGENVALUES);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers. 
  // Note: "essential" means that solution value is prescribed.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_BOTTOM, BDY_RIGHT, BDY_TOP, BDY_LEFT));

  // Enter Dirichlet boudnary values.
  BCValues bc_values;
  bc_values.add_zero(Hermes::vector<int>(BDY_BOTTOM, BDY_RIGHT, BDY_TOP, BDY_LEFT));

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);

  // Initialize the weak formulation for the left hand side i.e. H 
  WeakForm wf_left, wf_right;
  wf_left.add_matrix_form(callback(bilinear_form_left));
  wf_right.add_matrix_form(callback(bilinear_form_right));

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView sview_1("", new WinGeom(0, 0, 350, 250));
  sview_1.show_mesh(false);
  sview_1.fix_scale_width(60);
  ScalarView sview_2("", new WinGeom(360, 0, 350, 250));
  sview_2.show_mesh(false);
  sview_2.fix_scale_width(60);
  ScalarView sview_3("", new WinGeom(720, 0, 350, 250));
  sview_3.show_mesh(false);
  sview_3.fix_scale_width(60);
  ScalarView sview_4("", new WinGeom(0, 305, 350, 250));
  sview_4.show_mesh(false);
  sview_4.fix_scale_width(60);
  ScalarView sview_5("", new WinGeom(360, 305, 350, 250));
  sview_5.show_mesh(false);
  sview_5.fix_scale_width(60);
  ScalarView sview_6("", new WinGeom(720, 305, 350, 250));
  sview_6.show_mesh(false);
  sview_6.fix_scale_width(60);
  OrderView  oview("Polynomial orders", new WinGeom(1080, 0, 410, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  Solution sln[NUMBER_OF_EIGENVALUES], ref_sln[NUMBER_OF_EIGENVALUES];

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    info("Solving on reference mesh.");

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = construct_refined_space(&space);
    int ref_ndof = Space::get_num_dofs(ref_space);
    info("ref_ndof: %d.", ref_ndof);

    // Initialize matrices and matrix solver on referenc emesh.
    RCP<SparseMatrix> matrix_left = rcp(new CSCMatrix());
    RCP<SparseMatrix> matrix_right = rcp(new CSCMatrix());
    Solver* solver = create_linear_solver(matrix_solver, matrix_left.get());

    // Assemble the matrices on reference mesh.
    bool is_linear = true;
    DiscreteProblem* dp_left = new DiscreteProblem(&wf_left, ref_space, is_linear);
    dp_left->assemble(matrix_left.get());
    DiscreteProblem* dp_right = new DiscreteProblem(&wf_right, ref_space, is_linear);
    dp_right->assemble(matrix_right.get());

    // Time measurement.
    cpu_time.tick();

    // Time measurement.
    cpu_time.tick(HERMES_SKIP);

    EigenSolver es(matrix_left, matrix_right);
    info("Calling Pysparse...");
    es.solve(NUMBER_OF_EIGENVALUES, TARGET_VALUE, TOL, MAX_ITER);
    info("Pysparse finished.");
    es.print_eigenvalues();

    // Initializing solution vector, solution and ScalarView.
    double* ref_coeff_vec;
    //Solution sln[NUMBER_OF_EIGENVALUES], ref_sln[NUMBER_OF_EIGENVALUES];
    //ScalarView view("Solution", new WinGeom(0, 0, 440, 350));

    // Reading solution vectors from file and visualizing.
    double eigenval[NUMBER_OF_EIGENVALUES];
    int neig = es.get_n_eigs();
    if (neig != NUMBER_OF_EIGENVALUES) error("Mismatched number of eigenvectors in the eigensolver output file.");  
    for (int ieig = 0; ieig < NUMBER_OF_EIGENVALUES; ieig++) {
      eigenval[ieig] = es.get_eigenvalue(ieig);
      int n;
      es.get_eigenvector(ieig, &ref_coeff_vec, &n);
      // Convert coefficient vector into a Solution.
      Solution::vector_to_solution(ref_coeff_vec, ref_space, &(ref_sln[ieig]));

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution %d on coarse mesh.", ieig);
      OGProjection::project_global(&space, &(ref_sln[ieig]), &(sln[ieig]), matrix_solver);
    }  

    // FIXME: Below, the adaptivity is done for the last eigenvector only,
    // this needs to be changed to take into account all eigenvectors.

    // View the coarse mesh solution and polynomial orders.
    
    char title[100];
    if (NUMBER_OF_EIGENVALUES > 0) {
      sprintf(title, "Solution 0, val = %g", eigenval[0]);
      sview_1.set_title(title);
      sview_1.show(&(sln[0]));
    }
    if (NUMBER_OF_EIGENVALUES > 1) {
      sprintf(title, "Solution 1, val = %g", eigenval[1]);
      sview_2.set_title(title);
      sview_2.show(&(sln[1]));
    }
    if (NUMBER_OF_EIGENVALUES > 2) {
      sprintf(title, "Solution 2, val = %g", eigenval[2]);
      sview_3.set_title(title);
      sview_3.show(&(sln[2]));
    }
    if (NUMBER_OF_EIGENVALUES > 3) {
      sprintf(title, "Solution 3, val = %g", eigenval[3]);
      sview_4.set_title(title);
      sview_4.show(&(sln[3]));
    }
    if (NUMBER_OF_EIGENVALUES > 4) {
      sprintf(title, "Solution 4, val = %g", eigenval[4]);
      sview_5.set_title(title);
      sview_5.show(&(sln[4]));
    }
    if (NUMBER_OF_EIGENVALUES > 5) {
      sprintf(title, "Solution 5, val = %g", eigenval[5]);
      sview_6.set_title(title);
      sview_6.show(&(sln[5]));
    }
    oview.show(&space);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate.");
    Hermes::vector<Space *> spaces;
    for(int i = 0; i < NUMBER_OF_EIGENVALUES; i++) spaces.push_back(&space);
    Adapt* adaptivity = new Adapt(spaces);
 
    Hermes::vector<Solution *> slns;
    if (NUMBER_OF_EIGENVALUES > 0) slns.push_back(&sln[0]);
    if (NUMBER_OF_EIGENVALUES > 1) slns.push_back(&sln[1]);
    if (NUMBER_OF_EIGENVALUES > 2) slns.push_back(&sln[2]);
    if (NUMBER_OF_EIGENVALUES > 3) slns.push_back(&sln[3]);
    if (NUMBER_OF_EIGENVALUES > 4) slns.push_back(&sln[4]);
    if (NUMBER_OF_EIGENVALUES > 5) slns.push_back(&sln[5]);
    
    Hermes::vector<Solution *> ref_slns;
    if (NUMBER_OF_EIGENVALUES > 0) ref_slns.push_back(&ref_sln[0]);
    if (NUMBER_OF_EIGENVALUES > 1) ref_slns.push_back(&ref_sln[1]);
    if (NUMBER_OF_EIGENVALUES > 2) ref_slns.push_back(&ref_sln[2]);
    if (NUMBER_OF_EIGENVALUES > 3) ref_slns.push_back(&ref_sln[3]);
    if (NUMBER_OF_EIGENVALUES > 4) ref_slns.push_back(&ref_sln[4]);
    if (NUMBER_OF_EIGENVALUES > 5) ref_slns.push_back(&ref_sln[5]);
    Hermes::vector<double> component_errors;
    double err_est_rel = adaptivity->calc_err_est(slns, ref_slns, &component_errors) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d.", Space::get_num_dofs(&space), Space::get_num_dofs(ref_space));
    if (NUMBER_OF_EIGENVALUES > 0) info("err_est_rel[0]: %g%%", component_errors[0] * 100);
    if (NUMBER_OF_EIGENVALUES > 1) info("err_est_rel[1]: %g%%", component_errors[1] * 100);
    if (NUMBER_OF_EIGENVALUES > 2) info("err_est_rel[2]: %g%%", component_errors[2] * 100);
    if (NUMBER_OF_EIGENVALUES > 3) info("err_est_rel[3]: %g%%", component_errors[3] * 100);
    if (NUMBER_OF_EIGENVALUES > 4) info("err_est_rel[4]: %g%%", component_errors[4] * 100);
    if (NUMBER_OF_EIGENVALUES > 5) info("err_est_rel[5]: %g%%", component_errors[5] * 100);
   
    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else
    {
      info("Adapting coarse mesh.");
      Hermes::vector<RefinementSelectors::Selector *> selectors;
      for(int i = 0; i < NUMBER_OF_EIGENVALUES; i++)
        selectors.push_back(&selector);
      done = adaptivity->adapt(selectors, THRESHOLD, STRATEGY, MESH_REGULARITY);

      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
    if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;

    // Clean up.
    delete solver;
    delete adaptivity;
    if(done == false) delete ref_space->get_mesh();
    delete ref_space;
    delete dp_left;
    delete dp_right;
  }
  while (done == false);

  // Wait for all views to be closed.
  View::wait();
  return 0;
};

