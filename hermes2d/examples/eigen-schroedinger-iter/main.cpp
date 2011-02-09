#define HERMES_REPORT_ALL
#include "hermes2d.h"
#include <stdio.h>

using namespace RefinementSelectors;

//  This example shows how one can perform adaptivity to a selected eigenfunction
//  without calling the eigensolver again in each adaptivity step. The eigensolver 
//  is only called once at the beginning. 
//
//  PDE: -Laplace u + V*u = lambda_k u,
//  where lambda_0, lambda_1, ... are the eigenvalues.
//
//  Domain: Square (-pi/2, pi/2)^2 ... file domain_square.mesh,
//          L-Shape domain         ... file domain_lshape.mesh.
//
//  BC:  Homogeneous Dirichlet.
//
//  The following parameters can be changed:

int TARGET_EIGENFUNCTION = 1;                     // Desired eigenfunction: 1 for the first, 2 for the second, etc.

int ITERATIVE_METHOD = 2;                         // 1 = Newton, 2 = Picard.

int P_INIT = 2;                                   // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial mesh refinements.
const double THRESHOLD = 0.2;                     // This is a quantitative parameter of the adapt(...) function and
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
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                                  // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double ERR_STOP = 1e-1;                     // Stopping criterion for adaptivity (rel. error tolerance between the
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Pysparse parameters.
const double PYSPARSE_TARGET_VALUE = 2.0;         // PySparse parameter: Eigenvalues in the vicinity of this number will be computed. 
const double PYSPARSE_TOL = 1e-10;                // PySparse parameter: Error tolerance.
const int PYSPARSE_MAX_ITER = 1000;               // PySparse parameter: Maximum number of iterations.

// Parameters for the Newton's and Picard's method.
const double NEWTON_TOL = 1e-3;
const int NEWTON_MAX_ITER = 10;
const double PICARD_TOL = 1e-3;
const int PICARD_MAX_ITER = 50;

// Problem parameters.
double V(double x, double y) {
  return 0;
  //double r = sqrt(x*x + y*y);
  //return -1./(0.001 + r*r);
}

// Boundary markers.
const int BDY_MARKER = 1;

// Weak forms.
#include "forms.cpp"

// Extras.
#include "extras.cpp"

// Main function.
int main(int argc, char* argv[])
{
  info("Desired eigenfunction to calculate: %d.", TARGET_EIGENFUNCTION);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain_lshape.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_MARKER));

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_zero(Hermes::vector<int>(BDY_MARKER));

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);

  // Initialize the weak formulation for the left hand side, i.e., H.
  WeakForm wf_S, wf_M;
  wf_S.add_matrix_form(bilinear_form_S, bilinear_form_S_ord);
  wf_M.add_matrix_form(callback(bilinear_form_M));

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView sview("", new WinGeom(0, 0, 440, 350));
  sview.fix_scale_width(50);
  OrderView oview("", new WinGeom(450, 0, 410, 350));

  // DOF convergence graph.
  SimpleGraph graph_dof;

  // Initialize matrices and matrix solver.
  SparseMatrix* matrix_S = create_matrix(matrix_solver);
  SparseMatrix* matrix_M = create_matrix(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix_S);

  // Assemble matrices S and M on coarse mesh.
  info("Assembling matrices S and M on coarse mesh.");
  bool is_linear = true;
  DiscreteProblem dp_S(&wf_S, &space, is_linear);
  dp_S.assemble(matrix_S);
  DiscreteProblem dp_M(&wf_M, &space, is_linear);
  dp_M.assemble(matrix_M);

  // Write matrix_S in MatrixMarket format.
  write_matrix_mm("mat_S.mtx", matrix_S);

  // Write matrix_M in MatrixMarket format.
  write_matrix_mm("mat_M.mtx", matrix_M);

  // Call Python eigensolver. Solution will be written to "eivecs.dat".
  info("Calling Pysparse.");
  char call_cmd[255];
  sprintf(call_cmd, "python solveGenEigenFromMtx.py mat_S.mtx mat_M.mtx %g %d %g %d", 
	  PYSPARSE_TARGET_VALUE, TARGET_EIGENFUNCTION, PYSPARSE_TOL, PYSPARSE_MAX_ITER);
  system(call_cmd);
  info("Pysparse finished.");

  // Initialize eigenvector.
  double* coeff_vec = new double[ndof];

  // Read solution vectors from file and visualize it.
  double* eigenval =new double[TARGET_EIGENFUNCTION];
  FILE *file = fopen("eivecs.dat", "r");
  char line [64];                  // Maximum line size.
  fgets(line, sizeof line, file);  // ndof
  int n = atoi(line);            
  if (n != ndof) error("Mismatched ndof in the eigensolver output file.");  
  fgets(line, sizeof line, file);  // Number of eigenvectors in the file.
  int neig = atoi(line);
  if (neig != TARGET_EIGENFUNCTION) error("Mismatched number of eigenvectors in the eigensolver output file.");  
  for (int ieig = 0; ieig < neig; ieig++) {
    // Get next eigenvalue from the file
    fgets(line, sizeof line, file);  // eigenval
    eigenval[ieig] = atof(line);            
    // Get the corresponding eigenvector.
    for (int i = 0; i < ndof; i++) {  
      fgets(line, sizeof line, file);
      coeff_vec[i] = atof(line);
    }
    // Normalize the eigenvector.
    normalize((UMFPackMatrix*)matrix_M, coeff_vec, ndof);
  }  
  fclose(file);

  // Eigenvalue.
  double lambda = eigenval[TARGET_EIGENFUNCTION-1];
  info("Eigenvalue on coarse mesh: %g", lambda);
  info("Once more just to check: %g", calc_mass_product((UMFPackMatrix*)matrix_S, coeff_vec, ndof)
      / calc_mass_product((UMFPackMatrix*)matrix_M, coeff_vec, ndof));

  // Convert eigenvector into eigenfunction. After this, the 
  // eigenvector on the coarse mesh will not be needed anymore.
  Solution sln;
  Solution::vector_to_solution(coeff_vec, &space, &sln);
  delete [] coeff_vec;

  // Visualize the eigenfunction.
  info("Plotting initial eigenfunction on coarse mesh.");
  char title[100];
  sprintf(title, "Eigenfunction %d on initial mesh", neig);
  sview.set_title(title);
  sview.show_mesh(false);
  sview.show(&sln);
  sprintf(title, "Initial mesh");
  oview.set_title(title);
  oview.show(&space);
  //View::wait(HERMES_WAIT_KEYPRESS);

  /*** Begin adaptivity ***/

  // Adaptivity loop:
  Solution ref_sln;
  Space* ref_space = NULL;  
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    ref_space = construct_refined_space(&space);
    int ndof_ref = Space::get_num_dofs(ref_space);
    info("ndof: %d, ndof_ref: %d", ndof, ndof_ref);

    // Obtain initial approximation on new reference mesh.
    double* coeff_vec_ref = new double[ndof_ref];
    if (as == 1) {
      // Project the coarse mesh eigenfunction to the reference mesh.
      info("Projecting coarse mesh solution to reference mesh.");
      OGProjection::project_global(ref_space, &sln, coeff_vec_ref, matrix_solver);     
    }
    else {
      // Project the last reference mesh solution to the reference mesh.
      info("Projecting last reference mesh solution to new reference mesh.");
      OGProjection::project_global(ref_space, &ref_sln, coeff_vec_ref, matrix_solver);     
    }
    Solution::vector_to_solution(coeff_vec_ref, ref_space, &ref_sln);      

    // Initialize matrices and matrix solver on reference mesh.
    SparseMatrix* matrix_S_ref = create_matrix(matrix_solver);
    SparseMatrix* matrix_M_ref = create_matrix(matrix_solver);

    // Assemble matrices S and M on reference mesh.
    info("Assembling matrices S and M on reference mesh.");
    is_linear = true;
    DiscreteProblem dp_S_ref(&wf_S, ref_space, is_linear);
    matrix_S_ref->zero();
    dp_S_ref.assemble(matrix_S_ref);
    DiscreteProblem dp_M_ref(&wf_M, ref_space, is_linear);
    matrix_M_ref->zero();
    dp_M_ref.assemble(matrix_M_ref);

    // Calculate eigenvalue corresponding to the new reference solution.
    lambda = calc_mass_product((UMFPackMatrix*)matrix_S_ref, coeff_vec_ref, ndof_ref)
             / calc_mass_product((UMFPackMatrix*)matrix_M_ref, coeff_vec_ref, ndof_ref);
    info("Initial guess for eigenvalue on reference mesh: %.12f", lambda);

    if (ITERATIVE_METHOD == 1) {
      // Newton's method on the reference mesh.
      if(!solve_newton_eigen(ref_space, (UMFPackMatrix*)matrix_S_ref, (UMFPackMatrix*)matrix_M_ref, 
	  		     coeff_vec_ref, lambda, matrix_solver, NEWTON_TOL, NEWTON_MAX_ITER))
        error("Newton's method failed.");
    }
    else {
      // Picard's method on the reference mesh.
      if(!solve_picard_eigen(ref_space, (UMFPackMatrix*)matrix_S_ref, (UMFPackMatrix*)matrix_M_ref, 
	  		     coeff_vec_ref, lambda, matrix_solver, PICARD_TOL, PICARD_MAX_ITER))
        error("Picard's method failed.");
    }
    Solution::vector_to_solution(coeff_vec_ref, ref_space, &ref_sln);

    // Clean up.
    delete matrix_S_ref;
    delete matrix_M_ref;
    delete [] coeff_vec_ref;

    // Project reference solution to coarse mesh for error estimation.
    if (as > 1) {
      // Project reference solution to coarse mesh.
      info("Projecting reference solution to coarse mesh for error calculation.");
      OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 

      // Visualize the projection.
      info("Plotting projection of reference solution to new coarse mesh.");
      char title[100];
      sprintf(title, "Coarse mesh projection");
      sview.set_title(title);
      sview.show_mesh(false);
      sview.show(&sln);
      sprintf(title, "Coarse mesh, step %d", as);
      oview.set_title(title);
      oview.show(&space);
    }

    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    Adapt* adaptivity = new Adapt(&space);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
      Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_est_rel);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      
      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
    ndof = Space::get_num_dofs(&space);
    if (ndof >= NDOF_STOP) done = true;

    // Clean up.
    delete adaptivity;
    //delete ref_space->get_mesh();
    delete ref_space;

    // Wait for keypress.
    //View::wait(HERMES_WAIT_KEYPRESS);
  }
  while (done == false);

  // Wait for all views to be closed.
  info("Computation finished.");
  View::wait();

  return 0; 
};

