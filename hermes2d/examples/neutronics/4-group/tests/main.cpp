#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"

#include "definitions.h"
#include "problem_data.h"

// This test makes sure that example "neutronics/4-group" works correctly.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

const int INIT_REF_NUM = 3;                       // Number of initial uniform mesh refinements.
const int P_INIT_1 = 2,                           // Initial polynomial degree for approximation of group 1 fluxes.
          P_INIT_2 = 2,                           // Initial polynomial degree for approximation of group 2 fluxes.
          P_INIT_3 = 2,                           // Initial polynomial degree for approximation of group 3 fluxes.
          P_INIT_4 = 2;                           // Initial polynomial degree for approximation of group 4 fluxes.
const double ERROR_STOP = 1e-5;                   // Tolerance for the eigenvalue.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h)

// Initial eigenvalue approximation.
double k_eff = 1.0;         

// Element markers.
std::string reflector = "reflector";
std::string core = "core";

// Boundary markers.
std::string bdy_vacuum = "vacuum boundary";
std::string bdy_symmetry = "symmetry plane";

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// General testing class. It can test
//   > if the computed and expected results agree to within a specified tolerance (allowed_delta),
//   > if the computed results do not overshoot the expected results by more than allowed_delta,
//   > if the computed results do not undershoot the expected results by more than allowed_delta.
// When, for the first time, the expected behavior is not achieved, 'passed' is set to false and nothing more is tested.
// If used with a user-defined type, its class must overload all operators appearing in the intended test method.
template<typename T>
struct TestSubject
{
  T allowed_delta;
  bool passed;
  
  TestSubject(T allowed_delta) : allowed_delta(allowed_delta), passed(true) {};
  
  void test_equality(const T& computed, const T& expected) { 
    if(passed) passed = (computed <= expected + allowed_delta && computed >= expected - allowed_delta);
  }
  void test_overshoot(const T& computed, const T& expected) { 
    if(passed) passed = (computed <= expected + allowed_delta);
  }
  void test_undershoot(const T& computed, const T& expected) { 
    if(passed) passed = (computed >= expected - allowed_delta);
  }
};

// Structure that contains information about an extremum. It will be used
// for testing the equality of computed and expected maxima, therefore it
// overloads operators +,-,>=,<=..
struct Extremum
{
  double val, x, y;
  
  Extremum(double val, double x, double y) : val(val), x(x), y(y) {};
  
  inline Extremum operator+ (const Extremum &ex) const {
    return Extremum(val + ex.val, x + ex.x, y + ex.y);
  }
  inline Extremum operator- (const Extremum &ex) const {
    return Extremum(val - ex.val, x - ex.x, y - ex.y);
  }
  inline bool operator>= (const Extremum &ex) const {
    return (val >= ex.val && x >= ex.x && y >= ex.y);
  }
  inline bool operator<= (const Extremum &ex) const {
    return (val <= ex.val && x <= ex.x && y <= ex.y);
  }
  
  std::string str() { 
    char ret[50];
    sprintf(ret, "%lf, at (%lf,%lf)", val, x, y);
    return std::string(ret);
  }
};

// Calculates maximum of a given function, including its coordinates.
Extremum get_peak(MeshFunction *sln)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  Element* e;
  Mesh* mesh = sln->get_mesh();
  
  scalar peak = 0.0;
  double pos_x = 0.0;
  double pos_y = 0.0;
  
  for_all_active_elements(e, mesh)
  {
    update_limit_table(e->get_mode());
    sln->set_active_element(e);
    RefMap* ru = sln->get_refmap();
    int o = sln->get_fn_order() + ru->get_inv_ref_order();
    limit_order(o);
    sln->set_quad_order(o, H2D_FN_VAL);
    scalar *uval = sln->get_fn_values();
    int np = quad->get_num_points(o);
    double* x = ru->get_phys_x(o);
    double* y = ru->get_phys_y(o);
    
    for (int i = 0; i < np; i++)
      if (uval[i] > peak) {
        peak = uval[i];
        pos_x = x[i];
        pos_y = y[i];
      }
  }
  
  return Extremum(peak, pos_x, pos_y);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;
  
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load((std::string("../") + mesh_file).c_str(), &mesh);
  
  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  
  // Solution variables.
  Solution sln1, sln2, sln3, sln4;
  Hermes::vector<Solution*> solutions(&sln1, &sln2, &sln3, &sln4);
  
  // Define initial conditions.
  info("Setting initial conditions.");
  Solution iter1, iter2, iter3, iter4;
  iter1.set_const(&mesh, 1.00);
  iter2.set_const(&mesh, 1.00);
  iter3.set_const(&mesh, 1.00);
  iter4.set_const(&mesh, 1.00);
  Hermes::vector<MeshFunction*> iterates(&iter1, &iter2, &iter3, &iter4);
  
  // Create H1 spaces with default shapesets.
  H1Space space1(&mesh, P_INIT_1);
  H1Space space2(&mesh, P_INIT_2);
  H1Space space3(&mesh, P_INIT_3);
  H1Space space4(&mesh, P_INIT_4);
  Hermes::vector<Space*> spaces(&space1, &space2, &space3, &space4);
  
  int ndof = Space::get_num_dofs(spaces);
  info("ndof = %d.", ndof);
  
  // Load physical data of the problem for the 4 energy groups.
  MaterialPropertyMaps matprop(4);
  matprop.set_D(D);
  matprop.set_Sigma_r(Sr);
  matprop.set_Sigma_s(Ss);
  matprop.set_scattering_multigroup_structure(Ss_nnz);
  matprop.set_fission_multigroup_structure(chi_nnz);
  matprop.set_Sigma_a(Sa);
  matprop.set_Sigma_f(Sf);
  matprop.set_nu(nu);
  matprop.set_chi(chi);
  matprop.validate();
  
  std::cout << matprop;
  
  // Initialize the weak formulation.
  CustomWeakForm wf(matprop, iterates, k_eff, bdy_vacuum);
  
  // Initialize the FE problem.
  DiscreteProblem dp(&wf, spaces);
  
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  
  if (matrix_solver == SOLVER_AZTECOO) 
  {
    ((AztecOOSolver*) solver)->set_solver(iterative_method);
    ((AztecOOSolver*) solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }
  
  // Time measurement.
  TimePeriod cpu_time, solver_time;
  
  // Initial coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[ndof];
  
  // Force the Jacobian assembling in the first iteration.
  bool Jacobian_changed = true;
  
  // In the following iterations, Jacobian will not be changing; its LU factorization
  // may be reused.
  solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);
  
  // Main power iteration loop:
  int it = 1; bool done = false;
  do
  {
    info("------------ Power iteration %d:", it);
    
    info("Newton's method (matrix problem solved by %s).", MatrixSolverNames[matrix_solver].c_str());
    
    memset(coeff_vec, 0.0, ndof*sizeof(scalar)); //TODO: Why it doesn't work without zeroing coeff_vec in each iteration?
    
    solver_time.tick(HERMES_SKIP);      
    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, Jacobian_changed, 1e-8, 10, true)) 
      error("Newton's iteration failed.");
    solver_time.tick();
    
    Solution::vector_to_solutions(solver->get_solution(), spaces, solutions);
    
    // Compute eigenvalue.
    
    SourceFilter source(solutions, matprop);
    SourceFilter source_prev(iterates, matprop);
    
    double k_new = k_eff * (integrate(&source, core) / integrate(&source_prev, core));
    info("Largest eigenvalue: %.8g, rel. difference from previous it.: %g", k_new, fabs((k_eff - k_new) / k_new));
    
    // Stopping criterion.
    if (fabs((k_eff - k_new) / k_new) < ERROR_STOP) done = true;
    
    // Update eigenvalue.
    k_eff = k_new;
    wf.update_keff(k_eff);
    
    if (!done)
    {
      // Save solutions for the next iteration.
      iter1.copy(&sln1);    
      iter2.copy(&sln2);
      iter3.copy(&sln3);    
      iter4.copy(&sln4);
      
      // Don't need to reassemble the system matrix in further iterations,
      // only the rhs changes to reflect the progressively updated source.
      Jacobian_changed = false;
      
      it++;
    }
  }
  while (!done);
  
  delete [] coeff_vec;
  
  // Time measurement.
  cpu_time.tick();
  solver_time.tick(HERMES_SKIP);
  
  // Print timing information.
  verbose("Average solver time for one power iteration: %g s", solver_time.accumulated() / it);
  
  // Clean up.
  delete matrix;
  delete rhs;
  delete solver;

  // Print timing information.
  verbose("Total running time: %g s", cpu_time.accumulated());
  
  // Calculation results for testing.
  info("Number of iterations: %d", it);
  
  // Pointwise results.
  Extremum max1 = get_peak(&sln1); 
  Extremum max2 = get_peak(&sln2);
  Extremum max3 = get_peak(&sln3);
  Extremum max4 = get_peak(&sln4);
  
  info("Peak flux in group 1: %s", max1.str().c_str()); 
  info("Peak flux in group 2: %s", max2.str().c_str());
  info("Peak flux in group 3: %s", max3.str().c_str());
  info("Peak flux in group 4: %s", max4.str().c_str());
  
  // Integral results.
  info("Core eigenvalue: %lf", k_eff);
  
  TestSubject<int> num_iter(2);
  num_iter.test_overshoot(it, 48);
  
  TestSubject<Extremum> peak(Extremum(1e-2, 1e-2, 1e-2));
  peak.test_equality(max1, Extremum(1.020625,1.690169,5.497631));
  peak.test_equality(max2, Extremum(2.272566,1.807669,5.497631));
  peak.test_equality(max3, Extremum(0.337373,1.807669,5.497631));
  peak.test_equality(max4, Extremum(4.603256,1.149095,5.497631));
                                             
  TestSubject<double> eigenvalue(1e-5);      
  eigenvalue.test_equality(k_eff, 1.140726); 
  
  if (num_iter.passed && peak.passed && eigenvalue.passed) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    if (!peak.passed) printf("Peak flux test failed.\n");
    if (!eigenvalue.passed) printf("Eigenvalue test failed.\n");
    if (!num_iter.passed) printf("Number of iterations test failed.\n");
    return ERR_FAILURE;
  }
}
