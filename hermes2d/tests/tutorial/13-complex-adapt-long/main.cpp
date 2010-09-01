#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This is a long version of test example "13-complex-adapt-long": function solve_linear_adapt() is not used.
// This test makes sure that example 13-complex-adapt-long works correctly.

const int INIT_REF_NUM = 1;              // Number of initial uniform mesh refinements.
const int P_INIT = 1;                    // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.


// Problem parameters.
double mu_0 = 4.0*3.141592654E-7;
double J_wire = 5000000.0;
double freq = 5E3;
double omega = 2*3.141592654*freq;
double gamma_iron = 6E6;
double mu_iron = 1000*mu_0;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker==1) {return BC_NATURAL;}
  if (marker==2) {return BC_ESSENTIAL;}
  if (marker==3) {return BC_ESSENTIAL;}
  if (marker==4) {return BC_ESSENTIAL;}
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return cplx(0.0,0.0);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain2.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space* space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form_iron), H2D_SYM, 3);
  wf.add_matrix_form(callback(bilinear_form_wire), H2D_SYM, 2);
  wf.add_matrix_form(callback(bilinear_form_air), H2D_SYM, 1);
  wf.add_vector_form(callback(linear_form_wire), 2);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize matrix solver.
  bool is_complex = true;
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, get_num_dofs(space), mat, rhs, solver, is_complex);

  // Adaptivity loop.
  Solution *sln = new Solution();
  Solution *ref_sln = new Solution();
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    info("Solving on reference mesh.");

    // Construct globally refined reference mesh
    // and setup reference space.
    Mesh *ref_mesh = new Mesh();
    ref_mesh->copy(space->get_mesh());
    ref_mesh->refine_all_elements();
    Space* ref_space = space->dup(ref_mesh);
    int order_increase = 1;
    ref_space->copy_orders(space, order_increase);

    // Solve the reference problem.
    // The NULL pointer means that we do not want the resulting coefficient vector.
    solve_linear(ref_space, &wf, matrix_solver, ref_sln, NULL, is_complex);

    // Project the reference solution on the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    // NULL means that we do not want to know the resulting coefficient vector.
    project_global(space, H2D_H1_NORM, ref_sln, sln, NULL, is_complex); 

    // Calculate element errors.
    info("Calculating error (est).");
    Adapt hp(space, H2D_H1_NORM);
    hp.set_solutions(sln, ref_sln);
    hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL);
 
    // Calculate error estimate for each solution component.
    double err_est_abs = calc_abs_error(sln, ref_sln, H2D_H1_NORM);
    double norm_est = calc_norm(ref_sln, H2D_H1_NORM);
    double err_est_rel = err_est_abs / norm_est * 100.;

    // Report results.
    info("ndof: %d, ref_ndof: %d, err_est_rel_total: %g%%", 
         get_num_dofs(space), get_num_dofs(ref_space), err_est_rel);

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else {
      info("Adapting the coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (get_num_dofs(space) >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (done == false);

  int ndof = get_num_dofs(space);

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  printf("ndof allowed = %d\n", 650);
  printf("ndof actual = %d\n", ndof);
  if (ndof < 650) {      // ndofs was 625 atthe time this test was created
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
