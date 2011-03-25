#include "hermes2d.h"

// This test makes sure that example 08-system works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).
const int P_INIT = 8;                                      // Initial polynomial degree of all elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;           // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                           // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const std::string BDY_1 = "1", BDY_3 = "3";

// Problem parameters.
const double E  = 200e9;                                   // Young modulus (steel).
const double nu = 0.3;                                     // Poisson ratio.
const double rho = 8000.0;                                 // Density.
const double g1 = -9.81;                                   // Gravitational acceleration.
const double f0  = 0;                                      // External force in x-direction.
const double f1  = 8e4;                                    // External force in y-direction.

// Weak forms.
#include "../definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);

  // Perform uniform mesh refinement.
  mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst zero_disp(BDY_1, 0.0);
  EssentialBCs bcs(&zero_disp);

  // Create x- and y- displacement space using the default H1 shapeset.
  H1Space u1_space(&mesh, &bcs, P_INIT);
  H1Space u2_space(&mesh, &bcs, P_INIT);
  info("ndof = %d.", Space::get_num_dofs(Hermes::vector<Space *>(&u1_space, &u2_space)));

  // Initialize the weak formulation.
  CustomWeakFormLinearElasticity wf(E, nu, rho*g1, BDY_3, f0, f1);

  // Testing n_dof and correctness of solution vector
  // for p_init = 1, 2, ..., 10
  int success = 1;
  Solution xsln, ysln;
  for (int p_init = 1; p_init <= 10; p_init++) {
    printf("********* p_init = %d *********\n", p_init);
    u1_space.set_uniform_order(p_init);
    u2_space.set_uniform_order(p_init);

    // Initialize the FE problem.
    bool is_linear = true;
    DiscreteProblem dp(&wf, Hermes::vector<Space *>(&u1_space, &u2_space), is_linear);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Initialize the solutions.
    Solution u1_sln, u2_sln;

    // Assemble the stiffness matrix and right-hand side vector.
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    // Solve the linear system and if successful, obtain the solutions.
    info("Solving the matrix problem.");
    if(solver->solve())
      Solution::vector_to_solutions(solver->get_solution(), Hermes::vector<Space *>(&u1_space, &u2_space), 
                                    Hermes::vector<Solution *>(&u1_sln, &u2_sln));
    else
      error ("Matrix solver failed.\n");

    int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&u1_space, &u2_space));
    printf("ndof = %d\n", ndof);
    double sum = 0;
    for (int i=0; i < ndof; i++) sum += solver->get_solution()[i];
    printf("coefficient sum = %g\n", sum);

    // Actual test. The values of 'sum' depend on the
    // current shapeset. If you change the shapeset,
    // you need to correct these numbers.
    if (p_init == 1 && fabs(sum - 3.50185e-06) > 1e-3) success = 0;
    if (p_init == 2 && fabs(sum - 4.34916e-06) > 1e-3) success = 0;
    if (p_init == 3 && fabs(sum - 4.60553e-06) > 1e-3) success = 0;
    if (p_init == 4 && fabs(sum - 4.65616e-06) > 1e-3) success = 0;
    if (p_init == 5 && fabs(sum - 4.62893e-06) > 1e-3) success = 0;
    if (p_init == 6 && fabs(sum - 4.64336e-06) > 1e-3) success = 0;
    if (p_init == 7 && fabs(sum - 4.63724e-06) > 1e-3) success = 0;
    if (p_init == 8 && fabs(sum - 4.64491e-06) > 1e-3) success = 0;
    if (p_init == 9 && fabs(sum - 4.64582e-06) > 1e-3) success = 0;
    if (p_init == 10 && fabs(sum - 4.65028e-06) > 1e-3) success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

