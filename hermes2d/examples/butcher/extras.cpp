#include "hermes2d.h"
#include "../../src/weakform.h"
#include <string>

bool HERMES_RESIDUAL_AS_VECTOR = false;
bool rk_time_step(ButcherTable* bt, double time_step,
                  scalar* coeff_vec, DiscreteProblem* dp, MatrixSolverType matrix_solver, 
                  double newton_tol, int newton_max_iter, bool verbose = true, 
                  double newton_damping_coeff = 1.0, double newton_max_allowed_residual_norm = 1e6)
{
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Get number of stages.
  int num_stages = bt->get_size();

  // Space.
  Space* space = dp->get_space(0);

  // Get ndof.
  int ndof = space->get_num_dofs();

  // Extract mesh from space.
  Mesh* mesh = space->get_mesh();

  // Create num_stages spaces for stage solutions.
  Hermes::Tuple<Space*> stage_spaces;
  stage_spaces.push_back(dp->get_space(0));
  for (int i = 1; i < num_stages; i++) {
    stage_spaces.push_back(space->dup(mesh));
    stage_spaces[i]->copy_orders(space);
  }

  // Initialize stage solutions: One for each stage.
  Hermes::Tuple<MeshFunction*> stage_solutions;
  for (int i=0; i < num_stages; i++) {
    Solution* stage_sln = new Solution(mesh);
    stage_sln->set_zero(mesh);
    stage_solutions.push_back(stage_sln);
  }

  // Extract the weak formulation from the original DiscreteProblem.
  WeakForm* wf = dp->get_weak_formulation();
  if (wf->get_neq() != 1) error("wf->neq != 1 not implemented yet.");
  Hermes::Tuple<WeakForm::MatrixFormVol> mfvol = wf->get_mfvol();
  Hermes::Tuple<WeakForm::MatrixFormSurf> mfsurf = wf->get_mfsurf();
  Hermes::Tuple<WeakForm::VectorFormVol> vfvol = wf->get_vfvol();
  Hermes::Tuple<WeakForm::VectorFormSurf> vfsurf = wf->get_vfsurf();

  // TODO: these weak forms need to be duplicated and antered as
  // blocks of the stage weah formulation.






  // Initialize the stage weak formulation.
  WeakForm stage_wf(num_stages);
  for (int i=0; i < num_stages; i++) 
    for (int j=0; j < num_stages; j++) 
      stage_wf.add_matrix_form(i, j, callback(jac), HERMES_NONSYM, HERMES_ANY, stage_solutions[i]);
  for (int i=0; i < num_stages; i++) 
    stage_wf.add_vector_form(i, callback(res), HERMES_ANY, stage_solutions[i]);

  // Create a new DiscreteProblem for the stage slutions.
  bool is_linear = dp->get_is_linear();
  DiscreteProblem stage_dp(&stage_wf, stage_spaces, is_linear);

  // Stage vector of length num_stages * ndof, initialize with zeros.
  scalar* stage_vec = new scalar[num_stages*ndof];
  memset(stage_vec, 0, num_stages * ndof * sizeof(scalar));

  // Helper vector.
  scalar* vec = new scalar[ndof];

  // The Newton's loop.
  double residual_norm;
  int it = 1;
  while (true)
  {
    // Prepare external solution for each stage.
    for (int r = 0; r < num_stages; r++) {
      memset(vec, 0, ndof * sizeof(scalar));
      double increment;
      for (int i = 0; i < ndof; i++) {
        increment = 0;
        for (int s = 0; s < num_stages; s++) {
          increment += bt->get_A(r, s) * stage_vec[s*ndof + i]; 
        }
        vec[i] = coeff_vec[i] + time_step * increment;
      }
      Solution::vector_to_solution(vec, space, (Solution*)stage_solutions[r]);
    } 

    // Calculating weight coefficients for blocks in the 
    // Stage jacobian matrix. 
    Table block_weights(num_stages);
    for (int r = 0; r < num_stages; r++) {
      for (int s = 0; s < num_stages; s++) {
        block_weights.set_A(r, s, bt->get_A(r, s) * time_step);
      }
    }

    // Assemble the stage Jacobian matrix and residual vector.
    // Blocks that would be zeroed will not be assembled, and 
    // all assembled blocks will be weighted according to the 
    // Butcher's table.
    bool rhs_only = false;
    dp->assemble(NULL, matrix, rhs, rhs_only, &block_weights);

    // Add -1 to each diagonal element of the matrix.
    matrix->add_to_diagonal(-1);

    // Subtract stage_vec from rhs.
    for (int i = 0; i < num_stages*ndof; i++) rhs->add(i, -stage_vec[i]);

    // Multiply the residual vector with -1 since the matrix 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    rhs->change_sign();
    
    // Measure the residual norm.
    if (HERMES_RESIDUAL_AS_VECTOR) {
      // Calculate the l2-norm of residual vector.
      residual_norm = get_l2_norm(rhs);
    }
    else {
      // Translate the residual vector into a residual function (or multiple functions) 
      // in the corresponding finite element space(s) and measure their norm(s) there.
      // This is more meaningful since not all components in the coefficient vector 
      // have the same weight when translated into the finite element space.
      Hermes::Tuple<Solution*> residuals;
      Hermes::Tuple<bool> add_dir_lift;
      for (int i = 0; i < num_stages; i++) {
        residuals.push_back((Solution*)stage_solutions[i]);
        add_dir_lift.push_back(false);
      }
      Solution::vector_to_solutions(rhs, dp->get_spaces(), residuals, add_dir_lift);
      residual_norm = calc_norms(residuals);
    }

    // Info for the user.
    if (verbose) info("---- Newton iter %d, ndof %d, residual norm %g", it, ndof, residual_norm);

    // If maximum allowed residual norm is exceeded, fail.
    if (residual_norm > max_allowed_residual_norm) {
      if (verbose) {
        info("Current residual norm: %g", residual_norm);
        info("Maximum allowed residual norm: %g", max_allowed_residual_norm);
        info("Newton solve not successful, returning false.");
      }
      return false;
    }

    // If residual norm is within tolerance, or the maximum number 
    // of iteration has been reached, then quit.
    if ((residual_norm < newton_tol || it > newton_max_iter) && it > 1) break;

    // Solve the linear system.
    if(!solver->solve()) error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < num_stages*ndof; i++) stage_vec[i] += damping_coeff * solver->get_solution()[i];

    it++;
  }

  // If max number of iterations was exceeded, fail. 
  if (it >= newton_max_iter) {
    if (verbose) info("Maximum allowed number of Newton iterations exceeded, returning false.");
    // Delete helper vector.
    delete [] vec;
    return false;
  }

  // Calculate new coefficient vector using the stage vector and the Butcher's table.
  for (int i = 0; i < ndof; i++) {
    double increment = 0;
    for (int s = 0; s < num_stages; s++) {
      increment += bt->get_B(s) * stage_vec[s*ndof + i]; 
    }
    coeff_vec[i] += time_step * increment;
  } 

  // Clean up.
  delete [] vec;
  delete matrix;
  delete rhs;
  delete solver;
  for (int i=0; i<num_stages; i++) delete stage_spaces[i];
  
  return true;
}
