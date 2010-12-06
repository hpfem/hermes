// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "common.h"
#include "limit_order.h"
#include "discrete_problem.h"
#include "linear_problem.h"
#include "weakform.h"
#include "solver.h"
#include "space/space.h"
#include "precalc.h"
#include "refmap.h"
#include "solution.h"
#include "integrals_h1.h"
#include "views/view.h"
#include "views/vector_view.h"
#include "tuple.h"
#include "norm.h"


LinearProblem::LinearProblem() : DiscreteProblem() {};
LinearProblem::LinearProblem(WeakForm* wf_) : DiscreteProblem(wf_) {};
LinearProblem::LinearProblem(WeakForm* wf_, Space* s_) : DiscreteProblem(wf_, s_) {};
LinearProblem::LinearProblem(WeakForm* wf_, Hermes::Tuple<Space*> spaces_) : DiscreteProblem(wf_, spaces_) {};
LinearProblem::~LinearProblem() {};

void LinearProblem::assemble(Matrix* mat_ext, Vector* rhs_ext, bool rhsonly, bool is_complex)
{
  int ndof = this->get_num_dofs();
  if (ndof == 0) error("ndof == 0 in LinearProblem::assemble().");
  Vector* dir_ext = new AVector(ndof, is_complex);
  // The vector dir represents the contribution of the Dirichlet lift, 
  // and for linear problems it has to be subtracted from the right hand side.
  // The NULL stands for the initial coefficient vector that is not used.
  DiscreteProblem::assemble(NULL, mat_ext, dir_ext, rhs_ext, rhsonly);
  // FIXME: Do we really need to handle the real and complex cases separately?
  if (is_complex) for (int i=0; i < ndof; i++) rhs_ext->add(i, -dir_ext->get_cplx(i));
  else for (int i=0; i < ndof; i++) rhs_ext->add(i, -dir_ext->get(i));
  delete dir_ext;
}

// Solve a typical linear problem (without automatic adaptivity).
// Feel free to adjust this function for more advanced applications.
bool solve_linear(Hermes::Tuple<Space *> spaces, WeakForm* wf, MatrixSolverType matrix_solver, 
                  Hermes::Tuple<Solution *> solutions, Vector* coeff_vec, bool is_complex) 
{
  // Initialize the linear problem.
  LinearProblem lp(wf, spaces);
  int ndof = get_num_dofs(spaces);
  //info("ndof = %d", ndof);

  // Select matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;
  init_matrix_solver(matrix_solver, ndof, mat, rhs, solver, is_complex);

  // Assemble stiffness matrix and rhs.
  bool rhsonly = false;
  lp.assemble(mat, rhs, rhsonly, is_complex);

  //mat->print();

  // Solve the matrix problem.
  if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

  // Convert coefficient vector into a Solution.
  for (int i=0; i < solutions.size(); i++) {
    solutions[i]->set_coeff_vector(spaces[i], rhs);
  }
	
  // Copy the coefficient vector into coeff_vec.
  if (coeff_vec != NULL) {
    if (coeff_vec->get_size() != rhs->get_size()) error("Mismatched vector lengths in solve_linear().");
    if (is_complex)
      for (int i = 0; i < ndof; i++) coeff_vec->set(i, rhs->get_cplx(i));
    else
      for (int i = 0; i < ndof; i++) coeff_vec->set(i, rhs->get(i));
  }

  // Free memory.
  mat->free_data();
  rhs->free_data();
  //solver->free_data();  // FIXME: to be implemented. A default destructor is called for the time being.
  delete solver;

  return true;
}

// Solve a typical linear problem using automatic adaptivity.
// Feel free to adjust this function for more advanced applications.
bool solve_linear_adapt(Hermes::Tuple<Space *> spaces, WeakForm* wf, Vector* coeff_vec, 
                        MatrixSolverType matrix_solver, Hermes::Tuple<int> proj_norms, 
                        Hermes::Tuple<Solution *> slns, Hermes::Tuple<Solution *> ref_slns, 
                        Hermes::Tuple<WinGeom *> sln_win_geom, Hermes::Tuple<WinGeom *> mesh_win_geom, 
                        Hermes::Tuple<RefinementSelectors::Selector *> selectors, AdaptivityParamType* apt,
                        bool verbose, Hermes::Tuple<ExactSolution *> exact_slns, bool is_complex) 
{
  // sanity checks
  if (spaces.size() != selectors.size()) 
    error("There must be a refinement selector for each solution component in solve_linear_adapt().");
  if (spaces.size() != proj_norms.size()) 
    error("There must be a projection norm for each solution component in solve_linear_adapt().");

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity parameters.
  double err_stop = apt->err_stop; 
  int ndof_stop = apt->ndof_stop;
  double threshold = apt->threshold;
  int strategy = apt->strategy; 
  int mesh_regularity = apt->mesh_regularity;
  double to_be_processed = apt->to_be_processed;
  int total_error_flag = apt->total_error_flag;
  int elem_error_flag = apt->elem_error_flag;

  // Number of physical fields in the problem.
  int num_comps = spaces.size();

  // Number of degreeso of freedom 
  int ndof = get_num_dofs(spaces);

  // Number of exact solutions given.
  if (exact_slns.size() != 0 && exact_slns.size() != num_comps) 
    error("Number of exact solutions does not match number of equations.");
  bool is_exact_solution;
  if (exact_slns.size() == num_comps) is_exact_solution = true;
  else is_exact_solution = false;

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, ndof, mat, rhs, solver, is_complex);

  // Initialize views.
  ScalarView* s_view[H2D_MAX_COMPONENTS];
  VectorView* v_view[H2D_MAX_COMPONENTS];
  OrderView*  o_view[H2D_MAX_COMPONENTS];
  for (int i = 0; i < num_comps; i++) {
    char* title = (char*)malloc(100*sizeof(char));
    if (sln_win_geom != Hermes::Tuple<WinGeom *>() && sln_win_geom[i] != NULL) {
      if (num_comps == 1) sprintf(title, "Solution", i); 
      else sprintf(title, "Solution[%d]", i); 
      switch (proj_norms[i]) {
        case H2D_L2_NORM:    s_view[i] = new ScalarView(title, sln_win_geom[i]);
                             s_view[i]->show_mesh(false);
                             v_view[i] = NULL;
                             break;
        case H2D_H1_NORM:    s_view[i] = new ScalarView(title, sln_win_geom[i]);
                             s_view[i]->show_mesh(false);
                             v_view[i] = NULL;
                             break;
        case H2D_HCURL_NORM: s_view[i] = NULL;
                             v_view[i] = new VectorView(title, sln_win_geom[i]);
                             break;
        case H2D_HDIV_NORM:  s_view[i] = NULL;
		             v_view[i] = new VectorView(title, sln_win_geom[i]);
                             break;
      default: error("Unknown norm in solve_linear_adapt().");
      }
    }
    else {
      s_view[i] = NULL;
      v_view[i] = NULL;
    }
    if (mesh_win_geom != Hermes::Tuple<WinGeom *>() && mesh_win_geom[i] != NULL) {
      if (num_comps == 1) sprintf(title, "Mesh", i); 
      else sprintf(title, "Mesh[%d]", i); 
      o_view[i] = new OrderView(title, mesh_win_geom[i]);
    }
    else o_view[i] = NULL;
  }

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_exact, graph_cpu_exact;

  // Conversion from Hermes::Tuple<Solution *> to Hermes::Tuple<MeshFunction *>
  // so that project_global() below compiles. 
  Hermes::Tuple<MeshFunction *> ref_slns_mf;
  for (int i = 0; i < num_comps; i++) {
    ref_slns_mf.push_back((MeshFunction*)ref_slns[i]);
  }

  int as = 1; bool done = false;
  do
  {
    if (verbose) {
      info("---- Adaptivity step %d:", as);
      info("Solving on reference mesh.");
    }

    // Construct globally refined reference mesh(es)
    // and setup reference space(s).
    Hermes::Tuple<Space *> ref_spaces;
    Hermes::Tuple<Mesh *> ref_meshes;
    for (int i = 0; i < num_comps; i++) {
      ref_meshes.push_back(new Mesh());
      Mesh *ref_mesh = ref_meshes.back();
      ref_mesh->copy(spaces[i]->get_mesh());
      ref_mesh->refine_all_elements();
      ref_spaces.push_back(spaces[i]->dup(ref_mesh));
      int order_increase = 1;
      ref_spaces[i]->copy_orders(spaces[i], order_increase);
    }

    // Solve the reference problem.
    // The NULL pointer means that we do not want the resulting coefficient vector.
    solve_linear(ref_spaces, wf, matrix_solver, ref_slns, NULL, is_complex);

    // Project the reference solution on the coarse mesh.
    if (verbose) info("Projecting reference solution on coarse mesh.");
    project_global(spaces, proj_norms, ref_slns_mf, slns, coeff_vec, is_complex); 

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution(s).
    for (int i = 0; i < num_comps; i++) {
      if (proj_norms[i] == H2D_H1_NORM || proj_norms[i] == H2D_L2_NORM) {
        if (s_view[i] != NULL) s_view[i]->show(slns[i]);
      }
      if (proj_norms[i] == H2D_HCURL_NORM || proj_norms[i] == H2D_HDIV_NORM) {
        if (v_view[i] != NULL) v_view[i]->show(slns[i]);
      }
      if (o_view[i] != NULL) o_view[i]->show(spaces[i]);
    }

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors.
    if (verbose) info("Calculating error (est).");
    Adapt hp(spaces, proj_norms);
    // Pass special error forms if any.
    for (int k = 0; k < apt->error_form_val.size(); k++) {
      hp.set_error_form(apt->error_form_i[k], apt->error_form_j[k], 
                        apt->error_form_val[k], apt->error_form_ord[k]);
    }
    hp.set_solutions(slns, ref_slns);
    // Below, apt->total_error_flag = H2D_TOTAL_ERROR_REL or H2D_TOTAL_ERROR_ABS
    // and apt->elem_error_flag = H2D_ELEMENT_ERROR_REL or H2D_ELEMENT_ERROR_ABS
    hp.calc_elem_errors(total_error_flag | elem_error_flag);
 
    // Calculate error estimate for each solution component. Note, these can 
    // have different norms, such as L2, H1, Hcurl and Hdiv. 
    double err_est_abs[H2D_MAX_COMPONENTS];
    double norm_est[H2D_MAX_COMPONENTS];
    double err_est_abs_total = 0;
    double norm_est_total = 0;
    double err_est_rel_total;
    for (int i = 0; i < num_comps; i++) {
      err_est_abs[i] = calc_abs_error(slns[i], ref_slns[i], proj_norms[i]);
      norm_est[i] = calc_norm(ref_slns[i], proj_norms[i]);
      err_est_abs_total += err_est_abs[i] * err_est_abs[i];
      norm_est_total += norm_est[i] * norm_est[i];
    }
    err_est_abs_total = sqrt(err_est_abs_total);
    norm_est_total = sqrt(norm_est_total);
    err_est_rel_total = err_est_abs_total / norm_est_total * 100.;

    // Calculate exact error for each solution component.   
    double err_exact_abs[H2D_MAX_COMPONENTS];
    double norm_exact[H2D_MAX_COMPONENTS];
    double err_exact_abs_total = 0;
    double norm_exact_total = 0;
    double err_exact_rel_total;
    if (is_exact_solution == true) {
      for (int i = 0; i < num_comps; i++) {
        err_exact_abs[i] = calc_abs_error(slns[i], exact_slns[i], proj_norms[i]);
        norm_exact[i] = calc_norm(exact_slns[i], proj_norms[i]);
        err_exact_abs_total += err_exact_abs[i] * err_exact_abs[i];
        norm_exact_total += norm_exact[i] * norm_exact[i];
      }
      err_exact_abs_total = sqrt(err_exact_abs_total);
      norm_exact_total = sqrt(norm_exact_total);
      err_exact_rel_total = err_exact_abs_total / norm_exact_total * 100.;
    }

    // Report results.
    if (verbose) {
      if (num_comps == 1) {
        info("ndof: %d, ref_ndof: %d, err_est_rel_total: %g%%", 
             get_num_dofs(spaces), get_num_dofs(ref_spaces), err_est_rel_total);
        if (is_exact_solution == true) info("err_exact_rel_total: %g%%", err_exact_rel_total);
      }
      else {
        for (int i = 0; i < num_comps; i++) {
          info("ndof[%d]: %d, ref_ndof[%d]: %d, err_est_rel[%d]: %g%%", 
               i, spaces[i]->get_num_dofs(), i, ref_spaces[i]->get_num_dofs(),
               i, err_est_abs[i]/norm_est[i]*100);
          if (is_exact_solution == true) info("err_exact_rel[%d]: %g%%", 
                                              i, err_exact_abs[i]/norm_exact[i]*100);
        }
        info("ndof: %d, ref_ndof: %d, err_est_rel_total: %g%%", 
             get_num_dofs(spaces), get_num_dofs(ref_spaces), err_est_rel_total);
      }
    }

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(get_num_dofs(spaces), err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");
    if (is_exact_solution == true) {
      graph_dof_exact.add_values(get_num_dofs(spaces), err_exact_rel_total);
      graph_dof_exact.save("conv_dof_exact.dat");
      graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel_total);
      graph_cpu_exact.save("conv_cpu_exact.dat");
    }

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < err_stop) done = true;
    else {
      if (verbose) info("Adapting the coarse mesh.");
      done = hp.adapt(selectors, threshold, strategy, mesh_regularity, to_be_processed);

      if (get_num_dofs(spaces) >= ndof_stop) done = true;
    }

    // Free reference meshes and spaces.
    for (int i = 0; i < num_comps; i++) {
      ref_spaces[i]->free(); // This does not free the associated mesh, we must do it separately.
      ref_meshes[i]->free();
    }

    as++;
  }
  while (done == false);

  if (verbose) info("Total running time: %g s", cpu_time.accumulated());
	
  return true;
}


