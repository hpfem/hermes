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

#include "h2d_common.h"
#include "quadrature/quad_all.h"

template<typename Scalar> class MeshFunction;
class Transformable;

template<typename Scalar>
const std::string Hermes2D<Scalar>::get_quad_order_str(const int quad_order) {
  std::stringstream str;
  str << "(H:" << H2D_GET_H_ORDER(quad_order) << ";V:" << H2D_GET_V_ORDER(quad_order) << ")";
  return str.str();
}

template<typename Scalar>
int Hermes2D<Scalar>::make_edge_order(int mode, int edge, int encoded_order)
{
  assert(edge < 4);

  if (mode == HERMES_MODE_TRIANGLE || edge == 0 || edge == 2)
    return H2D_GET_H_ORDER(encoded_order);
  else
    return H2D_GET_V_ORDER(encoded_order);
}

template<typename Scalar>
double Hermes2D<Scalar>::get_l2_norm(Vector<Scalar>* vec) const 
{
  _F_
  Scalar val = 0;
  for (unsigned int i = 0; i < vec->length(); i++) {
    Scalar inc = vec->get(i);
    val = val + inc*conj(inc);
  }
  return sqrt(std::abs(val));
}

template<typename Scalar>
bool Hermes2D<Scalar>::solve_newton(Scalar* coeff_vec, DiscreteProblem<Scalar>* dp, Solver<Scalar>* solver, SparseMatrix<Scalar>* matrix,
                  Vector<Scalar>* rhs, double newton_tol, int newton_max_iter, bool verbose,
                  bool residual_as_function,
                  double damping_coeff, double max_allowed_residual_norm) const
{
  // Prepare solutions for measuring residual norm.
  int num_spaces = dp->get_spaces().size();
  Hermes::vector<Solution<Scalar>*> solutions;
  Hermes::vector<bool> dir_lift_false;
  for (int i=0; i < num_spaces; i++) {
    if (residual_as_function) solutions.push_back(new Solution<Scalar>());
    dir_lift_false.push_back(false);      // No Dirichlet lifts will be considered.
  }

  // The Newton's loop.
  double residual_norm;
  int it = 1;
  while (1)
  {
    // Obtain the number of degrees of freedom.
    int ndof = dp->get_num_dofs();

    // Assemble the Jacobian matrix and residual vector.
    dp->assemble(coeff_vec, matrix, rhs, false);

    // Multiply the residual vector with -1 since the matrix
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    rhs->change_sign();

    // Measure the residual norm.
    if (residual_as_function) {
      // Translate the residual vector into a residual function (or multiple functions)
      // in the corresponding finite element space(s) and measure their norm(s) there.
      // This is more meaningful than just measuring the l2-norm of the residual vector,
      // since in the FE space not all components in the residual vector have the same weight.
      // On the other hand, this is slower as it requires global norm calculation, and thus
      // numerical integration over the entire domain. Therefore this option is off by default.
      Solution<Scalar>::vector_to_solutions(rhs, dp->get_spaces(), solutions, dir_lift_false);
      residual_norm = calc_norms(solutions);
    }
    else {
      // Calculate the l2-norm of residual vector, this is the traditional way.
      residual_norm = get_l2_norm(rhs);
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
      for (unsigned int i = 0; i < solutions.size(); i++)
        delete solutions[i];
      return false;
    }

    // If residual norm is within tolerance, or the maximum number
    // of iteration has been reached, then quit.
    if ((residual_norm < newton_tol || it > newton_max_iter) && it > 1) break;

    // Solve the linear system.
    if(!solver->solve()) error ("Matrix<Scalar> solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof; i++) coeff_vec[i] += damping_coeff * solver->get_solution()[i];

    it++;
  }

  for (unsigned int i = 0; i < solutions.size(); i++)
    delete solutions[i];

  if (it >= newton_max_iter) {
    if (verbose) info("Maximum allowed number of Newton iterations exceeded, returning false.");
    return false;
  }

  return true;
}

// Perform Picard's iteration.
template<typename Scalar>
bool Hermes2D<Scalar>::solve_picard(WeakForm<Scalar>* wf, Space<Scalar>* space, Solution<Scalar>* sln_prev_iter,
                  MatrixSolverType matrix_solver, double picard_tol,
                  int picard_max_iter, bool verbose) const
{
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix<Scalar>* matrix = create_matrix<Scalar>(matrix_solver);
  Vector<Scalar>* rhs = create_vector<Scalar>(matrix_solver);
  Solver<Scalar>* solver = create_linear_solver<Scalar>(matrix_solver, matrix, rhs);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem<Scalar> dp(wf, space, is_linear);

  int iter_count = 0;
  while (true) {
    // Assemble the stiffness matrix and right-hand side.
    dp.assemble(matrix, rhs);

    // Solve the linear system and if successful, obtain the solution.
    Solution<Scalar> sln_new;
    if(solver->solve()) Solution<Scalar>::vector_to_solution(solver->get_solution(), space, &sln_new);
    else error ("Matrix<Scalar> solver failed.\n");

    double rel_error = calc_abs_error(sln_prev_iter, &sln_new, HERMES_H1_NORM)
                       / calc_norm(&sln_new, HERMES_H1_NORM) * 100;
    if (verbose) info("---- Picard iter %d, ndof %d, rel. error %g%%",
      iter_count+1, space->get_num_dofs(), rel_error);

    // Stopping criterion.
    if (rel_error < picard_tol) {
      sln_prev_iter->copy(&sln_new);
      delete matrix;
      delete rhs;
      delete solver;
      return true;
    }

    if (iter_count >= picard_max_iter) {
      delete matrix;
      delete rhs;
      delete solver;
      if (verbose) info("Maximum allowed number of Picard iterations exceeded, returning false.");
      return false;
    }

    // Saving solution for the next iteration;
    sln_prev_iter->copy(&sln_new);

    iter_count++;
  }
}


template<typename Scalar>
double Hermes2D<Scalar>::calc_abs_error(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, int norm_type) const
{
  // sanity checks
  if (sln1 == NULL) error("sln1 is NULL in calc_abs_error().");
  if (sln2 == NULL) error("sln2 is NULL in calc_abs_error().");

  Quad2D* quad = &g_quad_2d_std;
  sln1->set_quad_2d(quad);
  sln2->set_quad_2d(quad);

  Mesh* meshes[2] = { sln1->get_mesh(), sln2->get_mesh() };
  Transformable* tr[2] = { sln1, sln2 };
  Traverse trav;
  trav.begin(2, meshes, tr);

  double error = 0.0;
  Element** ee;
  while ((ee = trav.get_next_state(NULL, NULL)) != NULL)
  {
    update_limit_table(ee[0]->get_mode());

    RefMap* ru = sln1->get_refmap();
    RefMap* rv = sln2->get_refmap();
    switch (norm_type) {
      case HERMES_L2_NORM:
        error += error_fn_l2(sln1, sln2, ru, rv);
        break;
      case HERMES_H1_NORM:
        error += error_fn_h1(sln1, sln2, ru, rv);
        break;
      case HERMES_HCURL_NORM:
        error += error_fn_hc(sln1, sln2, ru, rv);
        break;
      case HERMES_HDIV_NORM:
        error += error_fn_hdiv(sln1, sln2, ru, rv);
        break;
      default: error("Unknown norm in calc_error().");
    }
  }
  trav.finish();
  return sqrt(error);
}

template<typename Scalar>
double Hermes2D<Scalar>::calc_norm(MeshFunction<Scalar>* sln, int norm_type) const
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);

  double norm = 0.0;
  Element* e;
  Mesh* mesh = sln->get_mesh();

  for_all_active_elements(e, mesh)
  {
    // set maximum integration order for use in integrals, see limit_order()
    update_limit_table(e->get_mode());

    sln->set_active_element(e);
    RefMap* ru = sln->get_refmap();

    switch (norm_type) {
      case HERMES_L2_NORM:
        norm += norm_fn_l2(sln, ru);
        break;
      case HERMES_H1_NORM:
        norm += norm_fn_h1(sln, ru);
        break;
      case HERMES_HCURL_NORM:
        norm += norm_fn_hc(sln, ru);
        break;
      case HERMES_HDIV_NORM:
        norm += norm_fn_hdiv(sln, ru);
        break;
      default: error("Unknown norm in calc_norm().");
    }
  }
  return sqrt(norm);
}

// Calculate norm of a (possibly vector-valued) solution.
// Take norm from spaces where these solutions belong.
template<typename Scalar>
double Hermes2D<Scalar>::calc_norms(Hermes::vector<Solution<Scalar>*> slns) const
{
  // Calculate norms for all solutions.
  Hermes::vector<double> norms;
  int n = slns.size();
  for (int i=0; i<n; i++) {
    switch (slns[i]->get_space_type()) {
      case HERMES_H1_SPACE: norms.push_back(calc_norm(slns[i], HERMES_H1_NORM)); break;
      case HERMES_HCURL_SPACE: norms.push_back(calc_norm(slns[i], HERMES_HCURL_NORM)); break;
      case HERMES_HDIV_SPACE: norms.push_back(calc_norm(slns[i], HERMES_HDIV_NORM)); break;
      case HERMES_L2_SPACE: norms.push_back(calc_norm(slns[i], HERMES_L2_NORM)); break;
      default: error("Internal in calc_norms(): unknown space type.");
    }
  }
  // Calculate the resulting norm.
  double result = 0;
  for (int i=0; i<n; i++) result += norms[i]*norms[i];
  return sqrt(result);
}


template<typename Scalar>
bool Hermes2D<Scalar>::calc_errors(Hermes::vector<Solution<Scalar>* > left, Hermes::vector<Solution<Scalar>*> right, Hermes::vector<double> & err_abs, Hermes::vector<double> & norm_vals,
                 double & err_abs_total, double & norm_total, double & err_rel_total, Hermes::vector<ProjNormType> norms) const
{
  bool default_norms = false;
  // Checks.
  if(left.size() != right.size())
    return false;
  if (norms != Hermes::vector<ProjNormType>())
  {
    if(left.size() != norms.size())
      return false;
  }
  else
    default_norms = true;

  // Zero the resulting Tuples.
  err_abs.clear();
  norm_vals.clear();

  // Zero the sums.
  err_abs_total = 0;
  norm_total = 0;
  err_rel_total = 0;

  // Calculation.
  for(unsigned int i = 0; i < left.size(); i++)
  {
    err_abs.push_back(calc_abs_error(left[i], right[i], default_norms ? HERMES_H1_NORM : norms[i]));
    norm_vals.push_back(calc_norm(right[i], default_norms ? HERMES_H1_NORM : norms[i]));
    err_abs_total += err_abs[i] * err_abs[i];
    norm_total += norm_vals[i] * norm_vals[i];
  }

  err_abs_total = sqrt(err_abs_total);
  norm_total = sqrt(norm_total);
  err_rel_total = err_abs_total / norm_total * 100.;

  // Everything went well, return appropriate flag.
  return true;
}


/// Calculates the absolute error between sln1 and sln2 using function fn
template<typename Scalar>
double Hermes2D<Scalar>::calc_abs_error(double (*fn)(MeshFunction<Scalar>*, MeshFunction<Scalar>*, RefMap*, RefMap*), MeshFunction<Scalar>* sln1,
                      MeshFunction<Scalar>* sln2) const
{
  // sanity checks
  if (fn == NULL) error("error norm function is NULL in calc_abs_error().");
  if (sln1 == NULL) error("sln1 is NULL in calc_abs_error().");
  if (sln2 == NULL) error("sln2 is NULL in calc_abs_error().");

  Quad2D* quad = &g_quad_2d_std;
  sln1->set_quad_2d(quad);
  sln2->set_quad_2d(quad);

  Mesh* meshes[2] = { sln1->get_mesh(), sln2->get_mesh() };
  Transformable* tr[2] = { sln1, sln2 };
  Traverse trav;
  trav.begin(2, meshes, tr);

  double error = 0.0;
  Element** ee;
  while ((ee = trav.get_next_state(NULL, NULL)) != NULL)
  {
    update_limit_table(ee[0]->get_mode());

    RefMap* ru = sln1->get_refmap();
    RefMap* rv = sln2->get_refmap();

    error += fn(sln1, sln2, ru, rv);
  }
  trav.finish();
  return sqrt(error);
}


/// Calculates the norm of sln using function fn
template<typename Scalar>
double Hermes2D<Scalar>::calc_norm(double (*fn)(MeshFunction<Scalar>*, RefMap*), MeshFunction<Scalar>* sln) const
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);

  double norm = 0.0;
  Element* e;
  Mesh* mesh = sln->get_mesh();

  for_all_active_elements(e, mesh)
  {
    // set maximum integration order for use in integrals, see limit_order()
    update_limit_table(e->get_mode());

    sln->set_active_element(e);
    RefMap* ru = sln->get_refmap();

    norm += fn(sln, ru);
  }
  return sqrt(norm);
}

template<typename Scalar>
double Hermes2D<Scalar>::calc_rel_error(MeshFunction<Scalar>* sln, MeshFunction<Scalar>* ref_sln, int norm_type) const
{
  double error = calc_abs_error(sln, ref_sln, norm_type);
  double norm = calc_norm(ref_sln, norm_type);

  return error/norm;
}

//// H1 space //////////////////////////////////////////////////////////////////////////////////////
// function used to calculate error in H1 norm
template<typename Scalar>
double Hermes2D<Scalar>::error_fn_h1(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv) const
{
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2*std::max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o);
  sln2->set_quad_order(o);

  Scalar *uval, *vval, *dudx, *dudy, *dvdx, *dvdy;
  uval = sln1->get_fn_values();
  vval = sln2->get_fn_values();
  sln1->get_dx_dy_values(dudx, dudy);
  sln2->get_dx_dy_values(dvdx, dvdy);

  double result = 0.0;
  h1_integrate_expression(sqr(uval[i] - vval[i]) +
                          sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]));
  return result;
}

// function used to calculate H1 norm of the solution
template<typename Scalar>
double Hermes2D<Scalar>::norm_fn_h1(MeshFunction<Scalar>* sln, RefMap* ru) const
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 * sln->get_fn_order() + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o);

  Scalar *uval, *dudx, *dudy;
  uval = sln->get_fn_values();
  sln->get_dx_dy_values(dudx, dudy);

  double result = 0.0;
  h1_integrate_expression(sqr(uval[i]) + sqr(dudx[i]) + sqr(dudy[i]));
  return result;
}

//// L2 space //////////////////////////////////////////////////////////////////////////////////////
// function used to calculate error in L2 norm
template<typename Scalar>
double Hermes2D<Scalar>::error_fn_l2(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv) const
{
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2*std::max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o, H2D_FN_VAL);
  sln2->set_quad_order(o, H2D_FN_VAL);

  Scalar *uval, *vval;
  uval = sln1->get_fn_values();
  vval = sln2->get_fn_values();

  double result = 0.0;
  h1_integrate_expression(sqr(uval[i] - vval[i]));
  return result;
}

// function used to calculate L2 norm of the solution
template<typename Scalar>
double Hermes2D<Scalar>::norm_fn_l2(MeshFunction<Scalar>* sln, RefMap* ru) const
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 *sln->get_fn_order() + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o, H2D_FN_VAL);

  Scalar* uval = sln->get_fn_values();

  double result = 0.0;
  h1_integrate_expression(sqr(uval[i]));
  return result;
}

//// Hcurl space ///////////////////////////////////////////////////////////////////////////////////
// function used to calculate error in Hcurl norm
template<typename Scalar>
double Hermes2D<Scalar>::error_fn_hc(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv) const
{
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2 * std::max(sln1->get_fn_order(), sln2->get_fn_order()) + 2 + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o);
  sln2->set_quad_order(o);


  Scalar *uval0 = sln1->get_fn_values(0), *uval1 = sln1->get_fn_values(1);
  Scalar *udx1  = sln1->get_dx_values(1), *udy0  = sln1->get_dy_values(0);
  Scalar *vval0 = sln2->get_fn_values(0), *vval1 = sln2->get_fn_values(1);
  Scalar *vdx1  = sln2->get_dx_values(1), *vdy0  = sln2->get_dy_values(0);

  double result = 0.0;
  h1_integrate_expression(sqr(uval0[i] - vval0[i]) + sqr(uval1[i] - vval1[i]) +
                          sqr((udx1[i] - udy0[i]) - (vdx1[i] - vdy0[i])));
  return result;
}

// function used to calculate Hcurl norm
template<typename Scalar>
double Hermes2D<Scalar>::norm_fn_hc(MeshFunction<Scalar>* sln, RefMap* ru) const
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 * sln->get_fn_order() + 2 + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o);

  Scalar *uval0 = sln->get_fn_values(0), *uval1 = sln->get_fn_values(1);
  Scalar *udx1  = sln->get_dx_values(1), *udy0  = sln->get_dy_values(0);

  double result = 0.0;
  h1_integrate_expression(sqr(uval0[i]) + sqr(uval1[i]) + sqr(udx1[i] - udy0[i]));
  return result;
}

// function used to calculate error in Hcurl norm
template<typename Scalar>
double Hermes2D<Scalar>::error_fn_hcl2(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv) const
{
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2 * std::max(sln1->get_fn_order(), sln2->get_fn_order()) + 2 + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o);
  sln2->set_quad_order(o);


  Scalar *uval0 = sln1->get_fn_values(0), *uval1 = sln1->get_fn_values(1);
  Scalar *vval0 = sln2->get_fn_values(0), *vval1 = sln2->get_fn_values(1);

  double result = 0.0;
  h1_integrate_expression(sqr(uval0[i] - vval0[i]) + sqr(uval1[i] - vval1[i]));
  return result;
}

// function used to calculate Hcurl norm
template<typename Scalar>
double Hermes2D<Scalar>::norm_fn_hcl2(MeshFunction<Scalar>* sln, RefMap* ru) const
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 * sln->get_fn_order() + 2 + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o);

  Scalar *uval0 = sln->get_fn_values(0), *uval1 = sln->get_fn_values(1);
  Scalar *udx1  = sln->get_dx_values(1), *udy0  = sln->get_dy_values(0);

  double result = 0.0;
  h1_integrate_expression(sqr(uval0[i]) + sqr(uval1[i]));
  return result;
}

//// Hdiv space ///////////////////////////////////////////////////////////////////////////////////
// function used to calculate error in Hcurl norm
template<typename Scalar>
double Hermes2D<Scalar>::error_fn_hdiv(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv) const
{
  error("error_fn_hdiv() not implemented yet.");

  // Hcurl code
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2 * std::max(sln1->get_fn_order(), sln2->get_fn_order()) + 2 + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o);
  sln2->set_quad_order(o);


  Scalar *uval0 = sln1->get_fn_values(0), *uval1 = sln1->get_fn_values(1);
  Scalar *udx1  = sln1->get_dx_values(1), *udy0  = sln1->get_dy_values(0);
  Scalar *vval0 = sln2->get_fn_values(0), *vval1 = sln2->get_fn_values(1);
  Scalar *vdx1  = sln2->get_dx_values(1), *vdy0  = sln2->get_dy_values(0);

  double result = 0.0;
  h1_integrate_expression(sqr(uval0[i] - vval0[i]) + sqr(uval1[i] - vval1[i]) +
                          sqr((udx1[i] - udy0[i]) - (vdx1[i] - vdy0[i])));
  return result;
}

// function used to calculate Hcurl norm
template<typename Scalar>
double Hermes2D<Scalar>::norm_fn_hdiv(MeshFunction<Scalar>* sln, RefMap* ru) const
{
  error("norm_fn_hdiv() not implemented yet.");

  // Hcurl code
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 * sln->get_fn_order() + 2 + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o);

  Scalar *uval0 = sln->get_fn_values(0), *uval1 = sln->get_fn_values(1);
  Scalar *udx1  = sln->get_dx_values(1), *udy0  = sln->get_dy_values(0);

  double result = 0.0;
  h1_integrate_expression(sqr(uval0[i]) + sqr(uval1[i]) + sqr(udx1[i] - udy0[i]));
  return result;
}


//// python support //////////////////////////////////////////////////////////////////////////////////
HERMES_API void throw_exception(char *text)
{
  throw std::runtime_error(text);
}

template class HERMES_API Hermes2D<double>;
template class HERMES_API Hermes2D<std::complex<double>>;