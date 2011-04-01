// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_COMMON_H_
#define __H2D_COMMON_H_

#include "../../hermes_common/common.h"
#include "../../hermes_common/matrix.h"

// H2D-specific error codes.
#define H2D_ERR_EDGE_INDEX_OUT_OF_RANGE         "Edge index out of range."

// get_elem_marker(e) 
#define get_elem_marker(e) wf->get_element_markers_conversion()->get_user_marker(e->elem_marker)

#define H2D_NUM_MODES 2 ///< A number of modes, see enum ElementMode2D.

// how many bits the order number takes
const int H2D_ORDER_BITS = 5;
const int H2D_ORDER_MASK = (1 << H2D_ORDER_BITS) - 1;

// macros for combining quad horizontal and vertical orders
#define H2D_MAKE_QUAD_ORDER(h_order, v_order) (((v_order) << H2D_ORDER_BITS) + (h_order))
#define H2D_GET_H_ORDER(order) ((order) & H2D_ORDER_MASK)
#define H2D_GET_V_ORDER(order) ((order) >> H2D_ORDER_BITS)

// Enabling second derivatives in weak forms. Turned off by default. Second
// derivatives are employed, among others, by stabilization methods for
// transport equations. For usage see the example linear-convection-diffusion.
#define H2D_SECOND_DERIVATIVES_ENABLED

/* Uncomment this line to disable internal mesh compatibility
   tests in Traverse:begin(). */
//#define H2D_DISABLE_MULTIMESH_TESTS

class MeshFunction;
class Solution;
enum ProjNormType;
class RefMap;
class DiscreteProblem;
class Space;
class WeakForm;

// Class for all global functions.
class HERMES_API Hermes2D {
public:
  // Error calculation in Hermes, useful for non-adaptive computations.
  double calc_rel_error(MeshFunction* sln1, MeshFunction* sln2, int norm_type) const; // Note: coarse mesh sln has to be first, then
                                                                                                  // ref_sln (because the abs. error is divided
                                                                                                  // by the norm of the latter).
  double calc_abs_error(MeshFunction* sln1, MeshFunction* sln2, int norm_type) const;
  double calc_norm(MeshFunction* sln, int norm_type) const;
  double calc_norms(Hermes::vector<Solution*> slns) const;         // Norms are determined from space_type in each Solution.

  // Function calculating errors between solutions in right and left vectors, returning all necessary parameters.
  // returns correct parameters only if the return value is true.
  // coarse mesh sln has to be first, then ref_sln.
  bool calc_errors(Hermes::vector<Solution* > left, Hermes::vector<Solution *> right,
                              Hermes::vector<double> & err_abs, Hermes::vector<double> & norm_vals,
                              double & err_abs_total, double & norm_total, double & err_rel_total,
                              Hermes::vector<ProjNormType> norms = Hermes::vector<ProjNormType>()) const;

  // Helper function.
  /// DEPRECATED.
  double calc_abs_error(double (*fn)(MeshFunction*, MeshFunction*, RefMap*, RefMap*),
                                          MeshFunction* sln1, MeshFunction* sln2) const;
  // Helper function.
  /// DEPRECATED.
  double calc_norm(double (*fn)(MeshFunction*, RefMap*), MeshFunction* sln) const;
  
  double error_fn_l2(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv) const;
  double norm_fn_l2(MeshFunction* sln, RefMap* ru) const;

  double error_fn_h1(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv) const;
  double norm_fn_h1(MeshFunction* sln, RefMap* ru) const;

  double error_fn_hc(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv) const;
  double norm_fn_hc(MeshFunction* sln, RefMap* ru) const;

  double error_fn_hcl2(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv) const;
  double norm_fn_hcl2(MeshFunction* sln, RefMap* ru) const;

  double error_fn_hdiv(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv) const;
  double norm_fn_hdiv(MeshFunction* sln, RefMap* ru) const;

  ///< Returns string representation of the quad order: used for debugging purposses.
  static const std::string get_quad_order_str(const int quad_order);
  ///< Returns the correct axial order for given edge.
  static int make_edge_order(int edge, int encoded_order, int mode);
  


  double get_l2_norm(Vector* vec) const;

  /// New interface, still in developement
  /// HERMES_API bool solve_newton(scalar* coeff_vec, DiscreteProblem* dp, Solver* solver, SparseMatrix* matrix,
  ///		                   Vector* rhs, double NEWTON_TOL, int NEWTON_MAX_ITER, bool verbose,
  ///                              unsigned int stop_condition = NEWTON_WATCH_RESIDUAL);
  bool solve_newton(scalar* coeff_vec, DiscreteProblem* dp, Solver* solver, SparseMatrix* matrix,
		    Vector* rhs, double NEWTON_TOL, int NEWTON_MAX_ITER, bool verbose = false,
                    bool residual_as_function = false,
                    double damping_coeff = 1.0, double max_allowed_residual_norm = 1e6) const;

  bool solve_picard(WeakForm* wf, Space* space, Solution* sln_prev_iter, 
                    MatrixSolverType matrix_solver, double picard_tol, 
                    int picard_max_iter, bool verbose) const;
};

#endif

