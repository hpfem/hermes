// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file hermes2d_common_defs.h
    \brief Common definitions for Hermes2D.
*/
#ifndef __H2D_COMMON_H_
#define __H2D_COMMON_H_

#include "../../hermes_common/include/hermes_common.h"

// H2D-specific error codes.
#define H2D_ERR_EDGE_INDEX_OUT_OF_RANGE         "Edge index out of range."

enum // node types
{
  HERMES_TYPE_VERTEX = 0,
  HERMES_TYPE_EDGE = 1
};

enum ElementMode2D {
	HERMES_MODE_TRIANGLE = 0,
	HERMES_MODE_QUAD = 1
};

class Ord2
{
  public:
    Ord2(int order_h, int order_v) : order_h(order_h), order_v(order_v) {};
    Ord2(int order) : order_h(order), order_v(order) {};
    int order_h;
    int order_v;
};

enum SpaceType {
  HERMES_H1_SPACE = 0,
  HERMES_HCURL_SPACE = 1,
  HERMES_HDIV_SPACE = 2,
  HERMES_L2_SPACE = 3,
  HERMES_INVALID_SPACE = -9999
};

#define H2D_MAX_ELEMENT_SONS 4 ///< A maximum number of sons of an element.

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

// Geometrical type of weak forms.
enum GeomType
{
  HERMES_PLANAR = 0,         // Planar problem.
  HERMES_AXISYM_X = 1,       // Axisymmetric problem where x-axis is the axis of symmetry.
  HERMES_AXISYM_Y = 2        // Axisymmetric problem where y-axis is the axis of symmetry.
};

// Enabling second derivatives in weak forms. Turned off by default. Second
// derivatives are employed, among others, by stabilization methods for
// transport equations. For usage see the example linear-convection-diffusion.
#define H2D_SECOND_DERIVATIVES_ENABLED

/* Uncomment this line to disable internal mesh compatibility
   tests in Traverse:begin(). */
//#define H2D_DISABLE_MULTIMESH_TESTS

template<typename Scalar> class MeshFunction;
template<typename Scalar> class Solution;
enum ProjNormType;
class RefMap;
template<typename Scalar> class DiscreteProblem;
template<typename Scalar> class Space;
template<typename Scalar> class WeakForm;
class Quad2D;
class Quad1DStd;
class Quad2DStd;

// Class for all global functions.
template<typename Scalar>
class HERMES_API Hermes2D {
public:
  // Error calculation in Hermes, useful for non-adaptive computations.
  double calc_rel_error(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, int norm_type) const; // Note: coarse mesh sln has to be first, then
                                                                                                  // ref_sln (because the abs. error is divided
                                                                                                  // by the norm of the latter).
  double calc_abs_error(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, int norm_type) const;
  double calc_norm(MeshFunction<Scalar>* sln, int norm_type) const;
  double calc_norms(Hermes::vector<Solution<Scalar>*> slns) const;         // Norms are determined from space_type in each Solution.

  // Function calculating errors between solutions in right and left vectors, returning all necessary parameters.
  // returns correct parameters only if the return value is true.
  // coarse mesh sln has to be first, then ref_sln.
  bool calc_errors(Hermes::vector<Solution<Scalar>* > left, Hermes::vector<Solution<Scalar>*> right,
                              Hermes::vector<double> & err_abs, Hermes::vector<double> & norm_vals,
                              double & err_abs_total, double & norm_total, double & err_rel_total,
                              Hermes::vector<ProjNormType> norms = Hermes::vector<ProjNormType>()) const;

  // Helper function.
  /// DEPRECATED.
  double calc_abs_error(double (*fn)(MeshFunction<Scalar>*, MeshFunction<Scalar>*, RefMap*, RefMap*),
                                          MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2) const;
  // Helper function.
  /// DEPRECATED.
  double calc_norm(double (*fn)(MeshFunction<Scalar>*, RefMap*), MeshFunction<Scalar>* sln) const;
  
  double error_fn_l2(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv) const;
  double norm_fn_l2(MeshFunction<Scalar>* sln, RefMap* ru) const;

  double error_fn_h1(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv) const;
  double norm_fn_h1(MeshFunction<Scalar>* sln, RefMap* ru) const;

  double error_fn_hc(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv) const;
  double norm_fn_hc(MeshFunction<Scalar>* sln, RefMap* ru) const;

  double error_fn_hcl2(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv) const;
  double norm_fn_hcl2(MeshFunction<Scalar>* sln, RefMap* ru) const;

  double error_fn_hdiv(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv) const;
  double norm_fn_hdiv(MeshFunction<Scalar>* sln, RefMap* ru) const;

  ///< Returns string representation of the quad order: used for debugging purposses.
  static const std::string get_quad_order_str(const int quad_order);
  ///< Returns the correct axial order for given edge.
  static int make_edge_order(int edge, int encoded_order, int mode);

  double get_l2_norm(Vector<Scalar>* vec) const;

  bool solve_newton(Scalar* coeff_vec, DiscreteProblem<Scalar>* dp, Hermes::Solvers::Solver<Scalar>* solver, SparseMatrix<Scalar>* matrix,
			       Vector<Scalar>* rhs, bool jacobian_changed = true, double newton_tol = 1e-8, 
                    int newton_max_iter = 100, bool verbose = false,
                    bool residual_as_function = false,
                    double damping_coeff = 1.0, double max_allowed_residual_norm = 1e6) const;

  bool solve_picard(WeakForm<Scalar>* wf, Space<Scalar>* space, Solution<Scalar>* sln_prev_iter, 
                    Hermes::MatrixSolverType matrix_solver, double tol = 1e-8, 
                    int max_iter = 100, bool verbose = false) const;
};
#endif

