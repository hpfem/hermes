/// This file is part of Hermes2D.
///
/// Hermes2D is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 2 of the License, or
/// (at your option) any later version.
///
/// Hermes2D is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY;without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with Hermes2D. If not, see <http:///www.gnu.org/licenses/>.

#ifndef __H2D_DISCRETE_PROBLEM_H
#define __H2D_DISCRETE_PROBLEM_H

#include "hermes_common.h"
#include "adapt/adapt.h"
#include "graph.h"
#include "forms.h"
#include "weakform/weakform.h"
#include "function/function.h"
#include "neighbor.h"
#include "refinement_selectors/selector.h"
#include "exceptions.h"

namespace Hermes
{
  namespace Hermes2D
  {
    class PrecalcShapeset;

    /// Discrete problem class.
    ///
    /// This class does assembling into external matrix / vector structures.
    ///
    template<typename Scalar>
    class HERMES_API DiscreteProblem : public DiscreteProblemInterface<Scalar>
    {
    public:
      /// Constructor for multiple components / equations.
      DiscreteProblem(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar> *> spaces);

      /// Constructor for one equation.
      DiscreteProblem(const WeakForm<Scalar>* wf, const Space<Scalar>* space);

      /// Non-parameterized constructor (currently used only in KellyTypeAdapt to gain access to NeighborSearch methods).
      DiscreteProblem();

      /// Destuctor.
      virtual ~DiscreteProblem();

      /// Assembling.
      /// General assembling procedure for nonlinear problems. coeff_vec is the
      /// previous Newton vector. If force_diagonal_block == true, then (zero) matrix
      /// antries are created in diagonal blocks even if corresponding matrix weak
      /// forms do not exist. This is useful if the matrix is later to be merged with
      /// a matrix that has nonzeros in these blocks. The Table serves for optional
      /// weighting of matrix blocks in systems. The parameter add_dir_lift decides
      /// whether Dirichlet lift will be added while coeff_vec is converted into
      /// Solutions.
      void assemble(Scalar* coeff_vec, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL,
        bool force_diagonal_blocks = false, Table* block_weights = NULL);

      /// Assembling.
      /// Without the matrix.
      void assemble(Scalar* coeff_vec, Vector<Scalar>* rhs = NULL,
        bool force_diagonal_blocks = false, Table* block_weights = NULL);

      /// Light version passing NULL for the coefficient vector. External solutions
      /// are initialized with zeros.
      void assemble(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL, bool force_diagonal_blocks = false,
        Table* block_weights = NULL);

      /// Light version passing NULL for the coefficient vector. External solutions
      /// are initialized with zeros.
      /// Without the matrix.
      void assemble(Vector<Scalar>* rhs = NULL, bool force_diagonal_blocks = false,
        Table* block_weights = NULL);

      void invalidate_matrix();

      /// Set this problem to Finite Volume.
      void set_fvm();

    protected:

      void init_assembling(Scalar* coeff_vec, PrecalcShapeset*** pss , PrecalcShapeset*** spss, RefMap*** refmaps, Solution<Scalar>*** u_ext, AsmList<Scalar>*** als, Hermes::vector<MeshFunction<Scalar>*>& ext_functions, MeshFunction<Scalar>*** ext, 
          Hermes::vector<MatrixFormVol<Scalar>*>* mfvol, Hermes::vector<MatrixFormSurf<Scalar>*>* mfsurf, Hermes::vector<VectorFormVol<Scalar>*>* vfvol, Hermes::vector<VectorFormSurf<Scalar>*>* vfsurf);

      void deinit_assembling(PrecalcShapeset*** pss , PrecalcShapeset*** spss, RefMap*** refmaps, Solution<Scalar>*** u_ext, AsmList<Scalar>*** als, Hermes::vector<MeshFunction<Scalar>*>& ext_functions, MeshFunction<Scalar>*** ext, 
          Hermes::vector<MatrixFormVol<Scalar>*>* mfvol, Hermes::vector<MatrixFormSurf<Scalar>*>* mfsurf, Hermes::vector<VectorFormVol<Scalar>*>* vfvol, Hermes::vector<VectorFormSurf<Scalar>*>* vfsurf);

      /// The form will be assembled.
      bool form_to_be_assembled(MatrixForm<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormVol<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormSurf<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorForm<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormVol<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormSurf<Scalar>* form, Traverse::State* current_state);

      // Return scaling coefficient.
      double block_scaling_coeff(MatrixForm<Scalar>* form);
      
      /// The current stage contains DG forms.
      void is_DG_stage();

      /// Get the number of unknowns.
      int get_num_dofs();

      /// Get info about presence of a matrix.
      bool is_matrix_free();

      // GET functions.
      /// Get pointer to n-th space.
      const Space<Scalar>* get_space(int n);

      /// Get the weak forms.
      const WeakForm<Scalar>* get_weak_formulation();

      /// Get all spaces as a Hermes::vector.
      Hermes::vector<const Space<Scalar>*> get_spaces();

      /// Preassembling.
      /// Precalculate matrix sparse structure.
      /// If force_diagonal_block == true, then (zero) matrix
      /// antries are created in diagonal blocks even if corresponding matrix weak
      /// forms do not exist. This is useful if the matrix is later to be merged with
      /// a matrix that has nonzeros in these blocks. The Table serves for optional
      /// weighting of matrix blocks in systems.
      void create_sparse_structure();
      void create_sparse_structure(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL);

      /// Initialize a state, returns a non-NULL Element.
      void init_state(PrecalcShapeset** current_pss, PrecalcShapeset** current_spss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state);
      void init_surface_state(AsmList<Scalar>** current_als, Traverse::State* current_state);

      /// Set the special handling of external functions of Runge-Kutta methods, including information how many spaces were there in the original problem.
      inline void set_RK(int original_spaces_count) { this->RungeKutta = true; RK_original_spaces_count = original_spaces_count; }

      /// Assemble one state.
      void assemble_one_state(PrecalcShapeset** current_pss, PrecalcShapeset** current_spss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state,
           MatrixFormVol<Scalar>** current_mfvol, MatrixFormSurf<Scalar>** current_mfsurf, VectorFormVol<Scalar>** current_vfvol, VectorFormSurf<Scalar>** current_vfsurf);
      
      /// Adjusts order to refmaps.
      void adjust_order_to_refmaps(Form<Scalar> *form, int& order, Hermes::Ord* o, RefMap** current_refmaps);

      /// Matrix volumetric forms - calculate the integration order.
      int calc_order_matrix_form(MatrixForm<Scalar>* mfv, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, Traverse::State* current_state);

      /// Matrix volumetric forms - assemble the form.
      void assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state);

      /// Vector volumetric forms - calculate the integration order.
      int calc_order_vector_form(VectorForm<Scalar>* mfv, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, Traverse::State* current_state);

      /// Vector volumetric forms - assemble the form.
      void assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state);

      /// \ingroup Helper methods inside {calc_order_*, assemble_*}
      /// Init geometry, jacobian * weights, return the number of integration points.
      int init_geometry_points(RefMap* reference_mapping, int order, Geom<double>*& geometry, double*& jacobian_x_weights);
      int init_surface_geometry_points(RefMap* reference_mapping, int& order, Traverse::State* current_state, Geom<double>*& geometry, double*& jacobian_x_weights);
      
      /// \ingroup Helper methods inside {calc_order_*, assemble_*}
      /// Calculates orders for external functions.
      void init_ext_orders(Form<Scalar> *form, Func<Hermes::Ord>** oi, ExtData<Hermes::Ord>* oext, Solution<Scalar>** current_u_ext, Traverse::State* current_state);
      /// \ingroup Helper methods inside {calc_order_*, assemble_*}
      /// Cleans up after init_ext_orders.
      void deinit_ext_orders(Form<Scalar> *form, Func<Hermes::Ord>** oi, ExtData<Hermes::Ord>* oext);

      /// \ingroup Helper methods inside {calc_order_*, assemble_*}
      /// Calculates external functions.
      void init_ext(Form<Scalar> *form, Func<Scalar>** u_ext, ExtData<Scalar>* ext, int order, Solution<Scalar>** current_u_ext, Traverse::State* current_state);
      /// \ingroup Helper methods inside {calc_order_*, assemble_*}
      /// Cleans up after init_ext.
      void deinit_ext(Form<Scalar> *form, Func<Scalar>** u_ext, ExtData<Scalar>* ext);
      
      /// Init function. Common code for the constructors.
      void init();

      /// Returns the matrix_buffer of the size n.
      Scalar** get_matrix_buffer(int n);

      /// Matrix structure as well as spaces and weak formulation is up-to-date.
      bool is_up_to_date();

      /// Minimum identifier of the meshes used in DG assembling in one stage.
      unsigned int min_dg_mesh_seq;

      /// Weak formulation.
      const WeakForm<Scalar>* wf;

      /// Seq number of the WeakForm.
      int wf_seq;

      /// Space instances for all equations in the system.
      Hermes::vector<const Space<Scalar>*> spaces;
      Hermes::vector<unsigned int> spaces_first_dofs;

      /// Seq numbers of Space instances in spaces.
      int* sp_seq;

      /// Number of DOFs of all Space instances in spaces.
      int ndof;

      /// Element usage flag: iempty[i] == true if the current state does not posses an active element in the i-th space.
      Hermes::vector<bool> isempty;

      /// Instance of the class Geom used in the calculation of integration order.
      Geom<Hermes::Ord> geom_ord;

      /// Fake weight used in the calculation of integration order.
      static double fake_wt;

      /// If the problem has only constant test functions, there is no need for order calculation,
      /// which saves time.
      bool is_fvm;

      Scalar** matrix_buffer;///< buffer for holding square matrix (during assembling)

      int matrix_buffer_dim;///< dimension of the matrix held by 'matrix_buffer'

      /// Matrix structure can be reused.
      /// If other conditions apply.
      bool have_matrix;

      /// There is a matrix form set on DG_INNER_EDGE area or not.
      bool DG_matrix_forms_present;

      /// There is a vector form set on DG_INNER_EDGE area or not.
      bool DG_vector_forms_present;

      /// Turn on Runge-Kutta specific handling of external functions.
      bool RungeKutta;

      /// Number of spaces in the original problem in a Runge-Kutta method.
      int RK_original_spaces_count;

      /// Storing assembling info.
      SparseMatrix<Scalar>* current_mat;
      Vector<Scalar>* current_rhs;
      bool current_force_diagonal_blocks;
      Table* current_block_weights;

      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class NewtonSolver;
      template<typename T> friend class PicardSolver;
      template<typename T> friend class RungeKutta;
    };
  }
}
#endif
