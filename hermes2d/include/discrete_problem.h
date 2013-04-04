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
#include "exceptions.h"
#include "mixins2d.h"
#include "discrete_problem/dg/discrete_problem_dg_assembler.h"
#include "discrete_problem/discrete_problem_cache.h"
#include "discrete_problem/discrete_problem_helpers.h"
#include "discrete_problem/discrete_problem_thread_assembler.h"
#include "discrete_problem/discrete_problem_integration_order_calculator.h"

namespace Hermes
{
  namespace Hermes2D
  {
    class PrecalcShapeset;
    /// @ingroup inner
    /// Discrete problem class.
    ///
    /// This class does assembling into external matrix / vector structures.
    ///
    template<typename Scalar>
    class HERMES_API DiscreteProblem : public DiscreteProblemInterface<Scalar>,
      public Hermes::Mixins::TimeMeasurable, 
      public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>, 
      public Hermes::Hermes2D::Mixins::StateQueryable,
      public Hermes::Hermes2D::Mixins::DiscreteProblemRungeKutta<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemWeakForm<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemCacheSettings,
      public Hermes::Hermes2D::Mixins::DiscreteProblemMatrixVector<Scalar>
    {
    public:
      /// Constructor for multiple components / equations.
      DiscreteProblem(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> > spaces);

      /// Constructor for one equation.
      DiscreteProblem(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar> space);

      /// Non-parameterized constructor.
      DiscreteProblem();

      /// Destuctor.
      virtual ~DiscreteProblem();

      /// Sets new spaces for the instance.
      virtual void set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> > spaces);
      virtual void set_space(SpaceSharedPtr<Scalar> space);

      /// Set the weak forms.
      void set_weak_formulation(WeakForm<Scalar>* wf);

      /// Get all spaces as a Hermes::vector.
      virtual Hermes::vector<SpaceSharedPtr<Scalar> > get_spaces() const;

      /// Get info about presence of a matrix.
      bool is_matrix_free() const;

      /// set time information for time-dependent problems.
      virtual void set_time(double time);
      virtual void set_time_step(double time_step);

      /// Assembling.
      void assemble(Scalar* coeff_vec, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL);

      /// Assembling.
      /// Without the matrix.
      void assemble(Scalar* coeff_vec, Vector<Scalar>* rhs = NULL);

      /// Light version passing NULL for the coefficient vector. External solutions
      /// are initialized with zeros.
      virtual void assemble(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL);

      /// Light version passing NULL for the coefficient vector. External solutions
      /// are initialized with zeros.
      /// Without the matrix.
      void assemble(Vector<Scalar>* rhs);

      

    protected:
      /// State querying helpers.
      virtual bool isOkay() const;
      virtual inline std::string getClassName() const { return "DiscreteProblem"; }

      void init_assembling(Scalar* coeff_vec, PrecalcShapeset*** pss, RefMap*** refmaps, Solution<Scalar>*** u_ext, AsmList<Scalar>*** als, AsmList<Scalar>**** alsSurface, WeakForm<Scalar>** weakforms, int num_threads);

      void deinit_assembling(PrecalcShapeset*** pss, RefMap*** refmaps, Solution<Scalar>*** u_ext, AsmList<Scalar>*** als, AsmList<Scalar>**** alsSurface, WeakForm<Scalar>** weakforms, int num_threads);

      /// The form will be assembled.
      bool form_to_be_assembled(MatrixForm<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormVol<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(MatrixFormSurf<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorForm<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormVol<Scalar>* form, Traverse::State* current_state);
      bool form_to_be_assembled(VectorFormSurf<Scalar>* form, Traverse::State* current_state);

      /// Preassembling.
      /// Precalculate matrix sparse structure.
      /// If force_diagonal_block == true, then (zero) matrix
      /// antries are created in diagonal blocks even if corresponding matrix weak
      /// forms do not exist. This is useful if the matrix is later to be merged with
      /// a matrix that has nonzeros in these blocks. The Table serves for optional
      /// weighting of matrix blocks in systems.
      void create_sparse_structure();
      void create_sparse_structure(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL);
      
      /// Assemble one state - needs recalculation?
      /// \return if one needs to recalculate, the method calculate_cache_records is called.
      typename DiscreteProblemCache<Scalar>::CacheRecord* get_state_cache(Traverse::State* state, PrecalcShapeset** current_pss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, AsmList<Scalar>*** current_alsSurface, WeakForm<Scalar>* current_wf, int& order);

      /// Assemble one state.
      void assemble_one_state(typename DiscreteProblemCache<Scalar>::CacheRecord* cache_record, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state, WeakForm<Scalar>* current_wf);

      /// Matrix volumetric forms - assemble the form.
      virtual void assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext,
        AsmList<Scalar>* current_als_i, AsmList<Scalar>* current_als_j, Traverse::State* current_state, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights);

      /// Vector volumetric forms - assemble the form.
      void assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext, 
        AsmList<Scalar>* current_als, Traverse::State* current_state, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights);

      /// Init function. Common code for the constructors.
      void init();

      /// Matrix structure as well as spaces and weak formulation is up-to-date.
      bool is_up_to_date() const;

      /// Space instances for all equations in the system.
      Hermes::vector<SpaceSharedPtr<Scalar> > spaces;
      int spaces_size;

      /// Seq numbers of Space instances in spaces.
      int* sp_seq;
      
      /// Internal.
      bool is_linear;

      /// Matrix structure can be reused.
      /// If other conditions apply.
      bool matrix_structure_reusable;
      
      /// DiscreteProblemMatrixVector methods.
      virtual void set_matrix(SparseMatrix<Scalar>* mat);
      virtual void set_rhs(Vector<Scalar>* rhs);

      /// Exception caught in a parallel region.
      std::exception* caughtException;

      /// The cache.
      DiscreteProblemCache<Scalar> cache;

      /// The DG assembler.
      DiscreteProblemDGAssembler<Scalar> dgAssembler;
      
      /// Assembly data.
      DiscreteProblemThreadAssembler<Scalar>* threadAssembler;

      /// Integration order calculator.
      DiscreteProblemIntegrationOrderCalculator<Scalar> integrationOrderCalculator;

      int num_threads_used;

      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class LinearSolver;
      template<typename T> friend class NewtonSolver;
      template<typename T> friend class PicardSolver;
      template<typename T> friend class RungeKutta;
    };
  }
}
#endif
