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
#include "weakform/weakform.h"
#include "function/function.h"
#include "exceptions.h"
#include "mixins2d.h"
#include "discrete_problem/discrete_problem_cache.h"
#include "discrete_problem/discrete_problem_helpers.h"
#include "discrete_problem/discrete_problem_thread_assembler.h"

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
      public Hermes::Mixins::IntegrableWithGlobalOrder,
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
      /// The inner-most version.
      void assemble(Solution<Scalar>** u_ext_sln, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs);

      /// Initialize states.
      void init_assembling(Traverse::State**& states, int& num_states, Solution<Scalar>** u_ext_sln);
      void deinit_assembling(Traverse::State** states, int num_states);


      /// set time information for time-dependent problems.
      virtual void set_time(double time);
      virtual void set_time_step(double time_step);

      /// Sets new spaces for the instance.
      virtual void set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> > spaces);
      virtual void set_space(SpaceSharedPtr<Scalar> space);

      /// Set the weak forms.
      void set_weak_formulation(WeakForm<Scalar>* wf);

      /// Get all spaces as a Hermes::vector.
      virtual Hermes::vector<SpaceSharedPtr<Scalar> > get_spaces() const;

    protected:
      /// State querying helpers.
      virtual bool isOkay() const;
      virtual inline std::string getClassName() const { return "DiscreteProblem"; }

      /// Preassembling.
      /// Precalculate matrix sparse structure.
      /// If force_diagonal_block == true, then (zero) matrix
      /// antries are created in diagonal blocks even if corresponding matrix weak
      /// forms do not exist. This is useful if the matrix is later to be merged with
      /// a matrix that has nonzeros in these blocks. The Table serves for optional
      /// weighting of matrix blocks in systems.
      void create_sparse_structure();
      void create_sparse_structure(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL);
      
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
      bool nonlinear;

      /// Matrix structure can be reused.
      /// If other conditions apply.
      bool matrix_structure_reusable;
      
      /// DiscreteProblemMatrixVector methods.
      virtual void set_matrix(SparseMatrix<Scalar>* mat);
      virtual void set_rhs(Vector<Scalar>* rhs);

      /// The cache.
      DiscreteProblemCache<Scalar> cache;
      
      /// Assembly data.
      DiscreteProblemThreadAssembler<Scalar>* threadAssembler;

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
