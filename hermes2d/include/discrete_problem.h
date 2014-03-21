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
    class HERMES_API DiscreteProblem : 
      public Hermes::Mixins::Loggable,
      public Hermes::Mixins::TimeMeasurable,
      public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>, 
      public Hermes::Mixins::StateQueryable,
      public Hermes::Hermes2D::Mixins::DiscreteProblemRungeKutta<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemWeakForm<Scalar>,
      public Hermes::Mixins::IntegrableWithGlobalOrder,
      public Hermes::Hermes2D::Mixins::DiscreteProblemMatrixVector<Scalar>,
      public Hermes::Hermes2D::Mixins::Parallel
    {
    public:
      /// Constructor for multiple components / equations.
      /// Making this DiscreteProblem linear does 2 things
      /// 1 - turns off initialization of previous iterations for nonlinear solvers.
      /// 2 - allows for assembling Dirichlet boundary conditions using a Dirichlet lift.
      /// \param[in] dirichlet_lift_accordingly If true, the appropriate settings for (linear / nonlinear)
      /// problem will be used (use Dirichlet lift iff the problem is linear). If false, the other way round.
      DiscreteProblem(WeakFormSharedPtr<Scalar> wf, SpaceSharedPtrVector<Scalar> spaces, bool linear = false, bool dirichlet_lift_accordingly = true);
      /// Constructor for one equation.
      /// Making this DiscreteProblem linear does 2 things
      /// 1 - turns off initialization of previous iterations for nonlinear solvers.
      /// 2 - allows for assembling Dirichlet boundary conditions using a Dirichlet lift.
      /// \param[in] dirichlet_lift_accordingly If true, the appropriate settings for (linear / nonlinear)
      /// problem will be used (use Dirichlet lift iff the problem is linear). If false, the other way round.
      DiscreteProblem(WeakFormSharedPtr<Scalar> wf, SpaceSharedPtr<Scalar> space, bool linear = false, bool dirichlet_lift_accordingly = true);
      /// Non-parameterized constructor.
      /// Making this DiscreteProblem linear does 2 things
      /// 1 - turns off initialization of previous iterations for nonlinear solvers.
      /// 2 - allows for assembling Dirichlet boundary conditions using a Dirichlet lift.
      /// \param[in] dirichlet_lift_accordingly If true, the appropriate settings for (linear / nonlinear)
      /// problem will be used (use Dirichlet lift iff the problem is linear). If false, the other way round.
      DiscreteProblem(bool linear = false, bool dirichlet_lift_accordingly = true);
      /// Destuctor.
      virtual ~DiscreteProblem();

      /// Assembling.
      bool assemble(Scalar* coeff_vec, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = nullptr);
      /// Assembling.
      /// Without the matrix.
      bool assemble(Scalar* coeff_vec, Vector<Scalar>* rhs = nullptr);
      /// Light version passing nullptr for the coefficient vector. External solutions
      /// are initialized with zeros.
      bool assemble(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = nullptr);
      /// Light version passing nullptr for the coefficient vector. External solutions
      /// are initialized with zeros.
      /// Without the matrix.
      bool assemble(Vector<Scalar>* rhs);
      /// Assembling.
      bool assemble(Solution<Scalar>** u_ext_sln, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs);

      /// set time information for time-dependent problems.
      void set_time(double time);
      void set_time_step(double time_step);

      /// Sets new_ spaces for the instance.
      void set_spaces(SpaceSharedPtrVector<Scalar> spaces);
      void set_space(SpaceSharedPtr<Scalar> space);

      /// Set the weak forms.
      void set_weak_formulation(WeakFormSharedPtr<Scalar> wf);

      /// Get all spaces as a std::vector.
      SpaceSharedPtrVector<Scalar> get_spaces();

      /// Experimental.
      typedef void(*reassembled_states_reuse_linear_system_fn)(Traverse::State**& states, unsigned int& num_states, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs);
      void set_reassembled_states_reuse_linear_system_fn(reassembled_states_reuse_linear_system_fn fn) {
        this->reassembled_states_reuse_linear_system = fn;
      }
      reassembled_states_reuse_linear_system_fn reassembled_states_reuse_linear_system;
      void set_reusable_DOFs(bool **reusable_DOFs, bool **reusable_Dirichlet)
      {
        for (int i = 0; i < this->num_threads_used; i++)
        {
          this->threadAssembler[i]->reusable_DOFs = reusable_DOFs;
          this->threadAssembler[i]->reusable_Dirichlet = reusable_Dirichlet;
        }
      }

      /// See Hermes::Mixins::Loggable.
      virtual void set_verbose_output(bool to_set);

    protected:
      /// Initialize states.
      void init_assembling(Traverse::State**& states, unsigned int& num_states, Solution<Scalar>** u_ext_sln, MeshSharedPtrVector& meshes);
      void deinit_assembling(Traverse::State** states, unsigned  int num_states);

      /// RungeKutta helpers.
      void set_RK(int original_spaces_count, bool force_diagonal_blocks = nullptr, Table* block_weights = nullptr);

      /// State querying helpers.
      bool isOkay() const;
      inline std::string getClassName() const { return "DiscreteProblem"; }

      /// Init function. Common code for the constructors.
      void init(bool linear, bool dirichlet_lift_accordingly);

      /// Space instances for all equations in the system.
      SpaceSharedPtrVector<Scalar> spaces;
      int spaces_size;

      /// Dirichlet lift rhs part.
      Vector<Scalar>* dirichlet_lift_rhs;

      /// Internal.
      bool nonlinear, add_dirichlet_lift;
      
      /// DiscreteProblemMatrixVector methods.
      bool set_matrix(SparseMatrix<Scalar>* mat);
      bool set_rhs(Vector<Scalar>* rhs);
      void invalidate_matrix();
      
      /// Assembly data.
      DiscreteProblemThreadAssembler<Scalar>** threadAssembler;

      /// Select the right things to assemble
      DiscreteProblemSelectiveAssembler<Scalar> selectiveAssembler;

      template<typename T> friend class Solver;
      template<typename T> friend class LinearSolver;
      template<typename T, typename S> friend class AdaptSolver;
      template<typename T> friend class NonlinearSolver;
      template<typename T> friend class NewtonSolver;
      template<typename T> friend class PicardSolver;
      template<typename T> friend class RungeKutta;
      template<typename T> friend class KellyTypeAdapt;
    };
  }
}
#endif
