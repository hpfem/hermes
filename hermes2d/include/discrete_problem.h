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
    class HERMES_API DiscreteProblem : 
      public Hermes::Mixins::TimeMeasurable, 
      public Hermes::Hermes2D::Mixins::SettableSpaces<Scalar>, 
      public Hermes::Hermes2D::Mixins::StateQueryable,
      public Hermes::Hermes2D::Mixins::DiscreteProblemRungeKutta<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemWeakForm<Scalar>,
      public Hermes::Hermes2D::Mixins::DiscreteProblemCacheSettings,
      public Hermes::Mixins::IntegrableWithGlobalOrder,
      public Hermes::Hermes2D::Mixins::DiscreteProblemMatrixVector<Scalar>,
      public Hermes::Hermes2D::Mixins::Parallel
    {
    public:
      /// Constructor for multiple components / equations.
      DiscreteProblem(WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      /// Constructor for one equation.
      DiscreteProblem(WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar>& space);
      /// Non-parameterized constructor.
      DiscreteProblem();
      /// Destuctor.
      virtual ~DiscreteProblem();

      /// Make this DiscreteProblem linear.
      /// Does 2 things
      /// 1 - turns off initialization of previous iterations for nonlinear solvers.
      /// 2 - allows for assembling Dirichlet boundary conditions using a Dirichlet lift.
      /// \param[in] dirichlet_lift_accordingly If true, the appropriate settings for (linear / nonlinear)
      /// problem will be used (use Dirichlet lift iff the problem is linear). If false, the other way round.
      void set_linear(bool to_set = true, bool dirichlet_lift_accordingly = true);

      /// Assembling.
      void assemble(Scalar* coeff_vec, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL);
      /// Assembling.
      /// Without the matrix.
      void assemble(Scalar* coeff_vec, Vector<Scalar>* rhs = NULL);
      /// Light version passing NULL for the coefficient vector. External solutions
      /// are initialized with zeros.
      void assemble(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs = NULL);
      /// Light version passing NULL for the coefficient vector. External solutions
      /// are initialized with zeros.
      /// Without the matrix.
      void assemble(Vector<Scalar>* rhs);
      /// Assembling.
      void assemble(Solution<Scalar>** u_ext_sln, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs);

      /// set time information for time-dependent problems.
      void set_time(double time);
      void set_time_step(double time_step);

      /// Sets new spaces for the instance.
      void set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spaces);
      void set_space(SpaceSharedPtr<Scalar>& space);

      /// Set the weak forms.
      void set_weak_formulation(WeakForm<Scalar>* wf);

      /// Get all spaces as a Hermes::vector.
      Hermes::vector<SpaceSharedPtr<Scalar> >& get_spaces();

      /// If the cache should not be used for any reason.
      void set_do_not_use_cache(bool to_set = true);
      
      /// Reports cache hits and misses.
      void set_report_cache_hits_and_misses(bool to_set = true);

    protected:
      /// Initialize states.
      void init_assembling(Traverse::State**& states, int& num_states, Solution<Scalar>** u_ext_sln, Hermes::vector<MeshSharedPtr>& meshes);
      void deinit_assembling(Traverse::State** states, int num_states);

      /// RungeKutta helpers.
      void set_RK(int original_spaces_count, bool force_diagonal_blocks = NULL, Table* block_weights = NULL);

      /// State querying helpers.
      bool isOkay() const;
      inline std::string getClassName() const { return "DiscreteProblem"; }

      /// Init function. Common code for the constructors.
      void init();

      /// Space instances for all equations in the system.
      Hermes::vector<SpaceSharedPtr<Scalar> > spaces;
      int spaces_size;
      
      /// Internal.
      bool nonlinear, add_dirichlet_lift;
      
      /// DiscreteProblemMatrixVector methods.
      void set_matrix(SparseMatrix<Scalar>* mat);
      void set_rhs(Vector<Scalar>* rhs);
      void invalidate_matrix();

      /// The cache.
      DiscreteProblemCache<Scalar> cache;
      
      /// Assembly data.
      DiscreteProblemThreadAssembler<Scalar>** threadAssembler;

      /// Select the right things to assemble
      DiscreteProblemSelectiveAssembler<Scalar> selectiveAssembler;

      template<typename T> friend class Solver;
      template<typename T> friend class LinearSolver;
      template<typename T> friend class NonlinearSolver;
      template<typename T> friend class NewtonSolver;
      template<typename T> friend class PicardSolver;
      template<typename T> friend class RungeKutta;
      template<typename T> friend class KellyTypeAdapt;
    };
  }
}
#endif
