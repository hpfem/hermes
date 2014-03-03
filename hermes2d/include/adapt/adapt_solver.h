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

#ifndef __H2D_ADAPT_SOLVER_H
#define __H2D_ADAPT_SOLVER_H

#include "adapt.h"
#include "error_calculator.h"
#include "../solver/linear_solver.h"
#include "../solver/newton_solver.h"
#include "../solver/picard_solver.h"
#include "../refinement_selectors/selector.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// A complete adaptivity solver class handling the matrix reuse.
    /// \ingroup g_adapt
    template<typename Scalar, typename SolverType>
    class HERMES_API AdaptSolver :
      public Hermes::Mixins::TimeMeasurable,
      public Hermes::Mixins::Loggable,
      public Hermes::Mixins::StateQueryable
    {
    public:
      /// Constructor.
      AdaptSolver(SpaceSharedPtrVector<Scalar> initial_spaces, WeakFormSharedPtr<Scalar> wf, ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* strategy, RefinementSelectors::SelectorVector<Scalar> selectors, double threshold);
      AdaptSolver(SpaceSharedPtr<Scalar> initial_space, WeakFormSharedPtr<Scalar> wf, ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* strategy, RefinementSelectors::Selector<Scalar>* selector, double threshold);

      /// Common code for the constructors.
      void init();

      /// Destruct this instance.
      ~AdaptSolver();

      /// The main method - solve.
      void solve();

      /// Get the solutions.
      MeshFunctionSharedPtrVector<Scalar> get_slns();

      /// Get i-th solution.
      MeshFunctionSharedPtr<Scalar> get_sln(int index);

      /// Get the solutions.
      MeshFunctionSharedPtrVector<Scalar> get_ref_slns();

      /// Get i-th solution.
      MeshFunctionSharedPtr<Scalar> get_ref_sln(int index);

      /// Setters.
      void set_initial_spaces(SpaceSharedPtrVector<Scalar>);
      void set_wf(WeakFormSharedPtr<Scalar>);
      void set_error_calculator(ErrorCalculator<Scalar>*);
      void set_strategy(AdaptivityStoppingCriterion<Scalar>*);
      void set_selectors(RefinementSelectors::SelectorVector<Scalar>);
      void set_threshold(double threshold);

      /// Getters.
      SpaceSharedPtrVector<Scalar> get_initial_spaces();
      WeakFormSharedPtr<Scalar> get_wf();
      ErrorCalculator<Scalar>* get_error_calculator();
      AdaptivityStoppingCriterion<Scalar>* get_strategy();
      RefinementSelectors::SelectorVector<Scalar> get_selectors();
      double get_threshold();

    private:
      
      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "AdaptSolver"; }

      /// Initialize and de-initialize data for one solution (one adaptivity loop)
      void init_solving();
      void deinit_solving();

      /// The threshold for the loop to stop - the quantity measured is the total error as measured by the provided instance of error calculator.
      double threshold;

      /// Information if the solve method is running.
      /// Used for banning descendants of this class to perform some actions 
      bool solve_method_running;

      /// Internal structures - Spaces.
      /// This class changes this instance during the solve() method: yes.
      /// Can user change this during adaptation: no.
      SpaceSharedPtrVector<Scalar> spaces;
      
      /// Internal structures - Weak form.
      /// This class changes this instance during the solve() method: no.
      /// Can user change this during adaptation: no [technically why not, but makes no sense].
      WeakFormSharedPtr<Scalar> wf;

      /// Internal structures - Error calculator.
      /// This class changes this instance during the solve() method: no.
      /// Can user change this during adaptation: yes [will be used from the following step onwards].
      ErrorCalculator<Scalar>* error_calculator;

      /// Internal structures - Stopping criterion for each refinement step.
      /// This class changes this instance during the solve() method: no.
      /// Can user change this during adaptation: [will be used from the following step onwards].
      AdaptivityStoppingCriterion<Scalar>* strategy;

      /// Internal structures - Stopping criterion for each refinement step.
      /// This class changes this instance during the solve() method: no.
      /// Can user change this during adaptation: [will be used from the following step onwards].
      RefinementSelectors::SelectorVector<Scalar> selectors;
      
      /// The solution being returned on demand.
      MeshFunctionSharedPtrVector<Scalar> ref_slns;
      MeshFunctionSharedPtrVector<Scalar> slns;

      /// Strictly private - Adapt instance.
      Adapt<Scalar>* adaptivity_internal;
      
      /// Strictly private - solver.
      SolverType* solver;

      /// Strictly private - counter.
      int adaptivity_step;
    };
  }
}
#endif