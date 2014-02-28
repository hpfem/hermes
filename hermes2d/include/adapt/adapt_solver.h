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
      AdaptSolver(SpaceSharedPtrVector<Scalar> initial_spaces, WeakFormSharedPtr<Scalar> wf, ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* strategy, RefinementSelectors::SelectorVector<Scalar> selectors);
      AdaptSolver(SpaceSharedPtr<Scalar> initial_space, WeakFormSharedPtr<Scalar> wf, ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* strategy, RefinementSelectors::Selector<Scalar> selector);

      /// Setters.
      void set_initial_spaces(SpaceSharedPtrVector<Scalar>);
      void set_wf(WeakFormSharedPtr<Scalar>);
      void set_error_calculator(ErrorCalculator<Scalar>*);
      void set_strategy(AdaptivityStoppingCriterion<Scalar>*);
      void set_adaptivity_internal(Adapt<Scalar>);
      void set_selectors(RefinementSelectors::SelectorVector<Scalar>);

      /// Getters.
      SpaceSharedPtrVector<Scalar> get_initial_spaces();
      WeakFormSharedPtr<Scalar> get_wf();
      ErrorCalculator<Scalar>* get_error_calculator();
      AdaptivityStoppingCriterion<Scalar>* get_strategy();
      Adapt<Scalar> get_adaptivity_internal();
      RefinementSelectors::SelectorVector<Scalar> get_selectors();


    private:
      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "AdaptSolver"; }

      SpaceSharedPtrVector<Scalar> initial_spaces;
      WeakFormSharedPtr<Scalar> wf;
      ErrorCalculator<Scalar>* error_calculator;
      AdaptivityStoppingCriterion<Scalar>* strategy;
      Adapt<Scalar> adaptivity_internal;
      RefinementSelectors::SelectorVector<Scalar> selectors;
      SolverType* solver;
    };
  }
}
#endif