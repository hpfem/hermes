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
    class HERMES_API AdaptSolverCriterion
    {
    public:
      AdaptSolverCriterion();
      virtual bool done(double error, int iteration) = 0;
    };

    class HERMES_API AdaptSolverCriterionErrorThreshold : public AdaptSolverCriterion
    {
    public:
      AdaptSolverCriterionErrorThreshold(double error_threshold);
      virtual bool done(double error, int iteration);
      double error_threshold;
    };

    class HERMES_API AdaptSolverCriterionFixed : public AdaptSolverCriterion
    {
    public:
      AdaptSolverCriterionFixed(int refinement_levels);
      virtual bool done(double error, int iteration);
      int refinement_levels;
    };

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
      AdaptSolver(SpaceSharedPtrVector<Scalar> initial_spaces, WeakFormSharedPtr<Scalar> wf, ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* stopping_criterion_single_step, RefinementSelectors::SelectorVector<Scalar> selectors, AdaptSolverCriterion* stopping_criterion_global);
      AdaptSolver(SpaceSharedPtr<Scalar> initial_space, WeakFormSharedPtr<Scalar> wf, ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* stopping_criterion_single_step, RefinementSelectors::Selector<Scalar>* selector, AdaptSolverCriterion* stopping_criterion_global);

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

      /// Switch visualization on / off.
      void switch_visualization(bool on_off);

      /// Setters.
      void set_initial_spaces(SpaceSharedPtrVector<Scalar>);
      void set_wf(WeakFormSharedPtr<Scalar>);
      void set_error_calculator(ErrorCalculator<Scalar>*);
      void set_stopping_criterion_single_step(AdaptivityStoppingCriterion<Scalar>*);
      void set_selectors(RefinementSelectors::SelectorVector<Scalar>);
      void set_stopping_criterion_global(AdaptSolverCriterion* stopping_criterion_global);

      /// Getters.
      SpaceSharedPtrVector<Scalar> get_initial_spaces();
      WeakFormSharedPtr<Scalar> get_wf();
      ErrorCalculator<Scalar>* get_error_calculator();
      AdaptivityStoppingCriterion<Scalar>* get_stopping_criterion_single_step();
      RefinementSelectors::SelectorVector<Scalar> get_selectors();
      AdaptSolverCriterion* get_stopping_criterion_global();

    private:
      
      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "AdaptSolver"; }

      /// Initialize data for one solution (one adaptivity loop)
      void init_solving();
      
      /// De-initialize data for one solution (one adaptivity loop)
      void deinit_solving();

      /// Fill the array element_ids_to_reassemble.
      /// This is in fact the main method responsible of any re-use logic:
      /// - not only it fills the array with elements that were directly changed in adaptivity (that is easy)
      /// - it has to also correctly identify all other elements that need to be reassembled (that is not easy).
      void mark_elements_to_reassemble();

      /// (Optional) visualization.
      void visualize(std::vector<SpaceSharedPtr<Scalar> >& ref_spaces);

      /// The stopping_criterion_global for the loop to stop - the quantity measured is the total error as measured by the provided instance of error calculator.
      AdaptSolverCriterion* stopping_criterion_global;

      /// Information if the solve method is running.
      /// Used for banning descendants of this class to perform some actions 
      bool solve_method_running;

      /// Internal structures - Spaces.
      /// This class changes this instance during the solve() method: yes.
      /// Can user change this during adaptation: no.
      SpaceSharedPtrVector<Scalar> spaces;

      /// This is to hold the ref_spaces for matrix reuse.
      SpaceSharedPtrVector<Scalar> ref_spaces;
      SpaceSharedPtrVector<Scalar> prev_ref_spaces;
      
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
      AdaptivityStoppingCriterion<Scalar>* stopping_criterion_single_step;

      /// Internal structures - Stopping criterion for each refinement step.
      /// This class changes this instance during the solve() method: no.
      /// Can user change this during adaptation: [will be used from the following step onwards].
      RefinementSelectors::SelectorVector<Scalar> selectors;
      
      /// The solution being returned on demand.
      MeshFunctionSharedPtrVector<Scalar> ref_slns;
      MeshFunctionSharedPtrVector<Scalar> slns;

      /// Views - used only if visualization is ON.
      std::vector<Views::ScalarView*> scalar_views;
      std::vector<Views::OrderView*> order_viewsRef;
      std::vector<Views::OrderView*> order_views;
      std::vector<Views::BaseView<Scalar>*> base_views;

      /// Strictly private - Adapt instance.
      Adapt<Scalar>* adaptivity_internal;
      
      /// Strictly private - solver.
      SolverType* solver;

      /// Strictly private - counter.
      int adaptivity_step;

      /// Strictly private - elements to reassemble.
      /// Internal data: std::pair: [0] - element id, [1] - component (for multimesh).
      std::set<std::pair<int, unsigned char> > elements_to_reassemble;
      std::set<std::pair<int, unsigned char> > DOFs_to_reassemble;

      CSCMatrix<Scalar>* prev_mat;
      Vector<Scalar>* prev_rhs;

      /// use Hermes views to display stuff.
      bool visualization;

      /// For info only.
      int total_elements_prev_spaces;
    };
  }
}
#endif