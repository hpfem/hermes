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
#include <unordered_set>

namespace Hermes
{
  namespace Hermes2D
  {
    class HERMES_API AdaptSolverCriterion
    {
    public:
      AdaptSolverCriterion();
      virtual bool done(double error, unsigned short iteration) = 0;
    };

    class HERMES_API AdaptSolverCriterionErrorThreshold : public AdaptSolverCriterion
    {
    public:
      AdaptSolverCriterionErrorThreshold(double error_threshold);
      virtual bool done(double error, unsigned short iteration);
      double error_threshold;
    };

    class HERMES_API AdaptSolverCriterionFixed : public AdaptSolverCriterion
    {
    public:
      AdaptSolverCriterionFixed(unsigned short refinement_levels);
      virtual bool done(double error, unsigned short iteration);
      unsigned short refinement_levels;
    };

    /// h-, p-, hp-adaptivity
    /// Also influences the selectors - i.e. for h-, or p- adaptivity 
    enum AdaptivityType
    {
      pAdaptivity,
      hAdaptivity,
      hpAdaptivity
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
      AdaptSolver(std::vector<SpaceSharedPtr<Scalar> > initial_spaces, WeakFormSharedPtr<Scalar> wf, ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* stopping_criterion_single_step, std::vector<RefinementSelectors::Selector<Scalar>*> selectors, AdaptSolverCriterion* stopping_criterion_global);
      AdaptSolver(SpaceSharedPtr<Scalar> initial_space, WeakFormSharedPtr<Scalar> wf, ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* stopping_criterion_single_step, RefinementSelectors::Selector<Scalar>* selector, AdaptSolverCriterion* stopping_criterion_global);

      /// Common code for the constructors.
      void init();

      /// Destruct this instance.
      ~AdaptSolver();

      /// The main method - solve.
      void solve(AdaptivityType adaptivityType);

      /// Get the solutions.
      std::vector<MeshFunctionSharedPtr<Scalar> > get_slns();

      /// Get i-th solution.
      MeshFunctionSharedPtr<Scalar> get_sln(int index);

      /// Get the solutions.
      std::vector<MeshFunctionSharedPtr<Scalar> > get_ref_slns();

      /// Get i-th solution.
      MeshFunctionSharedPtr<Scalar> get_ref_sln(int index);

      /// Switch visualization on / off.
      void switch_visualization(bool on_off);

      /// Add exact solutions for exact solver calculation.
      void set_exact_solution(MeshFunctionSharedPtr<Scalar> exact_sln);
      void set_exact_solutions(std::vector<MeshFunctionSharedPtr<Scalar> > exact_slns);

      /// Setters.
      void set_initial_spaces(std::vector<SpaceSharedPtr<Scalar> >);
      void set_wf(WeakFormSharedPtr<Scalar>);
      void set_error_calculator(ErrorCalculator<Scalar>*);
      void set_stopping_criterion_single_step(AdaptivityStoppingCriterion<Scalar>*);
      void set_selectors(std::vector<RefinementSelectors::Selector<Scalar>*>);
      void set_stopping_criterion_global(AdaptSolverCriterion* stopping_criterion_global);

      /// Getters.
      SolverType* get_solver();
      std::vector<SpaceSharedPtr<Scalar> > get_initial_spaces();
      WeakFormSharedPtr<Scalar> get_wf();
      ErrorCalculator<Scalar>* get_error_calculator();
      AdaptivityStoppingCriterion<Scalar>* get_stopping_criterion_single_step();
      std::vector<RefinementSelectors::Selector<Scalar>*> get_selectors();
      AdaptSolverCriterion* get_stopping_criterion_global();

      /// See Hermes::Mixins::Loggable.
      virtual void set_verbose_output(bool to_set);

      /// Views settings.
      static unsigned short view_size;
      static bool scalar_views_switch;
      static bool order_views_switch;
      static bool base_views_switch;
      static bool wait_on_show;

    private:

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "AdaptSolver"; }

      /// Initialize data for one solution (one adaptivity loop)
      void init_solving(AdaptivityType adaptivityType);

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
      std::vector<SpaceSharedPtr<Scalar> > spaces;

      /// This is to hold the ref_spaces for matrix reuse.
      std::vector<SpaceSharedPtr<Scalar> > ref_spaces;
      std::vector<SpaceSharedPtr<Scalar> > prev_ref_spaces;

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
      /// Can user change this during adaptation: yes [will be used from the following step onwards].
      AdaptivityStoppingCriterion<Scalar>* stopping_criterion_single_step;

      /// Internal structures - Stopping criterion for each refinement step.
      /// This class changes this instance during the solve() method: no.
      /// Can user change this during adaptation: yes [will be used from the following step onwards].
      std::vector<RefinementSelectors::Selector<Scalar>*> selectors;

      /// The solution being returned on demand.
      std::vector<MeshFunctionSharedPtr<Scalar> > ref_slns;
      std::vector<MeshFunctionSharedPtr<Scalar> > slns;
      std::vector<MeshFunctionSharedPtr<Scalar> > exact_slns;

      /// Views - used only if visualization is ON.
      std::vector<Views::ScalarView*> scalar_views;
      std::vector<Views::OrderView*> order_views;
      std::vector<Views::BaseView<Scalar>*> base_views;

      /// Strictly private - Adapt instance.
      Adapt<Scalar>* adaptivity_internal;

      /// Strictly private - solver.
      SolverType* solver;

      /// Strictly private - adaptivity steps counter.
      unsigned short adaptivity_step;

      /// Strictly private - elements to reassemble.
      /// Internal data: std::pair: [0] - element id, [1] - component (for multimesh).
      std::unordered_set<unsigned int> elements_to_reassemble[H2D_MAX_COMPONENTS];
      std::unordered_set<int> DOFs_to_reassemble[H2D_MAX_COMPONENTS];

      /// Previous algebraic structures for reusal - matrix.
      CSCMatrix<Scalar>* prev_mat;
      /// Previous algebraic structures for reusal - rhs.
      Vector<Scalar>* prev_rhs, *prev_dirichlet_lift_rhs;

      /// Use Hermes views to display stuff.
      bool visualization;

      /// For info only.
      unsigned int total_elements_prev_spaces;

      /// Utility simple selectors.
      std::vector<RefinementSelectors::Selector<Scalar>*> hOnlySelectors;
      std::vector<RefinementSelectors::Selector<Scalar>*> pOnlySelectors;

      /// Strictly private - size of the system (this dimension is what the size of all other data structures is compared to).
      unsigned char number_of_equations;
    };
  }
}
#endif