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

#include "adapt_solver.h"
#include "solver/linear_solver.h"
#include "solver/newton_solver.h"
#include "solver/picard_solver.h"
#include "ogprojection.h"
#include "views/scalar_view.h"
#include "views/base_view.h"
#include "views/order_view.h"

namespace Hermes
{
  namespace Hermes2D
  {
    AdaptSolverCriterion::AdaptSolverCriterion()
    {
    }

    AdaptSolverCriterionErrorThreshold::AdaptSolverCriterionErrorThreshold(double error_tolerance) : AdaptSolverCriterion()
    {
      this->error_threshold = error_threshold;
    }

    bool AdaptSolverCriterionErrorThreshold::done(double error, int iteration)
    {
      return error < this->error_threshold;
    }

    AdaptSolverCriterionFixed::AdaptSolverCriterionFixed(int refinement_levels) : AdaptSolverCriterion()
    {
      this->refinement_levels = refinement_levels;
    }

    bool AdaptSolverCriterionFixed::done(double error, int iteration)
    {
      return iteration >= this->refinement_levels;
    }

    template<typename Scalar, typename SolverType>
    AdaptSolver<Scalar, SolverType>::AdaptSolver(SpaceSharedPtrVector<Scalar> initial_spaces, WeakFormSharedPtr<Scalar> wf, ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* stopping_criterion_single_step, RefinementSelectors::SelectorVector<Scalar> selectors, AdaptSolverCriterion* stopping_criterion_global)
      : spaces(initial_spaces), wf(wf), error_calculator(error_calculator), stopping_criterion_single_step(stopping_criterion_single_step), selectors(selectors), stopping_criterion_global(stopping_criterion_global)
    {
      this->init();
    }

    template<typename Scalar, typename SolverType>
    AdaptSolver<Scalar, SolverType>::AdaptSolver(SpaceSharedPtr<Scalar> initial_space, WeakFormSharedPtr<Scalar> wf, ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* stopping_criterion_single_step, RefinementSelectors::Selector<Scalar>* selector, AdaptSolverCriterion* stopping_criterion_global)
      : wf(wf), stopping_criterion_single_step(stopping_criterion_single_step), error_calculator(error_calculator), stopping_criterion_global(stopping_criterion_global)
    {
      this->spaces.push_back(initial_space);
      this->selectors.push_back(selector);

      this->init();
    }

    static int currentAdaptSolverIteration;
    static int currentAdaptSolverRefinedElementsCount;
    static std::vector<int>* currentAdaptSolverRefinedElementsCount;

    void get_states_to_reassemble(Traverse::State**& states, int& num_states)
    {

    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::init()
    {
      this->solver = new SolverType(wf, spaces);
      this->solver->dp->refine_states = &get_states_to_reassemble;
      this->solver->output_matrix();
      this->solver->set_matrix_export_format(EXPORT_FORMAT_MATLAB_SIMPLE);
      this->adaptivity_internal = new Adapt<Scalar>(spaces, error_calculator, stopping_criterion_single_step);
      this->solve_method_running = false;
      this->visualization = false;
    }

    template<typename Scalar, typename SolverType>
    AdaptSolver<Scalar, SolverType>:: ~AdaptSolver()
    {
      delete this->solver;
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::switch_visualization(bool on_off)
    {
      this->visualization = on_off;
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::init_solving()
    {
      this->adaptivity_step = 1;
      this->solve_method_running = true;
      ref_slns.clear();
      slns.clear();
      for (int i = 0; i < this->spaces.size(); i++)
      {
        ref_slns.push_back(MeshFunctionSharedPtr<Scalar>(new Solution<Scalar>));
        slns.push_back(MeshFunctionSharedPtr<Scalar>(new Solution<Scalar>));
        if (this->visualization)
        {
          this->scalar_views.push_back(new Views::ScalarView("", new Views::WinGeom(i * 410, 10, 400, 300)));
          this->scalar_views.back()->set_title("Reference solution #%i", i);

          this->base_views.push_back(new Views::BaseView<Scalar>("", new Views::WinGeom(i * 410, 340, 400, 300)));
          this->base_views.back()->set_title("Reference space #%i - basis", i);

          this->order_views.push_back(new Views::OrderView("", new Views::WinGeom(i * 410, 670, 400, 300)));
          this->order_views.back()->set_title("Reference space #%i - orders", i);
        }
      }

      this->element_ids_to_reassemble.clear();
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::deinit_solving()
    {
      this->solve_method_running = false;

      if (this->visualization)
      for (int i = 0; i < this->spaces.size(); i++)
      {
        delete this->scalar_views[i];
        delete this->base_views[i];
        delete this->order_views[i];
      }
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::solve()
    {
      this->init_solving();
      do
      {
        this->info("AdaptSolver step %d:", this->adaptivity_step);

        // Construct globally refined reference meshes and setup reference spaces.
        std::vector<SpaceSharedPtr<Scalar> > ref_spaces;
        for (int i = 0; i < this->spaces.size(); i++)
        {
          Mesh::ReferenceMeshCreator ref_mesh_creator(spaces[i]->get_mesh());
          MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
          Space<Scalar>::ReferenceSpaceCreator u_ref_space_creator(spaces[i], ref_mesh);
          ref_spaces.push_back(u_ref_space_creator.create_ref_space());
        }

        // Set these to the solver.
        this->solver->set_spaces(ref_spaces);

        // Initialize reference problem.
        this->info("Solving on reference mesh.");


        // Perform solution.
        this->solver->solve();

        // Translate the resulting coefficient vector into the instance of Solution.
        Solution<Scalar>::vector_to_solutions(solver->get_sln_vector(), ref_spaces, ref_slns);

        if (this->visualization)
          this->visualize(ref_spaces);

        // Project the fine mesh solution onto the coarse mesh.
        this->info("Projecting reference solution on coarse mesh.");
        OGProjection<Scalar>::project_global(spaces, ref_slns, slns);

        // Calculate element errors.
        this->info("Calculating error estimate.");
        this->error_calculator->calculate_errors(slns, ref_slns, true);
        double error = this->error_calculator->get_total_error_squared() * 100;
        this->info("The error estimate: %f.", error);

        // If err_est too large, adapt the mesh.
        if (this->stopping_criterion_global->done(error, this->adaptivity_step))
        {
          this->deinit_solving();
          return;
        }
        else
        {
          this->info("Adapting coarse mesh.");
          this->adaptivity_internal->adapt(this->selectors);
          this->mark_elements_to_reassemble();
        }

        this->adaptivity_step++;
      } while (true);

      this->deinit_solving();
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::get_states_to_reassemble(Traverse::State**& states, int& num_states)
    {
      this->element_ids_to_reassemble.clear();
      for (int i = 0; i < this->adaptivity_internal->elements_to_refine_count; i++)
      {
        Element* current_element = current_state->e[this->elements_to_refine[i].comp];
        if (this->elements_to_refine[i].id == current_element->id || (current_element->parent && this->elements_to_refine[i].id == current_element->parent->id))
          info_array[current_state->e[this->num]->id] = this->elements_to_refine[i].split;
      }
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::mark_elements_to_reassemble()
    {
      this->element_ids_to_reassemble.clear();
      for (int i = 0; i < this->adaptivity_internal->elements_to_refine_count; i++)
      {
        Element* current_element = current_state->e[this->elements_to_refine[i].comp];
        if (this->elements_to_refine[i].id == current_element->id || (current_element->parent && this->elements_to_refine[i].id == current_element->parent->id))
          info_array[current_state->e[this->num]->id] = this->elements_to_refine[i].split;
      }
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::visualize(std::vector<SpaceSharedPtr<Scalar> >& ref_spaces)
    {
      for (int i = 0; i < this->spaces.size(); i++)
      {
        this->scalar_views[i]->show(this->ref_slns[i]);
        this->base_views[i]->show(ref_spaces[i]);
        this->order_views[i]->show(ref_spaces[i]);
      }

      Views::View::wait_for_keypress();
    }

    template<typename Scalar, typename SolverType>
    MeshFunctionSharedPtrVector<Scalar> AdaptSolver<Scalar, SolverType>::get_slns()
    {
      if (this->solve_method_running)
        throw Exceptions::Exception("AdaptSolver asked for solutions while it was running.");
      return this->slns;
    }

    template<typename Scalar, typename SolverType>
    MeshFunctionSharedPtr<Scalar> AdaptSolver<Scalar, SolverType>::get_sln(int index)
    {
      if (this->solve_method_running)
        throw Exceptions::Exception("AdaptSolver asked for a solution while it was running.");
      return this->slns[index];
    }

    template<typename Scalar, typename SolverType>
    MeshFunctionSharedPtrVector<Scalar> AdaptSolver<Scalar, SolverType>::get_ref_slns()
    {
      if (this->solve_method_running)
        throw Exceptions::Exception("AdaptSolver asked for solutions while it was running.");
      return this->ref_slns;
    }

    template<typename Scalar, typename SolverType>
    MeshFunctionSharedPtr<Scalar> AdaptSolver<Scalar, SolverType>::get_ref_sln(int index)
    {
      if (this->solve_method_running)
        throw Exceptions::Exception("AdaptSolver asked for a solution while it was running.");
      return this->ref_slns[index];
    }

    template<typename Scalar, typename SolverType>
    AdaptSolverCriterion* AdaptSolver<Scalar, SolverType>::get_stopping_criterion_global()
    {
      return this->stopping_criterion_global;
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::set_stopping_criterion_global(AdaptSolverCriterion* stopping_criterion_global)
    {
      this->stopping_criterion_global = stopping_criterion_global;
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::set_initial_spaces(SpaceSharedPtrVector<Scalar> spaces)
    {
      if (this->solve_method_running)
        throw Exceptions::Exception("AdaptSolver asked to change the initial spaces while it was running.");
      this->spaces = spaces;
      this->adaptivity_internal->set_spaces(this->spaces);
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::set_wf(WeakFormSharedPtr<Scalar> wf)
    {
      if (this->solve_method_running)
        throw Exceptions::Exception("AdaptSolver asked to change the weak formulation while it was running.");
      this->wf = wf;
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::set_error_calculator(ErrorCalculator<Scalar>* error_calculator)
    {
      this->error_calculator = error_calculator;
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::set_stopping_criterion_single_step(AdaptivityStoppingCriterion<Scalar>* stopping_criterion_single_step)
    {
      this->stopping_criterion_single_step = stopping_criterion_single_step;
      this->adaptivity_internal->set_strategy(stopping_criterion_single_step);
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::set_selectors(RefinementSelectors::SelectorVector<Scalar> selectors)
    {
      this->selectors = selectors;
    }

    template<typename Scalar, typename SolverType>
    SpaceSharedPtrVector<Scalar> AdaptSolver<Scalar, SolverType>::get_initial_spaces()
    {
      return this->spaces;
    }

    template<typename Scalar, typename SolverType>
    WeakFormSharedPtr<Scalar> AdaptSolver<Scalar, SolverType>::get_wf()
    {
      return this->wf;
    }

    template<typename Scalar, typename SolverType>
    ErrorCalculator<Scalar>* AdaptSolver<Scalar, SolverType>::get_error_calculator()
    {
      return this->error_calculator;
    }

    template<typename Scalar, typename SolverType>
    AdaptivityStoppingCriterion<Scalar>* AdaptSolver<Scalar, SolverType>::get_stopping_criterion_single_step()
    {
      return this->stopping_criterion_single_step;
    }

    template<typename Scalar, typename SolverType>
    RefinementSelectors::SelectorVector<Scalar> AdaptSolver<Scalar, SolverType>::get_selectors()
    {
      return this->selectors;
    }

    template<typename Scalar, typename SolverType>
    bool AdaptSolver<Scalar, SolverType>::isOkay() const
    {
      this->adaptivity_internal->check();
      for (int i = 0; i < this->spaces.size(); i++)
        this->spaces[i]->check();

      return true;
    }

    template HERMES_API class AdaptSolver<double, LinearSolver<double> >;
    template HERMES_API class AdaptSolver<std::complex<double>, LinearSolver<std::complex<double > > >;

    template HERMES_API class AdaptSolver<double, NewtonSolver<double> >;
    template HERMES_API class AdaptSolver<std::complex<double>, NewtonSolver<std::complex<double > > >;

    template HERMES_API class AdaptSolver<double, PicardSolver<double> >;
    template HERMES_API class AdaptSolver<std::complex<double>, PicardSolver<std::complex<double > > >;
  }
}