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

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::init()
    {
      this->solve_method_running = false;
      this->visualization = false;
      this->prev_mat = nullptr;
      this->prev_rhs = nullptr;
    }

    template<typename Scalar, typename SolverType>
    AdaptSolver<Scalar, SolverType>:: ~AdaptSolver()
    {
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::switch_visualization(bool on_off)
    {
      this->visualization = on_off;
    }

    template<typename Scalar>
    class StateReassemblyHelper
    {
    public:
      static int current_iteration;
      static std::vector<std::pair<int, int> >* current_elements_to_reassemble;
      static SpaceSharedPtrVector<Scalar>* current_ref_spaces;
      static SpaceSharedPtrVector<Scalar>* current_prev_ref_spaces;
      static CSCMatrix<Scalar>* prev_mat;
      static Scalar* prev_rhs;
    };

    template<typename Scalar>
    int StateReassemblyHelper<Scalar>::current_iteration;

    template<typename Scalar>
    std::vector<std::pair<int, int> >* StateReassemblyHelper<Scalar>::current_elements_to_reassemble;

    template<typename Scalar>
    SpaceSharedPtrVector<Scalar>* StateReassemblyHelper<Scalar>::current_ref_spaces;

    template<typename Scalar>
    SpaceSharedPtrVector<Scalar>* StateReassemblyHelper<Scalar>::current_prev_ref_spaces;

    template<typename Scalar>
    CSCMatrix<Scalar>* StateReassemblyHelper<Scalar>::prev_mat;

    template<typename Scalar>
    Scalar* StateReassemblyHelper<Scalar>::prev_rhs;

    template<typename Scalar>
    void get_states_to_reassemble(Traverse::State**& states, int& num_states, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs)
    {
      if (StateReassemblyHelper<Scalar>::current_iteration == 1)
        return;

      int spaces_size = StateReassemblyHelper<Scalar>::current_ref_spaces->size();
      int elements_to_reassemble_size = StateReassemblyHelper<Scalar>::current_elements_to_reassemble->size();
      std::vector<std::pair<int, int> > newSpace_elements_to_reassemble;

      // DoF to DOF map for the elements we can reuse. -1s are there to distinguish those we cannot.
      int prev_system_size = StateReassemblyHelper<Scalar>::prev_mat->get_size();
      assert(Space<Scalar>::get_num_dofs(*StateReassemblyHelper<Scalar>::current_prev_ref_spaces) == prev_system_size);
      int* DOF_to_DOF_map = new int[prev_system_size];
      for (int i = 0; i < prev_system_size; i++)
        DOF_to_DOF_map[i] = -1;

      // Utility assembly lists.
      AsmList<Scalar> al, al_prev;

      // First, identify the elements on the new reference spaces.
      Traverse trav(spaces_size * 2);
      MeshFunctionSharedPtrVector<Scalar> dummy_fns;
      for (int i = 0; i < spaces_size; i++)
        dummy_fns.push_back(new ZeroSolution<Scalar>(StateReassemblyHelper<Scalar>::current_prev_ref_spaces->at(i)->get_mesh()));
      for (int i = 0; i < spaces_size; i++)
        dummy_fns.push_back(new ZeroSolution<Scalar>(StateReassemblyHelper<Scalar>::current_ref_spaces->at(i)->get_mesh()));

      int num_states_local;
      Traverse::State** states_local = trav.get_states(dummy_fns, num_states_local);
      for (int local_state_i = 0; local_state_i < num_states_local; local_state_i++)
      {
        bool added = false;
        Traverse::State* current_local_state = states_local[local_state_i];
        for (int elem_i = 0; elem_i < elements_to_reassemble_size; elem_i++)
        {
          for (int space_i = 0; space_i < spaces_size; space_i++)
          {
            if (StateReassemblyHelper<Scalar>::current_elements_to_reassemble->at(elem_i).second == space_i)
            {
              if (StateReassemblyHelper<Scalar>::current_elements_to_reassemble->at(elem_i).first == current_local_state->e[space_i]->id)
              {
                newSpace_elements_to_reassemble.push_back(std::pair<int, int>(current_local_state->e[space_i + spaces_size]->id, space_i));
                added = true;
                break;
              }
            }
          }
          if (added)
            break;
        }
        if (!added)
        {
          for (int space_i = 0; space_i < spaces_size; space_i++)
          {
            StateReassemblyHelper<Scalar>::current_prev_ref_spaces->at(space_i)->get_element_assembly_list(current_local_state->e[space_i], &al_prev);
            StateReassemblyHelper<Scalar>::current_ref_spaces->at(space_i)->get_element_assembly_list(current_local_state->e[space_i + spaces_size], &al);
            for (short dof_i = 0; dof_i < al_prev.cnt; dof_i++)
            {
              assert(al_prev.idx[dof_i] == al.idx[dof_i]);
              DOF_to_DOF_map[al_prev.dof[dof_i]] = al.dof[dof_i];
            }
          }
        }
      }

      // Using the changed element on the new reference space, select the states from the original states that need to be recalculated.
      int newSpace_elements_to_reassemble_size = newSpace_elements_to_reassemble.size();
      Traverse::State** new_states = malloc_with_check<Traverse::State*>(num_states, true);
      int new_num_states = 0;
      for (int state_i = 0; state_i < num_states; state_i++)
      {
        for (int space_i = 0; space_i < spaces_size; space_i++)
        {
          bool added = false;
          for (int to_reassemble_i = 0; to_reassemble_i < newSpace_elements_to_reassemble_size; to_reassemble_i++)
          {
            if (space_i == newSpace_elements_to_reassemble.at(to_reassemble_i).second && states[state_i]->e[space_i]->id == newSpace_elements_to_reassemble.at(to_reassemble_i).first)
            {
              new_states[new_num_states++] = Traverse::State::clone(states[state_i]);
              added = true;
              break;
            }
            if (added)
              break;
          }
        }
      }
      free_with_check(states);
      new_states = realloc_with_check<Traverse::State*>(new_states, new_num_states);
      states = new_states;
      num_states = new_num_states;

      // Now we have to use the DOF to DOF map to fill in the necessary entries in the new matrix and rhs from the old ones.
      Scalar* Ax = StateReassemblyHelper<Scalar>::prev_mat->get_Ax();
      int* Ai = StateReassemblyHelper<Scalar>::prev_mat->get_Ai();
      int* Ap = StateReassemblyHelper<Scalar>::prev_mat->get_Ap();
      for (int i = 0; i < prev_system_size; i++)
      {
        if (DOF_to_DOF_map[i] != -1)
        {
          for (int j = 0; j < Ap[i + 1] - Ap[i]; j++)
          {
            if (DOF_to_DOF_map[Ai[Ap[i] + j]] != -1)
            {
              mat->add(DOF_to_DOF_map[i], DOF_to_DOF_map[Ai[Ap[i] + j]], Ax[Ap[i] + j]);
            }
          }
          rhs->add(DOF_to_DOF_map[i], StateReassemblyHelper<Scalar>::prev_rhs[i]);
        }
      }
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::init_solving()
    {
      this->adaptivity_step = 1;
      this->solve_method_running = true;
      ref_slns.clear();
      slns.clear();
      for(unsigned short i = 0; i < this->spaces.size(); i++)
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

      this->solver = new SolverType(wf, spaces);
      this->solver->output_matrix();
      this->solver->set_matrix_export_format(EXPORT_FORMAT_MATLAB_SIMPLE);
      this->solver->dp->set_reassembled_states_reuse_linear_system_fn(&get_states_to_reassemble<Scalar>);
      this->adaptivity_internal = new Adapt<Scalar>(spaces, error_calculator, stopping_criterion_single_step);

      this->elements_to_reassemble.clear();
      this->DOFs_to_reassemble.clear();
      this->ref_spaces.clear();
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::deinit_solving()
    {
      this->solve_method_running = false;

      if (this->visualization)
      for(unsigned short i = 0; i < this->spaces.size(); i++)
      {
        delete this->scalar_views[i];
        delete this->base_views[i];
        delete this->order_views[i];
      }

      delete this->adaptivity_internal;
      delete this->solver;
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::solve()
    {
      this->init_solving();
      do
      {
        this->info("AdaptSolver step %d:", this->adaptivity_step);

        // Construct globally refined reference meshes and setup reference spaces.
        prev_ref_spaces = ref_spaces;
        ref_spaces.clear();
        for(unsigned short i = 0; i < this->spaces.size(); i++)
        {
          Mesh::ReferenceMeshCreator ref_mesh_creator(spaces[i]->get_mesh());
          MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
          Space<Scalar>::ReferenceSpaceCreator u_ref_space_creator(spaces[i], ref_mesh);
          ref_spaces.push_back(u_ref_space_creator.create_ref_space());
        }

        // Set these to the solver.
        this->solver->set_spaces(ref_spaces);

        // Initialize the states handler.
        StateReassemblyHelper<Scalar>::current_iteration = this->adaptivity_step;
        StateReassemblyHelper<Scalar>::current_ref_spaces = &ref_spaces;
        StateReassemblyHelper<Scalar>::current_prev_ref_spaces = &prev_ref_spaces;
        StateReassemblyHelper<Scalar>::current_elements_to_reassemble = &this->elements_to_reassemble;

        // Perform solution.
        this->info("Solving on reference mesh.");
        this->solver->solve();

        // Update the stored (previous) linear system.
        if (this->prev_mat)
          delete this->prev_mat;
        this->prev_mat = (Hermes::Algebra::CSCMatrix<Scalar>*)(static_cast<Hermes::Solvers::MatrixSolver<Scalar>*>(this->solver))->get_jacobian()->duplicate();
        StateReassemblyHelper<Scalar>::prev_mat = this->prev_mat;
        if (this->prev_rhs)
          delete[] this->prev_rhs;
        this->prev_rhs = new Scalar[this->prev_mat->get_size()];
        (static_cast<Hermes::Solvers::MatrixSolver<Scalar>*>(this->solver))->get_residual()->extract(this->prev_rhs);
        StateReassemblyHelper<Scalar>::prev_rhs = this->prev_rhs;

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
    void AdaptSolver<Scalar, SolverType>::mark_elements_to_reassemble()
    {
      this->elements_to_reassemble.clear();
      this->DOFs_to_reassemble.clear();

      // Identify elements that changed.
      for (int i = 0; i < this->adaptivity_internal->elements_to_refine_count; i++)
      {
        ElementToRefine* element_to_refine = &this->adaptivity_internal->elements_to_refine[i];
        if (!element_to_refine->valid)
          continue;
        RefinementType refinement_type = element_to_refine->split;
        int component = element_to_refine->comp;
        Element* e = this->ref_spaces[component]->get_mesh()->get_element_fast(element_to_refine->id);
        for (int son_i = 0; son_i < H2D_MAX_ELEMENT_SONS; son_i++)
          this->elements_to_reassemble.push_back(std::pair<int, int>(e->sons[son_i]->id, component));
      }

      // Identify DOFs of the changed elements.
      AsmList<Scalar> al;
      for(unsigned short i = 0; i < this->elements_to_reassemble.size(); i++)
      {
        int component = this->elements_to_reassemble[i].second;
        Element* e = this->ref_spaces[component]->get_mesh()->get_element_fast(this->elements_to_reassemble[i].first);
        this->ref_spaces[component]->get_element_assembly_list(e, &al);
        for (int j = 0; j < al.cnt; j++)
          this->DOFs_to_reassemble.push_back(std::pair<int, int>(al.dof[j], component));
      }

      // Take a look at other elements if they share a DOF that changed.
      /// \todo This is ineffective, a more effective way is to employ an improved neighbor searching.
      for(unsigned short space_i = 0; space_i < this->ref_spaces.size(); space_i++)
      {
        for_all_active_elements_fast(this->ref_spaces[space_i]->get_mesh())
        {
          bool found = false;
          this->ref_spaces[space_i]->get_element_assembly_list(e, &al);
          for(unsigned short i = 0; i < this->DOFs_to_reassemble.size(); i++)
          {
            if (space_i == this->DOFs_to_reassemble[i].second)
            {
              for (int j = 0; j < al.cnt; j++)
              {
                if (al.dof[j] == this->DOFs_to_reassemble[i].first)
                {
                  this->elements_to_reassemble.push_back(std::pair<int, int>(e->id, space_i));
                  found = true;
                  break;
                }
              }
            }
            if (found)
              break;
          }
        }
      }
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::visualize(std::vector<SpaceSharedPtr<Scalar> >& ref_spaces)
    {
      for(unsigned short i = 0; i < this->spaces.size(); i++)
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
      for(unsigned short i = 0; i < this->spaces.size(); i++)
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