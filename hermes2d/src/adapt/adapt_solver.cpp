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
      static std::set<std::pair<int, unsigned char> >* current_elements_to_reassemble;
      static std::set<std::pair<int, unsigned char> >* current_DOFs_to_reassemble;
      static SpaceSharedPtrVector<Scalar>* current_ref_spaces;
      static SpaceSharedPtrVector<Scalar>* current_prev_ref_spaces;
      static CSCMatrix<Scalar>* prev_mat;
      static Scalar* prev_rhs;
      static bool** reusable_DOFs;
    };

    template<typename Scalar>
    int StateReassemblyHelper<Scalar>::current_iteration;

    template<typename Scalar>
    std::set<std::pair<int, unsigned char> >* StateReassemblyHelper<Scalar>::current_elements_to_reassemble;

    template<typename Scalar>
    std::set<std::pair<int, unsigned char> >* StateReassemblyHelper<Scalar>::current_DOFs_to_reassemble;

    template<typename Scalar>
    SpaceSharedPtrVector<Scalar>* StateReassemblyHelper<Scalar>::current_ref_spaces;

    template<typename Scalar>
    SpaceSharedPtrVector<Scalar>* StateReassemblyHelper<Scalar>::current_prev_ref_spaces;

    template<typename Scalar>
    CSCMatrix<Scalar>* StateReassemblyHelper<Scalar>::prev_mat;

    template<typename Scalar>
    Scalar* StateReassemblyHelper<Scalar>::prev_rhs;

    template<typename Scalar>
    bool** StateReassemblyHelper<Scalar>::reusable_DOFs;

    template<typename Scalar>
    void get_states_to_reassemble(Traverse::State**& states, int& num_states, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs)
    {
      if (StateReassemblyHelper<Scalar>::current_iteration == 1)
      {
        *StateReassemblyHelper<Scalar>::reusable_DOFs = nullptr;
        return;
      }

      // Create utility data
      // - time measurement
      Hermes::Mixins::TimeMeasurable cpu_time;
      // - size of spaces
      int spaces_size = StateReassemblyHelper<Scalar>::current_ref_spaces->size();
      // - shortcuts for spaces, to speed up indirect accesses
      Space<Scalar>** prev_ref_spaces = new Space<Scalar>*[spaces_size];
      Space<Scalar>** ref_spaces = new Space<Scalar>*[spaces_size];;
      // - number of elements of all new reference spaces combined
      int total_elements_new_spaces = 0;
      // - fill the above structures
      for (int i = 0; i < spaces_size; i++)
      {
        prev_ref_spaces[i] = StateReassemblyHelper<Scalar>::current_prev_ref_spaces->at(i).get();
        ref_spaces[i] = StateReassemblyHelper<Scalar>::current_ref_spaces->at(i).get();
        total_elements_new_spaces += ref_spaces[i]->get_mesh()->get_num_active_elements();
      }
      // - elements to reassemble on the new space
      // -- we could maybe directly used the resulting states (in which we are interested), but it would be complicated.
      std::set<std::pair<int, unsigned char> > newSpace_elements_to_reassemble;
      // - number of DOFs of the previous spaces.
      int prev_ref_system_size = StateReassemblyHelper<Scalar>::prev_mat->get_size();
      int ref_system_size = StateReassemblyHelper<Scalar>::prev_mat->get_size();
      // - DOF to DOF map of those DOFs we can reuse. -1s are there to distinguish those we cannot.
      int* DOF_to_DOF_map = malloc_with_check<int>(prev_ref_system_size);
      for (int i = 0; i < prev_ref_system_size; i++)
        DOF_to_DOF_map[i] = -1;
      // - DOF to reassemble.
      *StateReassemblyHelper<Scalar>::reusable_DOFs = calloc_with_check<bool>(ref_system_size, true);
      // - utility assembly lists.
      AsmList<Scalar> al, al_prev;
      // - dummy functions for traversing the previous and current reference spaces together
      MeshFunctionSharedPtrVector<Scalar> dummy_fns;
      for (int i = 0; i < spaces_size; i++)
        dummy_fns.push_back(new ZeroSolution<Scalar>(StateReassemblyHelper<Scalar>::current_prev_ref_spaces->at(i)->get_mesh()));
      for (int i = 0; i < spaces_size; i++)
        dummy_fns.push_back(new ZeroSolution<Scalar>(StateReassemblyHelper<Scalar>::current_ref_spaces->at(i)->get_mesh()));
      // - (for reporting) count of DOFs needed to be reassembled
      int DOF_to_reassemble_count = 0;

      // Start.
      Hermes::Mixins::Loggable::Static::info("\t   Handling Reusing matrix entries on the new Ref. Space:");
      cpu_time.tick();

      // Traverse the previous and current reference Spaces at once.
      Traverse trav(spaces_size * 2);
      int num_states_local;
      Traverse::State** states_local = trav.get_states(dummy_fns, num_states_local);
      for (int local_state_i = 0; local_state_i < num_states_local; local_state_i++)
      {
        Traverse::State* current_local_state = states_local[local_state_i];
        for (int space_i = 0; space_i < spaces_size; space_i++)
        {
          // Mark the appropriate elements
          Element* prev_ref_element = current_local_state->e[space_i];
          Element* ref_element = current_local_state->e[space_i + spaces_size];
          // If the element is changed on the previous ref space, mark the appropriate element on the new ref space as to-reassemble.
          if (StateReassemblyHelper<Scalar>::current_elements_to_reassemble->find(std::pair<int, unsigned char>(prev_ref_element->id, space_i)) != StateReassemblyHelper<Scalar>::current_elements_to_reassemble->end())
            newSpace_elements_to_reassemble.insert(std::pair<int, unsigned char>(ref_element->id, space_i));
          // Get assembly lists.
          prev_ref_spaces[space_i]->get_element_assembly_list(prev_ref_element, &al_prev);
          ref_spaces[space_i]->get_element_assembly_list(ref_element, &al);

          // Fun begins here - we need to figure out which DOFs on the new reference spaces can be reused and which cannot.
          unsigned short last_matched_index = 0;
          for (unsigned short i_al_prev = 0; i_al_prev < al_prev.cnt; i_al_prev++)
          {
            for (unsigned short j_al = last_matched_index; j_al < al.cnt; j_al++)
            {
              if (al.dof[j_al] < 0)
                continue;
              if (al_prev.idx[i_al_prev] == al.idx[j_al])
              {
                last_matched_index = j_al;
                if (StateReassemblyHelper<Scalar>::current_DOFs_to_reassemble->find(std::pair<int, unsigned char>(al_prev.dof[i_al_prev], space_i)) == StateReassemblyHelper<Scalar>::current_DOFs_to_reassemble->end())
                {
                  if (!*StateReassemblyHelper<Scalar>::reusable_DOFs[al.dof[j_al]])
                  {
                    DOF_to_reassemble_count++;
                    *StateReassemblyHelper<Scalar>::reusable_DOFs[al.dof[j_al]] = true;
                  }
                  DOF_to_DOF_map[al_prev.dof[i_al_prev]] = al.dof[j_al];
                }
              }
            }
          }
        }
      }

      cpu_time.tick();
      Hermes::Mixins::Loggable::Static::info("\t      No. of elements to reassemble: %i / %i total - %2.0f%%.", newSpace_elements_to_reassemble.size(), total_elements_new_spaces, ((float)newSpace_elements_to_reassemble.size() / (float)total_elements_new_spaces) * 100.);
      Hermes::Mixins::Loggable::Static::info("\t      No. of DOFs to reassemble: %i / %i total - %2.0f%%.", DOF_to_reassemble_count, ref_system_size, ((float)DOF_to_reassemble_count / (float)ref_system_size) * 100.);
      Hermes::Mixins::Loggable::Static::info("\t      Search for elements to reassemble: %4.3f s", cpu_time.last());
      cpu_time.tick();

      // Using the changed element on the new reference space, select the states from the original states that need to be recalculated.
      Traverse::State** new_states = malloc_with_check<Traverse::State*>(num_states, true);
      int new_num_states = 0;
      for (int state_i = 0; state_i < num_states; state_i++)
      {
        for (int space_i = 0; space_i < spaces_size; space_i++)
        {
          if (newSpace_elements_to_reassemble.find(std::pair<int, unsigned char>(states[state_i]->e[space_i]->id, space_i)) != newSpace_elements_to_reassemble.end())
          {
            new_states[new_num_states++] = Traverse::State::clone(states[state_i]);
            break;
          }
        }
      }
      free_with_check(states);
      new_states = realloc_with_check<Traverse::State*>(new_states, new_num_states);
      states = new_states;

      cpu_time.tick();
      Hermes::Mixins::Loggable::Static::info("\t      No. of states to reassemble: %i / %i total - %2.0f%%.", new_num_states, num_states, ((float)new_num_states / (float)num_states) * 100.);
      num_states = new_num_states;
      Hermes::Mixins::Loggable::Static::info("\t      Picking the new states: %4.3f s", cpu_time.last());
      cpu_time.tick();

      // Now we have to use the DOF to DOF map to fill in the necessary entries in the new matrix and rhs from the old ones.
      Scalar* Ax = StateReassemblyHelper<Scalar>::prev_mat->get_Ax();
      int* Ai = StateReassemblyHelper<Scalar>::prev_mat->get_Ai();
      int* Ap = StateReassemblyHelper<Scalar>::prev_mat->get_Ap();
      int total_entries = 0;
      int used_entries = 0;
      unsigned short current_row_entries;
      for (int i = 0; i < prev_ref_system_size; i++)
      {
        current_row_entries = Ap[i + 1] - Ap[i];
        total_entries += current_row_entries;
        if (DOF_to_DOF_map[i] != -1)
        {
          for (int j = 0; j < current_row_entries; j++)
          {
            if (DOF_to_DOF_map[Ai[Ap[i] + j]] != -1)
            {
              mat->add(DOF_to_DOF_map[i], DOF_to_DOF_map[Ai[Ap[i] + j]], Ax[Ap[i] + j]);
              used_entries++;
            }
          }
          rhs->add(DOF_to_DOF_map[i], StateReassemblyHelper<Scalar>::prev_rhs[i]);
        }
      }

      cpu_time.tick();
      Hermes::Mixins::Loggable::Static::info("\t      Reused linear system entries: %i / %i total - %2.0f%%.", used_entries, total_entries, ((float)used_entries / (float)total_entries) * 100.);
      Hermes::Mixins::Loggable::Static::info("\t      Copying the linear system: %4.3f s", cpu_time.last());
      free_with_check(DOF_to_DOF_map);
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
          this->scalar_views.back()->get_linearizer()->set_criterion(Views::LinearizerCriterionFixed(0));
          this->scalar_views.back()->set_title("Reference solution #%i", i);

          this->base_views.push_back(new Views::BaseView<Scalar>("", new Views::WinGeom(i * 410, 340, 400, 300)));
          this->base_views.back()->get_linearizer()->set_criterion(Views::LinearizerCriterionFixed(0));
          this->base_views.back()->set_title("Reference space #%i - basis", i);

          this->order_viewsRef.push_back(new Views::OrderView("", new Views::WinGeom(i * 410, 670, 400, 300)));
          this->order_viewsRef.back()->set_title("Reference space #%i - orders", i);

          this->order_views.push_back(new Views::OrderView("", new Views::WinGeom(i * 410, 1030, 400, 300)));
          this->order_views.back()->set_title("Coarse space #%i - orders", i);
        }
      }

      this->solver = new SolverType(wf, spaces);
      this->solver->dp->set_reassembled_states_reuse_linear_system_fn(&get_states_to_reassemble<Scalar>);
      this->adaptivity_internal = new Adapt<Scalar>(spaces, error_calculator, stopping_criterion_single_step);

      StateReassemblyHelper<Scalar>::reusable_DOFs = new bool*;
      this->solver->dp->set_reusable_DOFs(StateReassemblyHelper<Scalar>::reusable_DOFs);

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
        delete this->order_viewsRef[i];
        delete this->order_views[i];
      }

      delete this->adaptivity_internal;
      delete this->solver;
      delete StateReassemblyHelper<Scalar>::reusable_DOFs;
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::solve()
    {
      this->init_solving();
      do
      {
        this->info("\tAdaptSolver step %d:", this->adaptivity_step);

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
        this->info("\tSolving on reference mesh, %i DOFs.", Space<Scalar>::get_num_dofs(ref_spaces));
        this->solver->solve();

        // Free reusable DOFs data structures for this run.
        free_with_check<bool>(*StateReassemblyHelper<Scalar>::reusable_DOFs, true);

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
        this->info("\tProjecting reference solution on coarse mesh.");
        OGProjection<Scalar>::project_global(spaces, ref_slns, slns);

        // Calculate element errors.
        this->info("\tCalculating error estimate.");
        this->error_calculator->calculate_errors(slns, ref_slns, true);
        double error = this->error_calculator->get_total_error_squared() * 100;
        this->info("\tThe error estimate: %f.", error);

        // If err_est too large, adapt the mesh.
        if (this->stopping_criterion_global->done(error, this->adaptivity_step))
        {
          this->deinit_solving();
          return;
        }
        else
        {
          this->info("\tAdapting coarse mesh.");
          total_elements_prev_spaces = 0;
          for (int i = 0; i < this->spaces.size(); i++)
            total_elements_prev_spaces += this->spaces[i]->get_mesh()->get_num_active_elements();
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
      Hermes::Mixins::Loggable::Static::info("\t   Marking elements to reassemble on the already used Ref. Space:");
      this->tick();
      this->elements_to_reassemble.clear();
      this->DOFs_to_reassemble.clear();

      int total_elements_prev_ref_spaces = 0;
      for (int i = 0; i < this->spaces.size(); i++)
        total_elements_prev_ref_spaces += this->ref_spaces[i]->get_mesh()->get_num_active_elements();

      // Identify elements that changed.
      int valid_elements_to_refine_count = 0;
      for (int i = 0; i < this->adaptivity_internal->elements_to_refine_count; i++)
      {
        ElementToRefine* element_to_refine = &this->adaptivity_internal->elements_to_refine[i];
        if (!element_to_refine->valid)
          continue;
        valid_elements_to_refine_count++;
        RefinementType refinement_type = element_to_refine->split;
        int component = element_to_refine->comp;
        Element* e = this->ref_spaces[component]->get_mesh()->get_element_fast(element_to_refine->id);
        for (int son_i = 0; son_i < H2D_MAX_ELEMENT_SONS; son_i++)
          this->elements_to_reassemble.insert(std::pair<int, unsigned char>(e->sons[son_i]->id, component));
      }
      Hermes::Mixins::Loggable::Static::info("\t      No. of coarse mesh elements refined: %i / %i total - %2.0f%%.", valid_elements_to_refine_count, total_elements_prev_spaces, ((float)valid_elements_to_refine_count / (float)total_elements_prev_spaces) * 100.);
      Hermes::Mixins::Loggable::Static::info("\t      No. of fine mesh elements directly changed: %i / %i total - %2.0f%%.", this->elements_to_reassemble.size(), total_elements_prev_ref_spaces, ((float)this->elements_to_reassemble.size() / (float)total_elements_prev_ref_spaces) * 100.);
      this->tick();
      Hermes::Mixins::Loggable::Static::info("\t      Identify elements that changed: %4.3f s", this->last());
      this->tick();

      // Identify DOFs of the changed elements.
      AsmList<Scalar> al;
      for (std::set<std::pair<int, unsigned char> >::iterator it = elements_to_reassemble.begin(); it != elements_to_reassemble.end(); ++it)
      {
        int component = it->second;
        Element* e = this->ref_spaces[component]->get_mesh()->get_element_fast(it->first);
        this->ref_spaces[component]->get_element_assembly_list(e, &al);
        for (int j = 0; j < al.cnt; j++)
          this->DOFs_to_reassemble.insert(std::pair<int, unsigned char>(al.dof[j], component));
      }

      this->tick();
      Hermes::Mixins::Loggable::Static::info("\t      Identify DOFs that changed: %4.3f s", this->last());
      this->tick();

      // Take a look at other elements if they share a DOF that changed.
      /// \todo This is ineffective, a more effective way is to employ an improved neighbor searching.
      for(unsigned short space_i = 0; space_i < this->ref_spaces.size(); space_i++)
      {
        for_all_active_elements_fast(this->ref_spaces[space_i]->get_mesh())
        {
          this->ref_spaces[space_i]->get_element_assembly_list(e, &al);
          for (int j = 0; j < al.cnt; j++)
          {
            if (this->DOFs_to_reassemble.find(std::pair<int, unsigned char>(al.dof[j], space_i)) != this->DOFs_to_reassemble.end())
            {
              this->elements_to_reassemble.insert(std::pair<int, unsigned char>(e->id, space_i));
              break;
            }
          }
        }
      }

      Hermes::Mixins::Loggable::Static::info("\t      No. of all fine mesh elements that need to be reassembled: %i / %i total - %2.0f%%.", this->elements_to_reassemble.size(), total_elements_prev_ref_spaces, ((float)this->elements_to_reassemble.size() / (float)total_elements_prev_ref_spaces) * 100.);
      this->tick();
      Hermes::Mixins::Loggable::Static::info("\t      Looking for elements containing a changed DOF: %4.3f s", this->last());
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::visualize(std::vector<SpaceSharedPtr<Scalar> >& ref_spaces)
    {
      for(unsigned short i = 0; i < this->spaces.size(); i++)
      {
        this->scalar_views[i]->show(this->ref_slns[i]);
        this->base_views[i]->show(ref_spaces[i]);
        this->order_viewsRef[i]->show(ref_spaces[i]);
        this->order_views[i]->show(spaces[i]);
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