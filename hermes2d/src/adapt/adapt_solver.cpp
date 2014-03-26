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
#include "views/mesh_view.h"

namespace Hermes
{
  namespace Hermes2D
  {
    AdaptSolverCriterion::AdaptSolverCriterion()
    {
    }

    AdaptSolverCriterionErrorThreshold::AdaptSolverCriterionErrorThreshold(double error_tolerance) : AdaptSolverCriterion(), error_threshold(error_threshold)
    {
    }

    bool AdaptSolverCriterionErrorThreshold::done(double error, unsigned short iteration)
    {
      return error < this->error_threshold;
    }

    AdaptSolverCriterionFixed::AdaptSolverCriterionFixed(unsigned short refinement_levels) : AdaptSolverCriterion(), refinement_levels(refinement_levels)
    {
    }

    bool AdaptSolverCriterionFixed::done(double error, unsigned short iteration)
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
      this->number_of_equations = this->spaces.size();
      Helpers::check_length(this->selectors, this->number_of_equations);
      this->solve_method_running = false;
      this->visualization = false;
      this->prev_mat = nullptr;
      this->prev_rhs = nullptr;
      this->prev_dirichlet_lift_rhs = nullptr;
      this->solver = new SolverType(wf, spaces);
    }

    template<typename Scalar, typename SolverType>
    AdaptSolver<Scalar, SolverType>::~AdaptSolver()
    {
      if (this->prev_mat)
        delete this->prev_mat;

      if (this->prev_rhs)
        delete this->prev_rhs;

      if (this->prev_dirichlet_lift_rhs)
        delete this->prev_dirichlet_lift_rhs;

      delete this->solver;
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::set_verbose_output(bool to_set)
    {
      Loggable::set_verbose_output(to_set);
      if (this->error_calculator)
        this->error_calculator->set_verbose_output(to_set);
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
      static std::unordered_set<unsigned int>* current_elements_to_reassemble[H2D_MAX_COMPONENTS];
      static std::unordered_set<int>* current_DOFs_to_reassemble[H2D_MAX_COMPONENTS];
      static SpaceSharedPtrVector<Scalar>* current_ref_spaces;
      static SpaceSharedPtrVector<Scalar>* current_prev_ref_spaces;
      static unsigned char current_number_of_equations;
      static CSCMatrix<Scalar>* current_prev_mat;
      static Vector<Scalar>* current_prev_rhs;
      static Vector<Scalar>* current_prev_dirichlet_lift_rhs;
      
      static bool use_Dirichlet;
      static bool** reusable_DOFs;
      static bool** reusable_Dirichlet;
    };

    template<typename Scalar>
    int StateReassemblyHelper<Scalar>::current_iteration;

    template<typename Scalar>
    std::unordered_set<unsigned int>* StateReassemblyHelper<Scalar>::current_elements_to_reassemble[H2D_MAX_COMPONENTS];

    template<typename Scalar>
    std::unordered_set<int>* StateReassemblyHelper<Scalar>::current_DOFs_to_reassemble[H2D_MAX_COMPONENTS];

    template<typename Scalar>
    SpaceSharedPtrVector<Scalar>* StateReassemblyHelper<Scalar>::current_ref_spaces;

    template<typename Scalar>
    SpaceSharedPtrVector<Scalar>* StateReassemblyHelper<Scalar>::current_prev_ref_spaces;

    template<typename Scalar>
    CSCMatrix<Scalar>* StateReassemblyHelper<Scalar>::current_prev_mat;

    template<typename Scalar>
    Vector<Scalar>* StateReassemblyHelper<Scalar>::current_prev_rhs;

    template<typename Scalar>
    Vector<Scalar>* StateReassemblyHelper<Scalar>::current_prev_dirichlet_lift_rhs;

    template<typename Scalar>
    unsigned char StateReassemblyHelper<Scalar>::current_number_of_equations;

    template<typename Scalar>
    bool** StateReassemblyHelper<Scalar>::reusable_DOFs;

    template<typename Scalar>
    bool** StateReassemblyHelper<Scalar>::reusable_Dirichlet;

    template<typename Scalar>
    bool StateReassemblyHelper<Scalar>::use_Dirichlet;

    #define DEBUG_VIEWS

    template<typename Scalar>
    static bool compare(Scalar a, Scalar b);

    template<>
    static bool compare(std::complex<double> a, std::complex<double> b)
    {
      if (a.real() < b.real() && a.imag() < b.imag())
        return true;
      else
      return false;
    }

    template<>
    static bool compare(double a, double b)
    {
      return a < b;
    }

    template<typename Scalar>
    void get_states_to_reassemble(Traverse::State**& states, unsigned int& num_states, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, Vector<Scalar>* dirichlet_lift_rhs, Scalar*& coeff_vec)
    {
      int current_iteration = StateReassemblyHelper<Scalar>::current_iteration;
      if (StateReassemblyHelper<Scalar>::current_iteration == 1)
      {
        (*StateReassemblyHelper<Scalar>::reusable_DOFs) = nullptr;
        (*StateReassemblyHelper<Scalar>::reusable_Dirichlet) = nullptr;
        return;
      }

      // Create utility data
      // - time measurement
      Hermes::Mixins::TimeMeasurable cpu_time;
      // - size of spaces
      // - shortcuts for spaces, to speed up indirect accesses
      Space<Scalar>** prev_ref_spaces = new Space<Scalar>*[StateReassemblyHelper<Scalar>::current_number_of_equations];
      Space<Scalar>** ref_spaces = new Space<Scalar>*[StateReassemblyHelper<Scalar>::current_number_of_equations];
      // - number of elements of all new reference spaces combined
      int total_elements_new_spaces = 0;
      // - fill the above structures
      for (unsigned short i = 0; i < StateReassemblyHelper<Scalar>::current_number_of_equations; i++)
      {
        prev_ref_spaces[i] = StateReassemblyHelper<Scalar>::current_prev_ref_spaces->at(i).get();
        ref_spaces[i] = StateReassemblyHelper<Scalar>::current_ref_spaces->at(i).get();
        total_elements_new_spaces += ref_spaces[i]->get_mesh()->get_num_active_elements();
      }
      // - elements to reassemble on the new space
      // -- we could maybe directly used the resulting states (in which we are interested), but it would be complicated.
      std::unordered_set<unsigned int> newSpace_elements_to_reassemble[H2D_MAX_COMPONENTS];
      // - number of DOFs of the previous spaces.
      unsigned int prev_ref_system_size = StateReassemblyHelper<Scalar>::current_prev_mat->get_size();
      unsigned int ref_system_size = Space<Scalar>::get_num_dofs(*StateReassemblyHelper<Scalar>::current_ref_spaces);
      // - DOF to DOF map of those DOFs we can reuse. -1s are there to distinguish those we cannot.
      int* DOF_to_DOF_map = malloc_with_check<int>(prev_ref_system_size);
      for (unsigned int i = 0; i < prev_ref_system_size; i++)
        DOF_to_DOF_map[i] = -1;
      // - DOF to reassemble.
      (*StateReassemblyHelper<Scalar>::reusable_DOFs) = calloc_with_check<bool>(ref_system_size);
      (*StateReassemblyHelper<Scalar>::reusable_Dirichlet) = calloc_with_check<bool>(StateReassemblyHelper<Scalar>::current_number_of_equations);
      // - utility assembly lists.
      AsmList<Scalar> al, al_prev;
      // - dummy functions for traversing the previous and current reference spaces together
      MeshFunctionSharedPtrVector<Scalar> dummy_fns;
      for (unsigned short i = 0; i < StateReassemblyHelper<Scalar>::current_number_of_equations; i++)
        dummy_fns.push_back(new ZeroSolution<Scalar>(StateReassemblyHelper<Scalar>::current_prev_ref_spaces->at(i)->get_mesh()));
      for (unsigned short i = 0; i < StateReassemblyHelper<Scalar>::current_number_of_equations; i++)
        dummy_fns.push_back(new ZeroSolution<Scalar>(StateReassemblyHelper<Scalar>::current_ref_spaces->at(i)->get_mesh()));
      // - (for reporting) count of DOFs needed to be reassembled
      int reusable_DOFs_count = 0;
      
      // Start.
      Hermes::Mixins::Loggable::Static::info("\t   Handling Reusing matrix entries on the new Ref. Space:");
      cpu_time.tick();

      // Find out if Dirichlet is not directly changed.
      for (unsigned char space_i = 0; space_i < StateReassemblyHelper<Scalar>::current_number_of_equations; space_i++)
      {
        if (StateReassemblyHelper<Scalar>::current_DOFs_to_reassemble[space_i]->find(-1) == StateReassemblyHelper<Scalar>::current_DOFs_to_reassemble[space_i]->end())
          (*StateReassemblyHelper<Scalar>::reusable_Dirichlet)[space_i] = true;
      }

      // Traverse the previous and current reference Spaces at once.
      Traverse trav(StateReassemblyHelper<Scalar>::current_number_of_equations * 2);
      unsigned int num_states_local;
      Traverse::State** states_local = trav.get_states(dummy_fns, num_states_local);
      for (unsigned int local_state_i = 0; local_state_i < num_states_local; local_state_i++)
      {
        Traverse::State* current_local_state = states_local[local_state_i];
        for (unsigned char space_i = 0; space_i < StateReassemblyHelper<Scalar>::current_number_of_equations; space_i++)
        {
          // Mark the appropriate elements
          Element* prev_ref_element = current_local_state->e[space_i];
          Element* ref_element = current_local_state->e[space_i + StateReassemblyHelper<Scalar>::current_number_of_equations];
          // If the element is changed on the previous ref space, mark the appropriate element on the new ref space as to-reassemble.
          bool element_added = false;
          if (StateReassemblyHelper<Scalar>::current_elements_to_reassemble[space_i]->find(prev_ref_element->id) != StateReassemblyHelper<Scalar>::current_elements_to_reassemble[space_i]->end())
          {
            newSpace_elements_to_reassemble[space_i].insert(ref_element->id);
            element_added = true;
          }
          // Get assembly lists.
          prev_ref_spaces[space_i]->get_element_assembly_list(prev_ref_element, &al_prev);
          ref_spaces[space_i]->get_element_assembly_list(ref_element, &al);
          // Different assembly lists => reassemble element.
          if (!element_added && al_prev.cnt != al.cnt)
          {
            newSpace_elements_to_reassemble[space_i].insert(ref_element->id);
            element_added = true;
          }

          // Fun begins here - we need to figure out which DOFs on the new reference spaces can be reused and which cannot.
          unsigned short last_matched_index = 0;
          bool contains_Dirichlet = false, anything_changed = al_prev.cnt != al.cnt;
          for (unsigned short i_al_prev = 0; i_al_prev < al_prev.cnt; i_al_prev++)
          {
            bool reused = false;
            if (!element_added && al_prev.idx[i_al_prev] != al.idx[i_al_prev])
            {
              newSpace_elements_to_reassemble[space_i].insert(ref_element->id);
              element_added = true;
            }
            if (al_prev.dof[i_al_prev] < 0)
            {
              contains_Dirichlet = true;
              continue;
            }
            for (unsigned short j_al = last_matched_index; j_al < al.cnt; j_al++)
            {
              if (al.dof[j_al] < 0)
              {
                contains_Dirichlet = true;
                continue;
              }
              if (al_prev.idx[i_al_prev] == al.idx[j_al] && compare<Scalar>(-HermesEpsilon, al_prev.coef[i_al_prev] - al.coef[j_al]) && compare<Scalar>(al_prev.coef[i_al_prev] - al.coef[j_al], HermesEpsilon))
              {
                last_matched_index = j_al + 1;
                if (StateReassemblyHelper<Scalar>::current_DOFs_to_reassemble[space_i]->find(al_prev.dof[i_al_prev]) == StateReassemblyHelper<Scalar>::current_DOFs_to_reassemble[space_i]->end())
                {
                  if (!(*StateReassemblyHelper<Scalar>::reusable_DOFs)[al.dof[j_al]])
                    reusable_DOFs_count++;
                  (*StateReassemblyHelper<Scalar>::reusable_DOFs)[al.dof[j_al]] = true;
                  DOF_to_DOF_map[al_prev.dof[i_al_prev]] = al.dof[j_al];
                  reused = true;
                }
                break;
              }
            }
            if (!reused)
            {
              assert(DOF_to_DOF_map[al_prev.dof[i_al_prev]] == -1);
              anything_changed = true;
            }
          }
          // Dirichlet - indirectly changed.
          if (anything_changed && contains_Dirichlet)
            (*StateReassemblyHelper<Scalar>::reusable_Dirichlet)[space_i] = false;
        }
      }

      // Now we have to add all elements containing Dirichlet DOF (changed or not)
      if (StateReassemblyHelper<Scalar>::use_Dirichlet)
      {
        for (unsigned int local_state_i = 0; local_state_i < num_states_local; local_state_i++)
        {
          Traverse::State* current_local_state = states_local[local_state_i];
          for (unsigned char space_i = 0; space_i < StateReassemblyHelper<Scalar>::current_number_of_equations; space_i++)
          {
            // If Dirichlet is not necessary to be recalculated, go on.
            if ((*StateReassemblyHelper<Scalar>::reusable_Dirichlet)[space_i])
              continue;
            // Get current element.
            Element* ref_element = current_local_state->e[space_i + StateReassemblyHelper<Scalar>::current_number_of_equations];
            // If element is already being recalculated, continue.
            if (newSpace_elements_to_reassemble[space_i].find(ref_element->id) != newSpace_elements_to_reassemble[space_i].end())
              continue;

            // Mark the previous element.
            Element* prev_ref_element = current_local_state->e[space_i];

            // Get assembly lists.
            prev_ref_spaces[space_i]->get_element_assembly_list(prev_ref_element, &al_prev);
            ref_spaces[space_i]->get_element_assembly_list(ref_element, &al);

            // Since the element has not been added before we know the length of the prev and current assembly lists is the same.
            for (unsigned short al_i = 0; al_i < al_prev.cnt; al_i++)
            {
              if (al_prev.dof[al_i] < 0 || al.dof[al_i] < 0)
              {
                newSpace_elements_to_reassemble[space_i].insert(ref_element->id);
                break;
              }
            }
          }
        }
      }

      cpu_time.tick();
      Hermes::Mixins::Loggable::Static::info("\t      No. of DOFs to reassemble: %i / %i total - %2.0f%%.", ref_system_size - reusable_DOFs_count, ref_system_size, ((float)(ref_system_size - reusable_DOFs_count) / (float)ref_system_size) * 100.);
      Hermes::Mixins::Loggable::Static::info("\t      Search for elements to reassemble: %4.3f s", cpu_time.last());
      cpu_time.tick();

      // Using the changed element on the new reference space, select the states from the original states that need to be recalculated.
#ifdef DEBUG_VIEWS
      bool** reassembled = new bool*[StateReassemblyHelper<Scalar>::current_number_of_equations];
      for (unsigned char space_i = 0; space_i < StateReassemblyHelper<Scalar>::current_number_of_equations; space_i++)
        reassembled[space_i] = (bool*)calloc(ref_spaces[space_i]->get_mesh()->get_max_element_id() + 1, sizeof(bool));
#endif
      Traverse::State** new_states = malloc_with_check<Traverse::State*>(num_states, true);
      int new_num_states = 0;
      for (unsigned int state_i = 0; state_i < num_states; state_i++)
      {
        for (unsigned char space_i = 0; space_i < StateReassemblyHelper<Scalar>::current_number_of_equations; space_i++)
        {
          if (newSpace_elements_to_reassemble[space_i].find(states[state_i]->e[space_i]->id) != newSpace_elements_to_reassemble[space_i].end())
          {
            new_states[new_num_states++] = Traverse::State::clone(states[state_i]);
#ifdef DEBUG_VIEWS
            for (unsigned char space_i = 0; space_i < StateReassemblyHelper<Scalar>::current_number_of_equations; space_i++)
              reassembled[space_i][states[state_i]->e[space_i]->id] = true;
#endif
            break;
          }
        }
      }

      for (unsigned int i = 0; i < num_states; i++)
        delete states[i];
      free_with_check(states);
      new_states = realloc_with_check<Traverse::State*>(new_states, new_num_states);
      states = new_states;

      cpu_time.tick();
      Hermes::Mixins::Loggable::Static::info("\t      No. of states to reassemble: %i / %i total - %2.0f%%.", new_num_states, num_states, ((float)new_num_states / (float)num_states) * 100.);
      num_states = new_num_states;
      Hermes::Mixins::Loggable::Static::info("\t      Picking the new states: %4.3f s", cpu_time.last());
      cpu_time.tick();

      // Now we have to use the DOF to DOF map to fill in the necessary entries in the new matrix and rhs from the old ones.
      Scalar* Ax = StateReassemblyHelper<Scalar>::current_prev_mat->get_Ax();
      int* Ai = StateReassemblyHelper<Scalar>::current_prev_mat->get_Ai();
      int* Ap = StateReassemblyHelper<Scalar>::current_prev_mat->get_Ap();
      Scalar* new_coeff_vec = nullptr;
      if (coeff_vec)
        new_coeff_vec = calloc_with_check<Scalar>(ref_system_size, true);
      int total_entries = 0;
      int used_entries = 0;
      unsigned short current_row_entries;
      for (unsigned int i = 0; i < prev_ref_system_size; i++)
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
          if(rhs)
            rhs->add(DOF_to_DOF_map[i], StateReassemblyHelper<Scalar>::current_prev_rhs->get(i));
          if (coeff_vec)
            new_coeff_vec[DOF_to_DOF_map[i]] = coeff_vec[i];
        }
      }
      if (dirichlet_lift_rhs && StateReassemblyHelper<Scalar>::use_Dirichlet)
      {
        unsigned int running_count = 0;
        for (unsigned short space_i = 0; space_i < StateReassemblyHelper<Scalar>::current_number_of_equations; space_i++)
        {
          if ((*StateReassemblyHelper<Scalar>::reusable_Dirichlet)[space_i])
          {
            for (unsigned int dof_i = running_count; dof_i < running_count + prev_ref_spaces[space_i]->get_num_dofs(); dof_i++)
            {
              if (DOF_to_DOF_map[dof_i] != -1)
                dirichlet_lift_rhs->add(DOF_to_DOF_map[dof_i], StateReassemblyHelper<Scalar>::current_prev_dirichlet_lift_rhs->get(dof_i));
            }
          }
          running_count += prev_ref_spaces[space_i]->get_num_dofs();
        }
      }

      cpu_time.tick();
      Hermes::Mixins::Loggable::Static::info("\t      Reused matrix entries: %i / %i total - %2.0f%%.", used_entries, total_entries, ((float)used_entries / (float)total_entries) * 100.);
      Hermes::Mixins::Loggable::Static::info("\t      Copying the linear system: %4.3f s", cpu_time.last());

      mat->export_to_file("MatrixAfterReuse", "A", EXPORT_FORMAT_MATLAB_SIMPLE, "%16.16f");
      if(rhs)
        rhs->export_to_file("RhsAfterReuse", "b", EXPORT_FORMAT_PLAIN_ASCII, "%16.16f");
      if (dirichlet_lift_rhs && StateReassemblyHelper<Scalar>::use_Dirichlet)
        dirichlet_lift_rhs->export_to_file("DirLiftAfterReuse", "bd", EXPORT_FORMAT_PLAIN_ASCII, "%16.16f");
      
      if (coeff_vec)
      {
        FILE* coeff_vec_prev_file = fopen("CoeffVecPrev", "w");
        FILE* coeff_vec_file = fopen("CoeffVec", "w");
        for (unsigned int i = 0; i < prev_ref_system_size; i++)
          fprintf(coeff_vec_prev_file, "%i (-> %i): %f\n", i, DOF_to_DOF_map[i], coeff_vec[i]);
        for (unsigned int i = 0; i < ref_system_size; i++)
          fprintf(coeff_vec_file, "%i, %f\n", i, new_coeff_vec[i]);
        fclose(coeff_vec_file);
        fclose(coeff_vec_prev_file);
      }

      FILE* DOF_to_DOF_file = fopen("DOF_to_DOF", "w");
      for (unsigned int i = 0; i < prev_ref_system_size; i++)
        fprintf(DOF_to_DOF_file, "%i -> %i\n", i, DOF_to_DOF_map[i]);
      fclose(DOF_to_DOF_file);

      FILE* DOF_to_DOF_plus_1_file = fopen("DOF_to_DOF_plus1", "w");
      for (unsigned int i = 0; i < prev_ref_system_size; i++)
        fprintf(DOF_to_DOF_plus_1_file, "%i -> %i\n", i + 1, DOF_to_DOF_map[i] + 1);
      fclose(DOF_to_DOF_plus_1_file);

#ifdef DEBUG_VIEWS

      std::vector<Views::ScalarView*> reassembled_views;
      std::vector<Views::BaseView<Scalar>*> base_views_prev;
      std::vector<Views::BaseView<Scalar>*> base_views_curr;

      for (unsigned char space_i = 0; space_i < StateReassemblyHelper<Scalar>::current_number_of_equations; space_i++)
      {
        reassembled_views.push_back(new Views::ScalarView("", new Views::WinGeom(630 + space_i * 310, 630, 300, 300)));
        reassembled_views.back()->get_linearizer()->set_criterion(Views::LinearizerCriterionFixed(1));
        MeshFunctionSharedPtr<double> reassembled_states_function(new ExactSolutionConstantArray<double, bool>(ref_spaces[space_i]->get_mesh(), reassembled[space_i]));
        reassembled_views.back()->show(reassembled_states_function);
        base_views_prev.push_back(new Views::BaseView<Scalar>("", new Views::WinGeom(630, space_i * 310, 300, 300)));
        base_views_curr.push_back(new Views::BaseView<Scalar>("", new Views::WinGeom(950, space_i * 310, 300, 300)));
        base_views_prev.back()->get_linearizer()->set_criterion(Views::LinearizerCriterionFixed(3));
        base_views_prev.back()->show(prev_ref_spaces[space_i]);
        base_views_curr.back()->get_linearizer()->set_criterion(Views::LinearizerCriterionFixed(3));
        base_views_curr.back()->show(ref_spaces[space_i]);
      }

      Views::View::wait_for_keypress();
      for (unsigned char space_i = 0; space_i < StateReassemblyHelper<Scalar>::current_number_of_equations; space_i++)
      {
        ::free(reassembled[space_i]);
        delete reassembled_views[space_i];
        //delete base_views_prev[space_i];
        //delete base_views_curr[space_i];
      }
      delete[] reassembled;
#endif

      if (coeff_vec)
      {
        free_with_check(coeff_vec, true);
        coeff_vec = new_coeff_vec;
      }

      free_with_check(DOF_to_DOF_map);
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::init_solving(AdaptivityType adaptivityType)
    {
      this->adaptivity_step = 1;
      this->solve_method_running = true;
      ref_slns.clear();
      slns.clear();
      for (unsigned char i = 0; i < this->spaces.size(); i++)
      {
        ref_slns.push_back(MeshFunctionSharedPtr<Scalar>(new Solution<Scalar>));
        slns.push_back(MeshFunctionSharedPtr<Scalar>(new Solution<Scalar>));
        if (this->visualization)
        {
          this->scalar_views.push_back(new Views::ScalarView("", new Views::WinGeom(i * 320, 0, 300, 300)));
          this->scalar_views.back()->get_linearizer()->set_criterion(Views::LinearizerCriterionFixed(4));
          this->scalar_views.back()->set_title("Reference solution #%i", i);

          this->order_views.push_back(new Views::OrderView("", new Views::WinGeom(i * 320, 320, 300, 300)));
          this->order_views.back()->set_title("Reference space #%i - orders", i);

          this->base_views.push_back(new Views::BaseView<Scalar>("", new Views::WinGeom(i * 320, 640, 300, 300)));
          this->base_views.back()->get_linearizer()->set_criterion(Views::LinearizerCriterionFixed(4));
        }
      }

      this->solver->set_verbose_output(this->get_verbose_output());
      this->solver->output_matrix();
      this->solver->set_matrix_export_format(EXPORT_FORMAT_MATLAB_SIMPLE);

      this->solver->output_rhs();
      this->solver->set_rhs_export_format(EXPORT_FORMAT_PLAIN_ASCII);

      this->adaptivity_internal = new Adapt<Scalar>(spaces, error_calculator, stopping_criterion_single_step);
      this->adaptivity_internal->set_verbose_output(this->get_verbose_output());

      StateReassemblyHelper<Scalar>::reusable_DOFs = new bool*;
      StateReassemblyHelper<Scalar>::reusable_Dirichlet = new bool*;
      this->solver->dp->set_reusable_DOFs(StateReassemblyHelper<Scalar>::reusable_DOFs, StateReassemblyHelper<Scalar>::reusable_Dirichlet);
      StateReassemblyHelper<Scalar>::use_Dirichlet = this->solver->dp->add_dirichlet_lift;
      StateReassemblyHelper<Scalar>::current_number_of_equations = this->number_of_equations;

      this->ref_spaces.clear();

      if (adaptivityType == hAdaptivity)
      {
        for (unsigned char selector_i = 0; selector_i < this->spaces.size(); selector_i++)
          hOnlySelectors.push_back(new RefinementSelectors::HOnlySelector<Scalar>());
      }
      if (adaptivityType == pAdaptivity)
      {
        for (unsigned char selector_i = 0; selector_i < this->spaces.size(); selector_i++)
          pOnlySelectors.push_back(new RefinementSelectors::POnlySelector<Scalar>());
      }
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::deinit_solving()
    {
      this->solve_method_running = false;

      if (this->visualization)
      for (unsigned short i = 0; i < this->number_of_equations; i++)
      {
        delete this->scalar_views[i];
        delete this->order_views[i];
        delete this->base_views[i];
      }

      delete this->adaptivity_internal;
      delete StateReassemblyHelper<Scalar>::reusable_DOFs;
      delete StateReassemblyHelper<Scalar>::reusable_Dirichlet;

      for (unsigned char selector_i = 0; selector_i < this->hOnlySelectors.size(); selector_i++)
        delete hOnlySelectors[selector_i];
      hOnlySelectors.clear();
      for (unsigned char selector_i = 0; selector_i < this->pOnlySelectors.size(); selector_i++)
        delete pOnlySelectors[selector_i];
      pOnlySelectors.clear();
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::solve(AdaptivityType adaptivityType)
    {
      // Initialization - allocations etc.
      this->init_solving(adaptivityType);

      do
      {
        this->info("\tAdaptSolver step %d:", this->adaptivity_step);

        // Construct globally refined reference meshes and setup reference spaces.
        prev_ref_spaces = ref_spaces;
        ref_spaces.clear();
        for (unsigned char i = 0; i < number_of_equations; i++)
        {
          Mesh::ReferenceMeshCreator ref_mesh_creator(spaces[i]->get_mesh());
          MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
          Space<Scalar>::ReferenceSpaceCreator u_ref_space_creator(spaces[i], adaptivityType == pAdaptivity ? spaces[i]->get_mesh() : ref_mesh, adaptivityType == hAdaptivity ? 0 : 0);
          ref_spaces.push_back(u_ref_space_creator.create_ref_space());
        }

        // Set these to the solver.
        this->solver->set_spaces(ref_spaces);

        // Initialize the states handler.
        StateReassemblyHelper<Scalar>::current_iteration = this->adaptivity_step;
        StateReassemblyHelper<Scalar>::current_ref_spaces = &ref_spaces;
        StateReassemblyHelper<Scalar>::current_prev_ref_spaces = &prev_ref_spaces;
        for (unsigned char i = 0; i < number_of_equations; i++)
        {
          StateReassemblyHelper<Scalar>::current_elements_to_reassemble[i] = &this->elements_to_reassemble[i];
          StateReassemblyHelper<Scalar>::current_DOFs_to_reassemble[i] = &this->DOFs_to_reassemble[i];
        }

        // Perform solution.
        this->info("\tSolving on reference mesh, %i DOFs.", Space<Scalar>::get_num_dofs(ref_spaces));
        this->solver->dp->set_reassembled_states_reuse_linear_system_fn(&get_states_to_reassemble<Scalar>);
        if (adaptivity_step > 1 && dynamic_cast<NewtonSolver<Scalar>*>(this->solver))
        {
          Scalar* new_sln_vector = realloc_with_check<Scalar>(this->solver->sln_vector, Space<Scalar>::get_num_dofs(ref_spaces));
          this->solver->sln_vector = nullptr;
          memset(new_sln_vector + Space<Scalar>::get_num_dofs(prev_ref_spaces), 0, (Space<Scalar>::get_num_dofs(ref_spaces) - Space<Scalar>::get_num_dofs(prev_ref_spaces)) * sizeof(Scalar));
          this->solver->solve(new_sln_vector);
          free_with_check(new_sln_vector, true);
        }
        else
          this->solver->solve();

        if (this->solver->dp->add_dirichlet_lift)
          this->solver->dp->dirichlet_lift_rhs->export_to_file("DirLift", "bd", EXPORT_FORMAT_PLAIN_ASCII, "%16.16f");

        this->solver->get_linear_matrix_solver()->free();

        // Free reusable DOFs data structures for this run.
        free_with_check<bool>(*StateReassemblyHelper<Scalar>::reusable_DOFs);
        free_with_check<bool>(*StateReassemblyHelper<Scalar>::reusable_Dirichlet);

        // Update the stored (previous) linear system.
        if (this->prev_mat)
          delete this->prev_mat;
        this->prev_mat = (Hermes::Algebra::CSCMatrix<Scalar>*)(this->solver->get_jacobian()->duplicate());
        StateReassemblyHelper<Scalar>::current_prev_mat = this->prev_mat;

        if (this->prev_rhs)
          delete this->prev_rhs;
        this->prev_rhs = this->solver->get_residual()->duplicate();
        if (this->solver->dp->add_dirichlet_lift)
          this->prev_rhs->subtract_vector(this->solver->dp->dirichlet_lift_rhs);
        StateReassemblyHelper<Scalar>::current_prev_rhs = this->prev_rhs;

        if (this->solver->dp->add_dirichlet_lift)
        {
          if (this->prev_dirichlet_lift_rhs)
            delete this->prev_dirichlet_lift_rhs;
          this->prev_dirichlet_lift_rhs = this->solver->dp->dirichlet_lift_rhs->duplicate();
          StateReassemblyHelper<Scalar>::current_prev_dirichlet_lift_rhs = this->prev_dirichlet_lift_rhs;
        }

        // Translate the resulting coefficient vector into the instance of Solution.
        Solution<Scalar>::vector_to_solutions(solver->get_sln_vector(), ref_spaces, ref_slns);

        if (this->visualization)
          this->visualize(ref_spaces);

        // Project the fine mesh solution onto the coarse mesh.
        this->info("\tProjecting reference solution on coarse mesh.");
        OGProjection<Scalar>::project_global(spaces, ref_slns, slns);

        // Calculate element errors.
        if (!this->exact_slns.empty())
        {
          this->info("\tCalculating exact error.");
          this->error_calculator->calculate_errors(slns, exact_slns, false);
          double error = this->error_calculator->get_total_error_squared() * 100;
          this->info("\tExact error: %f.", error);
        }

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
          for (unsigned short i = 0; i < number_of_equations; i++)
            total_elements_prev_spaces += this->spaces[i]->get_mesh()->get_num_active_elements();
          // For hp-adaptivity, we use the provided selectors.
          if (adaptivityType == hpAdaptivity)
            this->adaptivity_internal->adapt(this->selectors);
          else if (adaptivityType == hAdaptivity)
            this->adaptivity_internal->adapt(hOnlySelectors);
          else
            this->adaptivity_internal->adapt(pOnlySelectors);

          this->mark_elements_to_reassemble();
        }

        this->adaptivity_step++;
        this->info("\n");
      } while (true);

      this->deinit_solving();
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::mark_elements_to_reassemble()
    {
      Hermes::Mixins::Loggable::Static::info("\t   Marking elements to reassemble on the already used Ref. Space:");
      this->tick();

      // Clear the arrays.
      for (unsigned char i = 0; i < this->spaces.size(); i++)
      {
        this->elements_to_reassemble[i].clear();
        this->DOFs_to_reassemble[i].clear();
      }

      int total_elements_prev_ref_spaces = 0;
      for (unsigned short i = 0; i < number_of_equations; i++)
        total_elements_prev_ref_spaces += this->ref_spaces[i]->get_mesh()->get_num_active_elements();

      // Identify elements that changed.
      int valid_elements_to_refine_count = 0;
      int ref_elements_changed_count = 0;
      for (int i = 0; i < this->adaptivity_internal->elements_to_refine_count; i++)
      {
        ElementToRefine* element_to_refine = &this->adaptivity_internal->elements_to_refine[i];
        if (!element_to_refine->valid)
          continue;
        valid_elements_to_refine_count++;
        RefinementType refinement_type = element_to_refine->split;
        int component = element_to_refine->comp;
        Element* e = this->ref_spaces[component]->get_mesh()->get_element_fast(element_to_refine->id);
        if (e->active)
        {
          this->elements_to_reassemble[component].insert(e->id);
          ref_elements_changed_count++;
        }
        else
        {
          for (int son_i = 0; son_i < H2D_MAX_ELEMENT_SONS; son_i++)
          {
            this->elements_to_reassemble[component].insert(e->sons[son_i]->id);
            ref_elements_changed_count++;
          }
        }
      }
      Hermes::Mixins::Loggable::Static::info("\t      No. of coarse mesh elements refined: %i / %i total - %2.0f%%.", valid_elements_to_refine_count, total_elements_prev_spaces, ((float)valid_elements_to_refine_count / (float)total_elements_prev_spaces) * 100.);
      Hermes::Mixins::Loggable::Static::info("\t      No. of fine mesh elements directly changed: %i / %i total - %2.0f%%.", ref_elements_changed_count, total_elements_prev_ref_spaces, ((float)ref_elements_changed_count / (float)total_elements_prev_ref_spaces) * 100.);
      this->tick();
      Hermes::Mixins::Loggable::Static::info("\t      Identify elements that changed: %4.3f s", this->last());
      this->tick();

      // Identify DOFs of the changed elements.
      AsmList<Scalar> al;
      for (unsigned char space_i = 0; space_i < this->number_of_equations; space_i++)
      {
        for (std::unordered_set<unsigned int>::iterator it = elements_to_reassemble[space_i].begin(); it != elements_to_reassemble[space_i].end(); ++it)
        {
          Element* e = this->ref_spaces[space_i]->get_mesh()->get_element_fast(*it);
          this->ref_spaces[space_i]->get_element_assembly_list(e, &al);
          for (int j = 0; j < al.cnt; j++)
            this->DOFs_to_reassemble[space_i].insert(al.dof[j]);
        }
      }

      this->tick();
      Hermes::Mixins::Loggable::Static::info("\t      Identify DOFs that changed: %4.3f s", this->last());
      this->tick();

      // Take a look at other elements if they share a DOF that changed.
      /// \todo This is ineffective, a more effective way is to employ an improved neighbor searching.
      ref_elements_changed_count = 0;
      for (unsigned char space_i = 0; space_i < this->number_of_equations; space_i++)
      {
        for_all_active_elements_fast(this->ref_spaces[space_i]->get_mesh())
        {
          this->ref_spaces[space_i]->get_element_assembly_list(e, &al);
          for (int j = 0; j < al.cnt; j++)
          {
            if (this->DOFs_to_reassemble[space_i].find(al.dof[j]) != this->DOFs_to_reassemble[space_i].end())
            {
              this->elements_to_reassemble[space_i].insert(e->id);
              ref_elements_changed_count++;
              break;
            }
          }
        }
      }

      Hermes::Mixins::Loggable::Static::info("\t      No. of all fine mesh elements that need to be reassembled: %i / %i total - %2.0f%%.", ref_elements_changed_count, total_elements_prev_ref_spaces, ((float)ref_elements_changed_count / (float)total_elements_prev_ref_spaces) * 100.);
      this->tick();
      Hermes::Mixins::Loggable::Static::info("\t      Looking for elements containing a changed DOF: %4.3f s", this->last());
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::visualize(std::vector<SpaceSharedPtr<Scalar> >& ref_spaces)
    {
      for (unsigned short i = 0; i < this->number_of_equations; i++)
      {
        this->scalar_views[i]->show(this->ref_slns[i]);
        this->order_views[i]->show(ref_spaces[i]);
        //this->base_views[i]->show(ref_spaces[i]);
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
    void AdaptSolver<Scalar, SolverType>::set_exact_solutions(MeshFunctionSharedPtrVector<Scalar> exact_slns)
    {
      if (this->solve_method_running)
        throw Exceptions::Exception("AdaptSolver asked to change the exact_slns while it was running.");
      this->exact_slns = exact_slns;
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::set_initial_spaces(SpaceSharedPtrVector<Scalar> spaces)
    {
      if (this->solve_method_running)
        throw Exceptions::Exception("AdaptSolver asked to change the initial spaces while it was running.");
      Helpers::check_length(this->spaces, spaces);
      this->spaces = spaces;

      this->adaptivity_internal->set_spaces(this->spaces);
      this->solver->set_spaces(this->spaces);
    }

    template<typename Scalar, typename SolverType>
    void AdaptSolver<Scalar, SolverType>::set_wf(WeakFormSharedPtr<Scalar> wf)
    {
      if (this->solve_method_running)
        throw Exceptions::Exception("AdaptSolver asked to change the weak formulation while it was running.");
      this->wf = wf;
      this->solver->set_weak_formulation(wf);
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
      Helpers::check_length(this->selectors, selectors);
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
    SolverType* AdaptSolver<Scalar, SolverType>::get_solver()
    {
      return this->solver;
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
      for (unsigned short i = 0; i < this->number_of_equations; i++)
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