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

#include "discrete_problem/dg/discrete_problem_dg_assembler.h"
#include "discrete_problem/discrete_problem_helpers.h"
#include "discrete_problem.h"
#include "function/exact_solution.h"
#include "mesh/traverse.h"
#include "space/space.h"
#include "function/solution.h"
#include "api2d.h"

using namespace Hermes::Algebra::DenseMatrixOperations;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    DiscreteProblem<Scalar>::DiscreteProblem(WeakForm<Scalar>* wf_, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces)
    {
      init();
      this->set_spaces(spaces);
      this->set_weak_formulation(wf_);
    }

    template<typename Scalar>
    DiscreteProblem<Scalar>::DiscreteProblem(WeakForm<Scalar>* wf_, SpaceSharedPtr<Scalar>& space)
    {
      init();
      this->set_space(space);
      this->set_weak_formulation(wf_);
    }

    template<typename Scalar>
    DiscreteProblem<Scalar>::DiscreteProblem()
    {
      init();
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init()
    {
      this->spaces_size = this->spaces.size();

      this->nonlinear = true;
      this->add_dirichlet_lift = false;

      // Local number of threads - to avoid calling it over and over again, and against faults caused by the
      // value being changed while assembling.
      this->threadAssembler = new DiscreteProblemThreadAssembler<Scalar>*[this->num_threads_used];
      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i] = new DiscreteProblemThreadAssembler<Scalar>(&this->selectiveAssembler);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_linear(bool to_set, bool dirichlet_lift_accordingly)
    {
      this->nonlinear = !to_set;
      if(dirichlet_lift_accordingly)
        this->add_dirichlet_lift = !this->nonlinear;
      else
        this->add_dirichlet_lift = this->nonlinear;
    }

    template<typename Scalar>
    DiscreteProblem<Scalar>::~DiscreteProblem()
    {
      
      for(int i = 0; i < this->num_threads_used; i++)
        delete this->threadAssembler[i];
      delete [] this->threadAssembler;
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::isOkay() const
    {
      if(!this->wf)
        return false;

      if(this->spaces_size == 0)
        return false;

      // Initial check of meshes and spaces.
      for(unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
        this->spaces[space_i]->check();

      for(unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
        if(!this->spaces[space_i]->is_up_to_date())
          throw Exceptions::Exception("Space is out of date, if you manually refine it, you have to call assign_dofs().");

      return true;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_time(double time)
    {
      Space<Scalar>::update_essential_bc_values(spaces, time);
      this->wf->set_current_time(time);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_time_step(double time_step)
    {
      this->wf->set_current_time_step(time_step);
    }

    template<typename Scalar>
    Hermes::vector<SpaceSharedPtr<Scalar> >& DiscreteProblem<Scalar>::get_spaces()
    {
      return this->spaces;
    }
    
    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_RK(int original_spaces_count, bool force_diagonal_blocks_, Table* block_weights_)
    {
      Mixins::DiscreteProblemRungeKutta<Scalar>::set_RK(original_spaces_count, force_diagonal_blocks_, block_weights_);

      this->selectiveAssembler.set_RK(original_spaces_count, force_diagonal_blocks_, block_weights_);

      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i]->set_RK(original_spaces_count, force_diagonal_blocks_, block_weights_);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_do_not_use_cache(bool to_set)
    {
      DiscreteProblemCacheSettings::set_do_not_use_cache(to_set);
      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i]->set_do_not_use_cache(to_set);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_report_cache_hits_and_misses(bool to_set)
    {
      DiscreteProblemCacheSettings::set_report_cache_hits_and_misses(to_set);
      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i]->set_report_cache_hits_and_misses(to_set);
    }
      
    template<typename Scalar>
    void DiscreteProblem<Scalar>::invalidate_matrix()
    {
      this->selectiveAssembler.matrix_structure_reusable = false;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_weak_formulation(WeakForm<Scalar>* wf_)
    {
      Mixins::DiscreteProblemWeakForm<Scalar>::set_weak_formulation(wf_);

      this->selectiveAssembler.set_weak_formulation(wf_);
      this->selectiveAssembler.matrix_structure_reusable = false;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_matrix(SparseMatrix<Scalar>* mat)
    {
      Mixins::DiscreteProblemMatrixVector<Scalar>::set_matrix(mat);
      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i]->set_matrix(mat);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_rhs(Vector<Scalar>* rhs)
    {
      Mixins::DiscreteProblemMatrixVector<Scalar>::set_rhs(rhs);
      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i]->set_rhs(rhs);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spacesToSet)
    {
      if(this->spaces_size != spacesToSet.size() && this->spaces_size > 0)
        throw Hermes::Exceptions::LengthException(0, spacesToSet.size(), this->spaces_size);

      for(unsigned int i = 0; i < spacesToSet.size(); i++)
      {
        if(!spacesToSet[i])
          throw Exceptions::NullException(0, i);
        spacesToSet[i]->check();
      }

      this->spaces_size = spacesToSet.size();
      this->spaces = spacesToSet;

      this->selectiveAssembler.set_spaces(spacesToSet);

      /// \todo TEMPORARY There is something wrong with caching vector shapesets.
      for(unsigned int i = 0; i < spacesToSet.size(); i++)
        if(spacesToSet[i]->get_shapeset()->get_num_components() > 1)
          this->set_do_not_use_cache();

      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i]->init_spaces(spaces);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_space(SpaceSharedPtr<Scalar>& space)
    {
      Hermes::vector<SpaceSharedPtr<Scalar> > spaces;
      spaces.push_back(space);
      this->set_spaces(spaces);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs)
    {
      assemble((Solution<Scalar>**)NULL, mat, rhs);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(Scalar* coeff_vec, Vector<Scalar>* rhs)
    {
      assemble(coeff_vec, NULL, rhs);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(Vector<Scalar>* rhs)
    {
      assemble((Solution<Scalar>**)NULL, NULL, rhs);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(Scalar* coeff_vec, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs)
    {
      Solution<Scalar>** u_ext_sln = NULL;

      if(this->nonlinear && coeff_vec)
      {
        u_ext_sln = new Solution<Scalar>*[spaces_size];
        int first_dof = 0;
        for(int i = 0; i < this->spaces_size; i++)
        {
          u_ext_sln[i] = new Solution<Scalar>(spaces[i]->get_mesh());
          Solution<Scalar>::vector_to_solution(coeff_vec, spaces[i], u_ext_sln[i], !this->rungeKutta, first_dof);
          first_dof += spaces[i]->get_num_dofs();
        }
      }

      assemble(u_ext_sln, mat, rhs);

      if(this->nonlinear && coeff_vec)
      {
        for(int i = 0; i < this->spaces_size; i++)
          delete u_ext_sln[i];
        delete [] u_ext_sln;
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_assembling(Traverse::State**& states, int& num_states, Solution<Scalar>** u_ext_sln, Hermes::vector<MeshSharedPtr >& meshes)
    {
      // Optionally zero cache hits and misses.
      if(this->report_cache_hits_and_misses)
        this->zero_cache_hits_and_misses();

      // Vector of meshes.
      for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
        meshes.push_back(spaces[space_i]->get_mesh());
      for(unsigned int ext_i = 0; ext_i < this->wf->ext.size(); ext_i++)
        meshes.push_back(this->wf->ext[ext_i]->get_mesh());
      for(unsigned int form_i = 0; form_i < this->wf->get_forms().size(); form_i++)
        for(unsigned int ext_i = 0; ext_i < this->wf->get_forms()[form_i]->ext.size(); ext_i++)
          if(this->wf->get_forms()[form_i]->ext[ext_i])
            meshes.push_back(this->wf->get_forms()[form_i]->ext[ext_i]->get_mesh());

      if(this->nonlinear)
      {
        for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
          meshes.push_back(spaces[space_i]->get_mesh());
      }
      
      // Important.
      // This must be here, because the weakforms may have changed since set_weak_formulation (where the following calls
      // used to be in development). And since the following clones the passed WeakForm, this has to be called
      // only after the weak forms are ready for calculation.
      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i]->set_weak_formulation(this->wf);

      Traverse trav(this->spaces_size);
      states = trav.get_states(meshes, num_states);

      // Init the caught parallel exception message.
      this->exceptionMessageCaughtInParallelBlock.clear();
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(Solution<Scalar>** u_ext_sln, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs)
    {
      // Check.
      this->check();

      // Set the matrices.
      this->set_matrix(mat);
      this->set_rhs(rhs);

      // Initialize cache.
      this->cache.init_assembling();

      // Initialize states && previous iterations.
      int num_states;
      Traverse::State** states;
      Hermes::vector<MeshSharedPtr> meshes;
      this->init_assembling(states, num_states, u_ext_sln, meshes);

      // Creating matrix sparse structure.
      // If there are no states, return.
      if(this->selectiveAssembler.prepare_sparse_structure(this->current_mat, this->current_rhs, this->spaces, states, num_states))
      {
        // Is this a DG assembling.
        bool is_DG = this->wf->is_DG();

  #pragma omp parallel num_threads(this->num_threads_used)
        {
          int thread_number = omp_get_thread_num();
          int start = (num_states / this->num_threads_used) * thread_number;
          int end = (num_states / this->num_threads_used) * (thread_number + 1);
          if(thread_number == this->num_threads_used - 1)
            end = num_states;

          try
          {
            this->threadAssembler[thread_number]->init_assembling(u_ext_sln, spaces, this->nonlinear, this->add_dirichlet_lift);

            DiscreteProblemDGAssembler<Scalar>* dgAssembler;
            if(is_DG)
              dgAssembler = new DiscreteProblemDGAssembler<Scalar>(this->threadAssembler[thread_number], this->spaces, meshes);

            for(int state_i = start; state_i < end; state_i++)
            {
              Traverse::State* current_state = states[state_i];

              this->threadAssembler[thread_number]->init_assembling_one_state(spaces, current_state);

              this->threadAssembler[thread_number]->handle_cache(spaces, &this->cache);

              this->threadAssembler[thread_number]->assemble_one_state();

              if(is_DG)
              {
                dgAssembler->init_assembling_one_state(current_state);
                dgAssembler->assemble_one_state();
                dgAssembler->deinit_assembling_one_state();
              }

              this->threadAssembler[thread_number]->deinit_assembling_one_state();
            }

            if(is_DG)
              delete dgAssembler;

            this->threadAssembler[thread_number]->deinit_assembling();
          }
          catch(std::exception& exception)
          {
            this->exceptionMessageCaughtInParallelBlock = exception.what();
          }
        }
      }

      // Deinitialize states && previous iterations.
      this->deinit_assembling(states, num_states);

      // Deinitialize cache.
      this->cache.free_unused();

      /// Finish the algebraic structures for solving.
      if(this->current_mat)
        this->current_mat->finish();
      if(this->current_rhs)
        this->current_rhs->finish();

      if(!this->exceptionMessageCaughtInParallelBlock.empty())
      {
        throw Hermes::Exceptions::Exception(this->exceptionMessageCaughtInParallelBlock.c_str());
        return;
      }

      Element* e;
      for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
      {
        for_all_active_elements(e, spaces[space_i]->get_mesh())
        {
          spaces[space_i]->edata[e->id].changed_in_last_adaptation = false;
          e->visited = false;
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::deinit_assembling(Traverse::State** states, int num_states)
    {
      for(int i = 0; i < this->num_threads_used; i++)
      {
        this->threadAssembler[i]->deinit_assembling();

        if(this->report_cache_hits_and_misses)
          this->add_cache_hits_and_misses(this->threadAssembler[i]);
      }

      for(int i = 0; i < num_states; i++)
        delete states[i];
      free(states);
    }

    template class HERMES_API DiscreteProblem<double>;
    template class HERMES_API DiscreteProblem<std::complex<double> >;
  }
}