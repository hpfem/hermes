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
    DiscreteProblem<Scalar>::DiscreteProblem(WeakForm<Scalar>* wf_, Hermes::vector<SpaceSharedPtr<Scalar> > spaces)
      : Hermes::Solvers::DiscreteProblemInterface<Scalar>()
    {
      init();
      this->set_spaces(spaces);
      this->set_weak_formulation(wf_);
    }

    template<typename Scalar>
    DiscreteProblem<Scalar>::DiscreteProblem(WeakForm<Scalar>* wf_, SpaceSharedPtr<Scalar> space)
      : Hermes::Solvers::DiscreteProblemInterface<Scalar>()
    {
      init();
      this->set_space(space);
      this->set_weak_formulation(wf_);
    }

    template<typename Scalar>
    DiscreteProblem<Scalar>::DiscreteProblem() : Hermes::Solvers::DiscreteProblemInterface<Scalar>()
    {
      init();
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init()
    {
      this->spaces_size = this->spaces.size();

      // Internal variables settings.
      sp_seq = spaces_size == 0 ? NULL : new int[wf->get_neq()];
      if(sp_seq)
        memset(sp_seq, -1, sizeof(int) * wf->get_neq());

      // Matrix<Scalar> related settings.
      matrix_structure_reusable = false;

      this->nonlinear = true;

      current_mat = NULL;
      current_rhs = NULL;

      this->do_not_use_cache = false;

      // Local number of threads - to avoid calling it over and over again, and against faults caused by the
      // value being changed while assembling.
      num_threads_used = Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads);
      this->threadAssembler = new DiscreteProblemThreadAssembler<Scalar>[num_threads_used];
    }

    template<typename Scalar>
    DiscreteProblem<Scalar>::~DiscreteProblem()
    {
      if(wf)
        memset(sp_seq, -1, sizeof(int) * wf->get_neq());

      if(sp_seq)
        delete [] sp_seq;

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
    Hermes::vector<SpaceSharedPtr<Scalar> > DiscreteProblem<Scalar>::get_spaces() const
    {
      return this->spaces;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_weak_formulation(WeakForm<Scalar>* wf_)
    {
      if(!wf_)
        throw Hermes::Exceptions::NullException(0);

      this->wf = wf_;

      this->matrix_structure_reusable = false;

      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i].set_weak_formulation(wf_);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_matrix(SparseMatrix<Scalar>* mat)
    {
      Mixins::DiscreteProblemMatrixVector<Scalar>::set_matrix(mat);
      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i].set_matrix(mat);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_rhs(Vector<Scalar>* rhs)
    {
      Mixins::DiscreteProblemMatrixVector<Scalar>::set_rhs(rhs);
      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i].set_rhs(rhs);
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::is_up_to_date() const
    {
      // check if we can reuse the matrix structure
      bool up_to_date = true;
      if(!matrix_structure_reusable)
        up_to_date = false;

      for (unsigned int i = 0; i < wf->get_neq(); i++)
      {
        if(spaces[i]->get_seq() != sp_seq[i])
        {
          up_to_date = false;
          break;
        }
      }

      return up_to_date;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> > spacesToSet)
    {
      if(this->spaces_size != spacesToSet.size() && this->spaces_size > 0)
        throw Hermes::Exceptions::LengthException(0, spacesToSet.size(), this->spaces_size);
      
      for(unsigned int i = 0; i < spacesToSet.size(); i++)
      {
        if(!spacesToSet[i])
          throw Exceptions::NullException(0, i);
        spacesToSet[i]->check();
      }

      if(this->spaces_size == 0)
      {
        // Internal variables settings.
        this->spaces_size = spacesToSet.size();
        sp_seq = new int[spaces_size];
        memset(sp_seq, -1, sizeof(int) * spaces_size);

        // Matrix<Scalar> related settings.
        matrix_structure_reusable = false;
      }
      else
      {
        for(unsigned int i = 0; i < spaces_size; i++)
          sp_seq[i] = spacesToSet[i]->get_seq();
      }

      this->spaces = spacesToSet;

      /// \todo TEMPORARY There is something wrong with caching vector shapesets.
      for(unsigned int i = 0; i < spacesToSet.size(); i++)
        if(spacesToSet[i]->get_shapeset()->get_num_components() > 1)
          this->set_do_not_use_cache();
      
      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i].init_spaces(spaces);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_space(SpaceSharedPtr<Scalar> space)
    {
      Hermes::vector<SpaceSharedPtr<Scalar> > spaces;
      spaces.push_back(space);
      this->set_spaces(spaces);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::create_sparse_structure(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs)
    {
      this->current_mat = mat;
      if(rhs)
        this->current_rhs = rhs;
      this->create_sparse_structure();
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::create_sparse_structure()
    {
      int ndof = Space<Scalar>::assign_dofs(spaces);

      if(is_up_to_date())
      {
        if(current_mat)
          current_mat->zero();
        if(current_rhs)
        {
          // If we use e.g. a new NewtonSolver (providing a new Vector) for this instance of DiscreteProblem that already assembled a system,
          // we end up with everything up_to_date, but unallocated Vector.
          if(current_rhs->length() == 0)
            current_rhs->alloc(ndof);
          else
            current_rhs->zero();
        }
        return;
      }

      // For DG, the sparse structure is different as we have to
      // account for over-edge calculations.
      bool is_DG = false;
      for(unsigned int i = 0; i < this->wf->mfsurf.size(); i++)
      {
        if(!this->wf->mfDG.empty())
        {
          is_DG = true;
          break;
        }
      }
      for(unsigned int i = 0; i < this->wf->vfsurf.size() && is_DG == false; i++)
      {
        if(!this->wf->vfDG.empty())
        {
          is_DG = true;
          break;
        }
      }

      if(current_mat)
      {
        // Spaces have changed: create the matrix from scratch.
        matrix_structure_reusable = true;
        current_mat->free();
        current_mat->prealloc(ndof);

        AsmList<Scalar>* al = new AsmList<Scalar>[wf->get_neq()];
        MeshSharedPtr* meshes = new MeshSharedPtr[wf->get_neq()];
        bool **blocks = wf->get_blocks(this->force_diagonal_blocks);

        // Init multi-mesh traversal.
        for (unsigned int i = 0; i < wf->get_neq(); i++)
          meshes[i] = spaces[i]->get_mesh();

        Traverse trav(true);
        trav.begin(wf->get_neq(), meshes);

        Traverse::State* current_state;
        // Loop through all elements.
        while (current_state = trav.get_next_state())
        {
          // Obtain assembly lists for the element at all spaces.
          /// \todo do not get the assembly list again if the element was not changed.
          for (unsigned int i = 0; i < wf->get_neq(); i++)
            if(current_state->e[i])
              spaces[i]->get_element_assembly_list(current_state->e[i], &(al[i]));

          if(is_DG)
          {
            // Number of edges ( =  number of vertices).
            int num_edges = current_state->e[0]->nvert;

            // Allocation an array of arrays of neighboring elements for every mesh x edge.
            Element **** neighbor_elems_arrays = new Element ***[wf->get_neq()];
            for(unsigned int i = 0; i < wf->get_neq(); i++)
              neighbor_elems_arrays[i] = new Element **[num_edges];

            // The same, only for number of elements
            int ** neighbor_elems_counts = new int *[wf->get_neq()];
            for(unsigned int i = 0; i < wf->get_neq(); i++)
              neighbor_elems_counts[i] = new int[num_edges];

            // Get the neighbors.
            for(unsigned int el = 0; el < wf->get_neq(); el++)
            {
              NeighborSearch<Scalar> ns(current_state->e[el], meshes[el]);

              // Ignoring errors (and doing nothing) in case the edge is a boundary one.
              ns.set_ignore_errors(true);

              for(int ed = 0; ed < num_edges; ed++)
              {
                ns.set_active_edge(ed);
                const Hermes::vector<Element *> *neighbors = ns.get_neighbors();

                neighbor_elems_counts[el][ed] = ns.get_num_neighbors();
                neighbor_elems_arrays[el][ed] = new Element *[neighbor_elems_counts[el][ed]];
                for(int neigh = 0; neigh < neighbor_elems_counts[el][ed]; neigh++)
                  neighbor_elems_arrays[el][ed][neigh] = (*neighbors)[neigh];
              }
            }

            // Pre-add into the stiffness matrix.
            for (unsigned int m = 0; m < wf->get_neq(); m++)
              for(unsigned int el = 0; el < wf->get_neq(); el++)
                for(int ed = 0; ed < num_edges; ed++)
                  for(int neigh = 0; neigh < neighbor_elems_counts[el][ed]; neigh++)
                    if((blocks[m][el] || blocks[el][m]) && current_state->e[m])
                    {
                      AsmList<Scalar>*am = &(al[m]);
                      AsmList<Scalar>*an = new AsmList<Scalar>;
                      spaces[el]->get_element_assembly_list(neighbor_elems_arrays[el][ed][neigh], an);

                      // pretend assembling of the element stiffness matrix
                      // register nonzero elements
                      for (unsigned int i = 0; i < am->cnt; i++)
                        if(am->dof[i] >= 0)
                          for (unsigned int j = 0; j < an->cnt; j++)
                            if(an->dof[j] >= 0)
                            {
                              if(blocks[m][el]) current_mat->pre_add_ij(am->dof[i], an->dof[j]);
                              if(blocks[el][m]) current_mat->pre_add_ij(an->dof[j], am->dof[i]);
                            }
                            delete an;
                    }

                    // Deallocation an array of arrays of neighboring elements
                    // for every mesh x edge.
                    for(unsigned int el = 0; el < wf->get_neq(); el++)
                    {
                      for(int ed = 0; ed < num_edges; ed++)
                        delete [] neighbor_elems_arrays[el][ed];
                      delete [] neighbor_elems_arrays[el];
                    }
                    delete [] neighbor_elems_arrays;

                    // The same, only for number of elements.
                    for(unsigned int el = 0; el < wf->get_neq(); el++)
                      delete [] neighbor_elems_counts[el];
                    delete [] neighbor_elems_counts;
          }

          // Go through all equation-blocks of the local stiffness matrix.
          for (unsigned int m = 0; m < wf->get_neq(); m++)
          {
            for (unsigned int n = 0; n < wf->get_neq(); n++)
            {
              if(blocks[m][n] && current_state->e[m] && current_state->e[n])
              {
                AsmList<Scalar>*am = &(al[m]);
                AsmList<Scalar>*an = &(al[n]);

                // Pretend assembling of the element stiffness matrix.
                for (unsigned int i = 0; i < am->cnt; i++)
                  if(am->dof[i] >= 0)
                    for (unsigned int j = 0; j < an->cnt; j++)
                      if(an->dof[j] >= 0)
                        current_mat->pre_add_ij(am->dof[i], an->dof[j]);
              }
            }
          }
        }

        trav.finish();
        delete [] al;
        delete [] meshes;
        delete [] blocks;

        current_mat->alloc();
      }

      // WARNING: unlike Matrix<Scalar>::alloc(), Vector<Scalar>::alloc(ndof) frees the memory occupied
      // by previous vector before allocating
      if(current_rhs)
        current_rhs->alloc(ndof);

      // save space seq numbers and weakform seq number, so we can detect their changes
      for (unsigned int i = 0; i < wf->get_neq(); i++)
        sp_seq[i] = spaces[i]->get_seq();
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

      if(this->nonlinear)
      {
        u_ext_sln = new Solution<Scalar>*[spaces_size];
        int first_dof = 0;
        for(int i = 0; i < this->spaces_size; i++)
        {
          u_ext_sln[i] = new Solution<Scalar>(spaces[i]->get_mesh());
          Solution<Scalar>::vector_to_solution(coeff_vec, spaces[i], u_ext_sln[i], !rungeKutta, first_dof);
          first_dof += spaces[i]->get_num_dofs();
        }
      }

      assemble(u_ext_sln, mat, rhs);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::init_assembling(Traverse::State**& states, int& num_states, Solution<Scalar>** u_ext_sln)
    {
      // Vector of meshes.
      Hermes::vector<MeshSharedPtr > meshes;
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

      Traverse trav_master(true);
      states = trav_master.get_states(meshes, num_states);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(Solution<Scalar>** u_ext_sln, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs)
    {
      // Check.
      this->check();

      // Set the matrices.
      this->set_matrix(mat);
      this->set_rhs(rhs);

      // Creating matrix sparse structure.
      create_sparse_structure();

      // Initialize cache.
      this->cache.init_assembling();

      // Initialize states && previous iterations.
      int num_states;
      Traverse::State** states;
      this->init_assembling(states, num_states, u_ext_sln);

      // Is this a DG assembling.
      bool is_DG = this->wf->is_DG();

#pragma omp parallel num_threads(num_threads_used)
      {
        int thread_number = omp_get_thread_num();
        int start = (num_states / num_threads_used) * thread_number;
        int end = (num_states / num_threads_used) * (thread_number + 1);
        if(thread_number == num_threads_used - 1)
          end = num_states;

        this->threadAssembler[thread_number].init_assembling(u_ext_sln, spaces, this->nonlinear);
        DiscreteProblemDGAssembler<Scalar>* dgAssembler;
        if(is_DG)
          dgAssembler = new DiscreteProblemDGAssembler<Scalar>(&this->threadAssembler[thread_number], this->spaces);

        for(int state_i = start; state_i < end; state_i++)
        {
          Traverse::State* current_state = states[state_i];

          this->threadAssembler[thread_number].init_assembling_one_state(spaces, current_state);
          
          this->threadAssembler[thread_number].handle_cache(spaces, &this->cache, this->do_not_use_cache);

          this->threadAssembler[thread_number].assemble_one_state();

          if(is_DG)
          {
            dgAssembler->init_assembling_one_state(current_state);
            dgAssembler->assemble_one_state();
            dgAssembler->deinit_assembling_one_state();
          }

          this->threadAssembler[thread_number].deinit_assembling_one_state();
        }

        if(is_DG)
          delete dgAssembler;
      }

      // Deinitialize states && previous iterations.
      this->deinit_assembling(states, num_states);

      // Deinitialize cache.
      this->cache.free_unused();

      /// Finish the algebraic structures for solving.
      if(current_mat)
        current_mat->finish();
      if(current_rhs)
        current_rhs->finish();

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
        this->threadAssembler[i].deinit_assembling();

      for(int i = 0; i < num_states; i++)
        delete states[i];
      free(states);
    }

    template class HERMES_API DiscreteProblem<double>;
    template class HERMES_API DiscreteProblem<std::complex<double> >;
  }
}