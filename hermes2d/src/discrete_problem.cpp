//#define DEBUG_DG_ASSEMBLING
//#define DEBUG_DG_ASSEMBLING_ELEMENT 44
//#define DEBUG_DG_ASSEMBLING_ISURF 3
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

#include "discrete_problem.h"
#include "function/exact_solution.h"
#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "global.h"
#include "integrals/h1.h"
#include "quadrature/limit_order.h"
#include "mesh/traverse.h"
#include "space/space.h"
#include "shapeset/precalc.h"
#include "mesh/refmap.h"
#include "function/solution.h"
#include "neighbor.h"
#include "api2d.h"

using namespace Hermes::Algebra::DenseMatrixOperations;

namespace Hermes
{
  namespace Hermes2D
  {

    template<typename Scalar>
    DiscreteProblem<Scalar>::DiscreteProblem(WeakForm<Scalar>* wf_, Hermes::vector<SpaceSharedPtr<Scalar> > spaces)
      : Hermes::Solvers::DiscreteProblemInterface<Scalar>(), dgAssembler(wf, spaces)
    {
      init();
      this->set_spaces(spaces);
      this->set_weak_formulation(wf_);
    }

    template<typename Scalar>
    DiscreteProblem<Scalar>::DiscreteProblem(WeakForm<Scalar>* wf_, SpaceSharedPtr<Scalar> space)
      : Hermes::Solvers::DiscreteProblemInterface<Scalar>(), dgAssembler(wf, space)
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

      this->is_linear = false;

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

      for(int i = 0; i < this->num_threads_used; i++)
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

      this->integrationOrderCalculator.set_weak_formulation(wf_);
      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i].set_weak_formulation(wf_);
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::is_matrix_free() const
    {
      return wf->is_matrix_free();
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_matrix(SparseMatrix<Scalar>* mat)
    {
      this->current_mat = mat;
      this->dgAssembler.set_matrix(mat);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::set_rhs(Vector<Scalar>* rhs)
    {
      this->current_rhs = rhs;
      this->dgAssembler.set_rhs(rhs);
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
      for(unsigned int i = 0; i < spacesToSet.size(); i++)
      {
        if(!spacesToSet[i])
          throw Exceptions::NullException(0, i);

        spacesToSet[i]->check();
      }

      if(this->spaces_size != spacesToSet.size() && this->spaces_size > 0)
        throw Hermes::Exceptions::LengthException(0, spacesToSet.size(), this->spaces_size);

      int originalSize = this->spaces.size();

      this->spaces = spacesToSet;

      this->dgAssembler.set_spaces(spacesToSet);

      /// \todo TEMPORARY There is something wrong with caching vector shapesets.
      for(unsigned int i = 0; i < spacesToSet.size(); i++)
        if(spacesToSet[i]->get_shapeset()->get_num_components() > 1)
          this->set_do_not_use_cache();

      matrix_structure_reusable = false;

      if(originalSize == 0)
      {
        // Internal variables settings.
        sp_seq = new int[spaces.size()];
        memset(sp_seq, -1, sizeof(int) * spaces.size());

        // Matrix<Scalar> related settings.
        matrix_structure_reusable = false;
      }
      else
      {
        for(unsigned int i = 1; i < spacesToSet.size(); i++)
        {
          sp_seq[i] = spacesToSet[i]->get_seq();
        }
      }

      this->spaces_size = this->spaces.size();

      for(int i = 0; i < this->num_threads_used; i++)
        this->threadAssembler[i].init_spaces(spaces)
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
      assemble(NULL, mat, rhs);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(Vector<Scalar>* rhs)
    {
      assemble(NULL, NULL, rhs);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(Scalar* coeff_vec, SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs)
    {
      // Check.
      this->check();

      // Important, sets the current caughtException to NULL.
      this->caughtException = NULL;

      this->set_matrix(mat);
      this->set_rhs(rhs);

      // Creating matrix sparse structure.
      create_sparse_structure();
      this->cache.init_assembling();

      if(this->is_linear)
      {
        Solution<Scalar>** u_ext_sln = new Solution<Scalar>*[spaces_size];
        int first_dof = 0;
        for(int i = 0; i < this->spaces_size; i++)
        {
          u_ext_sln[i] = new Solution<Scalar>(spaces[i]->get_mesh());
          Solution<Scalar>::vector_to_solution(coeff_vec, spaces[i], u_ext[i], !rungeKutta, first_dof);
          first_dof += spaces[i]->get_num_dofs();
        }

        for(int i = 0; i < this->num_threads_used; i++)
          this->threadAssembler[i].init_u_ext(spaces, u_ext_sln)
      }

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
      for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
        meshes.push_back(spaces[space_i]->get_mesh());

      Traverse trav_master(true);
      int num_states;
      Traverse::State** states = trav_master.get_states(meshes, num_states);

#pragma omp parallel num_threads(num_threads_used)
      {
        int thread_number = omp_get_thread_num();
        int start = (num_states / num_threads_used) * thread_number;
        int end = (num_states / num_threads_used) * (thread_number + 1);
        if(thread_number == num_threads_used - 1)
          end = num_states;

        this->threadAssembler[thread_number].init_assembling();

        int order;

        for(int state_i = start; state_i < end; state_i++)
        {
          if(this->caughtException)
            break;
          try
          {
            this->threadAssembler[thread_number].init_assembling_one_state(spaces, states[state_i]);

            typename DiscreteProblemCache<Scalar>::CacheRecord* cache_record;
            if(!this->do_not_use_cache)
              cache_record = this->get_state_cache(current_state, this->threadAssembler[thread_number], order);
            else
            {
              cache_record = new typename DiscreteProblemCache<Scalar>::CacheRecord;
              order = this->integrationOrderCalculator.calculate_order(spaces, current_state, current_refmaps, current_u_ext, current_weakform);
              cache_record->init(this->spaces, current_state, current_pss, current_refmaps, current_u_ext, current_als, current_als_surface, current_weakform, order);
            }

            assemble_one_state(cache_record, current_refmaps, current_u_ext, current_als, current_state, current_weakform);
            
            if(dgAssembler.DG_matrix_forms_present || dgAssembler.DG_vector_forms_present)
              this->dgAssembler.assemble_one_DG_state(current_pss, current_refmaps, current_u_ext, current_als, current_state, current_weakform->mfDG, current_weakform->vfDG, &fns[thread_number].front(), current_weakform);

            if(this->do_not_use_cache)
              delete cache_record;
          }
          catch(Hermes::Exceptions::Exception& e)
          {
            if(!this->caughtException)
              this->caughtException = e.clone();
          }
          catch(std::exception& e)
          {
            if(!this->caughtException)
              this->caughtException = new std::exception(e);
          }
        }
        
        if(this->dgAssembler.DG_matrix_forms_present || this->dgAssembler.DG_vector_forms_present)
         for (unsigned int j = 0; j < spaces_size; j++)
            delete current_spss[j];
        delete [] current_spss;

      }

      this->cache.free_unused();

      deinit_assembling(pss, refmaps, u_ext, als, alsSurface, weakforms, num_threads_used);

      for(int i = 0; i < num_states; i++)
        delete states[i];
      free(states);

      for(unsigned int i = 0; i < num_threads_used; i++)
      {
        fns[i].clear();
      }
      delete [] fns;

      /// \todo Should this be really here? Or in assemble()?
      if(current_mat)
        current_mat->finish();
      if(current_rhs)
        current_rhs->finish();

      if(dgAssembler.DG_matrix_forms_present || dgAssembler.DG_vector_forms_present)
      {
        Element* element_to_set_nonvisited;
        for(unsigned int mesh_i = 0; mesh_i < meshes.size(); mesh_i++)
          for_all_elements(element_to_set_nonvisited, meshes[mesh_i])
          element_to_set_nonvisited->visited = false;
      }

      Element* e;
      for(unsigned int space_i = 0; space_i < spaces.size(); space_i++)
      {
        for_all_active_elements(e, spaces[space_i]->get_mesh())
        spaces[space_i]->edata[e->id].changed_in_last_adaptation = false;
      }

      if(this->caughtException)
        throw *(this->caughtException);
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble(Scalar* coeff_vec, Vector<Scalar>* rhs)
    {
      assemble(coeff_vec, NULL, rhs);
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(MatrixForm<Scalar>* form, Traverse::State* current_state)
    {
      if(current_state->e[form->i] && current_state->e[form->j])
      {
        if(fabs(form->scaling_factor) < 1e-12)
          return false;

        // If a block scaling table is provided, and if the scaling coefficient
        // A_mn for this block is zero, then the form does not need to be assembled.
        if(this->block_weights)
          if(fabs(this->block_weights->get_A(form->i, form->j)) < 1e-12)
            return false;
        return true;
      }
      return false;
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(MatrixFormVol<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((MatrixForm<Scalar>*)form, current_state))
        return false;

      if(form->assembleEverywhere)
        return true;

      int this_marker = current_state->rep->marker;
      for (unsigned int ss = 0; ss < form->areas_internal.size(); ss++)
        if(form->areas_internal[ss] == this_marker)
          return true;

      return false;
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(MatrixFormSurf<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((MatrixForm<Scalar>*)form, current_state))
        return false;

      if(current_state->rep->en[current_state->isurf]->marker == 0)
        return false;

      if(form->assembleEverywhere)
        return true;

      int this_marker = current_state->rep->en[current_state->isurf]->marker;
      for (unsigned int ss = 0; ss < form->areas_internal.size(); ss++)
        if(form->areas_internal[ss] == this_marker)
          return true;

      return false;
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(VectorForm<Scalar>* form, Traverse::State* current_state)
    {
      if(!current_state->e[form->i])
        return false;
      if(fabs(form->scaling_factor) < 1e-12)
        return false;

      return true;
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(VectorFormVol<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((VectorForm<Scalar>*)form, current_state))
        return false;

      if(form->assembleEverywhere)
        return true;

      int this_marker = current_state->rep->marker;
      for (unsigned int ss = 0; ss < form->areas_internal.size(); ss++)
        if(form->areas_internal[ss] == this_marker)
          return true;

      return false;
    }

    template<typename Scalar>
    bool DiscreteProblem<Scalar>::form_to_be_assembled(VectorFormSurf<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((VectorForm<Scalar>*)form, current_state))
        return false;

      if(current_state->rep->en[current_state->isurf]->marker == 0)
        return false;

      if(form->assembleEverywhere)
        return true;

      int this_marker = current_state->rep->en[current_state->isurf]->marker;
      for (unsigned int ss = 0; ss < form->areas_internal.size(); ss++)
        if(form->areas_internal[ss] == this_marker)
          return true;

      return false;
    }

    template<typename Scalar>
    typename DiscreteProblemCache<Scalar>::CacheRecord* DiscreteProblem<Scalar>::get_state_cache(Traverse::State* current_state, PrecalcShapeset** current_pss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, AsmList<Scalar>*** current_alsSurface, WeakForm<Scalar>* current_wf, int& order)
    {
      typename DiscreteProblemCache<Scalar>::CacheRecord* cache_record = NULL;
      if(this->cache.get(current_state->rep, current_state->rep_subidx, current_state->rep_i, cache_record))
      {
        bool reinit = false;
        for(unsigned int i = 0; i < this->spaces_size; i++)
        {
          if(!current_state->e[i])
            continue;

          if(this->spaces[i]->edata[current_state->e[i]->id].changed_in_last_adaptation)
          {
            reinit = true;
            break;
          }

          if(cache_record->asmlistCnt[i] != current_als[i]->cnt)
          {
            reinit = true;
            break;
          }
          else
          {
            for(unsigned int idx_i = 0; idx_i < current_als[i]->cnt; idx_i++)
              if(current_als[i]->idx[idx_i] != cache_record->asmlistIdx[i][idx_i])
              {
                reinit = true;
                break;
              }
          }
        }
        if(reinit)
        {
          cache_record->free();
          order = this->integrationOrderCalculator.calculate_order(spaces, current_state, current_refmaps, current_u_ext, current_wf);
          cache_record->init(this->spaces, current_state, current_pss, current_refmaps, current_u_ext, current_als, current_alsSurface, current_wf, order);
        }
        else
          order = cache_record->order;
      }
      else
      {
        order = this->integrationOrderCalculator.calculate_order(spaces, current_state, current_refmaps, current_u_ext, current_wf);
        cache_record->init(this->spaces, current_state, current_pss, current_refmaps, current_u_ext, current_als, current_alsSurface, current_wf, order);
      }
      return cache_record;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_one_state(typename DiscreteProblemCache<Scalar>::CacheRecord* cache_record, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, 
      Traverse::State* current_state, WeakForm<Scalar>* current_wf)
    {
      // Representing space.
      int rep_space_i = -1;

      // Get necessary (volumetric) assembly lists.
      for(unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
      {
        if(current_state->e[space_i])
        {
          rep_space_i = space_i;
          current_refmaps[space_i]->set_active_element(current_state->e[space_i]);
          spaces[space_i]->get_element_assembly_list(current_state->e[space_i], current_als[space_i]);
        }
      }

      if(rep_space_i == -1)
        return;

      // Element-wise parameters for WeakForm.
      (const_cast<WeakForm<Scalar>*>(current_wf))->set_active_state(current_state->e);

      // Assembly lists for surface forms.
      AsmList<Scalar>** current_alsSurface = NULL;
      if(current_state->isBnd && (current_wf->mfsurf.size() > 0 || current_wf->vfsurf.size() > 0 || current_wf->mfDG.size() > 0 || current_wf->vfDG.size() > 0))
      {
        current_alsSurface = new AsmList<Scalar>*[this->spaces_size];
        for(unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
        {
          if(!current_state->e[space_i])
            continue;
          current_alsSurface[space_i] = new AsmList<Scalar>[current_state->rep->nvert];
          for (current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
            if(current_state->bnd[current_state->isurf])
              spaces[space_i]->get_boundary_assembly_list(current_state->e[space_i], current_state->isurf, &current_alsSurface[space_i][current_state->isurf]);
        }
      }

      if(this->caughtException)
        return;

      // Ext functions.
      // - order
      int order = cache_record->order;

      // - u_ext
      Func<Scalar>** u_ext = NULL;
      int prevNewtonSize = this->wf->get_neq();

      if(!this->is_linear)
      {
        u_ext = new Func<Scalar>*[prevNewtonSize];
        if(current_u_ext)
          for(int u_ext_i = 0; u_ext_i < prevNewtonSize; u_ext_i++)
            if(current_u_ext[u_ext_i])
              if(current_u_ext[u_ext_i]->get_active_element())
                u_ext[u_ext_i] = init_fn(current_u_ext[u_ext_i], order);
              else
                u_ext[u_ext_i] = NULL;
            else
              u_ext[u_ext_i] = NULL;
        else
          for(int u_ext_i = 0; u_ext_i < prevNewtonSize; u_ext_i++)
            u_ext[u_ext_i] = NULL;
      }

      // - ext
      int current_extCount = this->wf->ext.size();
      Func<Scalar>** ext = NULL;
      if(current_extCount > 0)
      {
        ext = new Func<Scalar>*[current_extCount];
        for(int ext_i = 0; ext_i < current_extCount; ext_i++)
          if(current_wf->ext[ext_i])
            if(current_wf->ext[ext_i]->get_active_element())
              ext[ext_i] = init_fn(current_wf->ext[ext_i].get(), order);
            else
              ext[ext_i] = NULL;
          else
            ext[ext_i] = NULL;
      }

      if(rungeKutta)
        for(int ext_i = 0; ext_i < this->RK_original_spaces_count; ext_i++)
          u_ext[ext_i]->add(ext[current_extCount - this->RK_original_spaces_count + ext_i]);

      if(current_mat)
      {
        for(int current_mfvol_i = 0; current_mfvol_i < wf->mfvol.size(); current_mfvol_i++)
        {
          MatrixFormVol<Scalar>* mfv = current_wf->mfvol[current_mfvol_i];

          if(!form_to_be_assembled(mfv, current_state))
            continue;

          int form_i = mfv->i;
          int form_j = mfv->j;

          assemble_matrix_form(mfv, 
            cache_record->order, 
            cache_record->fns[form_j], 
            cache_record->fns[form_i], 
            ext,
            u_ext,
            current_als[form_i], 
            current_als[form_j], 
            current_state, 
            cache_record->n_quadrature_points, 
            cache_record->geometry, 
            cache_record->jacobian_x_weights);
        }
      }
      if(current_rhs)
      {
        for(int current_vfvol_i = 0; current_vfvol_i < wf->vfvol.size(); current_vfvol_i++)
        {
          VectorFormVol<Scalar>* vfv = current_wf->vfvol[current_vfvol_i];

          if(!form_to_be_assembled(vfv, current_state))
            continue;

          int form_i = vfv->i;

          assemble_vector_form(vfv, 
            cache_record->order, 
            cache_record->fns[form_i], 
            ext,
            u_ext, 
            current_als[form_i], 
            current_state, 
            cache_record->n_quadrature_points,
            cache_record->geometry, 
            cache_record->jacobian_x_weights);
        }
      }

      // Cleanup - u_ext
      if(current_u_ext)
      {
        for(int u_ext_i = 0; u_ext_i < prevNewtonSize; u_ext_i++)
          if(current_u_ext[u_ext_i] && current_u_ext[u_ext_i]->get_active_element())
          {
            u_ext[u_ext_i]->free_fn();
            delete u_ext[u_ext_i];
          }
          delete [] u_ext;
      }

      // Cleanup - ext
      for(int ext_i = 0; ext_i < current_extCount; ext_i++)
      {
        if(current_wf->ext[ext_i] && current_wf->ext[ext_i]->get_active_element())
        {
          ext[ext_i]->free_fn();
          delete ext[ext_i];
        }
      }
      delete [] ext;

      // Assemble surface integrals now: loop through surfaces of the element.
      if(current_state->isBnd && (current_wf->mfsurf.size() > 0 || current_wf->vfsurf.size() > 0))
      {
        for (current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
        {
          if(!current_state->bnd[current_state->isurf])
            continue;

          // Edge-wise parameters for WeakForm.
          (const_cast<WeakForm<Scalar>*>(current_wf))->set_active_edge_state(current_state->e, current_state->isurf);

          // Ext functions.
          // - order
          int orderSurf = cache_record->orderSurface[current_state->isurf];

          // - u_ext
          int prevNewtonSize = this->wf->get_neq();
          Func<Scalar>** u_extSurf = NULL;
          if(!this->is_linear)
          {
            u_extSurf = new Func<Scalar>*[prevNewtonSize];
            if(current_u_ext)
              for(int u_ext_surf_i = 0; u_ext_surf_i < prevNewtonSize; u_ext_surf_i++)
                if(current_u_ext[u_ext_surf_i])
                  u_extSurf[u_ext_surf_i] = current_state->e[u_ext_surf_i] == NULL ? NULL : init_fn(current_u_ext[u_ext_surf_i], orderSurf);
                else
                  u_extSurf[u_ext_surf_i] = NULL;
            else
              for(int u_ext_surf_i = 0; u_ext_surf_i < prevNewtonSize; u_ext_surf_i++)
                u_extSurf[u_ext_surf_i] = NULL;
          }
          // - ext
          int current_extCount = this->wf->ext.size();
          Func<Scalar>** extSurf = new Func<Scalar>*[current_extCount];
          for(int ext_surf_i = 0; ext_surf_i < current_extCount; ext_surf_i++)
            if(current_wf->ext[ext_surf_i])
              extSurf[ext_surf_i] = current_state->e[ext_surf_i] == NULL ? NULL : init_fn(current_wf->ext[ext_surf_i].get(), orderSurf);
            else
              extSurf[ext_surf_i] = NULL;

          if(rungeKutta)
            for(int ext_surf_i = 0; ext_surf_i < this->RK_original_spaces_count; ext_surf_i++)
              u_extSurf[ext_surf_i]->add(extSurf[current_extCount - this->RK_original_spaces_count + ext_surf_i]);

          if(current_mat)
          {
            for(int current_mfsurf_i = 0; current_mfsurf_i < wf->mfsurf.size(); current_mfsurf_i++)
            {
              if(!form_to_be_assembled(current_wf->mfsurf[current_mfsurf_i], current_state))
                continue;

              int form_i = current_wf->mfsurf[current_mfsurf_i]->i;
              int form_j = current_wf->mfsurf[current_mfsurf_i]->j;

              assemble_matrix_form(current_wf->mfsurf[current_mfsurf_i], 
                cache_record->orderSurface[current_state->isurf], 
                cache_record->fnsSurface[current_state->isurf][form_j], 
                cache_record->fnsSurface[current_state->isurf][form_i], 
                extSurf, 
                u_extSurf,
                &current_alsSurface[form_i][current_state->isurf], 
                &current_alsSurface[form_j][current_state->isurf], 
                current_state, 
                cache_record->n_quadrature_pointsSurface[current_state->isurf], 
                cache_record->geometrySurface[current_state->isurf], 
                cache_record->jacobian_x_weightsSurface[current_state->isurf]);
            }
          }

          if(current_rhs)
          {
            for(int current_vfsurf_i = 0; current_vfsurf_i < wf->vfsurf.size(); current_vfsurf_i++)
            {
              if(!form_to_be_assembled(current_wf->vfsurf[current_vfsurf_i], current_state))
                continue;

              int form_i = current_wf->vfsurf[current_vfsurf_i]->i;

              assemble_vector_form(current_wf->vfsurf[current_vfsurf_i], 
                cache_record->orderSurface[current_state->isurf], 
                cache_record->fnsSurface[current_state->isurf][form_i], 
                extSurf, 
                u_extSurf, 
                &current_alsSurface[form_i][current_state->isurf], 
                current_state, 
                cache_record->n_quadrature_pointsSurface[current_state->isurf], 
                cache_record->geometrySurface[current_state->isurf], 
                cache_record->jacobian_x_weightsSurface[current_state->isurf]);
            }
          }

          if(current_u_ext)
          {
            for(int u_ext_surf_i = 0; u_ext_surf_i < prevNewtonSize; u_ext_surf_i++)
              if(current_u_ext[u_ext_surf_i] && current_state->e[u_ext_surf_i])
              {
                u_extSurf[u_ext_surf_i]->free_fn();
                delete u_extSurf[u_ext_surf_i];
              }
              delete [] u_extSurf;
          }

          for(int ext_surf_i = 0; ext_surf_i < current_extCount; ext_surf_i++)
            if(current_wf->ext[ext_surf_i] && current_state->e[ext_surf_i])
            {
              extSurf[ext_surf_i]->free_fn();
              delete extSurf[ext_surf_i];
            }
            delete [] extSurf;
        }

        for(unsigned int i = 0; i < this->spaces_size; i++)
          if(current_state->e[i])
            delete [] current_alsSurface[i];
      }

      if(current_alsSurface)
        delete [] current_alsSurface;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext,
      AsmList<Scalar>* current_als_i, AsmList<Scalar>* current_als_j, Traverse::State* current_state, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights)
    {
      bool surface_form = (dynamic_cast<MatrixFormVol<Scalar>*>(form) == NULL);

      double block_scaling_coefficient = this->block_scaling_coeff(form);

      bool tra = (form->i != form->j) && (form->sym != 0);
      bool sym = (form->i == form->j) && (form->sym == 1);

      // Assemble the local stiffness matrix for the form form.
      Scalar **local_stiffness_matrix = new_matrix<Scalar>(std::max(current_als_i->cnt, current_als_j->cnt));

      Func<Scalar>** local_ext = ext;
      // If the user supplied custom ext functions for this form.
      if(form->ext.size() > 0)
      {
        int local_ext_count = form->ext.size();
        local_ext = new Func<Scalar>*[local_ext_count];
        for(int ext_i = 0; ext_i < local_ext_count; ext_i++)
          if(form->ext[ext_i])
            local_ext[ext_i] = current_state->e[ext_i] == NULL ? NULL : init_fn(form->ext[ext_i].get(), order);
          else
            local_ext[ext_i] = NULL;
      }

      // Account for the previous time level solution previously inserted at the back of ext.
      if(rungeKutta)
        u_ext += form->u_ext_offset;

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_als_i->cnt; i++)
      {
        if(current_als_i->dof[i] < 0)
          continue;

        if((!tra || surface_form) && current_als_i->dof[i] < 0)
          continue;
        if(std::abs(current_als_i->coef[i]) < 1e-12)
          continue;
        if(!sym)
        {
          for (unsigned int j = 0; j < current_als_j->cnt; j++)
          {
            if(current_als_j->dof[j] >= 0)
            {
              // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
              if(std::abs(current_als_j->coef[j]) < 1e-12)
                continue;

              Func<double>* u = base_fns[j];
              Func<double>* v = test_fns[i];

              if(surface_form)
                local_stiffness_matrix[i][j] = 0.5 * block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];
              else
                local_stiffness_matrix[i][j] = block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];
            }
          }
        }
        // Symmetric block.
        else
        {
          for (unsigned int j = 0; j < current_als_j->cnt; j++)
          {
            if(j < i && current_als_j->dof[j] >= 0)
              continue;
            if(current_als_j->dof[j] >= 0)
            {
              // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
              if(std::abs(current_als_j->coef[j]) < 1e-12)
                continue;

              Func<double>* u = base_fns[j];
              Func<double>* v = test_fns[i];

              Scalar val = block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];

              local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
            }
          }
        }
      }

      // Insert the local stiffness matrix into the global one.
      current_mat->add(current_als_i->cnt, current_als_j->cnt, local_stiffness_matrix, current_als_i->dof, current_als_j->dof);

      // Insert also the off-diagonal (anti-)symmetric block, if required.
      if(tra)
      {
        if(form->sym < 0)
          chsgn(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt);
        transpose(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt);

        current_mat->add(current_als_j->cnt, current_als_i->cnt, local_stiffness_matrix, current_als_j->dof, current_als_i->dof);
      }

      if(form->ext.size() > 0)
      {
        for(int ext_i = 0; ext_i < form->ext.size(); ext_i++)
          if(form->ext[ext_i])
          {
            local_ext[ext_i]->free_fn();
            delete local_ext[ext_i];
          }
          delete [] local_ext;
      }

      if(rungeKutta)
        u_ext -= form->u_ext_offset;

      // Cleanup.
      delete [] local_stiffness_matrix;
    }

    template<typename Scalar>
    void DiscreteProblem<Scalar>::assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext, 
      AsmList<Scalar>* current_als_i, Traverse::State* current_state, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights)
    {
      bool surface_form = (dynamic_cast<VectorFormVol<Scalar>*>(form) == NULL);

      Func<Scalar>** local_ext = ext;

      // If the user supplied custom ext functions for this form.
      if(form->ext.size() > 0)
      {
        int local_ext_count = form->ext.size();
        local_ext = new Func<Scalar>*[local_ext_count];
        for(int ext_i = 0; ext_i < local_ext_count; ext_i++)
          if(form->ext[ext_i])
            local_ext[ext_i] = init_fn(form->ext[ext_i].get(), order);
          else
            local_ext[ext_i] = NULL;
      }

      // Account for the previous time level solution previously inserted at the back of ext.
      if(rungeKutta)
        u_ext += form->u_ext_offset;

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_als_i->cnt; i++)
      {
        if(current_als_i->dof[i] < 0)
          continue;

        // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
        if(std::abs(current_als_i->coef[i]) < 1e-12)
          continue;

        Func<double>* v = test_fns[i];

        Scalar val;
        if(surface_form)
          val = 0.5 * form->value(n_quadrature_points, jacobian_x_weights, u_ext, v, geometry, local_ext) * form->scaling_factor * current_als_i->coef[i];
        else
          val = form->value(n_quadrature_points, jacobian_x_weights, u_ext, v, geometry, local_ext) * form->scaling_factor * current_als_i->coef[i];

        current_rhs->add(current_als_i->dof[i], val);
      }

      if(form->ext.size() > 0)
      {
        for(int ext_i = 0; ext_i < form->ext.size(); ext_i++)
          if(form->ext[ext_i])
          {
            local_ext[ext_i]->free_fn();
            delete local_ext[ext_i];
          }
          delete [] local_ext;
      }

      if(rungeKutta)
        u_ext -= form->u_ext_offset;
    }
    
    template class HERMES_API DiscreteProblem<double>;
    template class HERMES_API DiscreteProblem<std::complex<double> >;
  }
}