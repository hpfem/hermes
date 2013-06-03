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

#include "discrete_problem/discrete_problem_selective_assembler.h"
#include "neighbor_search.h"

namespace Hermes
{
  namespace Hermes2D
  {

    template<typename Scalar>
    DiscreteProblemSelectiveAssembler<Scalar>::DiscreteProblemSelectiveAssembler()
    : sp_seq(NULL), 
      spaces_size(0),
      matrix_structure_reusable(false), 
      vector_structure_reusable(false)
    {
      for(int i = 0; i < 2; i++)
      {
        for(int j = 0; j < 2; j++)
        {
          state_reuse_kept[i][j] = NULL;
					markers_size[i][j] = 0;
        }
      }
    }

    template<typename Scalar>
    DiscreteProblemSelectiveAssembler<Scalar>::~DiscreteProblemSelectiveAssembler()
    {
      if(sp_seq)
        delete [] sp_seq;

      for(int i = 0; i < 2; i++)
      {
        for(int j = 0; j < 2; j++)
        {
          if(state_reuse_kept[i][j])
            ::free(state_reuse_kept[i][j]);
        }
      }
    }

    template<typename Scalar>
    bool DiscreteProblemSelectiveAssembler<Scalar>::prepare_sparse_structure(SparseMatrix<Scalar>* mat, Vector<Scalar>* rhs, Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Traverse::State**& states, int& num_states)
    {
      int ndof = Space<Scalar>::get_num_dofs(spaces);

      if(matrix_structure_reusable && mat)
      {
        Traverse::State** recalculated_states = (Traverse::State**)malloc(sizeof(Traverse::State*) * num_states);
        int num_recalculated_states = 0;
      
        // First check if any marker is reused.
        bool markers_to_reuse = false;
        for(int i = 0; i < this->markers_size[WeakForm<Scalar>::FormVol][WeakForm<Scalar>::MatrixForm]; i++)
          if(this->state_reuse_kept[WeakForm<Scalar>::FormVol][WeakForm<Scalar>::MatrixForm][i])
            markers_to_reuse = true;

        // Just zero the matrix if nothing can be reused.
        // Or if we also have to assemble the right-hand-side, we can not change the states.
        if(!markers_to_reuse || rhs)
          mat->zero();
        else
        {
          for(int state_i = 0; state_i < num_states; state_i++)
          {
            Traverse::State* current_state = states[state_i];
            int marker = current_state->rep->marker;

            if(this->state_reuse_kept[WeakForm<Scalar>::FormVol][WeakForm<Scalar>::MatrixForm][marker])
              delete current_state;
            else
            {
              recalculated_states[num_recalculated_states++] = current_state;

              for (unsigned int m = 0; m < spaces_size; m++)
              {
                if(current_state->e[m])
                {
                  AsmList<Scalar> am;
                  spaces[m]->get_element_assembly_list(current_state->e[m], &am);

                  // Pretend assembling of the element stiffness matrix.
                  for (unsigned int i = 0; i < am.cnt; i++)
                    mat->set_row_zero(am.dof[i]);
                }
              } 
            }
          }
        
          free(states);
          states = recalculated_states;
          num_states = num_recalculated_states;
          if(!num_states)
            return false;
        }
      }

      if(vector_structure_reusable && rhs)
      {
        if(rhs->length() == 0)
          rhs->alloc(ndof);
        else
          rhs->zero();
      }

      if(!matrix_structure_reusable && mat)
      {
        // Spaces have changed: create the matrix from scratch.
        matrix_structure_reusable = true;
        mat->free();
        mat->prealloc(ndof);

        AsmList<Scalar>* al = new AsmList<Scalar>[spaces_size];
        bool **blocks = this->wf->get_blocks(this->force_diagonal_blocks);
        
        // Loop through all elements.
        for(int state_i = 0; state_i < num_states; state_i++)
        {
          Traverse::State* current_state = states[state_i];

          // Obtain assembly lists for the element at all spaces.
          /// \todo do not get the assembly list again if the element was not changed.
          for (unsigned int i = 0; i < spaces_size; i++)
            if(current_state->e[i])
              spaces[i]->get_element_assembly_list(current_state->e[i], &(al[i]));

          if(this->wf->is_DG())
          {
            // Number of edges ( =  number of vertices).
            int num_edges = current_state->e[0]->nvert;

            // Allocation an array of arrays of neighboring elements for every mesh x edge.
            Element **** neighbor_elems_arrays = new Element ***[spaces_size];
            for(unsigned int i = 0; i < spaces_size; i++)
              neighbor_elems_arrays[i] = new Element **[num_edges];

            // The same, only for number of elements
            int ** neighbor_elems_counts = new int *[spaces_size];
            for(unsigned int i = 0; i < spaces_size; i++)
              neighbor_elems_counts[i] = new int[num_edges];

            // Get the neighbors.
            for(unsigned int el = 0; el < spaces_size; el++)
            {
              NeighborSearch<Scalar> ns(current_state->e[el], spaces[el]->get_mesh());

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
            for (unsigned int m = 0; m < spaces_size; m++)
            {
              for(unsigned int el = 0; el < spaces_size; el++)
              {
                for(int ed = 0; ed < num_edges; ed++)
                {
                  for(int neigh = 0; neigh < neighbor_elems_counts[el][ed]; neigh++)
                  {
                    if((blocks[m][el] || blocks[el][m]) && current_state->e[m])
                    {
                      AsmList<Scalar>*am = &(al[m]);
                      AsmList<Scalar>*an = new AsmList<Scalar>;
                      spaces[el]->get_element_assembly_list(neighbor_elems_arrays[el][ed][neigh], an);

                      // pretend assembling of the element stiffness matrix
                      // register nonzero elements
                      for (unsigned int i = 0; i < am->cnt; i++)
                      {
                        if(am->dof[i] >= 0)
                        {
                          for (unsigned int j = 0; j < an->cnt; j++)
                          {
                            if(an->dof[j] >= 0)
                            {
                              if(blocks[m][el]) mat->pre_add_ij(am->dof[i], an->dof[j]);
                              if(blocks[el][m]) mat->pre_add_ij(an->dof[j], am->dof[i]);
                            }
                          }
                        }
                      }
                      delete an;
                    }
                  }
                }
              }
            }

            // Deallocation an array of arrays of neighboring elements
            // for every mesh x edge.
            for(unsigned int el = 0; el < spaces_size; el++)
            {
              for(int ed = 0; ed < num_edges; ed++)
                delete [] neighbor_elems_arrays[el][ed];
              delete [] neighbor_elems_arrays[el];
            }
            delete [] neighbor_elems_arrays;

            // The same, only for number of elements.
            for(unsigned int el = 0; el < spaces_size; el++)
              delete [] neighbor_elems_counts[el];
            delete [] neighbor_elems_counts;
          }

          // Go through all equation-blocks of the local stiffness matrix.
          for (unsigned int m = 0; m < spaces_size; m++)
          {
            for (unsigned int n = 0; n < spaces_size; n++)
            {
              if(blocks[m][n] && current_state->e[m] && current_state->e[n])
              {
                AsmList<Scalar>*am = &(al[m]);
                AsmList<Scalar>*an = &(al[n]);

                // Pretend assembling of the element stiffness matrix.
                for (unsigned int i = 0; i < am->cnt; i++)
                {
                  if(am->dof[i] >= 0)
                    for (unsigned int j = 0; j < an->cnt; j++)
                      if(an->dof[j] >= 0)
                        mat->pre_add_ij(am->dof[i], an->dof[j]);
                }
              }
            }
          }
        }

        delete [] al;
        delete [] blocks;

        mat->alloc();
      }

      // WARNING: unlike Matrix<Scalar>::alloc(), Vector<Scalar>::alloc(ndof) frees the memory occupied
      // by previous vector before allocating
      if(!vector_structure_reusable && rhs)
      {
        vector_structure_reusable = true;
        rhs->alloc(ndof);
      }

      return true;
    }

    template<typename Scalar>
    void DiscreteProblemSelectiveAssembler<Scalar>::set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> >& spacesToSet)
    {
      if(!sp_seq)
      {
        // Internal variables settings.
        this->spaces_size = spacesToSet.size();
        sp_seq = new int[spaces_size];
        memset(sp_seq, -1, sizeof(int) * spaces_size);
      }
      else
      {
        for(unsigned int i = 0; i < spaces_size; i++)
        {
          int new_sp_seq = spacesToSet[i]->get_seq();

          if(new_sp_seq != sp_seq[i])
          {
            matrix_structure_reusable = false;
            vector_structure_reusable = false;
          }
          sp_seq[i] = new_sp_seq;
        }
      }

      int volume_markers_count = spacesToSet[0]->get_mesh()->get_element_markers_conversion().size();
      int surface_markers_count = spacesToSet[0]->get_mesh()->get_boundary_markers_conversion().size();

      this->alloc_recalculation_tables_spaces_settings(WeakForm<Scalar>::FormVol, WeakForm<Scalar>::MatrixForm, volume_markers_count);
      this->alloc_recalculation_tables_spaces_settings(WeakForm<Scalar>::FormVol, WeakForm<Scalar>::VectorForm, volume_markers_count);

      this->alloc_recalculation_tables_spaces_settings(WeakForm<Scalar>::FormSurf, WeakForm<Scalar>::MatrixForm, surface_markers_count);
      this->alloc_recalculation_tables_spaces_settings(WeakForm<Scalar>::FormSurf, WeakForm<Scalar>::VectorForm, surface_markers_count);
    }

    template<typename Scalar>
    void DiscreteProblemSelectiveAssembler<Scalar>::alloc_recalculation_tables_spaces_settings(typename WeakForm<Scalar>::FormIntegrationDimension dimension, typename WeakForm<Scalar>::FormEquationSide equation_side, int new_markers_count)
    {
      if(new_markers_count != this->markers_size[dimension][equation_side])
      {
        markers_size[dimension][equation_side] = new_markers_count;

        if(this->state_reuse_kept[dimension][equation_side])
          free(this->state_reuse_kept[dimension][equation_side]);
        
        this->state_reuse_kept[dimension][equation_side] = (bool*)calloc(markers_size[dimension][equation_side], sizeof(bool));
      }
      else
      {
        memset(this->state_reuse_kept[dimension][equation_side], 0, markers_size[dimension][equation_side] * sizeof(bool));
      }
    }

    template<typename Scalar>
    void DiscreteProblemSelectiveAssembler<Scalar>::set_weak_formulation(WeakForm<Scalar>* wf_)
    {
      Mixins::DiscreteProblemWeakForm<Scalar>::set_weak_formulation(wf_);

      this->matrix_structure_reusable = false;
      this->vector_structure_reusable = false;

      if(spaces_size == 0)
        return;

      this->alloc_recalculation_tables_weakform_settings(WeakForm<Scalar>::FormVol, WeakForm<Scalar>::MatrixForm);
      this->alloc_recalculation_tables_weakform_settings(WeakForm<Scalar>::FormVol, WeakForm<Scalar>::VectorForm);

      this->alloc_recalculation_tables_weakform_settings(WeakForm<Scalar>::FormSurf, WeakForm<Scalar>::MatrixForm);
      this->alloc_recalculation_tables_weakform_settings(WeakForm<Scalar>::FormSurf, WeakForm<Scalar>::VectorForm);
    }

    template<typename Scalar>
    void DiscreteProblemSelectiveAssembler<Scalar>::alloc_recalculation_tables_weakform_settings(typename WeakForm<Scalar>::FormIntegrationDimension dimension, typename WeakForm<Scalar>::FormEquationSide equation_side)
    {
      if(this->state_reuse_kept[dimension][equation_side])
        memset(this->state_reuse_kept[dimension][equation_side], 0, markers_size[dimension][equation_side] * sizeof(bool));
    }

    template<typename Scalar>
    bool DiscreteProblemSelectiveAssembler<Scalar>::form_to_be_assembled(MatrixForm<Scalar>* form, Traverse::State* current_state)
    {
      if(current_state->e[form->i] && current_state->e[form->j])
      {
        if(fabs(form->scaling_factor) < Hermes::epsilon)
          return false;

        // If a block scaling table is provided, and if the scaling coefficient
        // A_mn for this block is zero, then the form does not need to be assembled.
        if(this->block_weights)
          if(fabs(this->block_weights->get_A(form->i, form->j)) < Hermes::epsilon)
            return false;
        return true;
      }
      return false;
    }

    template<typename Scalar>
    bool DiscreteProblemSelectiveAssembler<Scalar>::form_to_be_assembled(MatrixFormVol<Scalar>* form, Traverse::State* current_state)
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
    bool DiscreteProblemSelectiveAssembler<Scalar>::form_to_be_assembled(MatrixFormSurf<Scalar>* form, Traverse::State* current_state)
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
    bool DiscreteProblemSelectiveAssembler<Scalar>::form_to_be_assembled(MatrixFormDG<Scalar>* form, Traverse::State* current_state)
    {
      return form_to_be_assembled((MatrixForm<Scalar>*)form, current_state);
    }

    template<typename Scalar>
    bool DiscreteProblemSelectiveAssembler<Scalar>::form_to_be_assembled(VectorForm<Scalar>* form, Traverse::State* current_state)
    {
      if(!current_state->e[form->i])
        return false;
      if(fabs(form->scaling_factor) < Hermes::epsilon)
        return false;

      return true;
    }

    template<typename Scalar>
    bool DiscreteProblemSelectiveAssembler<Scalar>::form_to_be_assembled(VectorFormVol<Scalar>* form, Traverse::State* current_state)
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
    bool DiscreteProblemSelectiveAssembler<Scalar>::form_to_be_assembled(VectorFormSurf<Scalar>* form, Traverse::State* current_state)
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
    bool DiscreteProblemSelectiveAssembler<Scalar>::form_to_be_assembled(VectorFormDG<Scalar>* form, Traverse::State* current_state)
    {
      return form_to_be_assembled((VectorForm<Scalar>*)form, current_state);
    }

    template class HERMES_API DiscreteProblemSelectiveAssembler<double>;
    template class HERMES_API DiscreteProblemSelectiveAssembler<std::complex<double> >;
  }
}