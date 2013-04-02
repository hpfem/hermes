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

#include "discrete_problem_linear.h"
#include <iostream>
#include <algorithm>
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
    DiscreteProblemLinear<Scalar>::DiscreteProblemLinear(const WeakForm<Scalar>* wf, Hermes::vector<SpaceSharedPtr<Scalar> > spaces) : DiscreteProblem<Scalar>(wf, spaces)
    {
      this->is_linear = true;
    }

    template<typename Scalar>
    DiscreteProblemLinear<Scalar>::DiscreteProblemLinear(const WeakForm<Scalar>* wf, SpaceSharedPtr<Scalar> space) : DiscreteProblem<Scalar>(wf, space)
    {
      this->is_linear = true;
    }

    template<typename Scalar>
    DiscreteProblemLinear<Scalar>::DiscreteProblemLinear() : DiscreteProblem<Scalar>()
    {
      this->is_linear = true;
    }

    template<typename Scalar>
    DiscreteProblemLinear<Scalar>::~DiscreteProblemLinear()
    {
    }

    template<typename Scalar>
    void DiscreteProblemLinear<Scalar>::assemble(SparseMatrix<Scalar>* mat,
      Vector<Scalar>* rhs,
      bool force_diagonal_blocks,
      Table* block_weights)
    {
      // Check.
      this->check();

      // Important, sets the current caughtException to NULL.
      this->caughtException = NULL;

      this->current_mat = mat;
      this->current_rhs = rhs;
      this->current_force_diagonal_blocks = force_diagonal_blocks;
      this->current_block_weights = block_weights;

      // Local number of threads - to avoid calling it over and over again, and against faults caused by the
      // value being changed while assembling.
      int num_threads_used = Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads);

      // Check that the block scaling table have proper dimension.
      if(block_weights)
        if(block_weights->get_size() != this->wf->get_neq())
          throw Exceptions::LengthException(6, block_weights->get_size(),this-> wf->get_neq());

      // Creating matrix sparse structure.
      this->create_sparse_structure();

      // Initial check of meshes and spaces.
      for(unsigned int ext_i = 0; ext_i < this->wf->get_ext().size(); ext_i++)
      {
        if(!this->wf->get_ext()[ext_i]->isOkay())
          throw Hermes::Exceptions::Exception("Ext function %d is not okay in assemble().", ext_i);
      }

      // Structures that cloning will be done into.
      RefMap*** refmaps = new RefMap**[num_threads_used];
      AsmList<Scalar>*** als = new AsmList<Scalar>**[num_threads_used];
      AsmList<Scalar>**** alsSurface = new AsmList<Scalar>***[num_threads_used];
      WeakForm<Scalar>** weakforms = new WeakForm<Scalar>*[num_threads_used];
      PrecalcShapeset*** pss = new PrecalcShapeset**[num_threads_used];

      // Fill these structures.
      this->init_assembling(NULL, pss, refmaps, NULL, als, alsSurface, weakforms, num_threads_used);

      // Vector of meshes.
      Hermes::vector<MeshSharedPtr > meshes;
      for(unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
        meshes.push_back(this->spaces[space_i]->get_mesh());
      for(unsigned int ext_i = 0; ext_i < this->wf->ext.size(); ext_i++)
        meshes.push_back(this->wf->ext[ext_i]->get_mesh());
      for(unsigned int form_i = 0; form_i < this->wf->get_forms().size(); form_i++)
        for(unsigned int ext_i = 0; ext_i < this->wf->get_forms()[form_i]->ext.size(); ext_i++)
          if(this->wf->get_forms()[form_i]->ext[ext_i])
            meshes.push_back(this->wf->get_forms()[form_i]->ext[ext_i]->get_mesh());

      Traverse trav_master(true);
      int num_states;
      Traverse::State** states = trav_master.get_states(meshes, num_states);

      Hermes::vector<Transformable *>* fns = new Hermes::vector<Transformable *>[num_threads_used];
      for(unsigned int i = 0; i < num_threads_used; i++)
      {
        for (unsigned j = 0; j < this->spaces_size; j++)
        {
          fns[i].push_back(pss[i][j]);
        }
        for (unsigned j = 0; j < this->wf->ext.size(); j++)
        {
          fns[i].push_back(weakforms[i]->ext[j].get());
          weakforms[i]->ext[j]->set_quad_2d(&g_quad_2d_std);
        }
        for(unsigned int form_i = 0; form_i < this->wf->get_forms().size(); form_i++)
        {
          for(unsigned int ext_i = 0; ext_i < this->wf->get_forms()[form_i]->ext.size(); ext_i++)
            if(this->wf->get_forms()[form_i]->ext[ext_i])
            {
              fns[i].push_back(weakforms[i]->get_forms()[form_i]->ext[ext_i].get());
              weakforms[i]->get_forms()[form_i]->ext[ext_i]->set_quad_2d(&g_quad_2d_std);
            }
        }
      }

#pragma omp parallel num_threads(num_threads_used)
      {
        int thread_number = omp_get_thread_num();
        int start = (num_states / num_threads_used) * thread_number;
        int end = (num_states / num_threads_used) * (thread_number + 1);
        if(thread_number == num_threads_used - 1)
          end = num_states;

        AsmList<Scalar>** current_als = als[thread_number];
        AsmList<Scalar>*** current_als_surface = alsSurface[thread_number];
        RefMap** current_refmaps = refmaps[thread_number];
        WeakForm<Scalar>* current_weakform = weakforms[thread_number];
        PrecalcShapeset** current_pss = pss[thread_number];

        PrecalcShapeset** current_spss = new PrecalcShapeset*[this->spaces_size];
        if(this->DG_matrix_forms_present || this->DG_vector_forms_present)
          for (unsigned int j = 0; j < this->spaces_size; j++)
            current_spss[j] = new PrecalcShapeset(current_pss[j]);

        int order;

        for(int state_i = start; state_i < end; state_i++)
        {
          if(this->caughtException)
            break;
          try
          {
            Traverse::State* current_state = states[state_i];

            for(int j = 0; j < fns[thread_number].size(); j++)
            {
              if(current_state->e[j])
              {
                fns[thread_number][j]->set_active_element(current_state->e[j]);
                fns[thread_number][j]->set_transform(current_state->sub_idx[j]);
              }
            }
            
            bool state_skipped = true;
            for(int j = 0; j < this->spaces_size; j++)
            {
              if(current_state->e[j])
              {
                this->spaces[j]->get_element_assembly_list(current_state->e[j], current_als[j]);
                if(this->DG_matrix_forms_present || this->DG_vector_forms_present)
                {
                  current_spss[j]->set_active_element(current_state->e[j]);
                  current_spss[j]->set_master_transform();
                }
                current_refmaps[j]->set_active_element(current_state->e[j]);
                current_refmaps[j]->force_transform(current_pss[j]->get_transform(), current_pss[j]->get_ctm());
                state_skipped = false;
              }
            }

            if(state_skipped)
              continue;

            typename DiscreteProblemCache<Scalar>::CacheRecord* cache_record;
            if(!this->do_not_use_cache)
              cache_record = this->get_state_cache(current_state, current_pss, current_refmaps, NULL, current_als, current_als_surface, current_weakform, order);
            else
            {
              cache_record = new typename DiscreteProblemCache<Scalar>::CacheRecord();
              order = this->calculate_order(current_state, current_refmaps, NULL, current_weakform);
              cache_record->init(this->spaces, current_state, current_pss, current_refmaps, NULL, current_als, current_als_surface, current_weakform, order);
            }

            this->assemble_one_state(cache_record, current_refmaps, NULL, current_als, current_state, current_weakform);

            if(this->DG_matrix_forms_present || this->DG_vector_forms_present)
              this->assemble_one_DG_state(current_pss, current_spss, current_refmaps, NULL, current_als, current_state, current_weakform->mfDG, current_weakform->vfDG, &fns[thread_number].front(), current_weakform);

            if(this->do_not_use_cache)
              delete cache_record;
          }
          catch(Hermes::Exceptions::Exception& e)
          {
            if(this->caughtException == NULL)
              this->caughtException = e.clone();
          }
          catch(std::exception& e)
          {
            if(this->caughtException == NULL)
              this->caughtException = new std::exception(e);
          }
        }

        if(this->DG_matrix_forms_present || this->DG_vector_forms_present)
         for (unsigned int j = 0; j < this->spaces_size; j++)
            delete current_spss[j];
        delete [] current_spss;
      }

      this->cache.free_unused();

      this->deinit_assembling(pss, refmaps, NULL, als, alsSurface, weakforms, num_threads_used);

      for(int i = 0; i < num_states; i++)
        delete states[i];
      free(states);

      for(unsigned int i = 0; i < num_threads_used; i++)
      {
        fns[i].clear();
      }
      delete [] fns;

      /// \todo Should this be really here? Or in assemble()?
      if(this->current_mat)
        this->current_mat->finish();
      if(this->current_rhs)
        this->current_rhs->finish();

      if(this->DG_matrix_forms_present || this->DG_vector_forms_present)
      {
        Element* element_to_set_nonvisited;
        for(unsigned int mesh_i = 0; mesh_i < meshes.size(); mesh_i++)
          for_all_elements(element_to_set_nonvisited, meshes[mesh_i])
          element_to_set_nonvisited->visited = false;
      }

      Element* e;
      for(unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
      {
        for_all_active_elements(e, this->spaces[space_i]->get_mesh())
          this->spaces[space_i]->edata[e->id].changed_in_last_adaptation = false;
      }

      if(this->caughtException)
        throw *(this->caughtException);
    }

    template<typename Scalar>
    void DiscreteProblemLinear<Scalar>::assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext,
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
            // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
            if(std::abs(current_als_j->coef[j]) < 1e-12)
              continue;

            Func<double>* u = base_fns[j];
            Func<double>* v = test_fns[i];

            if(current_als_j->dof[j] >= 0)
            {
              if(surface_form)
                local_stiffness_matrix[i][j] = 0.5 * block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];
              else
                local_stiffness_matrix[i][j] = block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];
            }
            else
            {
              {
                if(surface_form)
                  this->current_rhs->add(current_als_i->dof[i], -0.5 * block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i]);
                else
                  this->current_rhs->add(current_als_i->dof[i], -block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i]);
              }
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
            // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
            if(std::abs(current_als_j->coef[j]) < 1e-12)
              continue;

            Func<double>* u = base_fns[j];
            Func<double>* v = test_fns[i];

            Scalar val = block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];

            if(current_als_j->dof[j] >= 0)
              local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
            else
            {
              this->current_rhs->add(current_als_i->dof[i], -val);
            }
          }
        }
      }

      // Insert the local stiffness matrix into the global one.
      this->current_mat->add(current_als_i->cnt, current_als_j->cnt, local_stiffness_matrix, current_als_i->dof, current_als_j->dof);

      // Insert also the off-diagonal (anti-)symmetric block, if required.
      if(tra)
      {
        if(form->sym < 0)
          chsgn(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt);
        transpose(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt);

        this->current_mat->add(current_als_j->cnt, current_als_i->cnt, local_stiffness_matrix, current_als_j->dof, current_als_i->dof);

        // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
        for (unsigned int j = 0; j < current_als_i->cnt; j++)
          if(current_als_i->dof[j] < 0)
            for (unsigned int i = 0; i < current_als_j->cnt; i++)
              if(current_als_j->dof[i] >= 0)
                this->current_rhs->add(current_als_j->dof[i], -local_stiffness_matrix[i][j]);
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

      // Cleanup.
      delete [] local_stiffness_matrix;
    }

    template class HERMES_API DiscreteProblemLinear<double>;
    template class HERMES_API DiscreteProblemLinear<std::complex<double> >;
  }
}
