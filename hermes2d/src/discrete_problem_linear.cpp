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
    DiscreteProblemLinear<Scalar>::DiscreteProblemLinear(const WeakForm<Scalar>* wf, Hermes::vector<const Space<Scalar> *> spaces) : DiscreteProblem<Scalar>(wf, spaces)
    {
    }

    template<typename Scalar>
    DiscreteProblemLinear<Scalar>::DiscreteProblemLinear(const WeakForm<Scalar>* wf, const Space<Scalar>* space) : DiscreteProblem<Scalar>(wf, space)
    {
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
      // Important, sets the current caughtException to NULL.
      this->caughtException = NULL;

      this->current_mat = mat;
      this->current_rhs = rhs;
      this->current_force_diagonal_blocks = force_diagonal_blocks;
      this->current_block_weights = block_weights;

      // Check that the block scaling table have proper dimension.
      if(block_weights != NULL)
        if(block_weights->get_size() != this->wf->get_neq())
          throw Exceptions::LengthException(6, block_weights->get_size(),this-> wf->get_neq());

      // Creating matrix sparse structure.
      this->create_sparse_structure();

      Hermes::vector<MeshFunction<Scalar>*> ext_functions;
      for(unsigned int form_i = 0; form_i < this->wf->mfvol.size(); form_i++)
        for(unsigned int ext_i = 0; ext_i < this->wf->mfvol.at(form_i)->ext.size(); ext_i++)
          ext_functions.push_back(this->wf->mfvol.at(form_i)->ext[ext_i]);
      for(unsigned int form_i = 0; form_i < this->wf->mfsurf.size(); form_i++)
        for(unsigned int ext_i = 0; ext_i < this->wf->mfsurf.at(form_i)->ext.size(); ext_i++)
          ext_functions.push_back(this->wf->mfsurf.at(form_i)->ext[ext_i]);
      for(unsigned int form_i = 0; form_i < this->wf->vfvol.size(); form_i++)
        for(unsigned int ext_i = 0; ext_i < this->wf->vfvol.at(form_i)->ext.size(); ext_i++)
          ext_functions.push_back(this->wf->vfvol.at(form_i)->ext[ext_i]);
      for(unsigned int form_i = 0; form_i < this->wf->vfsurf.size(); form_i++)
        for(unsigned int ext_i = 0; ext_i < this->wf->vfsurf.at(form_i)->ext.size(); ext_i++)
          ext_functions.push_back(this->wf->vfsurf.at(form_i)->ext[ext_i]);

      // Structures that cloning will be done into.
      PrecalcShapeset*** pss = new PrecalcShapeset**[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
      PrecalcShapeset*** spss = new PrecalcShapeset**[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
      RefMap*** refmaps = new RefMap**[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
      Solution<Scalar>*** u_ext = new Solution<Scalar>**[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
      AsmList<Scalar>*** als = new AsmList<Scalar>**[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
      MeshFunction<Scalar>*** ext = new MeshFunction<Scalar>**[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
      Hermes::vector<MatrixFormVol<Scalar>*>* mfvol = new Hermes::vector<MatrixFormVol<Scalar>*>[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
      Hermes::vector<MatrixFormSurf<Scalar>*>* mfsurf = new Hermes::vector<MatrixFormSurf<Scalar>*>[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
      Hermes::vector<VectorFormVol<Scalar>*>* vfvol = new Hermes::vector<VectorFormVol<Scalar>*>[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
      Hermes::vector<VectorFormSurf<Scalar>*>* vfsurf = new Hermes::vector<VectorFormSurf<Scalar>*>[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];

      // Fill these structures.
      this->init_assembling(NULL, pss, spss, refmaps, u_ext, als, ext_functions, ext, mfvol, mfsurf, vfvol, vfsurf);

      // Vector of meshes.
      Hermes::vector<const Mesh*> meshes;
      for(unsigned int space_i = 0; space_i < this->spaces.size(); space_i++)
        meshes.push_back(this->spaces[space_i]->get_mesh());
      for (unsigned j = 0; j < ext_functions.size(); j++)
        meshes.push_back(ext_functions[j]->get_mesh());

      Traverse trav_master(true);
      unsigned int num_states = trav_master.get_num_states(meshes);

      trav_master.begin(meshes.size(), &(meshes.front()));

      Traverse* trav = new Traverse[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
      Hermes::vector<Transformable *>* fns = new Hermes::vector<Transformable *>[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
      for(unsigned int i = 0; i < Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads); i++)
      {
        for (unsigned j = 0; j < this->spaces.size(); j++)
          fns[i].push_back(pss[i][j]);
        for (unsigned j = 0; j < ext_functions.size(); j++)
        {
          fns[i].push_back(ext[i][j]);
          ext[i][j]->set_quad_2d(&g_quad_2d_std);
        }
        trav[i].begin(meshes.size(), &(meshes.front()), &(fns[i].front()));
        trav[i].stack = trav_master.stack;
      }

      int state_i;

      PrecalcShapeset** current_pss;
      PrecalcShapeset** current_spss;
      RefMap** current_refmaps;
      AsmList<Scalar>** current_als;

      MatrixFormVol<Scalar>** current_mfvol;
      MatrixFormSurf<Scalar>** current_mfsurf;
      VectorFormVol<Scalar>** current_vfvol;
      VectorFormSurf<Scalar>** current_vfsurf;

#define CHUNKSIZE 1
      int num_threads_used = Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads);
#pragma omp parallel shared(trav_master, mat, rhs) private(state_i, current_pss, current_spss, current_refmaps, current_als, current_mfvol, current_mfsurf, current_vfvol, current_vfsurf) num_threads(num_threads_used)
      {
#pragma omp for schedule(dynamic, CHUNKSIZE)
        for(state_i = 0; state_i < num_states; state_i++)
        {
          try
          {
            Traverse::State current_state;

  #pragma omp critical (get_next_state)
            current_state = trav[omp_get_thread_num()].get_next_state(&trav_master.top, &trav_master.id);

            current_pss = pss[omp_get_thread_num()];
            current_spss = spss[omp_get_thread_num()];
            current_refmaps = refmaps[omp_get_thread_num()];
            current_als = als[omp_get_thread_num()];

            current_mfvol = mfvol[omp_get_thread_num()].size() == 0 ? NULL : &(mfvol[omp_get_thread_num()].front());
            current_mfsurf = mfsurf[omp_get_thread_num()].size() == 0 ? NULL : &(mfsurf[omp_get_thread_num()].front());
            current_vfvol = vfvol[omp_get_thread_num()].size() == 0 ? NULL : &(vfvol[omp_get_thread_num()].front());
            current_vfsurf = vfsurf[omp_get_thread_num()].size() == 0 ? NULL : &(vfsurf[omp_get_thread_num()].front());

            // One state is a collection of (virtual) elements sharing
            // the same physical location on (possibly) different meshes.
            // This is then the same element of the virtual union mesh.
            // The proper sub-element mappings to all the functions of
            // this stage is supplied by the function Traverse::get_next_state()
            // called in the while loop.
            this->assemble_one_state(current_pss, current_spss, current_refmaps, NULL, current_als, &current_state, current_mfvol, current_mfsurf, current_vfvol, current_vfsurf);

            if(this->DG_matrix_forms_present || this->DG_vector_forms_present)
              this->assemble_one_DG_state(current_pss, current_spss, current_refmaps, NULL, current_als, &current_state, current_mfsurf, current_vfsurf, trav[omp_get_thread_num()].fn);
          }
          catch(Hermes::Exceptions::Exception& e)
          {
            if(this->caughtException == NULL)
              this->caughtException = e.clone();
          }
          catch(std::exception& e)
          {
            if(this->caughtException == NULL)
              this->caughtException = new Hermes::Exceptions::Exception(e.what());
          }
        }
      }

      deinit_assembling(pss, spss, refmaps, u_ext, als, ext_functions, ext, mfvol, mfsurf, vfvol, vfsurf);

      trav_master.finish();
      for(unsigned int i = 0; i < Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads); i++)
        trav[i].finish();

      for(unsigned int i = 0; i < Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads); i++)
      {
        fns[i].clear();
      }
      delete [] fns;
      delete [] trav;

      /// \todo Should this be really here? Or in assemble()?
      if(this->current_mat != NULL)
        this->current_mat->finish();
      if(this->current_rhs != NULL)
        this->current_rhs->finish();

      if(this->DG_matrix_forms_present || this->DG_vector_forms_present)
      {
        Element* element_to_set_nonvisited;
        for(unsigned int mesh_i = 0; mesh_i < meshes.size(); mesh_i++)
          for_all_elements(element_to_set_nonvisited, meshes[mesh_i])
          element_to_set_nonvisited->visited = false;
      }

      if(this->caughtException != NULL)
        throw *(this->caughtException);
    }

    template<typename Scalar>
    void DiscreteProblemLinear<Scalar>::assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, Traverse::State* current_state)
    {
      bool surface_form = (dynamic_cast<MatrixFormVol<Scalar>*>(form) == NULL);

      double block_scaling_coef = this->block_scaling_coeff(form);

      bool tra = (form->i != form->j) && (form->sym != 0);
      bool sym = (form->i == form->j) && (form->sym == 1);

      // Assemble the local stiffness matrix for the form form.
      Scalar **local_stiffness_matrix = new_matrix<Scalar>(std::max(current_als[form->i]->cnt, current_als[form->j]->cnt));

      // Init external functions.
      Func<Scalar>** u_ext = new Func<Scalar>*[this->RungeKutta ? this->RK_original_spaces_count : this->wf->get_neq() - form->u_ext_offset];
      ExtData<Scalar> ext;
      init_ext(form, u_ext, &ext, order, current_u_ext, current_state);

      // Add the previous time level solution previously inserted at the back of ext.
      if(this->RungeKutta)
        for(int ext_i = 0; ext_i < this->RK_original_spaces_count; ext_i++)
          u_ext[ext_i]->add(*ext.fn[form->ext.size() - this->RK_original_spaces_count + ext_i]);

      // Init geometry.
      int n_quadrature_points;
      Geom<double>* geometry = NULL;
      double* jacobian_x_weights = NULL;
      if(surface_form)
        n_quadrature_points = init_surface_geometry_points(current_refmaps[form->i], order, current_state, geometry, jacobian_x_weights);
      else
        n_quadrature_points = init_geometry_points(current_refmaps[form->i], order, geometry, jacobian_x_weights);

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_als[form->i]->cnt; i++)
      {
        if(current_als[form->i]->dof[i] < 0)
          continue;

        if((!tra || surface_form) && current_als[form->i]->dof[i] < 0)
          continue;
        if(std::abs(current_als[form->i]->coef[i]) < 1e-12)
          continue;
        if(!sym)
        {
          for (unsigned int j = 0; j < current_als[form->j]->cnt; j++)
          {
            // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
            if(std::abs(current_als[form->j]->coef[j]) < 1e-12)
              continue;

            Func<double>* u = base_fns[j];
            Func<double>* v = test_fns[i];

            if(current_als[form->j]->dof[j] >= 0)
            {
              if(surface_form)
                local_stiffness_matrix[i][j] = 0.5 * block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als[form->j]->coef[j] * current_als[form->i]->coef[i];
              else
                local_stiffness_matrix[i][j] = block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als[form->j]->coef[j] * current_als[form->i]->coef[i];
            }
            else
            {
#pragma omp critical (rhs)
              {
                if(surface_form)
                  this->current_rhs->add(current_als[form->i]->dof[i], -0.5 * block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als[form->j]->coef[j] * current_als[form->i]->coef[i]);
                else
                  this->current_rhs->add(current_als[form->i]->dof[i], -block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als[form->j]->coef[j] * current_als[form->i]->coef[i]);
              }
            }
          }
        }
        // Symmetric block.
        else
        {
          for (unsigned int j = 0; j < current_als[form->j]->cnt; j++)
          {
            if(j < i && current_als[form->j]->dof[j] >= 0)
              continue;
            // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
            if(std::abs(current_als[form->j]->coef[j]) < 1e-12)
              continue;

            Func<double>* u = base_fns[j];
            Func<double>* v = test_fns[i];

            Scalar val = block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als[form->j]->coef[j] * current_als[form->i]->coef[i];

            if(current_als[form->j]->dof[j] >= 0)
              local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
            else
            {
#pragma omp critical (mat)
              this->current_rhs->add(current_als[form->i]->dof[i], -val);
            }
          }
        }
      }

      // Insert the local stiffness matrix into the global one.
#pragma omp critical (mat)
      this->current_mat->add(current_als[form->i]->cnt, current_als[form->j]->cnt, local_stiffness_matrix, current_als[form->i]->dof, current_als[form->j]->dof);

      // Insert also the off-diagonal (anti-)symmetric block, if required.
      if(tra)
      {
        if(form->sym < 0)
          chsgn(local_stiffness_matrix, current_als[form->i]->cnt, current_als[form->j]->cnt);
        transpose(local_stiffness_matrix, current_als[form->i]->cnt, current_als[form->j]->cnt);
#pragma omp critical (mat)
        this->current_mat->add(current_als[form->j]->cnt, current_als[form->i]->cnt, local_stiffness_matrix, current_als[form->j]->dof, current_als[form->i]->dof);

        // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
        for (unsigned int j = 0; j < current_als[form->i]->cnt; j++)
          if(current_als[form->i]->dof[j] < 0)
            for (unsigned int i = 0; i < current_als[form->j]->cnt; i++)
              if(current_als[form->j]->dof[i] >= 0)
                this->current_rhs->add(current_als[form->j]->dof[i], -local_stiffness_matrix[i][j]);
      }

      // Cleanup.
      deinit_ext(form, u_ext, &ext);
      delete [] local_stiffness_matrix;
      delete [] jacobian_x_weights;
      geometry->free();
      delete geometry;
    }

    template<typename Scalar>
    void DiscreteProblemLinear<Scalar>::assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, Solution<Scalar>** current_u_ext, 
      AsmList<Scalar>* current_als_i, AsmList<Scalar>* current_als_j, Traverse::State* current_state, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights)
    {
      bool surface_form = (dynamic_cast<MatrixFormVol<Scalar>*>(form) == NULL);

      double block_scaling_coef = this->block_scaling_coeff(form);

      bool tra = (form->i != form->j) && (form->sym != 0);
      bool sym = (form->i == form->j) && (form->sym == 1);

      // Assemble the local stiffness matrix for the form form.
      Scalar **local_stiffness_matrix = new_matrix<Scalar>(std::max(current_als_i->cnt, current_als_j->cnt));

      // Init external functions.
      Func<Scalar>** u_ext = new Func<Scalar>*[this->RungeKutta ? this->RK_original_spaces_count : this->wf->get_neq() - form->u_ext_offset];
      ExtData<Scalar> ext;
      init_ext(form, u_ext, &ext, order, current_u_ext, current_state);

      // Add the previous time level solution previously inserted at the back of ext.
      if(this->RungeKutta)
        for(int ext_i = 0; ext_i < this->RK_original_spaces_count; ext_i++)
          u_ext[ext_i]->add(*ext.fn[form->ext.size() - this->RK_original_spaces_count + ext_i]);

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
                local_stiffness_matrix[i][j] = 0.5 * block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];
              else
                local_stiffness_matrix[i][j] = block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];
            }
            else
            {
#pragma omp critical (rhs)
              {
                if(surface_form)
                  this->current_rhs->add(current_als_i->dof[i], -0.5 * block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i]);
                else
                  this->current_rhs->add(current_als_i->dof[i], -block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i]);
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

            Scalar val = block_scaling_coeff(form) * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, &ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];

            if(current_als_j->dof[j] >= 0)
              local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
            else
            {
#pragma omp critical (mat)
              this->current_rhs->add(current_als_i->dof[i], -val);
            }
          }
        }
      }

      // Insert the local stiffness matrix into the global one.
#pragma omp critical (mat)
      this->current_mat->add(current_als_i->cnt, current_als_j->cnt, local_stiffness_matrix, current_als_i->dof, current_als_j->dof);

      // Insert also the off-diagonal (anti-)symmetric block, if required.
      if(tra)
      {
        if(form->sym < 0)
          chsgn(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt);
        transpose(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt);
#pragma omp critical (mat)
        this->current_mat->add(current_als_j->cnt, current_als_i->cnt, local_stiffness_matrix, current_als_j->dof, current_als_i->dof);

        // Linear problems only: Subtracting Dirichlet lift contribution from the RHS:
        for (unsigned int j = 0; j < current_als_i->cnt; j++)
          if(current_als_i->dof[j] < 0)
            for (unsigned int i = 0; i < current_als_j->cnt; i++)
              if(current_als_j->dof[i] >= 0)
                this->current_rhs->add(current_als_j->dof[i], -local_stiffness_matrix[i][j]);
      }

      // Cleanup.
      deinit_ext(form, u_ext, &ext);
      delete [] local_stiffness_matrix;
    }

    template class HERMES_API DiscreteProblemLinear<double>;
    template class HERMES_API DiscreteProblemLinear<std::complex<double> >;
  }
}