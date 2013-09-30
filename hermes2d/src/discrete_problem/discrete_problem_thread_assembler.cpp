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

#include "discrete_problem/discrete_problem_thread_assembler.h"
#include "discrete_problem/discrete_problem_selective_assembler.h"
#include "shapeset/precalc.h"
#include "function/solution.h"
#include "weakform/weakform.h"
#include "function/exact_solution.h"

namespace Hermes
{
  namespace Hermes2D
  {
#ifdef _DEBUG
    static int cache_searches = 0;
    static int cache_record_found = 0;
    static int cache_record_found_reinit = 0;
    static int cache_record_not_found = 0;
#endif

    template<typename Scalar>
    DiscreteProblemThreadAssembler<Scalar>::DiscreteProblemThreadAssembler(DiscreteProblemSelectiveAssembler<Scalar>* selectiveAssembler) : 
      pss(NULL), refmaps(NULL), u_ext(NULL), als(NULL), alsSurface(NULL), 
      selectiveAssembler(selectiveAssembler), integrationOrderCalculator(selectiveAssembler)
    {
    }

    template<typename Scalar>
    DiscreteProblemThreadAssembler<Scalar>::~DiscreteProblemThreadAssembler()
    {
      this->free();
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_spaces(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces)
    {
      this->free_spaces();

      this->spaces_size = spaces.size();

      pss = new PrecalcShapeset*[spaces_size];
      refmaps = new RefMap*[spaces_size];
      als = new AsmList<Scalar>*[spaces_size];
      alsSurface = new AsmList<Scalar>**[spaces_size];

      for (unsigned int j = 0; j < spaces_size; j++)
      {
        pss[j] = new PrecalcShapeset(spaces[j]->shapeset);
        refmaps[j] = new RefMap();
        refmaps[j]->set_quad_2d(&g_quad_2d_std);
        als[j] = new AsmList<Scalar>();
        alsSurface[j] = new AsmList<Scalar>*[H2D_MAX_NUMBER_EDGES];
        for(unsigned int k = 0; k < H2D_MAX_NUMBER_EDGES; k++)
          alsSurface[j][k] = new AsmList<Scalar>();
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::set_weak_formulation(WeakForm<Scalar>* wf_)
    {
      this->free_weak_formulation();

      this->wf = wf_->clone();
      this->wf->cloneMembers(wf_);
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_u_ext(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Solution<Scalar>** u_ext_sln)
    {
#ifdef _DEBUG
      assert(this->spaces_size == spaces.size() && this->pss);
#endif
      free_u_ext();
      u_ext = new Solution<Scalar>*[spaces_size];

      for (unsigned int j = 0; j < spaces_size; j++)
      {
        if(u_ext_sln)
        {
          u_ext[j] = new Solution<Scalar>(spaces[j]->get_mesh());
          u_ext[j]->copy(u_ext_sln[j]);
        }
        else
        {
          if(spaces[j]->get_shapeset()->get_num_components() == 1)
            u_ext[j] = new ZeroSolution<Scalar>(spaces[j]->get_mesh());
          else
            u_ext[j] = new ZeroSolutionVector<Scalar>(spaces[j]->get_mesh());
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_assembling(Solution<Scalar>** u_ext_sln, const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, bool nonlinear_, bool add_dirichlet_lift_)
    {
      if(this->report_cache_hits_and_misses)
        this->zero_cache_hits_and_misses();

      this->nonlinear = nonlinear_;
      this->add_dirichlet_lift = add_dirichlet_lift_;

      fns.clear();
      for (unsigned j = 0; j < this->spaces_size; j++)
      {
        fns.push_back(pss[j]);
        pss[j]->set_quad_2d(&g_quad_2d_std);
      }
      for (unsigned j = 0; j < this->wf->ext.size(); j++)
      {
        fns.push_back(this->wf->ext[j].get());
        this->wf->ext[j]->set_quad_2d(&g_quad_2d_std);
      }

      for(unsigned int form_i = 0; form_i < this->wf->get_forms().size(); form_i++)
      {
        Form<Scalar>* form = this->wf->get_forms()[form_i];
        for(unsigned int ext_i = 0; ext_i < form->ext.size(); ext_i++)
        {
          if(form->ext[ext_i])
          {
            fns.push_back(form->ext[ext_i].get());
            form->ext[ext_i]->set_quad_2d(&g_quad_2d_std);
          }
        }
      }
      if(this->nonlinear)
      {
        init_u_ext(spaces, u_ext_sln);
        for (unsigned j = 0; j < this->wf->get_neq(); j++)
        {
          fns.push_back(u_ext[j]);
          u_ext[j]->set_quad_2d(&g_quad_2d_std);
        }
      }

      this->wf->processFormMarkers(spaces);
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_assembling_one_state(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Traverse::State* current_state_)
    {
      current_state = current_state_;

      for(int j = 0; j < fns.size(); j++)
      {
        if(current_state->e[j])
        {
          fns[j]->set_active_element(current_state->e[j]);
          fns[j]->set_transform(current_state->sub_idx[j]);
        }
      }

      for(int j = 0; j < this->spaces_size; j++)
      {
        if(current_state->e[j])
        {
          spaces[j]->get_element_assembly_list(current_state->e[j], als[j]);
          refmaps[j]->set_active_element(current_state->e[j]);
          refmaps[j]->force_transform(pss[j]->get_transform(), pss[j]->get_ctm());
        }
      }

      if(current_state->isBnd && !(this->wf->mfsurf.empty() && this->wf->vfsurf.empty()))
      {
        for(int j = 0; j < this->spaces_size; j++)
        {
          if(current_state->e[j])
          {
            for(int k = 0; k < current_state->rep->nvert; k++)
            {
              if(current_state->bnd[k])
                spaces[j]->get_boundary_assembly_list(current_state->e[j], k, alsSurface[j][k]);
            }
          }
        }
      }
    }


    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::handle_cache(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, DiscreteProblemCache<Scalar>* cache)
    {
      if(this->do_not_use_cache)
      {
        this->current_cache_record = new typename DiscreteProblemCache<Scalar>::CacheRecord;
        int order = this->integrationOrderCalculator.calculate_order(spaces, current_state, refmaps, u_ext, this->wf);
        this->current_cache_record->init(current_state, pss, refmaps, u_ext, als, alsSurface, this->wf, order);
        return;
      }

      if(this->report_cache_hits_and_misses)
        cache_searches++;

      this->current_cache_record = NULL;
      if(cache->get(current_state->rep, current_state->rep_subidx, current_state->rep_i, this->current_cache_record))
      {
        if(this->report_cache_hits_and_misses)
          cache_record_found++;

        bool reinit = false;
        for(unsigned int i = 0; i < this->spaces_size; i++)
        {
          if(!current_state->e[i])
            continue;

          if(spaces[i]->edata[current_state->e[i]->id].changed_in_last_adaptation)
          {
            reinit = true;
            break;
          }

          if(this->current_cache_record->asmlistCnt[i] != als[i]->cnt)
          {
            reinit = true;
            break;
          }
          else
          {
            for(unsigned int idx_i = 0; idx_i < als[i]->cnt; idx_i++)
            {
              if(als[i]->idx[idx_i] != this->current_cache_record->asmlistIdx[i][idx_i])
              {
                reinit = true;
                break;
              }
            }
            if(reinit)
              break;
          }
        }
        if(reinit)
        {
          if(this->report_cache_hits_and_misses)
            cache_record_found_reinit++;

          this->current_cache_record->free();
          int order = this->integrationOrderCalculator.calculate_order(spaces, current_state, refmaps, u_ext, this->wf);
          this->current_cache_record->init(current_state, pss, refmaps, u_ext, als, alsSurface, this->wf, order);
        }
      }
      else
      {
        if(this->report_cache_hits_and_misses)
          cache_record_not_found++;
        int order = this->integrationOrderCalculator.calculate_order(spaces, current_state, refmaps, u_ext, this->wf);
        this->current_cache_record->init(current_state, pss, refmaps, u_ext, als, alsSurface, this->wf, order);
      }
    }

    template<typename Scalar>
    Func<Scalar>** DiscreteProblemThreadAssembler<Scalar>::init_u_ext_values(int order)
    {
      Func<Scalar>** u_ext_func = NULL;
      if(this->nonlinear)
      {
#ifdef _DEBUG
        assert(this->u_ext);
#endif
        u_ext_func = new Func<Scalar>*[spaces_size];

        for(int i = 0; i < spaces_size; i++)
        {
#ifdef _DEBUG
          assert(u_ext[i]);
#endif
          if(u_ext[i]->get_active_element())
            u_ext_func[i] = init_fn(u_ext[i], order);
          else
            u_ext_func[i] = NULL;
        }
      }

      return u_ext_func;
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_u_ext_values(Func<Scalar>** u_ext_func)
    {
      if(u_ext_func)
      {
        for(int i = 0; i < spaces_size; i++)
        {
          if(u_ext_func[i])
          {
            u_ext_func[i]->free_fn();
            delete u_ext_func[i];
          }
        }
        delete [] u_ext_func;
      }
    }

    template<typename Scalar>
    Func<Scalar>** DiscreteProblemThreadAssembler<Scalar>::init_ext_values(Hermes::vector<MeshFunctionSharedPtr<Scalar> >& ext, int order, Func<Scalar>** u_ext_func)
    {
      Func<Scalar>** ext_func = NULL;
      if(ext.size() > 0)
      {
        ext_func = new Func<Scalar>*[ext.size()];
        for(int ext_i = 0; ext_i < ext.size(); ext_i++)
          if(ext[ext_i])
            if(ext[ext_i]->get_active_element())
            {
              UExtFunction<Scalar>* u_ext_fn = dynamic_cast<UExtFunction<Scalar>*>(ext[ext_i].get());
              if(u_ext_fn)
                ext_func[ext_i] = init_fn(u_ext_fn, u_ext_func, this->spaces_size, order);
              else
                ext_func[ext_i] = init_fn(ext[ext_i].get(), order);
            }
            else
              ext_func[ext_i] = NULL;
          else
            ext_func[ext_i] = NULL;
      }

      if(this->rungeKutta)
      {
        for(int ext_i = 0; ext_i < ext.size(); ext_i++)
          u_ext_func[ext_i]->add(ext_func[ext.size() - this->RK_original_spaces_count + ext_i]);
      }

      return ext_func;
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_ext_values(Hermes::vector<MeshFunctionSharedPtr<Scalar> >& ext, Func<Scalar>** ext_func)
    {
      for(int ext_i = 0; ext_i < ext.size(); ext_i++)
      {
        if(ext[ext_i] && ext[ext_i]->get_active_element())
        {
          ext_func[ext_i]->free_fn();
          delete ext_func[ext_i];
        }
      }
      if(ext_func)
        delete [] ext_func;
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::assemble_one_state()
    {
      // init - u_ext_func
      Func<Scalar>** u_ext_func = this->init_u_ext_values(this->current_cache_record->order);

      // init - ext
      Func<Scalar>** ext_func = this->init_ext_values(this->wf->ext, this->current_cache_record->order, u_ext_func);

      if(this->current_mat || this->add_dirichlet_lift)
      {
        for(int current_mfvol_i = 0; current_mfvol_i < this->wf->mfvol.size(); current_mfvol_i++)
        {
          MatrixFormVol<Scalar>* mfv = this->wf->mfvol[current_mfvol_i];

          if(!selectiveAssembler->form_to_be_assembled(mfv, current_state))
            continue;

          int form_i = mfv->i;
          int form_j = mfv->j;

          this->assemble_matrix_form(mfv,
            this->current_cache_record->order,
            current_cache_record->fns[form_j], current_cache_record->fns[form_i], 
            ext_func, u_ext_func,
            als[form_i], als[form_j],
            current_cache_record->n_quadrature_points, current_cache_record->geometry, current_cache_record->jacobian_x_weights);
        }
      }
      if(this->current_rhs)
      {
        for(int current_vfvol_i = 0; current_vfvol_i < this->wf->vfvol.size(); current_vfvol_i++)
        {
          VectorFormVol<Scalar>* vfv = this->wf->vfvol[current_vfvol_i];

          if(!selectiveAssembler->form_to_be_assembled(vfv, current_state))
            continue;

          int form_i = vfv->i;

          this->assemble_vector_form(vfv, 
            this->current_cache_record->order, 
            current_cache_record->fns[form_i], 
            ext_func, u_ext_func, 
            als[form_i],
            current_cache_record->n_quadrature_points, current_cache_record->geometry, current_cache_record->jacobian_x_weights);
        }
      }

      // deinit - u_ext_func
      this->deinit_u_ext_values(u_ext_func);

      // deinit - ext
      this->deinit_ext_values(this->wf->ext, ext_func);

      // Assemble surface integrals now: loop through surfaces of the element.
      if(current_state->isBnd && (this->wf->mfsurf.size() > 0 || this->wf->vfsurf.size() > 0))
      {
        for (int isurf = 0; isurf < current_state->rep->nvert; isurf++)
        {
          if(!current_state->bnd[isurf])
            continue;

          current_state->isurf = isurf;

          // Edge-wise parameters for WeakForm.
          this->wf->set_active_edge_state(current_state->e, isurf);

          // Ext functions.
          // - order
          int orderSurf = current_cache_record->orderSurface[isurf];

          // init - u_ext_func
          Func<Scalar>** u_ext_funcSurf = this->init_u_ext_values(orderSurf);

          // init - ext
          Func<Scalar>** ext_funcSurf = this->init_ext_values(this->wf->ext, orderSurf, u_ext_funcSurf);

          if(this->current_mat || this->add_dirichlet_lift)
          {
            for(int current_mfsurf_i = 0; current_mfsurf_i < this->wf->mfsurf.size(); current_mfsurf_i++)
            {
              if(!selectiveAssembler->form_to_be_assembled(this->wf->mfsurf[current_mfsurf_i], current_state))
                continue;

              int form_i = this->wf->mfsurf[current_mfsurf_i]->i;
              int form_j = this->wf->mfsurf[current_mfsurf_i]->j;

              this->assemble_matrix_form(this->wf->mfsurf[current_mfsurf_i], 
                current_cache_record->orderSurface[isurf], 
                current_cache_record->fnsSurface[isurf][form_j], 
                current_cache_record->fnsSurface[isurf][form_i], 
                ext_funcSurf, u_ext_funcSurf,
                alsSurface[form_i][isurf], alsSurface[form_j][isurf], 
                current_cache_record->n_quadrature_pointsSurface[isurf], current_cache_record->geometrySurface[isurf], current_cache_record->jacobian_x_weightsSurface[isurf]);
            }
          }

          if(this->current_rhs)
          {
            for(int current_vfsurf_i = 0; current_vfsurf_i < this->wf->vfsurf.size(); current_vfsurf_i++)
            {
              if(!selectiveAssembler->form_to_be_assembled(this->wf->vfsurf[current_vfsurf_i], current_state))
                continue;

              int form_i = this->wf->vfsurf[current_vfsurf_i]->i;

              this->assemble_vector_form(this->wf->vfsurf[current_vfsurf_i], 
                current_cache_record->orderSurface[isurf], 
                current_cache_record->fnsSurface[isurf][form_i], 
                ext_funcSurf, u_ext_funcSurf, 
                alsSurface[form_i][isurf], 
                current_cache_record->n_quadrature_pointsSurface[isurf], current_cache_record->geometrySurface[isurf], current_cache_record->jacobian_x_weightsSurface[isurf]);
            }
          }

          // deinit - u_ext_func
          this->deinit_u_ext_values(u_ext_funcSurf);

          // deinit - ext
          this->deinit_ext_values(this->wf->ext, ext_funcSurf);
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext,
      AsmList<Scalar>* current_als_i, AsmList<Scalar>* current_als_j, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights)
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
        local_ext = this->init_ext_values(form->ext, order, u_ext);

      // Account for the previous time level solution previously inserted at the back of ext.
      if(this->rungeKutta)
        u_ext += form->u_ext_offset;

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_als_i->cnt; i++)
      {
        if(current_als_i->dof[i] < 0 || std::abs(current_als_i->coef[i]) < Hermes::epsilon)
          continue;

        if((!tra || surface_form) && current_als_i->dof[i] < 0)
          continue;

        for (unsigned int j = 0; j < current_als_j->cnt; j++)
        {
          // Skip symmetric values that do not contribute to Dirichlet lift.
          if(sym && j < i && current_als_j->dof[j] >= 0)
            continue;

          // Skip anything that does not contribute to Dirichlet in the case of just rhs assembling.
          if(current_als_j->dof[j] >= 0 && !this->current_mat)
            continue;

          if(std::abs(current_als_j->coef[j]) < Hermes::epsilon)
            continue;

          Func<double>* u = base_fns[j];
          Func<double>* v = test_fns[i];

          Scalar val = block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];

          if(current_als_j->dof[j] >= 0)
          {
            if(surface_form)
              local_stiffness_matrix[i][j] = 0.5 * val;
            else
              local_stiffness_matrix[i][j] = val;

            if(sym)
              local_stiffness_matrix[j][i] = local_stiffness_matrix[i][j];
          }
          else if(this->add_dirichlet_lift && this->current_rhs)
          {
            this->current_rhs->add(current_als_i->dof[i], -val);
          }
        }
      }

      // Insert the local stiffness matrix into the global one.
      if(this->current_mat)
        this->current_mat->add(current_als_i->cnt, current_als_j->cnt, local_stiffness_matrix, current_als_i->dof, current_als_j->dof);

      // Insert also the off-diagonal (anti-)symmetric block, if required.
      if(tra)
      {
        if(form->sym < 0)
          chsgn(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt);
        transpose(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt);

        if(this->current_mat)
          this->current_mat->add(current_als_j->cnt, current_als_i->cnt, local_stiffness_matrix, current_als_j->dof, current_als_i->dof);

        if(this->add_dirichlet_lift && this->current_rhs)
        {
          for (unsigned int j = 0; j < current_als_i->cnt; j++)
            if(current_als_i->dof[j] < 0)
              for (unsigned int i = 0; i < current_als_j->cnt; i++)
                if(current_als_j->dof[i] >= 0)
                  this->current_rhs->add(current_als_j->dof[i], -local_stiffness_matrix[i][j]);
        }
      }

      if(form->ext.size() > 0)
        this->deinit_ext_values(form->ext, local_ext);

      if(this->rungeKutta)
        u_ext -= form->u_ext_offset;

      // Cleanup.
      delete [] local_stiffness_matrix;
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext, 
      AsmList<Scalar>* current_als_i, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights)
    {
      bool surface_form = (dynamic_cast<VectorFormVol<Scalar>*>(form) == NULL);

      Func<Scalar>** local_ext = ext;
      if(form->ext.size() > 0)
        local_ext = this->init_ext_values(form->ext, order, u_ext);

      // Account for the previous time level solution previously inserted at the back of ext.
      if(this->rungeKutta)
        u_ext += form->u_ext_offset;

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_als_i->cnt; i++)
      {
        if(current_als_i->dof[i] < 0)
          continue;

        // Is this necessary, i.e. is there a coefficient smaller than Hermes::epsilon?
        if(std::abs(current_als_i->coef[i]) < Hermes::epsilon)
          continue;

        Func<double>* v = test_fns[i];

        Scalar val;
        if(surface_form)
          val = 0.5 * form->value(n_quadrature_points, jacobian_x_weights, u_ext, v, geometry, local_ext) * form->scaling_factor * current_als_i->coef[i];
        else
          val = form->value(n_quadrature_points, jacobian_x_weights, u_ext, v, geometry, local_ext) * form->scaling_factor * current_als_i->coef[i];

        this->current_rhs->add(current_als_i->dof[i], val);
      }

      if(form->ext.size() > 0)
        this->deinit_ext_values(form->ext, local_ext);

      if(this->rungeKutta)
        u_ext -= form->u_ext_offset;
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_assembling_one_state()
    {
      if(this->do_not_use_cache)
        delete this->current_cache_record;
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_assembling()
    {
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::free()
    {
      this->free_spaces();
      this->free_weak_formulation();
      this->free_u_ext();
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::free_spaces()
    {
      if(!this->pss)
        return;

      for (unsigned int j = 0; j < spaces_size; j++)
        delete pss[j];
      delete [] pss;
      pss = NULL;

      for (unsigned int j = 0; j < spaces_size; j++)
        delete refmaps[j];
      delete [] refmaps;

      for (unsigned int j = 0; j < spaces_size; j++)
        delete als[j];
      delete [] als;

      for (unsigned int j = 0; j < spaces_size; j++)
      {
        for (unsigned int k = 0; k < H2D_MAX_NUMBER_EDGES; k++)
          delete alsSurface[j][k];
        delete [] alsSurface[j];
      }
      delete [] alsSurface;
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::free_weak_formulation()
    {
      if(this->wf)
      {
        this->wf->free_ext();
        delete this->wf;
        this->wf = NULL;
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::free_u_ext()
    {
      if(u_ext)
      {
        for (unsigned int j = 0; j < spaces_size; j++)
          delete u_ext[j];
        delete [] u_ext;
      }
    }

    template class HERMES_API DiscreteProblemThreadAssembler<double>;
    template class HERMES_API DiscreteProblemThreadAssembler<std::complex<double> >;
  }
}