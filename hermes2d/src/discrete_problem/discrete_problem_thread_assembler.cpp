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
#include "discrete_problem/discrete_problem_form_assembler.h"
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
    DiscreteProblemThreadAssembler<Scalar>::DiscreteProblemThreadAssembler() : pss(NULL), refmaps(NULL), u_ext(NULL), als(NULL), alsSurface(NULL)
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
    void DiscreteProblemThreadAssembler<Scalar>::init_assembling(Solution<Scalar>** u_ext_sln, const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, bool nonlinear_)
    {
#ifdef _DEBUG
    cache_searches = 0;
    cache_record_found = 0;
    cache_record_found_reinit = 0;
    cache_record_not_found = 0;
#endif

      this->nonlinear = nonlinear_;

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
    void DiscreteProblemThreadAssembler<Scalar>::handle_cache(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, DiscreteProblemCache<Scalar>* cache, bool do_not_use_cache_)
    {

      this->do_not_use_cache = do_not_use_cache_;

      if(this->do_not_use_cache)
      {
        this->current_cache_record = new typename DiscreteProblemCache<Scalar>::CacheRecord;
        int order = this->integrationOrderCalculator.calculate_order(spaces, current_state, refmaps, u_ext, this->wf);
        this->current_cache_record->init(current_state, pss, refmaps, u_ext, als, alsSurface, this->wf, order);
        return;
      }
      
#ifdef _DEBUG
    cache_searches++;
#endif

      this->current_cache_record = NULL;
      if(cache->get(current_state->rep, current_state->rep_subidx, current_state->rep_i, this->current_cache_record))
      {
#ifdef _DEBUG
    cache_record_found++;
#endif

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
#ifdef _DEBUG
    cache_record_found_reinit++;
#endif
          this->current_cache_record->free();
          int order = this->integrationOrderCalculator.calculate_order(spaces, current_state, refmaps, u_ext, this->wf);
          this->current_cache_record->init(current_state, pss, refmaps, u_ext, als, alsSurface, this->wf, order);
        }
      }
      else
      {
#ifdef _DEBUG
    cache_record_not_found++;
#endif
        int order = this->integrationOrderCalculator.calculate_order(spaces, current_state, refmaps, u_ext, wf);
        this->current_cache_record->init(current_state, pss, refmaps, u_ext, als, alsSurface, wf, order);
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::assemble_one_state()
    {
      /// \todo
      DiscreteProblemFormAssembler<Scalar> formAssembler;
      formAssembler.set_matrix(this->current_mat);
      formAssembler.set_rhs(this->current_rhs);

      // - u_ext_func
      Func<Scalar>** u_ext_func = NULL;
      int prevNewtonSize = this->wf->get_neq();

      if(this->nonlinear)
      {
        u_ext_func = new Func<Scalar>*[prevNewtonSize];
        if(u_ext)
          for(int u_ext_func_i = 0; u_ext_func_i < prevNewtonSize; u_ext_func_i++)
            if(u_ext[u_ext_func_i])
              if(u_ext[u_ext_func_i]->get_active_element())
                u_ext_func[u_ext_func_i] = init_fn(u_ext[u_ext_func_i], this->current_cache_record->order);
              else
                u_ext_func[u_ext_func_i] = NULL;
            else
              u_ext_func[u_ext_func_i] = NULL;
        else
          for(int u_ext_func_i = 0; u_ext_func_i < prevNewtonSize; u_ext_func_i++)
            u_ext_func[u_ext_func_i] = NULL;
      }

      // - ext
      int current_extCount = this->wf->ext.size();
      Func<Scalar>** ext = NULL;
      if(current_extCount > 0)
      {
        ext = new Func<Scalar>*[current_extCount];
        for(int ext_i = 0; ext_i < current_extCount; ext_i++)
          if(this->wf->ext[ext_i])
            if(this->wf->ext[ext_i]->get_active_element())
              ext[ext_i] = init_fn(this->wf->ext[ext_i].get(), this->current_cache_record->order);
            else
              ext[ext_i] = NULL;
          else
            ext[ext_i] = NULL;
      }

      if(rungeKutta)
        for(int ext_i = 0; ext_i < this->RK_original_spaces_count; ext_i++)
          u_ext_func[ext_i]->add(ext[current_extCount - this->RK_original_spaces_count + ext_i]);

      if(current_mat)
      {
        for(int current_mfvol_i = 0; current_mfvol_i < wf->mfvol.size(); current_mfvol_i++)
        {
          MatrixFormVol<Scalar>* mfv = this->wf->mfvol[current_mfvol_i];

          if(!formAssembler.form_to_be_assembled(mfv, current_state))
            continue;

          int form_i = mfv->i;
          int form_j = mfv->j;

          formAssembler.assemble_matrix_form(mfv, 
            this->current_cache_record->order, 
            current_cache_record->fns[form_j], 
            current_cache_record->fns[form_i], 
            ext,
            u_ext_func,
            als[form_i], 
            als[form_j], 
            current_state, 
            current_cache_record->n_quadrature_points, 
            current_cache_record->geometry, 
            current_cache_record->jacobian_x_weights);
        }
      }
      if(current_rhs)
      {
        for(int current_vfvol_i = 0; current_vfvol_i < wf->vfvol.size(); current_vfvol_i++)
        {
          VectorFormVol<Scalar>* vfv = this->wf->vfvol[current_vfvol_i];

          if(!formAssembler.form_to_be_assembled(vfv, current_state))
            continue;

          int form_i = vfv->i;

          formAssembler.assemble_vector_form(vfv, 
            this->current_cache_record->order, 
            current_cache_record->fns[form_i], 
            ext,
            u_ext_func, 
            als[form_i], 
            current_state, 
            current_cache_record->n_quadrature_points,
            current_cache_record->geometry, 
            current_cache_record->jacobian_x_weights);
        }
      }

      // Cleanup - u_ext_func
      if(u_ext)
      {
        for(int u_ext_func_i = 0; u_ext_func_i < prevNewtonSize; u_ext_func_i++)
          if(u_ext[u_ext_func_i] && u_ext[u_ext_func_i]->get_active_element())
          {
            u_ext_func[u_ext_func_i]->free_fn();
            delete u_ext_func[u_ext_func_i];
          }
          delete [] u_ext_func;
      }

      // Cleanup - ext
      for(int ext_i = 0; ext_i < current_extCount; ext_i++)
      {
        if(this->wf->ext[ext_i] && this->wf->ext[ext_i]->get_active_element())
        {
          ext[ext_i]->free_fn();
          delete ext[ext_i];
        }
      }
      delete [] ext;

      // Assemble surface integrals now: loop through surfaces of the element.
      if(current_state->isBnd && (this->wf->mfsurf.size() > 0 || this->wf->vfsurf.size() > 0))
      {
        for (current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
        {
          if(!current_state->bnd[current_state->isurf])
            continue;

          // Edge-wise parameters for WeakForm.
          (const_cast<WeakForm<Scalar>*>(this->wf))->set_active_edge_state(current_state->e, current_state->isurf);

          // Ext functions.
          // - order
          int orderSurf = current_cache_record->orderSurface[current_state->isurf];

          // - u_ext_func
          int prevNewtonSize = this->wf->get_neq();
          Func<Scalar>** u_ext_funcSurf = NULL;
          if(this->nonlinear)
          {
            u_ext_funcSurf = new Func<Scalar>*[prevNewtonSize];
            if(u_ext)
              for(int u_ext_func_surf_i = 0; u_ext_func_surf_i < prevNewtonSize; u_ext_func_surf_i++)
                if(u_ext[u_ext_func_surf_i])
                  u_ext_funcSurf[u_ext_func_surf_i] = current_state->e[u_ext_func_surf_i] == NULL ? NULL : init_fn(u_ext[u_ext_func_surf_i], orderSurf);
                else
                  u_ext_funcSurf[u_ext_func_surf_i] = NULL;
            else
              for(int u_ext_func_surf_i = 0; u_ext_func_surf_i < prevNewtonSize; u_ext_func_surf_i++)
                u_ext_funcSurf[u_ext_func_surf_i] = NULL;
          }
          // - ext
          int current_extCount = this->wf->ext.size();
          Func<Scalar>** extSurf = new Func<Scalar>*[current_extCount];
          for(int ext_surf_i = 0; ext_surf_i < current_extCount; ext_surf_i++)
            if(this->wf->ext[ext_surf_i])
              extSurf[ext_surf_i] = current_state->e[ext_surf_i] == NULL ? NULL : init_fn(this->wf->ext[ext_surf_i].get(), orderSurf);
            else
              extSurf[ext_surf_i] = NULL;

          if(rungeKutta)
            for(int ext_surf_i = 0; ext_surf_i < this->RK_original_spaces_count; ext_surf_i++)
              u_ext_funcSurf[ext_surf_i]->add(extSurf[current_extCount - this->RK_original_spaces_count + ext_surf_i]);

          if(current_mat)
          {
            for(int current_mfsurf_i = 0; current_mfsurf_i < wf->mfsurf.size(); current_mfsurf_i++)
            {
              if(!formAssembler.form_to_be_assembled(this->wf->mfsurf[current_mfsurf_i], current_state))
                continue;

              int form_i = this->wf->mfsurf[current_mfsurf_i]->i;
              int form_j = this->wf->mfsurf[current_mfsurf_i]->j;

              formAssembler.assemble_matrix_form(this->wf->mfsurf[current_mfsurf_i], 
                current_cache_record->orderSurface[current_state->isurf], 
                current_cache_record->fnsSurface[current_state->isurf][form_j], 
                current_cache_record->fnsSurface[current_state->isurf][form_i], 
                extSurf, 
                u_ext_funcSurf,
                alsSurface[form_i][current_state->isurf], 
                alsSurface[form_j][current_state->isurf], 
                current_state, 
                current_cache_record->n_quadrature_pointsSurface[current_state->isurf], 
                current_cache_record->geometrySurface[current_state->isurf], 
                current_cache_record->jacobian_x_weightsSurface[current_state->isurf]);
            }
          }

          if(current_rhs)
          {
            for(int current_vfsurf_i = 0; current_vfsurf_i < wf->vfsurf.size(); current_vfsurf_i++)
            {
              if(!formAssembler.form_to_be_assembled(this->wf->vfsurf[current_vfsurf_i], current_state))
                continue;

              int form_i = this->wf->vfsurf[current_vfsurf_i]->i;

              formAssembler.assemble_vector_form(this->wf->vfsurf[current_vfsurf_i], 
                current_cache_record->orderSurface[current_state->isurf], 
                current_cache_record->fnsSurface[current_state->isurf][form_i], 
                extSurf, 
                u_ext_funcSurf, 
                alsSurface[form_i][current_state->isurf], 
                current_state, 
                current_cache_record->n_quadrature_pointsSurface[current_state->isurf], 
                current_cache_record->geometrySurface[current_state->isurf], 
                current_cache_record->jacobian_x_weightsSurface[current_state->isurf]);
            }
          }

          if(u_ext)
          {
            for(int u_ext_func_surf_i = 0; u_ext_func_surf_i < prevNewtonSize; u_ext_func_surf_i++)
              if(u_ext[u_ext_func_surf_i] && current_state->e[u_ext_func_surf_i])
              {
                u_ext_funcSurf[u_ext_func_surf_i]->free_fn();
                delete u_ext_funcSurf[u_ext_func_surf_i];
              }
              delete [] u_ext_funcSurf;
          }

          for(int ext_surf_i = 0; ext_surf_i < current_extCount; ext_surf_i++)
            if(this->wf->ext[ext_surf_i] && current_state->e[ext_surf_i])
            {
              extSurf[ext_surf_i]->free_fn();
              delete extSurf[ext_surf_i];
            }
            delete [] extSurf;
        }
      }
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
#ifdef _DEBUG
    if(!this->do_not_use_cache)
    {
      std::cout << std::endl;
      std::cout << "cache_searches  " << cache_searches << std::endl;
      std::cout << "cache_record_found  " << cache_record_found << std::endl;
      std::cout << "cache_record_found_reinit  " << cache_record_found_reinit << std::endl;
      std::cout << "cache_record_not_found  " << cache_record_not_found << std::endl;
      std::cout << std::endl;
    }
#endif
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