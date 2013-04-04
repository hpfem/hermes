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
#include "shapeset/precalc.h"
#include "function/solution.h"
#include "weakform/weakform.h"
#include "function/exact_solution.h"

namespace Hermes
{
  namespace Hermes2D
  {

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
    void DiscreteProblemThreadAssembler<Scalar>::init_assembling()
    {
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
      for (unsigned j = 0; j < this->wf->get_neq(); j++)
      {
        fns.push_back(u_ext[j]);
        u_ext[j]->set_quad_2d(&g_quad_2d_std);
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_assembling_one_state(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Traverse::State* current_state)
    {
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
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_spaces(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces)
    {
      if(this->spaces_size == spaces.size() && this->pss)
        return;

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
      //this->wf->processFormMarkers(this->spaces);
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

      if(u_ext)
      {
        for (unsigned int j = 0; j < spaces_size; j++)
          delete u_ext[j];
        delete [] u_ext;
      }

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