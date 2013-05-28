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
#include "global.h"
#include "quadrature/limit_order.h"
#include "mesh/traverse.h"
#include "space/space.h"
#include "shapeset/precalc.h"
#include "mesh/refmap.h"
#include "function/solution.h"
#include "neighbor_search.h"
#include "api2d.h"
#include <algorithm>

using namespace Hermes::Algebra::DenseMatrixOperations;

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    DiscreteProblemIntegrationOrderCalculator<Scalar>::DiscreteProblemIntegrationOrderCalculator(DiscreteProblemSelectiveAssembler<Scalar>* selectiveAssembler) : selectiveAssembler(selectiveAssembler)
    {
    }

    template<typename Scalar>
    int DiscreteProblemIntegrationOrderCalculator<Scalar>::calculate_order(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Traverse::State* current_state, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, WeakForm<Scalar>* current_wf)
    {
      // Order calculation.
      int order = current_wf->global_integration_order_set ? current_wf->global_integration_order : 0;
      if(order == 0)
      {
        for(int current_mfvol_i = 0; current_mfvol_i < current_wf->mfvol.size(); current_mfvol_i++)
        {
          MatrixFormVol<Scalar>* current_mfvol = current_wf->mfvol[current_mfvol_i];
          if(!selectiveAssembler->form_to_be_assembled(current_mfvol, current_state))
            continue;
          current_mfvol->wf = current_wf;
          int orderTemp = calc_order_matrix_form(spaces, current_mfvol, current_refmaps, current_u_ext, current_state);
          if(order < orderTemp)
            order = orderTemp;
        }

        for(int current_vfvol_i = 0; current_vfvol_i < current_wf->vfvol.size(); current_vfvol_i++)
        {
          VectorFormVol<Scalar>* current_vfvol = current_wf->vfvol[current_vfvol_i];
          if(!selectiveAssembler->form_to_be_assembled(current_vfvol, current_state))
            continue;
          current_vfvol->wf = current_wf;
          int orderTemp = calc_order_vector_form(spaces, current_vfvol, current_refmaps, current_u_ext, current_state);
          if(order < orderTemp)
            order = orderTemp;
        }

        // Surface forms.
        if(current_state->isBnd && (current_wf->mfsurf.size() > 0 || current_wf->vfsurf.size() > 0))
        {
          for (current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
          {
            if(!current_state->bnd[current_state->isurf])
              continue;
            for(int current_mfsurf_i = 0; current_mfsurf_i < current_wf->mfsurf.size(); current_mfsurf_i++)
            {
              MatrixFormSurf<Scalar>* current_mfsurf = current_wf->mfsurf[current_mfsurf_i];
              if(!selectiveAssembler->form_to_be_assembled(current_mfsurf, current_state))
                continue;
              current_mfsurf->wf = current_wf;
              int orderTemp = calc_order_matrix_form(spaces, current_mfsurf, current_refmaps, current_u_ext, current_state);
              if(order < orderTemp)
                order = orderTemp;
            }

            for(int current_vfsurf_i = 0; current_vfsurf_i < current_wf->vfsurf.size(); current_vfsurf_i++)
            {
              VectorFormSurf<Scalar>* current_vfsurf = current_wf->vfsurf[current_vfsurf_i];
              if(!selectiveAssembler->form_to_be_assembled(current_vfsurf, current_state))
                continue;

              current_vfsurf->wf = current_wf;
              int orderTemp = calc_order_vector_form(spaces, current_vfsurf, current_refmaps, current_u_ext, current_state);
              if(order < orderTemp)
                order = orderTemp;
            }
          }
        }
      }

      return order;
    } 

    template<typename Scalar>
    int DiscreteProblemIntegrationOrderCalculator<Scalar>::calc_order_matrix_form(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, MatrixForm<Scalar> *form, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, Traverse::State* current_state)
    {
      int order;

      // order of solutions from the previous Newton iteration etc..
      Func<Hermes::Ord>** u_ext_ord = current_u_ext == NULL ? NULL : new Func<Hermes::Ord>*[this->rungeKutta ? this->RK_original_spaces_count : form->wf->get_neq() - form->u_ext_offset];
      Func<Hermes::Ord>** ext_ord = NULL;
      int ext_size = std::max(form->ext.size(), form->wf->ext.size());
      if(ext_size > 0)
        ext_ord = new Func<Hermes::Ord>*[ext_size];
      init_ext_orders(form, u_ext_ord, ext_ord, current_u_ext, current_state);

      // Order of shape functions.
      int max_order_j = spaces[form->j]->get_element_order(current_state->e[form->j]->id);
      int max_order_i = spaces[form->i]->get_element_order(current_state->e[form->i]->id);
      if(H2D_GET_V_ORDER(max_order_i) > H2D_GET_H_ORDER(max_order_i))
        max_order_i = H2D_GET_V_ORDER(max_order_i);
      else
        max_order_i = H2D_GET_H_ORDER(max_order_i);
      if(H2D_GET_V_ORDER(max_order_j) > H2D_GET_H_ORDER(max_order_j))
        max_order_j = H2D_GET_V_ORDER(max_order_j);
      else
        max_order_j = H2D_GET_H_ORDER(max_order_j);

      for (unsigned int k = 0; k < current_state->rep->nvert; k++)
      {
        int eo = spaces[form->i]->get_edge_order(current_state->e[form->i], k);
        if(eo > max_order_i)
          max_order_i = eo;
        eo = spaces[form->j]->get_edge_order(current_state->e[form->j], k);
        if(eo > max_order_j)
          max_order_j = eo;
      }

      Func<Hermes::Ord>* ou = init_fn_ord(max_order_j + (spaces[form->j]->get_shapeset()->get_num_components() > 1 ? 1 : 0));
      Func<Hermes::Ord>* ov = init_fn_ord(max_order_i + (spaces[form->i]->get_shapeset()->get_num_components() > 1 ? 1 : 0));

      // Total order of the vector form.
      double fake_wt = 1.0;
      Geom<Hermes::Ord> *tmp = init_geom_ord();
      Hermes::Ord o = form->ord(1, &fake_wt, u_ext_ord, ou, ov, tmp, ext_ord);
      delete tmp;

      adjust_order_to_refmaps(form, order, &o, current_refmaps);

      // Cleanup.
      deinit_ext_orders(form, u_ext_ord, ext_ord);
      ou->free_ord();
      delete ou;
      ov->free_ord();
      delete ov;

      return order;
    }

    template<typename Scalar>
    int DiscreteProblemIntegrationOrderCalculator<Scalar>::calc_order_vector_form(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, VectorForm<Scalar> *form, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, Traverse::State* current_state)
    {
      int order;

      // order of solutions from the previous Newton iteration etc..
      Func<Hermes::Ord>** u_ext_ord = current_u_ext == NULL ? NULL : new Func<Hermes::Ord>*[this->rungeKutta ? this->RK_original_spaces_count : form->wf->get_neq() - form->u_ext_offset];
      Func<Hermes::Ord>** ext_ord = NULL;
      int ext_size = std::max(form->ext.size(), form->wf->ext.size());
      if(ext_size > 0)
        ext_ord = new Func<Hermes::Ord>*[ext_size];
      init_ext_orders(form, u_ext_ord, ext_ord, current_u_ext, current_state);

      // Order of shape functions.
      int max_order_i = spaces[form->i]->get_element_order(current_state->e[form->i]->id);
      if(H2D_GET_V_ORDER(max_order_i) > H2D_GET_H_ORDER(max_order_i))
        max_order_i = H2D_GET_V_ORDER(max_order_i);
      else
        max_order_i = H2D_GET_H_ORDER(max_order_i);

      for (unsigned int k = 0; k < current_state->rep->nvert; k++)
      {
        int eo = spaces[form->i]->get_edge_order(current_state->e[form->i], k);
        if(eo > max_order_i)
          max_order_i = eo;
      }
      Func<Hermes::Ord>* ov = init_fn_ord(max_order_i + (spaces[form->i]->get_shapeset()->get_num_components() > 1 ? 1 : 0));

      // Total order of the vector form.
      double fake_wt = 1.0;
      Geom<Hermes::Ord> *tmp = init_geom_ord();
      Hermes::Ord o = form->ord(1, &fake_wt, u_ext_ord, ov, tmp, ext_ord);
      delete tmp;

      adjust_order_to_refmaps(form, order, &o, current_refmaps);

      // Cleanup.
      deinit_ext_orders(form, u_ext_ord, ext_ord);
      ov->free_ord();
      delete ov;

      return order;
    }

    template<typename Scalar>
    void DiscreteProblemIntegrationOrderCalculator<Scalar>::init_ext_orders(Form<Scalar> *form, Func<Hermes::Ord>** oi, Func<Hermes::Ord>** oext, Solution<Scalar>** current_u_ext, Traverse::State* current_state)
    {
      unsigned int prev_size = this->rungeKutta ? this->RK_original_spaces_count : form->wf->get_neq() - form->u_ext_offset;
      bool surface_form = (current_state->isurf > -1);

      if(current_u_ext)
        for(int i = 0; i < prev_size; i++)
          if(current_u_ext[i + form->u_ext_offset])
            if(surface_form)
              oi[i] = init_fn_ord(current_u_ext[i + form->u_ext_offset]->get_edge_fn_order(current_state->isurf) + (current_u_ext[i + form->u_ext_offset]->get_num_components() > 1 ? 1 : 0));
            else
              oi[i] = init_fn_ord(current_u_ext[i + form->u_ext_offset]->get_fn_order() + (current_u_ext[i + form->u_ext_offset]->get_num_components() > 1 ? 1 : 0));
          else
            oi[i] = init_fn_ord(0);

      if(form->ext.size() > 0)
      {
        for (int i = 0; i < form->ext.size(); i++)
          if(surface_form)
            oext[i] = init_fn_ord(form->ext[i]->get_edge_fn_order(current_state->isurf) + (form->ext[i]->get_num_components() > 1 ? 1 : 0));
          else
            oext[i] = init_fn_ord(form->ext[i]->get_fn_order() + (form->ext[i]->get_num_components() > 1 ? 1 : 0));
      }

      else
      {
        for (int i = 0; i < form->wf->ext.size(); i++)
          if(surface_form)
            oext[i] = init_fn_ord(form->wf->ext[i]->get_edge_fn_order(current_state->isurf) + (form->wf->ext[i]->get_num_components() > 1 ? 1 : 0));
          else
            oext[i] = init_fn_ord(form->wf->ext[i]->get_fn_order() + (form->wf->ext[i]->get_num_components() > 1 ? 1 : 0));
      }
    }

    template<typename Scalar>
    void DiscreteProblemIntegrationOrderCalculator<Scalar>::deinit_ext_orders(Form<Scalar> *form, Func<Hermes::Ord>** oi, Func<Hermes::Ord>** oext)
    {
      unsigned int prev_size = oi ? (this->rungeKutta ? this->RK_original_spaces_count : form->wf->get_neq() - form->u_ext_offset) : 0;
      for(int i = 0; i < prev_size; i++)
      {
        oi[i]->free_ord();
        delete oi[i];
      }
      if(oi)
        delete [] oi;

      if(oext)
      {
        if(form->ext.size() > 0)
          for (int i = 0; i < form->ext.size(); i++)
          {
            oext[i]->free_ord();
            delete oext[i];
          }
        else
          for (int i = 0; i < form->wf->ext.size(); i++)
          {
            oext[i]->free_ord();
            delete oext[i];
          }

          delete [] oext;
      }
    }

    template<typename Scalar>
    void DiscreteProblemIntegrationOrderCalculator<Scalar>::adjust_order_to_refmaps(Form<Scalar> *form, int& order, Hermes::Ord* o, RefMap** current_refmaps)
    {
      // Increase due to reference map.
      int coordinate = form->i;
      order = current_refmaps[coordinate]->get_inv_ref_order();
      order += o->get_order();
      limit_order(order, current_refmaps[coordinate]->get_active_element()->get_mode());
    }

    template<typename Scalar>
    int DiscreteProblemIntegrationOrderCalculator<Scalar>::calc_order_dg_matrix_form(MatrixFormDG<Scalar>* mfDG, PrecalcShapeset* fu, PrecalcShapeset* fv, RefMap* ru, bool neighbor_supp_u, bool neighbor_supp_v, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_u)
    {
      throw Hermes::Exceptions::Exception("DiscreteProblemIntegrationOrderCalculator<Scalar>::calc_order_dg_matrix_form not implemented.");
      return -1;
    }

    template<typename Scalar>
    int DiscreteProblemIntegrationOrderCalculator<Scalar>::calc_order_dg_vector_form(VectorFormDG<Scalar>* vfDG, Hermes::vector<Solution<Scalar> > u_ext,
      PrecalcShapeset* fv, RefMap* ru, LightArray<NeighborSearch<Scalar>*>& neighbor_searches, int neighbor_index_v)
    {
      throw Hermes::Exceptions::Exception("DiscreteProblemIntegrationOrderCalculator<Scalar>::calc_order_dg_vector_form not implemented.");
      return -1;
    }

    template class HERMES_API DiscreteProblemIntegrationOrderCalculator<double>;
    template class HERMES_API DiscreteProblemIntegrationOrderCalculator<std::complex<double> >;
  }
}