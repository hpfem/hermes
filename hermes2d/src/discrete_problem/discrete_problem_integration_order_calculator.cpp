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
#include "discrete_problem/dg/discrete_problem_dg_assembler.h"
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
    /// Geometry instance for order calculation.
    static GeomVol<Hermes::Ord> geom_order_vol;
    static GeomSurf<Hermes::Ord> geom_order_surf;
    /// "Fake" integration weight for order calculation.
    double wt_order = 1.0;

    Func<Hermes::Ord> func_order[g_max_quad + 1] =
    {
      Func<Hermes::Ord>(0),
      Func<Hermes::Ord>(1),
      Func<Hermes::Ord>(2),
      Func<Hermes::Ord>(3),
      Func<Hermes::Ord>(4),
      Func<Hermes::Ord>(5),
      Func<Hermes::Ord>(6),
      Func<Hermes::Ord>(7),
      Func<Hermes::Ord>(8),
      Func<Hermes::Ord>(9),
      Func<Hermes::Ord>(10),
      Func<Hermes::Ord>(11),
      Func<Hermes::Ord>(12),
      Func<Hermes::Ord>(13),
      Func<Hermes::Ord>(14),
      Func<Hermes::Ord>(15),
      Func<Hermes::Ord>(16),
      Func<Hermes::Ord>(17),
      Func<Hermes::Ord>(18),
      Func<Hermes::Ord>(19),
      Func<Hermes::Ord>(20),
      Func<Hermes::Ord>(21),
      Func<Hermes::Ord>(22),
      Func<Hermes::Ord>(23),
      Func<Hermes::Ord>(24)
    };

    template<typename Scalar>
    DiscreteProblemIntegrationOrderCalculator<Scalar>::DiscreteProblemIntegrationOrderCalculator(DiscreteProblemSelectiveAssembler<Scalar>* selectiveAssembler) :
      selectiveAssembler(selectiveAssembler),
      current_state(nullptr),
      u_ext(nullptr)
    {
    }

    template<typename Scalar>
    int DiscreteProblemIntegrationOrderCalculator<Scalar>::calculate_order(const SpaceSharedPtrVector<Scalar>& spaces, RefMap** current_refmaps, WeakFormSharedPtr<Scalar> current_wf)
    {
      // Order set to constant.
      if (current_wf->global_integration_order_set)
        return current_wf->global_integration_order;

      // Order calculation.
      int order = 0;

      // init - u_ext_func
      Func<Hermes::Ord>** u_ext_func = this->init_u_ext_orders();
      // init - ext
      Func<Hermes::Ord>** ext_func = this->init_ext_orders(current_wf->ext, current_wf->u_ext_fn, u_ext_func);

      for(unsigned short current_mfvol_i = 0; current_mfvol_i < current_wf->mfvol.size(); current_mfvol_i++)
      {
        MatrixFormVol<Scalar>* current_mfvol = current_wf->mfvol[current_mfvol_i];
        if (!selectiveAssembler->form_to_be_assembled(current_mfvol, current_state))
          continue;
        current_mfvol->wf = current_wf.get();
        int orderTemp = calc_order_matrix_form(spaces, current_mfvol, current_refmaps, ext_func, u_ext_func);
        if (order < orderTemp)
          order = orderTemp;
      }

      for(unsigned short current_vfvol_i = 0; current_vfvol_i < current_wf->vfvol.size(); current_vfvol_i++)
      {
        VectorFormVol<Scalar>* current_vfvol = current_wf->vfvol[current_vfvol_i];
        if (!selectiveAssembler->form_to_be_assembled(current_vfvol, current_state))
          continue;
        current_vfvol->wf = current_wf.get();
        int orderTemp = calc_order_vector_form(spaces, current_vfvol, current_refmaps, ext_func, u_ext_func);
        if (order < orderTemp)
          order = orderTemp;
      }

      // Surface forms.
      if (current_state->isBnd && (current_wf->mfsurf.size() > 0 || current_wf->vfsurf.size() > 0))
      {
        for (current_state->isurf = 0; current_state->isurf < current_state->rep->nvert; current_state->isurf++)
        {
          if (!current_state->bnd[current_state->isurf])
            continue;

          // init - u_ext_func
          Func<Hermes::Ord>** u_ext_funcSurf = this->init_u_ext_orders();

          // init - ext
          Func<Hermes::Ord>** ext_funcSurf = this->init_ext_orders(current_wf->ext, current_wf->u_ext_fn, u_ext_func);

          for(unsigned short current_mfsurf_i = 0; current_mfsurf_i < current_wf->mfsurf.size(); current_mfsurf_i++)
          {
            MatrixFormSurf<Scalar>* current_mfsurf = current_wf->mfsurf[current_mfsurf_i];
            if (!selectiveAssembler->form_to_be_assembled(current_mfsurf, current_state))
              continue;
            current_mfsurf->wf = current_wf.get();
            int orderTemp = calc_order_matrix_form(spaces, current_mfsurf, current_refmaps, ext_funcSurf, u_ext_funcSurf);
            if (order < orderTemp)
              order = orderTemp;
          }

          for(unsigned short current_vfsurf_i = 0; current_vfsurf_i < current_wf->vfsurf.size(); current_vfsurf_i++)
          {
            VectorFormSurf<Scalar>* current_vfsurf = current_wf->vfsurf[current_vfsurf_i];
            if (!selectiveAssembler->form_to_be_assembled(current_vfsurf, current_state))
              continue;

            current_vfsurf->wf = current_wf.get();
            int orderTemp = calc_order_vector_form(spaces, current_vfsurf, current_refmaps, ext_funcSurf, u_ext_funcSurf);
            if (order < orderTemp)
              order = orderTemp;
          }

          // deinit - u_ext_func
          this->deinit_u_ext_orders(u_ext_funcSurf);

          // deinit - ext
          this->deinit_ext_orders(ext_funcSurf);
        }
      }

      // deinit - u_ext_func
      this->deinit_u_ext_orders(u_ext_func);

      // deinit - ext
      this->deinit_ext_orders(ext_func);

      return order;
    }

    template<typename Scalar>
    template<typename MatrixFormType>
    int DiscreteProblemIntegrationOrderCalculator<Scalar>::calc_order_matrix_form(const SpaceSharedPtrVector<Scalar>& spaces, MatrixFormType *form, RefMap** current_refmaps, Func<Hermes::Ord>** ext, Func<Hermes::Ord>** u_ext)
    {
      int order;

      Func<Hermes::Ord>** local_ext = ext;
      // If the user supplied custom ext functions for this form.
      if (form->ext.size() > 0)
        local_ext = this->init_ext_orders(form->ext, (form->u_ext_fn.size() > 0 ? form->u_ext_fn : form->wf->u_ext_fn), u_ext);

      // Order of shape functions.
      int max_order_j = spaces[form->j]->get_element_order(current_state->e[form->j]->id);
      int max_order_i = spaces[form->i]->get_element_order(current_state->e[form->i]->id);
      if (H2D_GET_V_ORDER(max_order_i) > H2D_GET_H_ORDER(max_order_i))
        max_order_i = H2D_GET_V_ORDER(max_order_i);
      else
        max_order_i = H2D_GET_H_ORDER(max_order_i);
      if (H2D_GET_V_ORDER(max_order_j) > H2D_GET_H_ORDER(max_order_j))
        max_order_j = H2D_GET_V_ORDER(max_order_j);
      else
        max_order_j = H2D_GET_H_ORDER(max_order_j);

      for (unsigned int k = 0; k < current_state->rep->nvert; k++)
      {
        int eo = spaces[form->i]->get_edge_order(current_state->e[form->i], k);
        if (eo > max_order_i)
          max_order_i = eo;
        eo = spaces[form->j]->get_edge_order(current_state->e[form->j], k);
        if (eo > max_order_j)
          max_order_j = eo;
      }

      Func<Hermes::Ord>* ou = &func_order[max_order_j + (spaces[form->j]->shapeset->num_components > 1 ? 1 : 0)];
      Func<Hermes::Ord>* ov = &func_order[max_order_i + (spaces[form->i]->shapeset->num_components > 1 ? 1 : 0)];

      // Total order of the vector form.
      Hermes::Ord o;
      if(dynamic_cast<MatrixFormVol<Scalar>*>(form))
        o = (dynamic_cast<MatrixFormVol<Scalar>*>(form))->ord(1, &wt_order, u_ext, ou, ov, &geom_order_vol, local_ext);
      else
        o = (dynamic_cast<MatrixFormSurf<Scalar>*>(form))->ord(1, &wt_order, u_ext, ou, ov, &geom_order_surf, local_ext);

      adjust_order_to_refmaps(form, order, &o, current_refmaps);

      // Cleanup.
      if (form->ext.size() > 0)
        this->deinit_ext_orders(local_ext);

      return order;
    }

    template<typename Scalar>
    template<typename VectorFormType>
    int DiscreteProblemIntegrationOrderCalculator<Scalar>::calc_order_vector_form(const SpaceSharedPtrVector<Scalar>& spaces, VectorFormType *form, RefMap** current_refmaps, Func<Hermes::Ord>** ext, Func<Hermes::Ord>** u_ext)
    {
      int order;

      Func<Hermes::Ord>** local_ext = ext;
      // If the user supplied custom ext functions for this form.
      if (form->ext.size() > 0)
        local_ext = this->init_ext_orders(form->ext, (form->u_ext_fn.size() > 0 ? form->u_ext_fn : form->wf->u_ext_fn), u_ext);

      // Order of shape functions.
      int max_order_i = spaces[form->i]->get_element_order(current_state->e[form->i]->id);
      if (H2D_GET_V_ORDER(max_order_i) > H2D_GET_H_ORDER(max_order_i))
        max_order_i = H2D_GET_V_ORDER(max_order_i);
      else
        max_order_i = H2D_GET_H_ORDER(max_order_i);

      for (unsigned int k = 0; k < current_state->rep->nvert; k++)
      {
        int eo = spaces[form->i]->get_edge_order(current_state->e[form->i], k);
        if (eo > max_order_i)
          max_order_i = eo;
      }
      Func<Hermes::Ord>* ov = &func_order[max_order_i + (spaces[form->i]->shapeset->num_components > 1 ? 1 : 0)];

      // Total order of the vector form.
      Hermes::Ord o;
      if (dynamic_cast<VectorFormVol<Scalar>*>(form))
        o = (dynamic_cast<VectorFormVol<Scalar>*>(form))->ord(1, &wt_order, u_ext, ov, &geom_order_vol, local_ext);
      else
        o = (dynamic_cast<VectorFormSurf<Scalar>*>(form))->ord(1, &wt_order, u_ext, ov, &geom_order_surf, local_ext);

      adjust_order_to_refmaps(form, order, &o, current_refmaps);

      // Cleanup.
      if (form->ext.size() > 0)
        this->deinit_ext_orders(local_ext);

      return order;
    }

    template<typename Scalar>
    Func<Hermes::Ord>** DiscreteProblemIntegrationOrderCalculator<Scalar>::init_u_ext_orders()
    {
      Func<Hermes::Ord>** u_ext_func = nullptr;
      bool surface_form = (this->current_state->isurf > -1);
      if (this->u_ext)
      {
        u_ext_func = new Func<Hermes::Ord>*[this->selectiveAssembler->spaces_size];

        for (int i = 0; i < this->selectiveAssembler->spaces_size; i++)
        {
          assert(this->u_ext[i]);

          if (this->u_ext[i]->get_active_element())
          {
            if (surface_form)
              u_ext_func[i] = &func_order[this->u_ext[i]->get_edge_fn_order(this->current_state->isurf) + (this->u_ext[i]->get_num_components() > 1 ? 1 : 0)];
            else
              u_ext_func[i] = &func_order[this->u_ext[i]->get_fn_order() + (this->u_ext[i]->get_num_components() > 1 ? 1 : 0)];
          }
          else
            u_ext_func[i] = &func_order[0];
        }
      }

      return u_ext_func;
    }

    template<typename Scalar>
    void DiscreteProblemIntegrationOrderCalculator<Scalar>::deinit_u_ext_orders(Func<Hermes::Ord>** u_ext_func)
    {
      free_with_check(u_ext_func);
    }


    template<typename Scalar>
    Func<Hermes::Ord>** DiscreteProblemIntegrationOrderCalculator<Scalar>::init_ext_orders(MeshFunctionSharedPtrVector<Scalar>& ext, std::vector<UExtFunctionSharedPtr<Scalar> >& u_ext_fns, Func<Hermes::Ord>** u_ext_func)
    {
      int ext_size = ext.size();
      int u_ext_fns_size = u_ext_fns.size();

      Func<Hermes::Ord>** ext_func = nullptr;

      bool surface_form = (this->current_state->isurf > -1);

      if (ext_size > 0 || u_ext_fns_size > 0)
      {
        ext_func = malloc_with_check<Func<Hermes::Ord>*>(ext_size + u_ext_fns_size);
        for(unsigned short ext_i = 0; ext_i < ext.size(); ext_i++)
        {
          if (ext[ext_i])
          {
            if (ext[ext_i]->get_active_element())
            {
              if (surface_form)
                ext_func[u_ext_fns_size + ext_i] = &func_order[ext[ext_i]->get_edge_fn_order(this->current_state->isurf) + (ext[ext_i]->get_num_components() > 1 ? 1 : 0)];
              else
                ext_func[u_ext_fns_size + ext_i] = &func_order[ext[ext_i]->get_fn_order() + (ext[ext_i]->get_num_components() > 1 ? 1 : 0)];
            }
            else
              ext_func[u_ext_fns_size + ext_i] = nullptr;
          }
          else
            ext_func[u_ext_fns_size + ext_i] = nullptr;
        }

        for (int ext_i = 0; ext_i < u_ext_fns_size; ext_i++)
        {
          if (u_ext_fns[ext_i])
          {
            ext_func[ext_i] = &func_order[0];
            u_ext_fns[ext_i]->ord(ext_func + u_ext_fns_size, u_ext_func, ext_func[ext_i]);
          }
          else
            ext_func[ext_i] = nullptr;
        }
      }

      return ext_func;
    }

    template<typename Scalar>
    void DiscreteProblemIntegrationOrderCalculator<Scalar>::deinit_ext_orders(Func<Hermes::Ord>** ext_func)
    {
      free_with_check(ext_func);
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
    DiscontinuousFunc<Ord>* DiscreteProblemIntegrationOrderCalculator<Scalar>::init_ext_fn_ord(NeighborSearch<Scalar>* ns, MeshFunctionSharedPtr<Scalar> fu)
    {
      int inc = (fu->get_num_components() == 2) ? 1 : 0;
      int central_order = fu->get_edge_fn_order(ns->active_edge) + inc;
      int neighbor_order = fu->get_edge_fn_order(ns->neighbor_edge.local_num_of_edge) + inc;
      return new DiscontinuousFunc<Ord>(&func_order[central_order], &func_order[neighbor_order]);
    }

    template<typename Scalar>
    DiscontinuousFunc<Hermes::Ord>** DiscreteProblemIntegrationOrderCalculator<Scalar>::init_ext_fns_ord(MeshFunctionSharedPtrVector<Scalar> &ext,
      NeighborSearch<Scalar>** neighbor_searches)
    {
      DiscontinuousFunc<Ord>** fake_ext_fns = new DiscontinuousFunc<Ord>*[ext.size()];
      for (unsigned int j = 0; j < ext.size(); j++)
        fake_ext_fns[j] = init_ext_fn_ord(DiscreteProblemDGAssembler<Scalar>::get_neighbor_search_ext(this->selectiveAssembler->get_weak_formulation(), neighbor_searches, j), ext[j]);

      return fake_ext_fns;
    }

    template<typename Scalar>
    template<typename FormType>
    void DiscreteProblemIntegrationOrderCalculator<Scalar>::deinit_ext_fns_ord(Form<Scalar> *form, FormType** oi, FormType** oext)
    {
      unsigned int prev_size = oi ? (this->rungeKutta ? this->RK_original_spaces_count : form->wf->get_neq() - form->u_ext_offset) : 0;

      if (oi)
      {
        for (int i = 0; i < prev_size; i++)
          delete oi[i];
        delete[] oi;
      }

      int ext_size = form->ext.size() ? form->ext.size() : form->wf->ext.size();
      int u_ext_fn_size = form->u_ext_fn.size() ? form->u_ext_fn.size() : form->wf->u_ext_fn.size();
      if (oext)
      {
        for (int i = 0; i < ext_size + u_ext_fn_size; i++)
          delete oext[i];
        delete[] oext;
      }
    }

    template<typename Scalar>
    int DiscreteProblemIntegrationOrderCalculator<Scalar>::calc_order_dg_matrix_form(const SpaceSharedPtrVector<Scalar> spaces, Traverse::State* current_state, MatrixFormDG<Scalar>* mfDG, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, bool neighbor_supp_u, bool neighbor_supp_v, NeighborSearch<Scalar>** neighbor_searches)
    {
      NeighborSearch<Scalar>* nbs_u = neighbor_searches[mfDG->j];

      unsigned short prev_size = this->rungeKutta ? this->RK_original_spaces_count : mfDG->wf->get_neq() - mfDG->u_ext_offset;

      // Order to return.
      int order = 0;

      DiscontinuousFunc<Hermes::Ord>** u_ext_ord = current_u_ext == nullptr ? nullptr : new DiscontinuousFunc<Hermes::Ord>*[this->rungeKutta ? this->RK_original_spaces_count : mfDG->wf->get_neq() - mfDG->u_ext_offset];

      if (current_u_ext)
      for (unsigned short i = 0; i < prev_size; i++)
      if (current_u_ext[i + mfDG->u_ext_offset])
        u_ext_ord[i] = init_ext_fn_ord(nbs_u, current_u_ext[i + mfDG->u_ext_offset]);
      else
        u_ext_ord[i] = new DiscontinuousFunc<Ord>(&func_order[0], false, false);

      // Order of additional external functions.
      DiscontinuousFunc<Ord>** ext_ord = nullptr;
      MeshFunctionSharedPtrVector<Scalar> ext_ord_fns = mfDG->ext.size() ? mfDG->ext : mfDG->wf->ext;
      if (ext_ord_fns.size() > 0)
        ext_ord = init_ext_fns_ord(ext_ord_fns, neighbor_searches);

      // Order of shape functions.
      int max_order_j = spaces[mfDG->j]->get_element_order(current_state->e[mfDG->j]->id);
      int max_order_i = spaces[mfDG->i]->get_element_order(current_state->e[mfDG->i]->id);
      if (H2D_GET_V_ORDER(max_order_i) > H2D_GET_H_ORDER(max_order_i))
        max_order_i = H2D_GET_V_ORDER(max_order_i);
      else
        max_order_i = H2D_GET_H_ORDER(max_order_i);
      if (H2D_GET_V_ORDER(max_order_j) > H2D_GET_H_ORDER(max_order_j))
        max_order_j = H2D_GET_V_ORDER(max_order_j);
      else
        max_order_j = H2D_GET_H_ORDER(max_order_j);

      // Order of shape functions.
      DiscontinuousFunc<Ord>* ou = new DiscontinuousFunc<Ord>(&func_order[max_order_j], neighbor_supp_u);
      DiscontinuousFunc<Ord>* ov = new DiscontinuousFunc<Ord>(&func_order[max_order_i], neighbor_supp_v);

      // Order of geometric attributes (eg. for multiplication of a solution with coordinates, normals, etc.).

      // Total order of the matrix form.
      Ord o = mfDG->ord(1, &wt_order, u_ext_ord, ou, ov, &geom_order_surf, ext_ord);

      adjust_order_to_refmaps(mfDG, order, &o, current_refmaps);

      // Cleanup.
      deinit_ext_fns_ord(mfDG, u_ext_ord, ext_ord);

      delete ou;
      delete ov;

      return order;
    }

    template<typename Scalar>
    int DiscreteProblemIntegrationOrderCalculator<Scalar>::calc_order_dg_vector_form(const SpaceSharedPtrVector<Scalar> spaces, Traverse::State* current_state, VectorFormDG<Scalar>* vfDG, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, bool neighbor_supp_v, NeighborSearch<Scalar>** neighbor_searches)
    {
      NeighborSearch<Scalar>* nbs_u = neighbor_searches[vfDG->i];

      unsigned short prev_size = this->rungeKutta ? this->RK_original_spaces_count : vfDG->wf->get_neq() - vfDG->u_ext_offset;

      // Order to return.
      int order = 0;

      DiscontinuousFunc<Hermes::Ord>** u_ext_ord = current_u_ext == nullptr ? nullptr : new DiscontinuousFunc<Hermes::Ord>*[this->rungeKutta ? this->RK_original_spaces_count : vfDG->wf->get_neq() - vfDG->u_ext_offset];

      if (current_u_ext)
      for (unsigned short i = 0; i < prev_size; i++)
      if (current_u_ext[i + vfDG->u_ext_offset])
        u_ext_ord[i] = init_ext_fn_ord(nbs_u, current_u_ext[i + vfDG->u_ext_offset]);
      else
        u_ext_ord[i] = new DiscontinuousFunc<Ord>(&func_order[0], false, false);

      // Order of additional external functions.
      DiscontinuousFunc<Ord>** ext_ord = nullptr;
      MeshFunctionSharedPtrVector<Scalar> ext_ord_fns = vfDG->ext.size() ? vfDG->ext : vfDG->wf->ext;
      if (ext_ord_fns.size() > 0)
        ext_ord = init_ext_fns_ord(ext_ord_fns, neighbor_searches);

      // Order of shape functions.
      int max_order_i = spaces[vfDG->i]->get_element_order(current_state->e[vfDG->i]->id);
      if (H2D_GET_V_ORDER(max_order_i) > H2D_GET_H_ORDER(max_order_i))
        max_order_i = H2D_GET_V_ORDER(max_order_i);
      else
        max_order_i = H2D_GET_H_ORDER(max_order_i);

      // Order of shape functions.
      DiscontinuousFunc<Ord>* ov = new DiscontinuousFunc<Ord>(&func_order[max_order_i], neighbor_supp_v);

      // Total order of the matrix form.
      Ord o = vfDG->ord(1, &wt_order, u_ext_ord, ov, &geom_order_surf, ext_ord);

      adjust_order_to_refmaps(vfDG, order, &o, current_refmaps);

      // Cleanup.
      deinit_ext_fns_ord(vfDG, u_ext_ord, ext_ord);
      delete ov;

      return order;
    }

    template class HERMES_API DiscreteProblemIntegrationOrderCalculator<double>;
    template class HERMES_API DiscreteProblemIntegrationOrderCalculator<std::complex<double> >;
  }
}