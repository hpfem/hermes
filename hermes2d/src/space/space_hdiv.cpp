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

#include "hermes2d_common_defs.h"
#include "space_hdiv.h"
#include "matrix.h"
#include "quad_all.h"
#include "shapeset/shapeset_hd_all.h"
#include "boundary_conditions/essential_boundary_conditions.h"
namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    double** HdivSpace<Scalar>::hdiv_proj_mat = NULL;
    template<typename Scalar>
    double*  HdivSpace<Scalar>::hdiv_chol_p   = NULL;
    template<typename Scalar>
    int      HdivSpace<Scalar>::hdiv_proj_ref = 0;

    template<typename Scalar>
    void HdivSpace<Scalar>::init(Shapeset* shapeset, Ord2 p_init)
    {
      if (shapeset == NULL)
      {
        this->shapeset = new HdivShapeset;
        this->own_shapeset = true;
      }
      if (this->shapeset->get_num_components() < 2) error("HdivSpace requires a vector shapeset.");

      if (!hdiv_proj_ref++)
      {
        this->precalculate_projection_matrix(0, hdiv_proj_mat, hdiv_chol_p);
      }

      this->proj_mat = hdiv_proj_mat;
      this->chol_p   = hdiv_chol_p;

      // set uniform poly order in elements
      if (p_init.order_h < 0 || p_init.order_v < 0) error("P_INIT must be >= 0 in an Hdiv space.");
      else this->set_uniform_order_internal(p_init, HERMES_ANY_INT);

      // enumerate basis functions
      this->assign_dofs();
    }

    template<typename Scalar>
    HdivSpace<Scalar>::HdivSpace(Mesh* mesh, EssentialBCs<Scalar>* essential_bcs, int p_init, Shapeset* shapeset)
      : Space<Scalar>(mesh, shapeset, essential_bcs, Ord2(p_init, p_init))
    {
      _F_;
      init(shapeset, Ord2(p_init, p_init));
    }

    template<typename Scalar>
    HdivSpace<Scalar>::HdivSpace(Mesh* mesh, int p_init, Shapeset* shapeset)
      : Space<Scalar>(mesh, shapeset, NULL, Ord2(p_init, p_init))
    {
      _F_;
      init(shapeset, Ord2(p_init, p_init));
    }

    template<typename Scalar>
    HdivSpace<Scalar>::~HdivSpace()
    {
      if (!--hdiv_proj_ref)
      {
        delete [] hdiv_proj_mat;
        delete [] hdiv_chol_p;
      }
      if (this->own_shapeset)
        delete this->shapeset;
    }


    template<typename Scalar>
    Space<Scalar>* HdivSpace<Scalar>::dup(Mesh* mesh, int order_increase) const
    {
      // FIXME
      // HdivSpace<Scalar>* space = new HdivSpace(mesh, essential_bcs, 0, this->shapeset);
      // space->copy_callbacks(this);
      // space->copy_orders(this, order_increase);
      // return space;
      return NULL;
    }

    template<typename Scalar>
    void HdivSpace<Scalar>::set_shapeset(Shapeset *shapeset)
    {
      if(shapeset->get_id() < 30 && shapeset->get_id() > 19)
      {
        this->shapeset = shapeset;
        this->own_shapeset = false;
      }
      else
        error("Wrong shapeset type in HdivSpace<Scalar>::set_shapeset()");
    }

    //// dof assignment ////////////////////////////////////////////////////////////////////////////////

    template<typename Scalar>
    void HdivSpace<Scalar>::assign_edge_dofs()
    {
      Node* en;
      for_all_edge_nodes(en, this->mesh)
      {
        if (en->ref > 1 || en->bnd || this->mesh->peek_vertex_node(en->p1, en->p2) != NULL)
        {
          int ndofs = this->get_edge_order_internal(en) + 1;
          this->ndata[en->id].n = ndofs;

          if (en->bnd)
            if(this->essential_bcs != NULL)
              if(this->essential_bcs->get_boundary_condition(this->mesh->get_boundary_markers_conversion().get_user_marker(en->marker)) != NULL)
                this->ndata[en->id].dof = this->H2D_CONSTRAINED_DOF;
              else 
              {
                this->ndata[en->id].dof = this->next_dof;
                this->next_dof += ndofs * this->stride;
              }
            else 
            {
              this->ndata[en->id].dof = this->next_dof;
              this->next_dof += ndofs * this->stride;
            }
          else 
          {
            this->ndata[en->id].dof = this->next_dof;
            this->next_dof += ndofs * this->stride;
          }
        }
        else
        {
          this->ndata[en->id].n = -1;
        }
      }
    }


    template<typename Scalar>
    void HdivSpace<Scalar>::assign_bubble_dofs()
    {
      Element* e;
      for_all_active_elements(e, this->mesh)
      {
        this->shapeset->set_mode(e->get_mode());
        typename Space<Scalar>::ElementData* ed = &this->edata[e->id];
        ed->bdof = this->next_dof;
        ed->n = this->shapeset->get_num_bubbles(ed->order);
        this->next_dof += ed->n * this->stride;
      }
    }



    //// assembly lists ////////////////////////////////////////////////////////////////////////////////

    template<typename Scalar>
    void HdivSpace<Scalar>::get_boundary_assembly_list_internal(Element* e, int surf_num, AsmList<Scalar>* al)
    {
      Node* en = e->en[surf_num];
      typename Space<Scalar>::NodeData* nd = &this->ndata[en->id];

      if (nd->n >= 0) // unconstrained
      {
        if (nd->dof >= 0)
        {
          int ori = (e->vn[surf_num]->id < e->vn[e->next_vert(surf_num)]->id) ? 0 : 1;
          for (int j = 0, dof = nd->dof; j < nd->n; j++, dof += this->stride)
            al->add_triplet(this->shapeset->get_edge_index(surf_num, ori, j), dof, 1.0);
        }
        else
        {
          for (int j = 0; j < nd->n; j++)
            al->add_triplet(this->shapeset->get_edge_index(surf_num, 0, j), -1, nd->edge_bc_proj[j]);
        }
      }
      else // constrained
      {
        int part = nd->part;
        int ori = part < 0 ? 1 : 0;
        if (part < 0) part ^= ~0;

        nd = &this->ndata[nd->base->id]; // ccc
        for (int j = 0, dof = nd->dof; j < nd->n; j++, dof += this->stride)
          al->add_triplet(this->shapeset->get_constrained_edge_index(surf_num, j, ori, part), dof, 1.0);
      }
    }


    template<typename Scalar>
    void HdivSpace<Scalar>::get_bubble_assembly_list(Element* e, AsmList<Scalar>* al)
    {
      typename Space<Scalar>::ElementData* ed = &this->edata[e->id];
      if (!ed->n) return;

      int* indices = this->shapeset->get_bubble_indices(ed->order);
      for (int i = 0, dof = ed->bdof; i < ed->n; i++, dof += this->stride)
        al->add_triplet(*indices++, dof, 1.0);
    }


    //// BC stuff //////////////////////////////////////////////////////////////////////////////////////

    template<typename Scalar>
    Scalar* HdivSpace<Scalar>::get_bc_projection(SurfPos* surf_pos, int order)
    {
      assert(order >= 0);
      Scalar* proj = new Scalar[order + 1];

      Quad1DStd quad1d;
      Scalar* rhs = proj;
      int mo = quad1d.get_max_order();
      double2* pt = quad1d.get_points(mo);

      Node* vn1 = this->mesh->get_node(surf_pos->v1);
      Node* vn2 = this->mesh->get_node(surf_pos->v2);
      double el = sqrt(sqr(vn1->x - vn2->x) + sqr(vn1->y - vn2->y));
      el *= 0.5 * (surf_pos->hi - surf_pos->lo);

      // get boundary values at integration points, construct rhs
      for (int i = 0; i <= order; i++)
      {
        rhs[i] = 0.0;
        int ii = this->shapeset->get_edge_index(0, 0, i);
        for (int j = 0; j < quad1d.get_num_points(mo); j++)
        {
          double t = (pt[j][0] + 1) * 0.5, s = 1.0 - t;
          surf_pos->t = surf_pos->lo * s + surf_pos->hi * t;

          // If the BC on this part of the boundary is constant.
          EssentialBoundaryCondition<Scalar> *bc = this->essential_bcs->get_boundary_condition(this->mesh->get_boundary_markers_conversion().get_user_marker(surf_pos->marker));

          if (bc->get_value_type() == EssentialBoundaryCondition<Scalar>::BC_CONST)
          {
            rhs[i] += pt[j][1] * this->shapeset->get_fn_value(ii, pt[j][0], -1.0, 1)
              * bc->value_const * el;
          }
          // If the BC is not constant.
          else if (bc->get_value_type() == EssentialBoundaryCondition<Scalar>::BC_FUNCTION)
          {
            // Find out the (x,y) coordinate.
            double x, y, n_x, n_y, t_x, t_y;
            Nurbs* nurbs = surf_pos->base->is_curved() ? surf_pos->base->cm->nurbs[surf_pos->surf_num] : NULL;
            CurvMap::nurbs_edge(surf_pos->base, nurbs, surf_pos->surf_num, 2.0*surf_pos->t - 1.0, x, y, n_x, n_y, t_x, t_y);
            // Calculate.
            rhs[i] += pt[j][1] * this->shapeset->get_fn_value(ii, pt[j][0], -1.0, 1)
              * bc->value(x, y, n_x, n_y, t_x, t_y) * el;
          }
        }
      }

      // solve the system using a precalculated Cholesky decomposed projection matrix
      cholsl(this->proj_mat, order + 1, this->chol_p, rhs, rhs);

      for (int i = 0; i < order+1; i++)
        proj[i] = 0.0;

      return proj;
    }



    //// hanging nodes /////////////////////////////////////////////////////////////////////////////////

    static Node* get_mid_edge_vertex_node(Element* e, int i, int j)
    {
      if (e->is_triangle()) return e->sons[3]->vn[e->prev_vert(i)];
      else if (e->sons[2] == NULL) return i == 1 ? e->sons[0]->vn[2] : i == 3 ? e->sons[0]->vn[3] : NULL;
      else if (e->sons[0] == NULL) return i == 0 ? e->sons[2]->vn[1] : i == 2 ? e->sons[2]->vn[2] : NULL;
      else return e->sons[i]->vn[j];
    }

    template<typename Scalar>
    void HdivSpace<Scalar>::update_constrained_nodes(Element* e, EdgeInfo* ei0, EdgeInfo* ei1, EdgeInfo* ei2, EdgeInfo* ei3)
    {
      int j;
      EdgeInfo* ei[4] = { ei0, ei1, ei2, ei3 };
      typename Space<Scalar>::NodeData* nd;

      // on non-refined elements all we have to do is update edge nodes lying on constrained edges
      if (e->active)
      {
        for (unsigned int i = 0; i < e->nvert; i++)
        {
          if (ei[i] != NULL)
          {
            nd = &this->ndata[e->en[i]->id];
            nd->base = ei[i]->node;
            nd->part = ei[i]->part;
            if (ei[i]->ori) nd->part ^= ~0;
          }
        }
      }
      // the element has sons - update mid-edge constrained vertex nodes
      else
      {

        // create new edge infos where we don't have them yet
        EdgeInfo ei_data[4];
        for (unsigned int i = 0; i < e->nvert; i++)
        {
          if (ei[i] == NULL)
          {
            j = e->next_vert(i);
            Node* mid_vn = get_mid_edge_vertex_node(e, i, j);
            if (mid_vn != NULL && mid_vn->is_constrained_vertex())
            {
              Node* mid_en = this->mesh->peek_edge_node(e->vn[i]->id, e->vn[j]->id);
              if (mid_en != NULL)
              {
                ei[i] = ei_data + i;
                ei[i]->node = mid_en;
                ei[i]->part = -1;
                ei[i]->lo = -1.0;
                ei[i]->hi =  1.0;
                ei[i]->ori = (e->vn[i]->id < e->vn[j]->id) ? 0 : 1;
              }
            }
          }
        }

        // create edge infos for half-edges
        EdgeInfo  half_ei_data[4][2];
        EdgeInfo* half_ei[4][2];
        for (unsigned int i = 0; i < e->nvert; i++)
        {
          if (ei[i] == NULL)
          {
            half_ei[i][0] = half_ei[i][1] = NULL;
          }
          else
          {
            EdgeInfo* h0 = half_ei[i][0] = half_ei_data[i];
            EdgeInfo* h1 = half_ei[i][1] = half_ei_data[i] + 1;

            h0->node = h1->node = ei[i]->node;
            h0->part = (ei[i]->part + 1) * 2;
            h1->part = h0->part + 1;
            h0->hi = h1->lo = (ei[i]->lo + ei[i]->hi) / 2;
            h0->lo = ei[i]->lo;
            h1->hi = ei[i]->hi;
            h1->ori = h0->ori = ei[i]->ori;
          }
        }

        // recur to sons
        if (e->is_triangle())
        {
          update_constrained_nodes(e->sons[0], half_ei[0][0], NULL, half_ei[2][1], NULL);
          update_constrained_nodes(e->sons[1], half_ei[0][1], half_ei[1][0], NULL, NULL);
          update_constrained_nodes(e->sons[2], NULL, half_ei[1][1], half_ei[2][0], NULL);
          update_constrained_nodes(e->sons[3], NULL, NULL, NULL, NULL);
        }
        else if (e->sons[2] == NULL) // 'horizontally' split quad
        {
          update_constrained_nodes(e->sons[0], ei[0], half_ei[1][0], NULL, half_ei[3][1]);
          update_constrained_nodes(e->sons[1], NULL, half_ei[1][1], ei[2], half_ei[3][0]);
        }
        else if (e->sons[0] == NULL) // 'vertically' split quad
        {
          update_constrained_nodes(e->sons[2], half_ei[0][0], NULL, half_ei[2][1], ei[3]);
          update_constrained_nodes(e->sons[3], half_ei[0][1], ei[1], half_ei[2][0], NULL);
        }
        else // fully split quad
        {
          update_constrained_nodes(e->sons[0], half_ei[0][0], NULL, NULL, half_ei[3][1]);
          update_constrained_nodes(e->sons[1], half_ei[0][1], half_ei[1][0], NULL, NULL);
          update_constrained_nodes(e->sons[2], NULL, half_ei[1][1], half_ei[2][0], NULL);
          update_constrained_nodes(e->sons[3], NULL, NULL, half_ei[2][1], half_ei[3][0]);
        }
      }
    }


    template<typename Scalar>
    void HdivSpace<Scalar>::update_constraints()
    {
      Element* e;
      for_all_base_elements(e, this->mesh)
        update_constrained_nodes(e, NULL, NULL, NULL, NULL);
    }

    template HERMES_API class HdivSpace<double>;
    template HERMES_API class HdivSpace<std::complex<double> >;
  }
}