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

#include "../h2d_common.h"
#include "space_hcurl.h"
#include "../../../hermes_common/matrix.h"
#include "../quad_all.h"
#include "../shapeset/shapeset_hc_all.h"


double** HcurlSpace::hcurl_proj_mat = NULL;
double*  HcurlSpace::hcurl_chol_p   = NULL;
int      HcurlSpace::hcurl_proj_ref = 0;

void HcurlSpace::init(Shapeset* shapeset, Ord2 p_init)
{
  if (shapeset == NULL)
  {
    this->shapeset = new HcurlShapeset;
    own_shapeset = true;
  }
  if (this->shapeset->get_num_components() < 2) error("HcurlSpace requires a vector shapeset.");

  if (!hcurl_proj_ref++)
  {
    precalculate_projection_matrix(0, hcurl_proj_mat, hcurl_chol_p);
  }

  proj_mat = hcurl_proj_mat;
  chol_p   = hcurl_chol_p;

  // set uniform poly order in elements
  if (p_init.order_h < 0 || p_init.order_v < 0) error("P_INIT must be >= 0 in an Hcurl space.");
  else this->set_uniform_order_internal(p_init);

  // enumerate basis functions
  this->assign_dofs();
}

HcurlSpace::HcurlSpace(Mesh* mesh, BCTypes* bc_types, int p_init, Shapeset* shapeset) 
  : Space(mesh, shapeset, bc_types, (BCValues*) NULL, Ord2(p_init, p_init))
{
  init(shapeset, Ord2(p_init, p_init));
}

HcurlSpace::HcurlSpace(Mesh* mesh, BCTypes*  bc_types, Ord2 p_init, Shapeset* shapeset)
  : Space(mesh, shapeset, bc_types, (BCValues*) NULL, p_init)
{
  init(shapeset, p_init);
}

HcurlSpace::HcurlSpace(Mesh* mesh, BCTypes* bc_types, BCValues* bc_values, int p_init, Shapeset* shapeset) 
  : Space(mesh, shapeset, bc_types, bc_values, Ord2(p_init, p_init))
{
  init(shapeset, Ord2(p_init, p_init));
}

HcurlSpace::HcurlSpace(Mesh* mesh, BCTypes*  bc_types, BCValues* bc_values, Ord2 p_init, Shapeset* shapeset)
  : Space(mesh, shapeset, bc_types, bc_values, p_init)
{
  init(shapeset, p_init);
}

// All the following constructors are DEPRECATED!
HcurlSpace::HcurlSpace(Mesh* mesh, BCTypes* bc_types,
                 scalar (*bc_value_callback_by_coord)(int, double, double), int p_init, 
                 Shapeset* shapeset) : Space(mesh, shapeset, bc_types, bc_value_callback_by_coord, Ord2(p_init, p_init))
{
  if (shapeset == NULL)
  {
    this->shapeset = new HcurlShapeset;
    own_shapeset = true;
  }
  if (this->shapeset->get_num_components() < 2) error("HcurlSpace requires a vector shapeset.");

  if (!hcurl_proj_ref++)
  {
    precalculate_projection_matrix(0, hcurl_proj_mat, hcurl_chol_p);
  }

  proj_mat = hcurl_proj_mat;
  chol_p   = hcurl_chol_p;

  // set uniform poly order in elements
  if (p_init < 0) error("P_INIT must be >= 0 in an Hcurl space.");
  else this->set_uniform_order_internal(Ord2(p_init, p_init));

  // enumerate basis functions
  this->assign_dofs();
}

HcurlSpace::HcurlSpace(Mesh* mesh, BCTypes*  bc_types,
		       scalar (*bc_value_callback_by_coord)(int, double, double), Ord2 p_init,
                       Shapeset* shapeset)
          : Space(mesh, shapeset, bc_types, bc_value_callback_by_coord, p_init)
{
  if (shapeset == NULL)
  {
    this->shapeset = new HcurlShapeset;
    own_shapeset = true;
  }
  if (this->shapeset->get_num_components() < 2) error("HcurlSpace requires a vector shapeset.");

  if (!hcurl_proj_ref++)
  {
    precalculate_projection_matrix(0, hcurl_proj_mat, hcurl_chol_p);
  }

  proj_mat = hcurl_proj_mat;
  chol_p   = hcurl_chol_p;

  // set uniform poly order in elements
  if (p_init.order_h < 0 || p_init.order_v < 0) error("P_INIT must be >= 0 in an Hcurl space.");
  else this->set_uniform_order_internal(p_init);

  // enumerate basis functions
  this->assign_dofs();
}

HcurlSpace::HcurlSpace(Mesh* mesh, BCType (*bc_type_callback)(int), 
                 scalar (*bc_value_callback_by_coord)(int, double, double), int p_init, 
                 Shapeset* shapeset) : Space(mesh, shapeset, bc_type_callback, bc_value_callback_by_coord, Ord2(p_init, p_init))
{
  if (shapeset == NULL)
  {
    this->shapeset = new HcurlShapeset;
    own_shapeset = true;
  }
  if (this->shapeset->get_num_components() < 2) error("HcurlSpace requires a vector shapeset.");

  if (!hcurl_proj_ref++)
  {
    precalculate_projection_matrix(0, hcurl_proj_mat, hcurl_chol_p);
  }

  proj_mat = hcurl_proj_mat;
  chol_p   = hcurl_chol_p;

  // set uniform poly order in elements
  if (p_init < 0) error("P_INIT must be >= 0 in an Hcurl space.");
  else this->set_uniform_order_internal(Ord2(p_init, p_init));

  // enumerate basis functions
  this->assign_dofs();
}

HcurlSpace::HcurlSpace(Mesh* mesh, BCType (*bc_type_callback)(int), 
		       scalar (*bc_value_callback_by_coord)(int, double, double), Ord2 p_init,
                       Shapeset* shapeset)
          : Space(mesh, shapeset, bc_type_callback, bc_value_callback_by_coord, p_init)
{
  if (shapeset == NULL)
  {
    this->shapeset = new HcurlShapeset;
    own_shapeset = true;
  }
  if (this->shapeset->get_num_components() < 2) error("HcurlSpace requires a vector shapeset.");

  if (!hcurl_proj_ref++)
  {
    precalculate_projection_matrix(0, hcurl_proj_mat, hcurl_chol_p);
  }

  proj_mat = hcurl_proj_mat;
  chol_p   = hcurl_chol_p;

  // set uniform poly order in elements
  if (p_init.order_h < 0 || p_init.order_v < 0) error("P_INIT must be >= 0 in an Hcurl space.");
  else this->set_uniform_order_internal(p_init);

  // enumerate basis functions
  this->assign_dofs();
}

HcurlSpace::~HcurlSpace()
{
  if (!--hcurl_proj_ref)
  {
    delete [] hcurl_proj_mat;
    delete [] hcurl_chol_p;
  }
  if (own_shapeset)
    delete this->shapeset;
}


Space* HcurlSpace::dup(Mesh* mesh) const
{
 // FIXME:
  HcurlSpace* space = new HcurlSpace(mesh, this->bc_types,
          this->bc_value_callback_by_coord, Ord2(0,0), this->shapeset);
  space->copy_callbacks(this);
  return space;
}

void HcurlSpace::set_shapeset(Shapeset *shapeset)
{
  if(shapeset->get_id() < 20 && shapeset->get_id() > 9)
  {
    this->shapeset = shapeset;
    own_shapeset = false;
  }
  else
    error("Wrong shapeset type in HcurlSpace::set_shapeset()");
}

// Sets element order and updates enumeration of dofs. Intended for 
// the user.
void HcurlSpace::set_element_order(int id, int order)
{
  set_element_order_internal(id, order);

  // since space changed, call assign_dofs()
  this->assign_dofs();
}

// Sets element order without updating the enumeration of dofs. For internal use.
void HcurlSpace::set_element_order_internal(int id, int order)
{
  assert_msg(mesh->get_element(id)->is_quad() || H2D_GET_V_ORDER(order) == 0, "Element #%d is triangle but vertical order is not zero", id);
  if (id < 0 || id >= mesh->get_max_element_id())
    error("Invalid element id.");
  H2D_CHECK_ORDER(order);

  resize_tables();
  edata[id].order = order;
  seq++;
}

//// dof assignment ////////////////////////////////////////////////////////////////////////////////

void HcurlSpace::assign_edge_dofs()
{
  Node* en;
  for_all_edge_nodes(en, mesh)
  {
    if (en->ref > 1 || en->bnd || mesh->peek_vertex_node(en->p1, en->p2) != NULL)
    {
      int ndofs = get_edge_order_internal(en) + 1;
      ndata[en->id].n = ndofs;

      if (en->bnd && this->bc_types->get_type(en->marker) == BC_ESSENTIAL)
      {
        ndata[en->id].dof = -1;
      }
      else
      {
        ndata[en->id].dof = next_dof;
        next_dof += ndofs * stride;
      }
    }
    else
    {
      ndata[en->id].n = -1;
    }
  }
}


void HcurlSpace::assign_bubble_dofs()
{
  Element* e;
  for_all_active_elements(e, mesh)
  {
    shapeset->set_mode(e->get_mode());
    ElementData* ed = &edata[e->id];
    ed->bdof = next_dof;
    ed->n = shapeset->get_num_bubbles(ed->order);
    next_dof += ed->n * stride;
  }
}

//// assembly lists ////////////////////////////////////////////////////////////////////////////////

void HcurlSpace::get_boundary_assembly_list_internal(Element* e, int surf_num, AsmList* al)
{
  Node* en = e->en[surf_num];
  NodeData* nd = &ndata[en->id];

  if (nd->n >= 0) // unconstrained
  {
    if (nd->dof >= 0)
    {
      int ori = (e->vn[surf_num]->id < e->vn[e->next_vert(surf_num)]->id) ? 0 : 1;
      for (int j = 0, dof = nd->dof; j < nd->n; j++, dof += stride)
        al->add_triplet(shapeset->get_edge_index(surf_num, ori, j), dof, 1.0);
    }
    else
    {
      for (int j = 0; j < nd->n; j++)
        al->add_triplet(shapeset->get_edge_index(surf_num, 0, j), -1, nd->edge_bc_proj[j]);
    }
  }
  else // constrained
  {
    int part = nd->part;
    int ori = part < 0 ? 1 : 0;
    if (part < 0) part ^= ~0;

    nd = &ndata[nd->base->id]; // ccc
    for (int j = 0, dof = nd->dof; j < nd->n; j++, dof += stride)
      al->add_triplet(shapeset->get_constrained_edge_index(surf_num, j, ori, part), dof, 1.0);
  }
}

//// BC stuff //////////////////////////////////////////////////////////////////////////////////////

scalar* HcurlSpace::get_bc_projection(SurfPos* surf_pos, int order)
{
  _F_
  assert(order >= 0);
  scalar* proj = new scalar[order + 1];

  Quad1DStd quad1d;
  scalar* rhs = proj;
  int mo = quad1d.get_max_order();
  double2* pt = quad1d.get_points(mo);

  Node* vn1 = mesh->get_node(surf_pos->v1);
  Node* vn2 = mesh->get_node(surf_pos->v2);
  double el = sqrt(sqr(vn1->x - vn2->x) + sqr(vn1->y - vn2->y));
  el *= 0.5 * (surf_pos->hi - surf_pos->lo);

  // get boundary values at integration points, construct rhs
  for (int i = 0; i <= order; i++)
  {
    rhs[i] = 0.0;
    int ii = shapeset->get_edge_index(0, 0, i);
    for (int j = 0; j < quad1d.get_num_points(mo); j++)
    {
      double t = (pt[j][0] + 1) * 0.5, s = 1.0 - t;
      surf_pos->t = surf_pos->lo * s + surf_pos->hi * t;
      // If the "old" callbacks are used.
      if(surf_pos->space->bc_value_callback_by_coord != NULL)
        rhs[i] += pt[j][1] * shapeset->get_fn_value(ii, pt[j][0], -1.0, 0)
                         * (bc_value_callback_by_edge(surf_pos)) * el;
      // If BCValues class is used.
      else {
        // If the BC on this part of the boundary is constant.
        if(bc_values->is_const(surf_pos->marker))
          rhs[i] += pt[j][1] * shapeset->get_fn_value(ii, pt[j][0], -1.0, 0)
          * (bc_values->calculate(surf_pos->marker)) * el;
        // If the BC is not constant.
        else {
          // Find out the (x,y) coordinate.
          double x, y;
          Nurbs* nurbs = surf_pos->base->is_curved() ? surf_pos->base->cm->nurbs[surf_pos->surf_num] : NULL;
          nurbs_edge(surf_pos->base, nurbs, surf_pos->surf_num, 2.0*surf_pos->t - 1.0, x, y);
          // Calculate.
          rhs[i] += pt[j][1] * shapeset->get_fn_value(ii, pt[j][0], -1.0, 0)
            * (bc_values->calculate(surf_pos->marker, x, y)) * el;
        }
      }
    }
  }

  // solve the system using a precalculated Cholesky decomposed projection matrix
  cholsl(proj_mat, order + 1, chol_p, rhs, rhs);
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

void HcurlSpace::update_constrained_nodes(Element* e, EdgeInfo* ei0, EdgeInfo* ei1, EdgeInfo* ei2, EdgeInfo* ei3)
{
  int j;
  EdgeInfo* ei[4] = { ei0, ei1, ei2, ei3 };
  NodeData* nd;

  // on non-refined elements all we have to do is update edge nodes lying on constrained edges
  if (e->active)
  {
    for (unsigned int i = 0; i < e->nvert; i++)
    {
      if (ei[i] != NULL)
      {
        nd = &ndata[e->en[i]->id];
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
          Node* mid_en = mesh->peek_edge_node(e->vn[i]->id, e->vn[j]->id);
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


void HcurlSpace::update_constraints()
{
  Element* e;
  for_all_base_elements(e, mesh)
    update_constrained_nodes(e, NULL, NULL, NULL, NULL);
}
