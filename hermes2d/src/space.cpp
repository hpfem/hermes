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

#include "common.h"
#include "space.h"
#include "matrix_old.h"
#include "auto_local_array.h"

Space::Space(Mesh* mesh, Shapeset* shapeset, BCType (*bc_type_callback)(int), 
             scalar (*bc_value_callback_by_coord)(int, double, double), int p_init)
     : mesh(mesh), shapeset(shapeset)
{
  if (mesh == NULL) error("Space must be initialized with an existing mesh.");
  this->default_tri_order = -1;
  this->default_quad_order = -1;
  this->ndata = NULL;
  this->edata = NULL;
  this->nsize = esize = 0;
  this->ndata_allocated = 0;
  this->mesh_seq = -1;
  this->seq = 0;
  this->was_assigned = false;
  this->ndof = 0;

  this->set_bc_types_init(bc_type_callback);
  this->set_essential_bc_values(bc_value_callback_by_coord);
  this->set_essential_bc_values((scalar (*)(EdgePos*)) NULL);
}

Space::~Space()
{
  free();
}

void Space::free()
{
  free_extra_data();
  if (nsize) { ::free(ndata); ndata=NULL; }
  if (esize) { ::free(edata); edata=NULL; }
}

//// element orders ///////////////////////////////////////////////////////////////////////////////

void Space::resize_tables()
{
  if ((nsize < mesh->get_max_node_id()) || (ndata == NULL))
  {
    //HACK: definition of allocated size and the result number of elements
    nsize = mesh->get_max_node_id();
    if ((nsize > ndata_allocated) || (ndata == NULL))
    {
      int prev_allocated = ndata_allocated;
      if (ndata_allocated == 0)
        ndata_allocated = 1024;
      while (ndata_allocated < nsize)
        ndata_allocated = ndata_allocated * 3 / 2;
      ndata = (NodeData*)realloc(ndata, ndata_allocated * sizeof(NodeData));
      for(int i = prev_allocated; i < ndata_allocated; i++)
        ndata[i].edge_bc_proj = NULL;
    }
  }

  if ((esize < mesh->get_max_element_id()) || (edata == NULL))
  {
    int oldsize = esize;
    if (!esize) esize = 1024;
    while (esize < mesh->get_max_element_id()) esize = esize * 3 / 2;
    edata = (ElementData*) realloc(edata, sizeof(ElementData) * esize);
    for (int i = oldsize; i < esize; i++)
      edata[i].order = -1;
  }
}


void Space::H2D_CHECK_ORDER(int order)
{
  if (H2D_GET_H_ORDER(order) < 0 || H2D_GET_V_ORDER(order) < 0)
    error("Order cannot be negative.");
  if (H2D_GET_H_ORDER(order) > 10 || H2D_GET_V_ORDER(order) > 10)
    error("Order = %d, maximum is 10.", order);
}

// if the user calls this, then the enumeration of dof 
// is updated
void Space::set_element_order(int id, int order)
{
  set_element_order_internal(id, order);

  // since space changed, enumerate basis functions
  this->assign_dofs();
}

// just sets the element order without enumerating dof
void Space::set_element_order_internal(int id, int order)
{
  //NOTE: We need to take into account that L2 and Hcurl may use zero orders. The latter has its own version of this method, however.
  assert_msg(mesh->get_element(id)->is_triangle() || get_type() == 3 || H2D_GET_V_ORDER(order) != 0, "Element #%d is quad but given vertical order is zero", id);
  assert_msg(mesh->get_element(id)->is_quad() || H2D_GET_V_ORDER(order) == 0, "Element #%d is triangle but vertical is not zero", id);
  if (id < 0 || id >= mesh->get_max_element_id())
    error("Invalid element id.");
  H2D_CHECK_ORDER(order);

  resize_tables();
  if (mesh->get_element(id)->is_quad() && get_type() != 3 && H2D_GET_V_ORDER(order) == 0) 
     order = H2D_MAKE_QUAD_ORDER(order, order);
  edata[id].order = order;
  seq++;
}


int Space::get_element_order(int id) const
{
  // sanity checks (for internal purposes)
  if (this->mesh == NULL) error("NULL Mesh pointer detected in Space::get_element_order().");
  if(edata == NULL) error("NULL edata detected in Space::get_element_order().");
  if (id >= esize) {
    warn("Element index %d in Space::get_element_order() while maximum is %d.", id, esize);
    error("Wring element index in Space::get_element_order().");
  }
  return edata[id].order;
}


void Space::set_uniform_order(int order, int marker)
{
  set_uniform_order_internal(order, marker);

  // since space changed, enumerate basis functions
  this->assign_dofs();
}
  
void Space::set_uniform_order_internal(int order, int marker)
{
  resize_tables();
  H2D_CHECK_ORDER(order);
  int quad_order = H2D_MAKE_QUAD_ORDER(order, order);

  Element* e;
  for_all_active_elements(e, mesh)
  {
    if (marker == H2D_ANY || e->marker == marker)
    {
      ElementData* ed = &edata[e->id];
      if (e->is_triangle())
        ed->order = order;
      else
        ed->order = quad_order;
    }
  }
  seq++;
}

void Space::set_element_orders(int* elem_orders_)
{
  resize_tables();
  
  Element* e;
  int counter = 0;
  for_all_elements(e, mesh)
  {
    H2D_CHECK_ORDER(elem_orders_[counter]);
    ElementData* ed = &edata[e->id];
    if (e->is_triangle())
      ed->order = elem_orders_[counter];
    else
      ed->order = H2D_MAKE_QUAD_ORDER(elem_orders_[counter], elem_orders_[counter]);
    counter++;
  }
}

void Space::set_default_order(int tri_order, int quad_order)
{
  if (quad_order == -1) quad_order = H2D_MAKE_QUAD_ORDER(tri_order, tri_order);
  default_tri_order = tri_order;
  default_quad_order = quad_order;
}


void Space::copy_orders_recurrent(Element* e, int order)
{
  if (e->active)
    edata[e->id].order = order;
  else
    for (int i = 0; i < 4; i++)
      if (e->sons[i] != NULL)
        copy_orders_recurrent(e->sons[i], order);
}


void Space::copy_orders(Space* space, int inc)
{
  Element* e;
  resize_tables();
  for_all_active_elements(e, space->get_mesh())
  {
    int oo = space->get_element_order(e->id);
    if (oo < 0) error("Source space has an uninitialized order (element id = %d)", e->id);

    int mo = shapeset->get_max_order();
    int lower_limit = (get_type() == 3 || get_type() == 1) ? 0 : 1; // L2 and Hcurl may use zero orders.
    int ho = std::max(lower_limit, std::min(H2D_GET_H_ORDER(oo) + inc, mo));
    int vo = std::max(lower_limit, std::min(H2D_GET_V_ORDER(oo) + inc, mo));
    oo = e->is_triangle() ? ho : H2D_MAKE_QUAD_ORDER(ho, vo);

    H2D_CHECK_ORDER(oo);
    copy_orders_recurrent(mesh->get_element/*sic!*/(e->id), oo);
  }
  seq++;

  // since space changed, enumerate basis functions
  this->assign_dofs();
}


int Space::get_edge_order(Element* e, int edge)
{
  Node* en = e->en[edge];
  if (en->id >= nsize || edge >= (int)e->nvert) return 0;

  if (ndata[en->id].n == -1)
    return get_edge_order_internal(ndata[en->id].base); // constrained node
  else
    return get_edge_order_internal(en);
}


int Space::get_edge_order_internal(Node* en)
{
  assert(en->type == H2D_TYPE_EDGE);
  Element** e = en->elem;
  int o1 = 1000, o2 = 1000;
  assert(e[0] != NULL || e[1] != NULL);

  if (e[0] != NULL)
  {
    if (e[0]->is_triangle() || en == e[0]->en[0] || en == e[0]->en[2])
      o1 = H2D_GET_H_ORDER(edata[e[0]->id].order);
    else
      o1 = H2D_GET_V_ORDER(edata[e[0]->id].order);
  }

  if (e[1] != NULL)
  {
    if (e[1]->is_triangle() || en == e[1]->en[0] || en == e[1]->en[2])
      o2 = H2D_GET_H_ORDER(edata[e[1]->id].order);
    else
      o2 = H2D_GET_V_ORDER(edata[e[1]->id].order);
  }

  if (o1 == 0) return o2 == 1000 ? 0 : o2;
  if (o2 == 0) return o1 == 1000 ? 0 : o1;
  return std::min(o1, o2);
}


void Space::set_mesh(Mesh* mesh)
{
  if (this->mesh == mesh) return;
  free();
  this->mesh = mesh;
  seq++;

  // since space changed, enumerate basis functions
  this->assign_dofs();
}


void Space::propagate_zero_orders(Element* e)
{
  warn_if(get_element_order(e->id) != 0, "zeroing order of an element ID:%d, original order (H:%d; V:%d)", e->id, H2D_GET_H_ORDER(get_element_order(e->id)), H2D_GET_V_ORDER(get_element_order(e->id)));
  set_element_order_internal(e->id, 0);
  if (!e->active)
    for (int i = 0; i < 4; i++)
      if (e->sons[i] != NULL)
        propagate_zero_orders(e->sons[i]);
}


void Space::distribute_orders(Mesh* mesh, int* parents)
{
  int num = mesh->get_max_element_id();
  AUTOLA_OR(int, orders, num+1);
  Element* e;
  for_all_active_elements(e, mesh)
  {
    int p = get_element_order(parents[e->id]);
    if (e->is_triangle() && (H2D_GET_V_ORDER(p) != 0))
      p = std::max(H2D_GET_H_ORDER(p), H2D_GET_V_ORDER(p));
    orders[e->id] = p;
  }
  for_all_active_elements(e, mesh)
    set_element_order_internal(e->id, orders[e->id]);

}


//// dof assignment ////////////////////////////////////////////////////////////////////////////////

int Space::assign_dofs(int first_dof, int stride)
{
  if (first_dof < 0) error("Invalid first_dof.");
  if (stride < 1)    error("Invalid stride.");

  resize_tables();

  Element* e;
  /** \todo Find out whether the following code this is crucial.
   *  If uncommented, this enforces 0 order for all sons if the base element has 0 order.
   *  In this case, an element with 0 order means an element which is left out from solution. */
  //for_all_base_elements(e, mesh)
  //  if (get_element_order(e->id) == 0)
  //    propagate_zero_orders(e);

  //check validity of orders
  for_all_active_elements(e, mesh) {
    if (e->id >= esize || edata[e->id].order < 0) {
      printf("e->id = %d\n", e->id);
      printf("esize = %d\n", esize);
      printf("edata[%d].order = %d\n", e->id, edata[e->id].order);
      error("Uninitialized element order.");
    }
  }

  this->first_dof = next_dof = first_dof;
  this->stride = stride;

  reset_dof_assignment();
  assign_vertex_dofs();
  assign_edge_dofs();
  assign_bubble_dofs();

  free_extra_data();
  update_essential_bc_values();
  update_constraints();
  post_assign();

  mesh_seq = mesh->get_seq();
  was_assigned = true;
  this->ndof = (next_dof - first_dof) / stride;

  return this->ndof;
}

void Space::reset_dof_assignment() {
  // First assume that all vertex nodes are part of a natural BC. the member NodeData::n
  // is misused for this purpose, since it stores nothing at this point. Also assume
  // that all DOFs are unassigned.
  int i, j;
  for (i = 0; i < mesh->get_max_node_id(); i++)
  {
    ndata[i].n = BC_NATURAL;
    ndata[i].dof = H2D_UNASSIGNED_DOF;
  }

  // next go through all boundary edge nodes constituting an essential BC and mark their
  // neighboring vertex nodes also as essential
  Element* e;
  for_all_active_elements(e, mesh)
  {
    for (unsigned int i = 0; i < e->nvert; i++)
    {
      if (e->en[i]->bnd && bc_type_callback(e->en[i]->marker) == BC_ESSENTIAL)
      {
        j = e->next_vert(i);
        ndata[e->vn[i]->id].n = BC_ESSENTIAL;
        ndata[e->vn[j]->id].n = BC_ESSENTIAL;
      }
    }
  }
}

//// assembly lists ///////////////////////////////////////////////////////////////////////////////

void AsmList::enlarge()
{
  cap = !cap ? 256 : cap * 2;
  idx = (int*) realloc(idx, sizeof(int) * cap);
  dof = (int*) realloc(dof, sizeof(int) * cap);
  coef = (scalar*) realloc(coef, sizeof(scalar) * cap);
}


void Space::get_element_assembly_list(Element* e, AsmList* al)
{
  // some checks
  if (e->id >= esize || edata[e->id].order < 0)
    error("Uninitialized element order (id = #%d).", e->id);
  if (!is_up_to_date())
    error("The space is out of date. You need to update it with assign_dofs()"
          " any time the mesh changes.");

  // add vertex, edge and bubble functions to the assembly list
  al->clear();
  shapeset->set_mode(e->get_mode());
  for (unsigned int i = 0; i < e->nvert; i++)
    get_vertex_assembly_list(e, i, al);
  for (unsigned int i = 0; i < e->nvert; i++)
    get_edge_assembly_list_internal(e, i, al);
  get_bubble_assembly_list(e, al);
}


void Space::get_edge_assembly_list(Element* e, int edge, AsmList* al)
{
  al->clear();
  shapeset->set_mode(e->get_mode());
  get_vertex_assembly_list(e, edge, al);
  get_vertex_assembly_list(e, e->next_vert(edge), al);
  get_edge_assembly_list_internal(e, edge, al);
}


void Space::get_bubble_assembly_list(Element* e, AsmList* al)
{
  ElementData* ed = &edata[e->id];

  if (!ed->n) return;

  int* indices = shapeset->get_bubble_indices(ed->order);
  for (int i = 0, dof = ed->bdof; i < ed->n; i++, dof += stride, indices++)
    al->add_triplet(*indices, dof, 1.0);
}

//// BC stuff /////////////////////////////////////////////////////////////////////////////////////

static BCType default_bc_type(int marker)
{
  return BC_NATURAL;
}

static scalar default_bc_value_by_coord(int marker, double x, double y)
{
  return 0;
}

scalar default_bc_value_by_edge(EdgePos* ep)
{
  double x, y;
  Nurbs* nurbs = ep->base->is_curved() ? ep->base->cm->nurbs[ep->edge] : NULL;
  nurbs_edge(ep->base, nurbs, ep->edge, 2.0*ep->t - 1.0, x, y);
  return ep->space->bc_value_callback_by_coord(ep->marker, x, y);
}


void Space::set_bc_types(BCType (*bc_type_callback)(int))
{
  if (bc_type_callback == NULL) bc_type_callback = default_bc_type;
  this->bc_type_callback = bc_type_callback;
  seq++;

  // since space changed, enumerate basis functions
  this->assign_dofs();
}

void Space::set_bc_types_init(BCType (*bc_type_callback)(int))
{
  if (bc_type_callback == NULL) bc_type_callback = default_bc_type;
  this->bc_type_callback = bc_type_callback;
  seq++;
}



void Space::set_essential_bc_values(scalar (*bc_value_callback_by_coord)(int, double, double))
{
  if (bc_value_callback_by_coord == NULL) bc_value_callback_by_coord = default_bc_value_by_coord;
  this->bc_value_callback_by_coord = bc_value_callback_by_coord;
  seq++;
}

void Space::set_essential_bc_values(scalar (*bc_value_callback_by_edge)(EdgePos*))
{
  if (bc_value_callback_by_edge == NULL) bc_value_callback_by_edge = default_bc_value_by_edge;
  this->bc_value_callback_by_edge = bc_value_callback_by_edge;
  seq++;
}


void Space::copy_callbacks(const Space* space)
{
  bc_type_callback = space->bc_type_callback;
  bc_value_callback_by_coord = space->bc_value_callback_by_coord;
  bc_value_callback_by_edge  = space->bc_value_callback_by_edge;
}


void Space::precalculate_projection_matrix(int nv, double**& mat, double*& p)
{
  int n = shapeset->get_max_order() + 1 - nv;
  mat = new_matrix<double>(n, n);
  int component = get_type() == 2 ? 1 : 0;

  Quad1DStd quad1d;
  //shapeset->set_mode(H2D_MODE_TRIANGLE);
  shapeset->set_mode(H2D_MODE_QUAD);
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      int o = i + j + 4;
      double2* pt = quad1d.get_points(o);
      int ii = shapeset->get_edge_index(0, 0, i + nv);
      int ij = shapeset->get_edge_index(0, 0, j + nv);
      double val = 0.0;
      for (int k = 0; k < quad1d.get_num_points(o); k++)
      {
        val += pt[k][1] * shapeset->get_fn_value(ii, pt[k][0], -1.0, component)
                        * shapeset->get_fn_value(ij, pt[k][0], -1.0, component);
      }
      mat[i][j] = val;
    }
  }

  p = new double[n];
  choldc(mat, n, p);
}


void Space::update_edge_bc(Element* e, EdgePos* ep)
{
  if (e->active)
  {
    Node* en = e->en[ep->edge];
    NodeData* nd = &ndata[en->id];
    nd->edge_bc_proj = NULL;

    if (nd->dof != H2D_UNASSIGNED_DOF && en->bnd && bc_type_callback(en->marker) == BC_ESSENTIAL)
    {
      int order = get_edge_order_internal(en);
      ep->marker = en->marker;
      nd->edge_bc_proj = get_bc_projection(ep, order);
      extra_data.push_back(nd->edge_bc_proj);

      int i = ep->edge, j = e->next_vert(i);
      ndata[e->vn[i]->id].vertex_bc_coef = nd->edge_bc_proj + 0;
      ndata[e->vn[j]->id].vertex_bc_coef = nd->edge_bc_proj + 1;
    }
  }
  else
  {
    int son1, son2;
    if (mesh->get_edge_sons(e, ep->edge, son1, son2) == 2)
    {
      double mid = (ep->lo + ep->hi) * 0.5, tmp = ep->hi;
      ep->hi = mid;
      update_edge_bc(e->sons[son1], ep);
      ep->lo = mid; ep->hi = tmp;
      update_edge_bc(e->sons[son2], ep);
    }
    else
      update_edge_bc(e->sons[son1], ep);
  }
}


void Space::update_essential_bc_values()
{
  Element* e;
  for_all_base_elements(e, mesh)
  {
    for (unsigned int i = 0; i < e->nvert; i++)
    {
      int j = e->next_vert(i);
      if (e->vn[i]->bnd && e->vn[j]->bnd)
      {
        EdgePos ep = { e->vn[i]->id, e->vn[j]->id, 0, i, 0.0, 0.0, 1.0, e, this, NULL, NULL };
        update_edge_bc(e, &ep);
      }
    }
  }
}


void Space::free_extra_data()
{
  for (unsigned int i = 0; i < extra_data.size(); i++)
    delete [] (scalar*) extra_data[i];
  extra_data.clear();
}

/*void Space::dump_node_info()
{
  Node* n;
  for_all_nodes(n, mesh)
  {
    NodeData* nd = &ndata[n->id];
    if (n->type == H2D_TYPE_VERTEX)
    {
      printf("vert node id=%d ref=%d bnd=%d x=%g y=%g dof=%d n=%d ",
             n->id, n->ref, n->bnd, n->x, n->y, nd->dof, nd->n);
      if (nd->dof < 0)
        printf("coef=%g", nd->vertex_bc_coef[0]);
    }
    else
    {
      printf("edge node id=%d ref=%d bnd=%d marker=%d p1=%d p2=%d dof=%d n=%d ",
             n->id, n->ref, n->bnd, n->marker, n->p1, n->p2, nd->dof, nd->n);
      if (nd->dof < 0)
        for (int i = 0; i < nd->n; i++)
          printf("proj[%d]=%g ", i, nd->edge_bc_proj[i+2]);
    }
    printf("\n");
  }
}*/

// new way of enumerating degrees of freedom
H2D_API int assign_dofs(Tuple<Space*> spaces) 
{
  int n = spaces.size();
  // assigning dofs to each space
  int ndof = 0;  
  for (int i = 0; i < n; i++) {
    ndof += spaces[i]->assign_dofs(ndof);
  }

  return ndof;
}

// updating time-dependent essential BC
H2D_API void update_essential_bc_values(Tuple<Space*> spaces) {
  int n = spaces.size();
  for (int i = 0; i < n; i++) {
    spaces[i]->update_essential_bc_values();
  }
}

H2D_API void update_essential_bc_values(Space *s) {
  return update_essential_bc_values(Tuple<Space*>(s));
}

