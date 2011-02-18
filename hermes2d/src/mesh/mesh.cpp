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
#include "mesh.h"
#include "h2d_reader.h"


//// nodes, element ////////////////////////////////////////////////////////////////////////////////

void Node::ref_element(Element* e)
{
  if (type == HERMES_TYPE_EDGE)
  {
    // store the element pointer in a free slot of 'elem'
    if (elem[0] == NULL) elem[0] = e;
    else if (elem[1] == NULL) elem[1] = e;
    else assert_msg(false, "No free slot 'elem'");
  }
  ref++;
}


void Node::unref_element(HashTable* ht, Element* e)
{
  if (type == HERMES_TYPE_VERTEX)
  {
    if (!--ref) ht->remove_vertex_node(id);
  }
  else
  {
    // remove the element from the array 'elem'
    if (elem[0] == e) elem[0] = NULL;
    else if (elem[1] == e) elem[1] = NULL;

    if (!--ref) ht->remove_edge_node(id);
  }
}


void Element::ref_all_nodes()
{
  for (unsigned int i = 0; i < nvert; i++)
  {
    vn[i]->ref_element();
    en[i]->ref_element(this);
  }
}


void Element::unref_all_nodes(HashTable* ht)
{
  for (unsigned int i = 0; i < nvert; i++)
  {
    vn[i]->unref_element(ht);
    en[i]->unref_element(ht, this);
  }
}


Element* Element::get_neighbor(int ie) const
{
  Element** elem = en[ie]->elem;
  if (elem[0] == this) return elem[1];
  if (elem[1] == this) return elem[0];
  assert(0);
}


double Element::get_area() const
{
  double ax, ay, bx, by;
  ax = vn[1]->x - vn[0]->x;
  ay = vn[1]->y - vn[0]->y;
  bx = vn[2]->x - vn[0]->x;
  by = vn[2]->y - vn[0]->y;

  double area = 0.5*(ax*by - ay*bx);
  if (is_triangle()) return area;

  ax = bx; ay = by;
  bx = vn[3]->x - vn[0]->x;
  by = vn[3]->y - vn[0]->y;

  return area + 0.5*(ax*by - ay*bx);
}


double Element::get_diameter() const
{
  double max, l;
  if (is_triangle())
  {
    max = 0.0;
    for (int i = 0; i < 3; i++)
    {
      int j = next_vert(i);
      l = sqr(vn[i]->x - vn[j]->x) + sqr(vn[i]->y - vn[j]->y);
      if (l > max) max = l;
    }
  }
  else
  {
    max = sqr(vn[0]->x - vn[2]->x) + sqr(vn[0]->y - vn[2]->y);
    l   = sqr(vn[1]->x - vn[3]->x) + sqr(vn[1]->y - vn[3]->y);
    if (l > max) max = l;
  }
  return sqrt(max);
}


//// mesh //////////////////////////////////////////////////////////////////////////////////////////

unsigned g_mesh_seq = 0;

Mesh::Mesh() : HashTable()
{
  nbase = nactive = ntopvert = ninitial = 0;
  seq = g_mesh_seq++;
  this->markers_conversion = new Mesh::MarkersConversion;
}


Element* Mesh::get_element(int id) const
{
  if (id < 0 || id >= elements.get_size())
    error("Invalid element ID %d, current range: [0; %d]", id, elements.get_size());
  return &(elements[id]);
}


int Mesh::get_edge_sons(Element* e, int edge, int& son1, int& son2)
{
  assert(!e->active);

  if (!e->is_triangle())
  {
    if (e->sons[2] == NULL) // horz quad
    {
      if (edge == 0 || edge == 2) { son1 = edge >> 1;   return 1; }
              else if (edge == 1) { son1 = 0; son2 = 1; return 2; }
                             else { son1 = 1; son2 = 0; return 2; }
    }
    else if (e->sons[0] == NULL) // vert quad
    {
      if (edge == 1 || edge == 3) { son1 = (edge == 1) ? 3 : 2; return 1; }
              else if (edge == 0) { son1 = 2; son2 = 3; return 2; }
                             else { son1 = 3; son2 = 2; return 2; }
    }
  }

  // triangle or 4-son quad
  son1 = edge;
  son2 = e->next_vert(edge);
  return 2;
}


//// low-level refinement //////////////////////////////////////////////////////////////////////////

/*  node and son numbering on a triangle:

                    vn[2]                                       vn[2]

                      *                                           *
                     / \                                         / \
                    /   \                                       /   \
                   /     \                                     /     \
                  /       \                                   / son[2]\
                 /         \                                 /_________\
        en[2]   /           \   en[1]                 vn[0] *           * vn[1]
               *             *                       vn[1]  *-----------*  vn[0]
              /               \                     vn[2] *  \         /  * vn[2]
             /                 \                         / \  \ son[3]/  / \
            /                   \                       /   \  \     /  /   \
           /                     \                     /     \  \   /  /     \
          /                       \                   / son[0]\  \ /  /son[1] \
         /                         \                 /         \  *  /         \
        *-------------*-------------*               *-----------*   *-----------*
                                               vn[0]      vn[1] vn[2] vn[0]      vn[1]
    vn[0]           en[0]           vn[1]


   node and son numbering on a quad:          refinement '0':

    vn[3]           en[2]           vn[2]       vn[3]        vn[2] vn[3]        vn[2]

        *-------------*-------------*               *------------* *------------*
        |                           |               |            | |            |
        |                           |               |            | |            |
        |                           |               |   son[3]   | |   son[2]   |
        |                           |               |            | |            |
        |                           |               |       vn[1]| |vn[0]       |
        |                           |         vn[0] *------------* *------------* vn[1]
 en[3]  *                           *  en[1]  vn[3] *------------* *------------* vn[2]
        |                           |               |       vn[2]| |vn[3]       |
        |                           |               |            | |            |
        |                           |               |   son[0]   | |   son[1]   |
        |                           |               |            | |            |
        |                           |               |            | |            |
        |                           |               *------------* *------------*
        *-------------*-------------*
                                                vn[0]        vn[1] vn[0]        vn[1]
    vn[0]           en[0]           vn[1]


  refinement '1':                             refinement '2':

    vn[3]                           vn[2]       vn[3]        vn[2] vn[3]        vn[2]

        *---------------------------*               *------------* *------------*
        |                           |               |            | |            |
        |                           |               |            | |            |
        |          son[1]           |               |            | |            |
        |                           |               |            | |            |
        |                           |               |            | |            |
  vn[0] *---------------------------* vn[1]         |            | |            |
  vn[3] *---------------------------* vn[2]         |   son[2]   | |   son[3]   |
        |                           |               |            | |            |
        |                           |               |            | |            |
        |          son[0]           |               |            | |            |
        |                           |               |            | |            |
        |                           |               |            | |            |
        *---------------------------*               *------------* *------------*

    vn[0]                           vn[1]       vn[0]        vn[1] vn[0]        vn[1]

*/

Element* Mesh::create_triangle(int marker, Node* v0, Node* v1, Node* v2, CurvMap* cm)
{
  // create a new element
  Element* e = elements.add();
  e->active = 1;
  e->marker = marker;
  e->userdata = 0;
  e->nvert = 3;
  e->iro_cache = -1;
  e->cm = cm;
  e->parent = NULL;
  e->visited = false;

  // set vertex and edge node pointers
  e->vn[0] = v0;
  e->vn[1] = v1;
  e->vn[2] = v2;
  e->en[0] = get_edge_node(v0->id, v1->id);
  e->en[1] = get_edge_node(v1->id, v2->id);
  e->en[2] = get_edge_node(v2->id, v0->id);

  // register in the nodes
  e->ref_all_nodes();

  return e;
}


Element* Mesh::create_quad(int marker, Node* v0, Node* v1, Node* v2, Node* v3, CurvMap* cm)
{
  // create a new element
  Element* e = elements.add();
  e->active = 1;
  e->marker = marker;
  e->userdata = 0;
  e->nvert = 4;
  e->iro_cache = -1;
  e->cm = cm;
  e->parent = NULL;
  e->visited = false;

  // set vertex and edge node pointers
  e->vn[0] = v0;
  e->vn[1] = v1;
  e->vn[2] = v2;
  e->vn[3] = v3;
  e->en[0] = get_edge_node(v0->id, v1->id);
  e->en[1] = get_edge_node(v1->id, v2->id);
  e->en[2] = get_edge_node(v2->id, v3->id);
  e->en[3] = get_edge_node(v3->id, v0->id);

  // register in the nodes
  e->ref_all_nodes();

  return e;
}


static CurvMap* create_son_curv_map(Element* e, int son)
{
  // if the top three bits of part are nonzero, we would overflow
  // -- make the element non-curvilinear
  if (e->cm->part & 0xe000000000000000ULL) return NULL;

  // if the parent element is already almost straight-edged,
  // the son will be even more straight-edged
  if (e->iro_cache == 0) return NULL;

  CurvMap* cm = new CurvMap;
  if (e->cm->toplevel == false)
  {
    cm->parent = e->cm->parent;
    cm->part = (e->cm->part << 3) + son + 1;
  }
  else
  {
    cm->parent = e;
    cm->part = (son + 1);
  }
  cm->toplevel = false;
  cm->order = 4;

  return cm;
}

void Mesh::refine_triangle_to_triangles(Element* e)
{
  // remember the markers of the edge nodes
  int bnd[3] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd    };
  int mrk[3] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker };

  // obtain three mid-edge vertex nodes
  Node* x0 = get_vertex_node(e->vn[0]->id, e->vn[1]->id);
  Node* x1 = get_vertex_node(e->vn[1]->id, e->vn[2]->id);
  Node* x2 = get_vertex_node(e->vn[2]->id, e->vn[0]->id);

  CurvMap* cm[4];
  memset(cm, 0, sizeof(cm));

  // adjust mid-edge coordinates if this is a curved element
  if (e->is_curved())
  {
    double2 pt[3] = { { 0.0,-1.0 }, { 0.0, 0.0 }, { -1.0, 0.0 } };
    e->cm->get_mid_edge_points(e, pt, 3);
    x0->x = pt[0][0]; x0->y = pt[0][1];
    x1->x = pt[1][0]; x1->y = pt[1][1];
    x2->x = pt[2][0]; x2->y = pt[2][1];

    // create CurvMaps for sons (pointer to parent element, part)
    for (int i = 0; i < 4; i++)
      cm[i] = create_son_curv_map(e, i);
  }

  // create the four sons
  Element* sons[4];
  sons[0] = create_triangle(e->marker, e->vn[0], x0, x2, cm[0]);
  sons[1] = create_triangle(e->marker, x0, e->vn[1], x1, cm[1]);
  sons[2] = create_triangle(e->marker, x2, x1, e->vn[2], cm[2]);
  sons[3] = create_triangle(e->marker, x1, x2, x0, cm[3]);

  // update coefficients of curved reference mapping
  for (int i = 0; i < 4; i++)
    if (sons[i]->is_curved())
      sons[i]->cm->update_refmap_coeffs(sons[i]);

  // deactivate this element and unregister from its nodes
  e->active = 0;
  nactive += 3;
  e->unref_all_nodes(this);
  // now the original edge nodes may no longer exist...

  // set correct boundary status and markers for the new nodes
  sons[0]->en[0]->bnd = bnd[0];  sons[0]->en[0]->marker = mrk[0];
  sons[0]->en[2]->bnd = bnd[2];  sons[0]->en[2]->marker = mrk[2];
  sons[1]->en[0]->bnd = bnd[0];  sons[1]->en[0]->marker = mrk[0];
  sons[1]->en[1]->bnd = bnd[1];  sons[1]->en[1]->marker = mrk[1];
  sons[2]->en[1]->bnd = bnd[1];  sons[2]->en[1]->marker = mrk[1];
  sons[2]->en[2]->bnd = bnd[2];  sons[2]->en[2]->marker = mrk[2];
  sons[3]->vn[0]->bnd = bnd[1];
  sons[3]->vn[1]->bnd = bnd[2];
  sons[3]->vn[2]->bnd = bnd[0];

  //set pointers to parent element for sons
  for(int i = 0; i < 4; i++)
	  if(sons[i] != NULL)
		  sons[i]->parent = e;

  // copy son pointers (could not have been done earlier because of the union)
  memcpy(e->sons, sons, 4 * sizeof(Element*));
}

void Mesh::refine_quad(Element* e, int refinement)
{
  int i, j;
  Element* sons[4];

  // remember the markers of the edge nodes
  int bnd[4] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd,    e->en[3]->bnd    };
  int mrk[4] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker, e->en[3]->marker };

  // deactivate this element and unregister from its nodes
  e->active = false;
  nactive--;
  e->unref_all_nodes(this);
  // now the original edge nodes may no longer exist...

  CurvMap* cm[4];
  memset(cm, 0, sizeof(cm));

  // default refinement: one quad to four quads
  if (refinement == 0)
  {
    // obtain four mid-edge vertex nodes and one mid-element vetex node
    Node* x0 = get_vertex_node(e->vn[0]->id, e->vn[1]->id);
    Node* x1 = get_vertex_node(e->vn[1]->id, e->vn[2]->id);
    Node* x2 = get_vertex_node(e->vn[2]->id, e->vn[3]->id);
    Node* x3 = get_vertex_node(e->vn[3]->id, e->vn[0]->id);
    Node* mid = get_vertex_node(x0->id, x2->id);

    // adjust mid-edge coordinates if this is a curved element
    if (e->is_curved())
    {
      double2 pt[5] = { { 0.0,-1.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }, { -1.0, 0.0 }, { 0.0, 0.0 } };
      e->cm->get_mid_edge_points(e, pt, 5);
      x0->x = pt[0][0];  x0->y = pt[0][1];
      x1->x = pt[1][0];  x1->y = pt[1][1];
      x2->x = pt[2][0];  x2->y = pt[2][1];
      x3->x = pt[3][0];  x3->y = pt[3][1];
      mid->x = pt[4][0]; mid->y = pt[4][1];

      // create CurvMaps for sons (pointer to parent element, part)
      for (i = 0; i < 4; i++)
        cm[i] = create_son_curv_map(e, i);
    }

    // create the four sons
    sons[0] = create_quad(e->marker, e->vn[0], x0, mid, x3, cm[0]);
    sons[1] = create_quad(e->marker, x0, e->vn[1], x1, mid, cm[1]);
    sons[2] = create_quad(e->marker, mid, x1, e->vn[2], x2, cm[2]);
    sons[3] = create_quad(e->marker, x3, mid, x2, e->vn[3], cm[3]);
    nactive += 4;

    // set correct boundary markers for the new edge nodes
    for (i = 0; i < 4; i++)
    {
      j = (i > 0) ? i-1 : 3;
      sons[i]->en[j]->bnd = bnd[j];  sons[i]->en[j]->marker = mrk[j];
      sons[i]->en[i]->bnd = bnd[i];  sons[i]->en[i]->marker = mrk[i];
      sons[i]->vn[j]->bnd = bnd[j];
    }
  }
  // refinement '1': one quad to two 'horizontal' quads
  else if (refinement == 1)
  {
    Node* x1 = get_vertex_node(e->vn[1]->id, e->vn[2]->id);
    Node* x3 = get_vertex_node(e->vn[3]->id, e->vn[0]->id);

    // adjust mid-edge coordinates if this is a curved element
    if (e->is_curved())
    {
      double2 pt[2] = { { 1.0, 0.0 }, { -1.0, 0.0 } };
      e->cm->get_mid_edge_points(e, pt, 2);
      x1->x = pt[0][0];  x1->y = pt[0][1];
      x3->x = pt[1][0];  x3->y = pt[1][1];

      // create CurvMaps for sons (pointer to parent element, part)
      for (i = 0; i < 2; i++)
        cm[i] = create_son_curv_map(e, i + 4);
    }

    sons[0] = create_quad(e->marker, e->vn[0], e->vn[1], x1, x3, cm[0]);
    sons[1] = create_quad(e->marker, x3, x1, e->vn[2], e->vn[3], cm[1]);
    sons[2] = sons[3] = NULL;
    nactive += 2;

    sons[0]->en[0]->bnd = bnd[0];  sons[0]->en[0]->marker = mrk[0];
    sons[0]->en[1]->bnd = bnd[1];  sons[0]->en[1]->marker = mrk[1];
    sons[0]->en[3]->bnd = bnd[3];  sons[0]->en[3]->marker = mrk[3];
    sons[1]->en[1]->bnd = bnd[1];  sons[1]->en[1]->marker = mrk[1];
    sons[1]->en[2]->bnd = bnd[2];  sons[1]->en[2]->marker = mrk[2];
    sons[1]->en[3]->bnd = bnd[3];  sons[1]->en[3]->marker = mrk[3];
    sons[0]->vn[2]->bnd = bnd[1];
    sons[0]->vn[3]->bnd = bnd[3];
  }
  // refinement '2': one quad to two 'vertical' quads
  else if (refinement == 2)
  {
    Node* x0 = get_vertex_node(e->vn[0]->id, e->vn[1]->id);
    Node* x2 = get_vertex_node(e->vn[2]->id, e->vn[3]->id);

    // adjust mid-edge coordinates if this is a curved element
    if (e->is_curved())
    {
      double2 pt[2] = { { 0.0, -1.0 }, { 0.0, 1.0 } };
      e->cm->get_mid_edge_points(e, pt, 2);
      x0->x = pt[0][0];  x0->y = pt[0][1];
      x2->x = pt[1][0];  x2->y = pt[1][1];

      // create CurvMaps for sons (pointer to parent element, part)
      for (i = 0; i < 2; i++)
        cm[i] = create_son_curv_map(e, i + 6);
    }

    sons[0] = sons[1] = NULL;
    sons[2] = create_quad(e->marker, e->vn[0], x0, x2, e->vn[3], cm[0]);
    sons[3] = create_quad(e->marker, x0, e->vn[1], e->vn[2], x2, cm[1]);
    nactive += 2;

    sons[2]->en[0]->bnd = bnd[0];  sons[2]->en[0]->marker = mrk[0];
    sons[2]->en[2]->bnd = bnd[2];  sons[2]->en[2]->marker = mrk[2];
    sons[2]->en[3]->bnd = bnd[3];  sons[2]->en[3]->marker = mrk[3];
    sons[3]->en[0]->bnd = bnd[0];  sons[3]->en[0]->marker = mrk[0];
    sons[3]->en[1]->bnd = bnd[1];  sons[3]->en[1]->marker = mrk[1];
    sons[3]->en[2]->bnd = bnd[2];  sons[3]->en[2]->marker = mrk[2];
    sons[2]->vn[1]->bnd = bnd[0];
    sons[2]->vn[2]->bnd = bnd[2];
  }
  else assert(0);

  // update coefficients of curved reference mapping
  for (i = 0; i < 4; i++)
    if (sons[i] != NULL && sons[i]->cm != NULL)
      sons[i]->cm->update_refmap_coeffs(sons[i]);

  // optimization: iro never gets worse
  if (e->iro_cache == 0)
    for (i = 0; i < 4; i++)
      if (sons[i] != NULL)
        sons[i]->iro_cache = 0;

  //set pointers to parent element for sons
  for(int i = 0; i < 4; i++)
	  if(sons[i] != NULL)
		  sons[i]->parent = e;

  // copy son pointers (could not have been done earlier because of the union)
  memcpy(e->sons, sons, sizeof(sons));
}


void Mesh::unrefine_element_internal(Element* e)
{
  assert(!e->active);
  unsigned int i;
  int s1, s2;

  // obtain markers and bnds from son elements
  int mrk[4], bnd[4];
  for (unsigned i = 0; i < e->nvert; i++)
  {
    get_edge_sons(e, i, s1, s2);
    assert(e->sons[s1]->active);
    mrk[i] = e->sons[s1]->en[i]->marker;
    bnd[i] = e->sons[s1]->en[i]->bnd;
  }

  // remove all sons
  for (i = 0; i < 4; i++)
  {
    Element* son = e->sons[i];
    if (son != NULL)
    {
      son->unref_all_nodes(this);
      if (son->cm != NULL) delete son->cm;
      elements.remove(son->id);
      nactive--;
    }
  }

  // recreate edge nodes
  for (i = 0; i < e->nvert; i++)
    e->en[i] = get_edge_node(e->vn[i]->id, e->vn[e->next_vert(i)]->id);

  e->ref_all_nodes();
  e->active = 1;
  nactive++;

  // restore edge node markers and bnds
  for (i = 0; i < e->nvert; i++)
  {
    e->en[i]->marker = mrk[i];
    e->en[i]->bnd = bnd[i];
  }
}


//// high-level element refinement /////////////////////////////////////////////////////////////////

void Mesh::refine_element(int id, int refinement)
{
  Element* e = get_element(id);
  if (!e->used) error("Invalid element id number.");
  if (!e->active) error("Attempt to refine element #%d which has been refined already.", e->id);

  if (e->is_triangle()) {
    if (refinement == 3) refine_triangle_to_quads(e);
    else refine_triangle_to_triangles(e);
  }
  else
    refine_quad(e, refinement);

  seq = g_mesh_seq++;
}


void Mesh::refine_all_elements(int refinement)
{
  Element* e;
  elements.set_append_only(true);
  for_all_active_elements(e, this)
    refine_element(e->id, refinement);
  elements.set_append_only(false);
}


void Mesh::refine_by_criterion(int (*criterion)(Element*), int depth)
{
  Element* e;
  elements.set_append_only(true);
  for (int r, i = 0; i < depth; i++) {
    for_all_active_elements(e, this) {
      if ((r = criterion(e)) >= 0) refine_element(e->id, r);
    }
  }
  elements.set_append_only(false);
}

static int rtv_id;

static int rtv_criterion(Element* e)
{
  for (unsigned int i = 0; i < e->nvert; i++)
    if (e->vn[i]->id == rtv_id)
      return 0;
  return -1;
}

void Mesh::refine_towards_vertex(int vertex_id, int depth)
{
  rtv_id = vertex_id;
  refine_by_criterion(rtv_criterion, depth);
}

static int rtb_marker;
static bool rtb_aniso;
static bool rtb_tria_to_quad;
static char* rtb_vert;

static int rtb_criterion(Element* e)
{
  unsigned int i;
  for (i = 0; i < e->nvert; i++) {
    if (e->en[i]->marker == rtb_marker || rtb_vert[e->vn[i]->id]) {
      break;
    }
  }

  if (i >= e->nvert) return -1;
  // triangle should be split into 3 quads
  if (e->is_triangle() && rtb_tria_to_quad) return 3;
  // triangle should be split into 4 triangles or quad should
  // be split into 4 quads
  if (e->is_triangle() || !rtb_aniso) return 0;

  // quads - anisotropic case 1
  if ((e->en[0]->marker == rtb_marker && !rtb_vert[e->vn[2]->id] && !rtb_vert[e->vn[3]->id]) ||
      (e->en[2]->marker == rtb_marker && !rtb_vert[e->vn[0]->id] && !rtb_vert[e->vn[1]->id]) ||
      (e->en[0]->marker == rtb_marker && e->en[2]->marker == rtb_marker &&
       e->en[1]->marker != rtb_marker && e->en[3]->marker != rtb_marker)) return 1;

  // quads - anisotropic case 2
  if ((e->en[1]->marker == rtb_marker && !rtb_vert[e->vn[3]->id] && !rtb_vert[e->vn[0]->id]) ||
      (e->en[3]->marker == rtb_marker && !rtb_vert[e->vn[1]->id] && !rtb_vert[e->vn[2]->id]) ||
      (e->en[1]->marker == rtb_marker && e->en[3]->marker == rtb_marker &&
       e->en[0]->marker != rtb_marker && e->en[2]->marker != rtb_marker)) return 2;

  return 0;
}

void Mesh::refine_towards_boundary(int marker, int depth, bool aniso, bool tria_to_quad)
{
  rtb_marker = marker;
  rtb_aniso  = aniso;
  rtb_tria_to_quad = tria_to_quad;

  for (int i = 0; i < depth; i++)
  {
    int size = get_max_node_id()+1;
    rtb_vert = new char[size];
    memset(rtb_vert, 0, sizeof(char) * size);

    Element* e;
    for_all_active_elements(e, this)
      for (unsigned int j = 0; j < e->nvert; j++) {
        if (e->en[j]->marker == marker) {
          rtb_vert[e->vn[j]->id] = rtb_vert[e->vn[e->next_vert(j)]->id] = 1;
        }
      }

    refine_by_criterion(rtb_criterion, 1);
    delete [] rtb_vert;
  }
}


void Mesh::refine_towards_boundary(std::string marker, int depth, bool aniso)
{
  if(this->markers_conversion == NULL)
        error("MarkersConversion class has to be used if string boundary/area markers are to be used.");
  this->refine_towards_boundary(this->markers_conversion->get_internal_boundary_marker(marker), depth, aniso);
}


void Mesh::unrefine_element(int id)
{
  Element* e = get_element(id);
  if (!e->used) error("Invalid element id number.");
  if (e->active) return;

  for (int i = 0; i < 4; i++)
    if (e->sons[i] != NULL)
      unrefine_element(e->sons[i]->id);

  unrefine_element_internal(e);
  seq = g_mesh_seq++;
}


void Mesh::unrefine_all_elements(bool keep_initial_refinements)
{
  // find inactive elements with active sons
  std::vector<int> list;
  Element* e;
  for_all_inactive_elements(e, this)
  {
    bool found = true;
    for (unsigned int i = 0; i < 4; i++)
      if (e->sons[i] != NULL && (!e->sons[i]->active ||
          (keep_initial_refinements && e->sons[i]->id < ninitial))  )
        { found = false; break; }

    if (found) list.push_back(e->id);
  }

  // unrefine the found elements
  for (unsigned int i = 0; i < list.size(); i++)
    unrefine_element(list[i]);
}

/// Returns a NURBS curve with reversed control points and inverted knot vector.
/// Used for curved edges inside a mesh, where two mirror Nurbs have to be created
/// for the adjacent elements
///
Nurbs* Mesh::reverse_nurbs(Nurbs* nurbs)
{
  Nurbs* rev = new Nurbs;
  *rev = *nurbs;
  rev->twin = true;

  rev->pt = new double3[nurbs->np];
  for (int i = 0; i < nurbs->np; i++)
  {
    rev->pt[nurbs->np-1 - i][0] = nurbs->pt[i][0];
    rev->pt[nurbs->np-1 - i][1] = nurbs->pt[i][1];
    rev->pt[nurbs->np-1 - i][2] = nurbs->pt[i][2];
  }

  rev->kv = new double[nurbs->nk];
  for (int i = 0; i < nurbs->nk; i++)
    rev->kv[i] = nurbs->kv[i];
  for (int i = nurbs->degree + 1; i < nurbs->nk - nurbs->degree - 1; i++)
    rev->kv[nurbs->nk-1 - i] = 1.0 - nurbs->kv[i];

  rev->arc = nurbs->arc;
  rev->angle = -nurbs->angle;
  return rev;
}

// computing vector length
double vector_length(double a_1, double a_2)
{
  return sqrt(sqr(a_1) + sqr(a_2));
}

// checking whether the points p, q, r lie on the same line
bool same_line(double p_1, double p_2, double q_1, double q_2, double r_1, double r_2)
{
  double pq_1 = q_1 - p_1, pq_2 = q_2 - p_2, pr_1 = r_1 - p_1, pr_2 = r_2 - p_2;
  double length_pq = vector_length(pq_1, pq_2);
  double length_pr = vector_length(pr_1, pr_2);
  double sin_angle = (pq_1*pr_2 - pq_2*pr_1)/(length_pq*length_pr);
  if(fabs(sin_angle) < 1e-8) return true;
  else return false;
}

// checking whether the angle of vectors 'a' and 'b' is between zero and Pi
bool is_convex(double a_1, double a_2, double b_1, double b_2)
{
  if(a_1*b_2 - a_2*b_1 > 0) return true;
  else return false;
}

void check_triangle(int i, Node *&v0, Node *&v1, Node *&v2)
{
  // checking that all edges have nonzero length
  double
    length_1 = vector_length(v1->x - v0->x, v1->y - v0->y),
    length_2 = vector_length(v2->x - v1->x, v2->y - v1->y),
    length_3 = vector_length(v0->x - v2->x, v0->y - v2->y);
  if(length_1 < 1e-14 || length_2 < 1e-14 || length_3 < 1e-14)
    error("Edge of triangular element #%d has length less than 1e-14.", i);

  // checking that vertices do not lie on the same line
  if(same_line(v0->x, v0->y, v1->x, v1->y, v2->x, v2->y))
    error("Triangular element #%d: all vertices lie on the same line.", i);

  // checking positive orientation. If not positive, swapping vertices
  if (!is_convex(v1->x - v0->x, v1->y - v0->y, v2->x - v0->x, v2->y - v0->y)) {
    warn("Triangular element #%d not positively oriented, swapping vertices.", i);
    std::swap(v1, v2);
  }
}

void check_quad(int i, Node *&v0, Node *&v1, Node *&v2, Node *&v3)
{
  // checking that all edges have nonzero length
  double
    length_1 = vector_length(v1->x - v0->x, v1->y - v0->y),
    length_2 = vector_length(v2->x - v1->x, v2->y - v1->y),
    length_3 = vector_length(v3->x - v2->x, v3->y - v2->y),
    length_4 = vector_length(v0->x - v3->x, v0->y - v3->y);
  if(length_1 < 1e-14 || length_2 < 1e-14 || length_3 < 1e-14 || length_4 < 1e-14)
    error("Edge of quad element #%d has length less than 1e-14.", i);

  // checking that both diagonals have nonzero length
  double
    diag_1 = vector_length(v2->x - v0->x, v2->y - v0->y),
    diag_2 = vector_length(v3->x - v1->x, v3->y - v1->y);
  if(diag_1 < 1e-14 || diag_2 < 1e-14)
    error("Diagonal of quad element #%d has length less than 1e-14.", i);

  // checking that vertices v0, v1, v2 do not lie on the same line
  if(same_line(v0->x, v0->y, v1->x, v1->y, v2->x, v2->y))
    error("Quad element #%d: vertices v0, v1, v2 lie on the same line.", i);
  // checking that vertices v0, v1, v3 do not lie on the same line
  if(same_line(v0->x, v0->y, v1->x, v1->y, v3->x, v3->y))
    error("Quad element #%d: vertices v0, v1, v3 lie on the same line.", i);
  // checking that vertices v0, v2, v3 do not lie on the same line
  if(same_line(v0->x, v0->y, v2->x, v2->y, v3->x, v3->y))
    error("Quad element #%d: vertices v0, v2, v3 lie on the same line.", i);
  // checking that vertices v1, v2, v3 do not lie on the same line
  if(same_line(v1->x, v1->y, v2->x, v2->y, v3->x, v3->y))
    error("Quad element #%d: vertices v1, v2, v3 lie on the same line.", i);

  // checking that vertex v1 lies on the right of the diagonal v2-v0
  int vertex_1_ok = is_convex(v1->x - v0->x, v1->y - v0->y, v2->x - v0->x, v2->y - v0->y);
  if(!vertex_1_ok) error("Vertex v1 of quad element #%d does not lie on the right of the diagonal v2-v0.", i);
  // checking that vertex v3 lies on the left of the diagonal v2-v0
  int vertex_3_ok = is_convex(v2->x - v0->x, v2->y - v0->y, v3->x - v0->x, v3->y - v0->y);
  if(!vertex_3_ok) error("Vertex v3 of quad element #%d does not lie on the left of the diagonal v2-v0.", i);
  // checking that vertex v2 lies on the right of the diagonal v3-v1
  int vertex_2_ok = is_convex(v2->x - v1->x, v2->y - v1->y, v3->x - v1->x, v3->y - v1->y);
  if(!vertex_2_ok) error("Vertex v2 of quad element #%d does not lie on the right of the diagonal v3-v1.", i);
  // checking that vertex v0 lies on the left of the diagonal v3-v1
  int vertex_0_ok = is_convex(v3->x - v1->x, v3->y - v1->y, v0->x - v1->x, v0->y - v1->y);
  if(!vertex_0_ok) error("Vertex v0 of quad element #%d does not lie on the left of the diagonal v2-v1.", i);
}


//// mesh::create //////////////////////////////////////////////////////////////////////////////////

void Mesh::create(int nv, double2* verts, int nt, int4* tris,
                  int nq, int5* quads, int nm, int3* mark)
{
  //printf("Calling Mesh::free() in Mesh::create().\n");
  free();

  // initialize hash table
  int size = 16;
  while (size < 2*nv) size *= 2;
  HashTable::init(size);

  // create vertex nodes
  for (int i = 0; i < nv; i++)
  {
    Node* node = nodes.add();
    assert(node->id == i);
    node->ref = TOP_LEVEL_REF;
    node->type = HERMES_TYPE_VERTEX;
    node->bnd = 0;
    node->p1 = node->p2 = -1;
    node->next_hash = NULL;
    node->x = verts[i][0];
    node->y = verts[i][1];
  }
  ntopvert = nv;

  // create triangles
  for (int i = 0; i < nt; i++)
    create_triangle(tris[i][3], &nodes[tris[i][0]], &nodes[tris[i][1]], &nodes[tris[i][2]], NULL);

  // create quads
  for (int i = 0; i < nq; i++)
    create_quad(quads[i][4], &nodes[quads[i][0]], &nodes[quads[i][1]], &nodes[quads[i][2]], &nodes[quads[i][3]], NULL);

  // set boundary markers
  for (int i = 0; i < nm; i++)
  {
    Node* en = peek_edge_node(mark[i][0], mark[i][1]);
    if (en == NULL) error("Boundary data error (edge does not exist)");
    en->marker = mark[i][2];

    if (en->marker > 0)
    {
      nodes[mark[i][0]].bnd = 1;
      nodes[mark[i][1]].bnd = 1;
      en->bnd = 1;
    }
  }

  nbase = nactive = ninitial = nt + nq;
  seq = g_mesh_seq++;
}

bool Mesh::rescale(double x_ref, double y_ref) 
{
  // Sanity checks.
  if (fabs(x_ref) < 1e-10) return false;
  if (fabs(y_ref) < 1e-10) return false;

  // If curvilinear, the mesh cannot be rescaled.
  bool curved = false;
  Element* e;
  for_all_elements(e, this) {
    if (e->cm != NULL) {
      curved = true;
      break;
    }
  }
  if (curved == true) return false;

  // Go through all vertices and rescale coordinates.
  Node* n;
  for_all_vertex_nodes(n, this) {
    n->x /= x_ref;
    n->y /= y_ref;
  }

  return true;
}

//// mesh copy /////////////////////////////////////////////////////////////////////////////////////

void Mesh::copy(const Mesh* mesh)
{
  unsigned int i;

  //printf("Calling Mesh::free() in Mesh::copy().\n");
  free();

  // copy nodes and elements
  HashTable::copy(mesh);
  elements.copy(mesh->elements);

  Element* e;
  for_all_elements(e, this)
  {
    // update vertex node pointers
    for (i = 0; i < e->nvert; i++)
      e->vn[i] = &nodes[e->vn[i]->id];

    if (e->active)
    {
      // update edge node pointers
      for (i = 0; i < e->nvert; i++)
        e->en[i] = &nodes[e->en[i]->id];
    }
    else
    {
      // update son pointers
      for (i = 0; i < 4; i++)
        if (e->sons[i] != NULL)
          e->sons[i] = &elements[e->sons[i]->id];
    }

    // copy CurvMap, update its parent
    if (e->cm != NULL)
    {
      e->cm = new CurvMap(e->cm);
      if (!e->cm->toplevel)
        e->cm->parent = &elements[e->cm->parent->id];
    }

    //update parent pointer
    if(e->parent != NULL)
      e->parent = &elements[e->parent->id];
  }

  // update element pointers in edge nodes
  Node* node;
  for_all_edge_nodes(node, this)
    for (i = 0; i < 2; i++)
      if (node->elem[i] != NULL)
        node->elem[i] = &elements[node->elem[i]->id];

  nbase = mesh->nbase;
  nactive = mesh->nactive;
  ntopvert = mesh->ntopvert;
  ninitial = mesh->ninitial;
  seq = mesh->seq;
  markers_conversion = mesh->markers_conversion;
}


Node* Mesh::get_base_edge_node(Element* base, int edge)
{
  while (!base->active) // we need to go down to an active element
  {
    int son1, son2;
    get_edge_sons(base, edge, son1, son2);
    base = base->sons[son1];
  }
  return base->en[edge];
}


void Mesh::copy_base(Mesh* mesh)
{
  //printf("Calling Mesh::free() in Mesh::copy_base().\n");
  free();
  HashTable::init();

  // copy top-level vertex nodes
  for (int i = 0; i < mesh->get_max_node_id(); i++)
  {
    Node* node = &(mesh->nodes[i]);
    if (node->ref < TOP_LEVEL_REF) break;
    Node* newnode = nodes.add();
    assert(newnode->id == i && node->type == HERMES_TYPE_VERTEX);
    memcpy(newnode, node, sizeof(Node));
    newnode->ref = TOP_LEVEL_REF;
  }

  // copy base elements
  Element* e;
  for_all_base_elements(e, mesh)
  {
    Element* enew;
    Node *v0 = &nodes[e->vn[0]->id], *v1 = &nodes[e->vn[1]->id], *v2 = &nodes[e->vn[2]->id];
    if (e->is_triangle())
      enew = create_triangle(e->marker, v0, v1, v2, NULL);
    else
      enew = create_quad(e->marker, v0, v1, v2, &nodes[e->vn[3]->id], NULL);

    // copy edge markers
    for (unsigned int j = 0; j < e->nvert; j++)
    {
      Node* en = get_base_edge_node(e, j);
      enew->en[j]->bnd = en->bnd; // copy bnd data from the active el.
      enew->en[j]->marker = en->marker;
    }

    enew->userdata = e->userdata;
    if (e->is_curved())
      enew->cm = new CurvMap(e->cm);
  }

  nbase = nactive = ninitial = mesh->nbase;
  ntopvert = mesh->ntopvert;
  seq = g_mesh_seq++;
}


//// free //////////////////////////////////////////////////////////////////////////////////////////

void Mesh::free()
{
  Element* e;
  for_all_elements(e, this)
    if (e->cm != NULL)
    {
      delete e->cm;
      e->cm = NULL; // fixme!!!
    }

  elements.free();
  HashTable::free();
}

void Mesh::copy_converted(Mesh* mesh)
{
  //printf("Calling Mesh::free() in Mesh::copy_converted().\n");
  free();
  HashTable::copy(mesh);
  // clear refernce for all nodes
  for(int i = 0; i < nodes.get_size(); i++)
  {
    Node& node = nodes[i];
    if (node.type == HERMES_TYPE_EDGE) { //process only edge nodes
      for(int k = 0; k < 2; k++)
        node.elem[k] = NULL;
    }
  }

  // copy active elements
  Element* e;

  for_all_active_elements(e, mesh)
  {
    Element* enew;
    Node *v0 = &nodes[e->vn[0]->id], *v1 = &nodes[e->vn[1]->id], *v2 = &nodes[e->vn[2]->id];
    Node *e0 = &nodes[e->en[0]->id], *e1 = &nodes[e->en[1]->id], *e2 = &nodes[e->en[2]->id];
    if (e->is_triangle())
    {
      // create a new element
      enew = elements.add();
      enew->active = 1;
      enew->marker = e->marker;
      enew->userdata = 0;
      enew->nvert = 3;
      enew->iro_cache = -1;
      enew->cm = e->cm;

      // set vertex and edge node pointers
      enew->vn[0] = v0;
      enew->vn[1] = v1;
      enew->vn[2] = v2;
      enew->en[0] = e0;
      enew->en[1] = e1;
      enew->en[2] = e2;
    }
    else
    {
      // create a new element
      Node *v3 = &nodes[e->vn[3]->id];
      Node *e3 = &nodes[e->en[3]->id];
      enew = elements.add();
      enew->active = 1;
      enew->marker = e->marker;
      enew->userdata = 0;
      enew->nvert = 4;
      enew->iro_cache = -1;
      enew->cm = e->cm;
      enew->parent = NULL;
      enew->visited = false;

      // set vertex and edge node pointers
      enew->vn[0] = v0;
      enew->vn[1] = v1;
      enew->vn[2] = v2;
      enew->vn[3] = v3;
      enew->en[0] = e0;
      enew->en[1] = e1;
      enew->en[2] = e2;
      enew->en[3] = e3;
    }

    // copy edge markers
    for (unsigned int j = 0; j < e->nvert; j++)
    {
      Node* en = get_base_edge_node(e, j);
      enew->en[j]->bnd = en->bnd;
      enew->en[j]->marker = en->marker;
    }

    enew->userdata = e->userdata;
    if (e->is_curved())
      enew->cm = new CurvMap(e->cm);
  }

  nbase = nactive = ninitial = mesh->nactive;
  ntopvert = mesh->ntopvert = get_num_nodes();
  seq = g_mesh_seq++;
}

////convert a triangle element into three quadrilateral elements///////

void Mesh::convert_triangles_to_quads()
{
  Element* e;
  Mesh tmp;

  elements.set_append_only(true);
  for_all_active_elements(e, this)
    refine_element_to_quads(e->id);
  elements.set_append_only(false);

  tmp.copy_converted(this);
  for (int i = 0; i < tmp.ntopvert; i++)
  {
     if (tmp.nodes[i].type == 1)
     {
       tmp.nodes[i].y = 0.0;
     }
  }
  char* mesh_tmp = NULL;
  mesh_tmp = tmpnam(NULL);
  tmp.save(mesh_tmp);
  Mesh mesh_tmp_for_convert;
  H2DReader loader_mesh_tmp_for_convert;
  loader_mesh_tmp_for_convert.load(mesh_tmp, &tmp);
  remove(mesh_tmp);
  copy(&tmp);
}

////convert a quad element into two triangle elements///////
void Mesh::convert_quads_to_triangles()
{
  Element* e;
  Mesh tmp;

  elements.set_append_only(true);
  for_all_active_elements(e, this)
    refine_element_to_triangles(e->id);
  elements.set_append_only(false);

  tmp.copy_converted(this);
  for (int i = 0; i < tmp.ntopvert; i++)
  {
     if (tmp.nodes[i].type == 1)
     {
       tmp.nodes[i].y = 0.0;
     }
  }
  char* mesh_tmp = NULL;
  mesh_tmp = tmpnam(NULL);
  tmp.save(mesh_tmp);
  Mesh mesh_tmp_for_convert;
  H2DReader loader_mesh_tmp_for_convert;
  loader_mesh_tmp_for_convert.load(mesh_tmp, &tmp);
  remove(mesh_tmp);
  copy(&tmp);
}

void Mesh::refine_triangle_to_quads(Element* e)
{
  // remember the markers of the edge nodes
  int bnd[3] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd    };
  int mrk[3] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker };

  // deactivate this element and unregister from its nodes
  e->active = false;
  nactive--;
  e->unref_all_nodes(this);

  // obtain three mid-edge and one gravity vertex nodes
  Node* x0 = get_vertex_node(e->vn[0]->id, e->vn[1]->id);
  Node* x1 = get_vertex_node(e->vn[1]->id, e->vn[2]->id);
  Node* x2 = get_vertex_node(e->vn[2]->id, e->vn[0]->id);
  Node* mid = get_vertex_node(x0->id, e->vn[1]->id);

  mid->x = (nodes[x0->id].x + nodes[x1->id].x + nodes[x2->id].x)/3;
  mid->y = (nodes[x0->id].y + nodes[x1->id].y + nodes[x2->id].y)/3;

  // adjust mid-edge and gravity coordinates if this is a curved element
  if (e->is_curved())
  {
    double2 pt[4] = { { 0.0,-1.0 }, { 0.0, 0.0 },{ -1.0, 0.0 }, { -0.33333333, -0.33333333 } };
    e->cm->get_mid_edge_points(e, pt, 4);
    x0->x = pt[0][0]; x0->y = pt[0][1];
    x1->x = pt[1][0]; x1->y = pt[1][1];
    x2->x = pt[2][0]; x2->y = pt[2][1];
    mid->x = pt[3][0]; mid->y = pt[3][1];
   }

  double angle2;
  int idx;
  CurvMap* cm[3];
  memset(cm, 0, sizeof(cm));

  // create CurvMaps for sons if this is a curved element
  if (e->is_curved())
  {
    for (idx = 0; idx < 2; idx++)
    {
      if (e->cm->nurbs[idx] != NULL)
      {
        cm[idx] = new CurvMap;
        memset(cm[idx], 0, sizeof(CurvMap));
        cm[idx+1] = new CurvMap;
        memset(cm[idx+1], 0, sizeof(CurvMap));
      }
    }

    idx=0;
    if (e->cm->nurbs[idx] != NULL)
    {
      angle2 = e->cm->nurbs[idx]->angle/2;
      Node* node_temp = get_vertex_node(e->vn[idx%3]->id, e->vn[(idx+1)%3]->id);

      for (int k = 0; k < 2; k++)
      {
        int p1, p2;
        int idx2;

        if (k == 0)
        {
          p1 = e->vn[(idx)%3]->id;
          p2 = node_temp->id;
          if (idx == 0) idx2 = 0;
          if (idx == 1) idx2 = 1;
          if (idx == 2) continue;
        }
        else if (k == 1)
        {
          p1 = node_temp->id;
          p2 = e->vn[(idx+1)%3]->id;
          idx = (idx+1)%3;
          if (idx == 0) continue;
          if (idx == 1) idx2 = 0;
          if (idx == 2) idx2 = 0;
        }

        Nurbs* nurbs = new Nurbs;
        bool cricle = true;

        nurbs->arc = cricle;
        nurbs->degree = 2;

        int inner = 1, outer;
        inner = 1;
        nurbs->np = inner + 2;
        nurbs->pt = new double3[nurbs->np];

        nurbs->pt[0][0] = nodes[p1].x;
        nurbs->pt[0][1] = nodes[p1].y;
        nurbs->pt[0][2] = 1.0;

        nurbs->pt[inner+1][0] = nodes[p2].x;
        nurbs->pt[inner+1][1] = nodes[p2].y;
        nurbs->pt[inner+1][2] = 1.0;

        double angle = angle2;
        double a = (180.0 - angle) / 180.0 * M_PI;
        nurbs->angle = angle;

        // generate one control point
        double x = 1.0 / tan(a * 0.5);
        nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
        nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
        nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

        int i;
        inner = 0;
        nurbs->nk = nurbs->degree + nurbs->np + 1;
        outer = nurbs->nk - inner;

        // knot vector is completed by 0.0 on the left and by 1.0 on the right
        nurbs->kv = new double[nurbs->nk];

        for (i = 0; i < outer/2; i++)
          nurbs->kv[i] = 0.0;
        for (i = outer/2 + inner; i < nurbs->nk; i++)
          nurbs->kv[i] = 1.0;
        nurbs->ref = 0;

        cm[idx]->toplevel = 1;
        cm[idx]->order = 4;
        cm[idx]->nurbs[idx2] = nurbs;
        nurbs->ref++;
      }
    }

    idx = 1;
    if (e->cm->nurbs[idx] != NULL)
    {
      angle2 = e->cm->nurbs[idx]->angle/2;
      Node* node_temp = get_vertex_node(e->vn[idx%3]->id, e->vn[(idx+1)%3]->id);
      for (int k = 0; k < 2; k++)
      {
        int p1, p2;
        int idx2;
        if (k == 0)
        {
          p1 = e->vn[(idx)%3]->id;
          p2 = node_temp->id;
          if (idx == 0) idx2 = 0;
          if (idx == 1) idx2 = 1;
          if (idx == 2) continue;
        }
        else if (k == 1)
        {
          p1 = node_temp->id;
          p2 = e->vn[(idx+1)%3]->id;
          idx = (idx+1)%3;
          if (idx == 0) continue;
          if (idx == 1) idx2 = 0;
          if (idx == 2) idx2 = 0;
        }

        Nurbs* nurbs = new Nurbs;
        bool cricle = true;

        nurbs->arc = cricle;
        nurbs->degree = 2;
        int inner = 1, outer;
        inner = 1;
        nurbs->np = inner + 2;
        nurbs->pt = new double3[nurbs->np];

        nurbs->pt[0][0] = nodes[p1].x;
        nurbs->pt[0][1] = nodes[p1].y;
        nurbs->pt[0][2] = 1.0;

        nurbs->pt[inner+1][0] = nodes[p2].x;
        nurbs->pt[inner+1][1] = nodes[p2].y;
        nurbs->pt[inner+1][2] = 1.0;

        double angle = angle2;
        double a = (180.0 - angle) / 180.0 * M_PI;
        nurbs->angle = angle;

        // generate one control point
        double x = 1.0 / tan(a * 0.5);
        nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
        nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
        nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

        int i;
        inner = 0;
        nurbs->nk = nurbs->degree + nurbs->np + 1;
        outer = nurbs->nk - inner;

        // knot vector is completed by 0.0 on the left and by 1.0 on the right
        nurbs->kv = new double[nurbs->nk];
        for (i = 0; i < outer/2; i++)
          nurbs->kv[i] = 0.0;
        for (i = outer/2 + inner; i < nurbs->nk; i++)
         nurbs->kv[i] = 1.0;
        nurbs->ref = 0;

        cm[idx]->toplevel = 1;
        cm[idx]->order = 4;
        cm[idx]->nurbs[idx2] = nurbs;
        nurbs->ref++;
      }
    }
  }

  // create the four sons
  Element* sons[4];
  sons[0] = create_quad(e->marker, e->vn[0], x0, mid, x2, cm[0]);
  sons[1] = create_quad(e->marker, x0, e->vn[1], x1, mid, cm[1]);
  sons[2] = create_quad(e->marker, x1, e->vn[2], x2, mid, cm[2]);
  sons[3] = NULL;

  // update coefficients of curved reference mapping
  for (int i = 0; i < 3; i++)
  {
    if (sons[i]->is_curved())
    {
      sons[i]->cm->update_refmap_coeffs(sons[i]);
    }
  }
  nactive += 3;
  // now the original edge nodes may no longer exist...
  // set correct boundary status and markers for the new nodes
  sons[0]->en[0]->bnd = bnd[0];  sons[0]->en[0]->marker = mrk[0];
  sons[0]->en[3]->bnd = bnd[2];  sons[0]->en[3]->marker = mrk[2];
  sons[0]->vn[1]->bnd = bnd[0];

  sons[1]->en[0]->bnd = bnd[0];  sons[1]->en[0]->marker = mrk[0];
  sons[1]->en[1]->bnd = bnd[1];  sons[1]->en[1]->marker = mrk[1];
  sons[1]->vn[2]->bnd = bnd[1];

  sons[2]->en[0]->bnd = bnd[1];  sons[2]->en[0]->marker = mrk[1];
  sons[2]->en[1]->bnd = bnd[2];  sons[2]->en[1]->marker = mrk[2];
  sons[2]->vn[2]->bnd = bnd[2];

  //set pointers to parent element for sons
  for(int i = 0; i < 4; i++)
	  if(sons[i] != NULL)
		  sons[i]->parent = e;

  // copy son pointers (could not have been done earlier because of the union)
  memcpy(e->sons, sons, 4 * sizeof(Element*));
}


void Mesh::refine_element_to_quads(int id)
{
  Element* e = get_element(id);
  if (!e->used) error("Invalid element id number.");
  if (!e->active) error("Attempt to refine element #%d which has been refined already.", e->id);

  if (e->is_triangle())
    refine_triangle_to_quads(e);
  else
    return;

  seq = g_mesh_seq++;
}


void Mesh::refine_quad_to_triangles(Element* e)
{
  // remember the markers of the edge nodes
  int bnd[4] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd,    e->en[3]->bnd };
  int mrk[4] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker, e->en[3]->marker };

  // deactivate this element and unregister from its nodes
  e->active = false;
  nactive--;
  e->unref_all_nodes(this);

  bool bcheck = true;  ///< if bcheck is true, it is default add a new edge between
                       ///<  vn[0] and vn[2]
  double length_x_0_2 = (e->vn[0]->x - e->vn[2]->x)*(e->vn[0]->x - e->vn[2]->x);
  double length_x_1_3 = (e->vn[1]->x - e->vn[3]->x)*(e->vn[1]->x - e->vn[3]->x);

  double length_y_0_2 = (e->vn[0]->y - e->vn[2]->y)*(e->vn[0]->y - e->vn[2]->y);
  double length_y_1_3 = (e->vn[1]->y - e->vn[3]->y)*(e->vn[1]->y - e->vn[3]->y);

  if ((length_x_0_2 + length_y_0_2) > (length_x_1_3 + length_y_1_3))
  {
    bcheck = false;
  }

  double angle2;
  unsigned int idx;
  CurvMap* cm[2];
  memset(cm, 0, sizeof(cm));

  // create CurvMaps for sons if this is a curved element
  if (e->is_curved())
  {
    int i_case2 = 0;
    if (bcheck == true)
    {
      if ((e->cm->nurbs[0] != NULL) || (e->cm->nurbs[1] != NULL))
      {
        cm[0] = new CurvMap;
        memset(cm[0], 0, sizeof(CurvMap));
      }
      if ((e->cm->nurbs[2] != NULL) || (e->cm->nurbs[3] != NULL))
      {
        cm[1] = new CurvMap;
        memset(cm[1], 0, sizeof(CurvMap));
      }
    }
    else if (bcheck == false)
    {
      if ((e->cm->nurbs[1] != NULL) || (e->cm->nurbs[2] != NULL))
      {
        cm[0] = new CurvMap;
        memset(cm[0], 0, sizeof(CurvMap));
      }
      if ((e->cm->nurbs[3] != NULL) || (e->cm->nurbs[0] != NULL))
      {
        cm[1] = new CurvMap;
        memset(cm[1], 0, sizeof(CurvMap));
      }
      i_case2 = 1; //switch to the shorter diagonal
    }

    for (unsigned int k = 0; k < 2; k++)
    {
      for (idx = 0 + 2*k; idx < 2 + 2*k; idx++)
      {
        if (e->cm->nurbs[(idx + i_case2)%4] != NULL)
        {
          angle2 = e->cm->nurbs[(idx + i_case2)%4]->angle;

          int p1, p2;
          unsigned int idx2 = idx;

          p1 = e->vn[(idx + i_case2)%4]->id;
          p2 = e->vn[(idx + i_case2 + 1)%4]->id;  //node_temp->id;

          Nurbs* nurbs = new Nurbs;
          bool cricle = true;

          nurbs->arc = cricle;
          nurbs->degree = 2;

          int inner = 1, outer;
          inner = 1;
          nurbs->np = inner + 2;
          nurbs->pt = new double3[nurbs->np];

          nurbs->pt[0][0] = nodes[p1].x;
          nurbs->pt[0][1] = nodes[p1].y;
          nurbs->pt[0][2] = 1.0;

          nurbs->pt[inner+1][0] = nodes[p2].x;
          nurbs->pt[inner+1][1] = nodes[p2].y;
          nurbs->pt[inner+1][2] = 1.0;

          double angle = angle2;
          double a = (180.0 - angle) / 180.0 * M_PI;
          nurbs->angle = angle;

          // generate one control point
          double x = 1.0 / tan(a * 0.5);
          nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
          nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
          nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

          int i;
          inner = 0;
          nurbs->nk = nurbs->degree + nurbs->np + 1;
          outer = nurbs->nk - inner;

          // knot vector is completed by 0.0 on the left and by 1.0 on the right
          nurbs->kv = new double[nurbs->nk];

          for (i = 0; i < outer/2; i++)
            nurbs->kv[i] = 0.0;
          for (i = outer/2 + inner; i < nurbs->nk; i++)
            nurbs->kv[i] = 1.0;
          nurbs->ref = 0;

          cm[k]->toplevel = 1;
          cm[k]->order = 4;
          cm[k]->nurbs[idx%2] = nurbs;
          nurbs->ref++;
        }
      }
    }
  }

  // create the four sons
  Element* sons[4];
  if (bcheck == true)
  {
    sons[0] = create_triangle(e->marker, e->vn[0], e->vn[1], e->vn[2], cm[0]);
    sons[1] = create_triangle(e->marker, e->vn[2], e->vn[3], e->vn[0], cm[1]);
    sons[2] = NULL;
    sons[3] = NULL;
  }
  else
  {
    sons[0] = create_triangle(e->marker, e->vn[1], e->vn[2], e->vn[3], cm[0]);
    sons[1] = create_triangle(e->marker, e->vn[3], e->vn[0], e->vn[1], cm[1]);
    sons[2] = NULL;
    sons[3] = NULL;
  }

  // update coefficients of curved reference mapping
  for (int i = 0; i < 2; i++)
  {
    if (sons[i]->is_curved())
    {
      sons[i]->cm->update_refmap_coeffs(sons[i]);
    }
  }
  nactive += 2;
  // now the original edge nodes may no longer exist...
  // set correct boundary status and markers for the new nodes
  if (bcheck == true)
  {
    sons[0]->en[0]->bnd = bnd[0];  sons[0]->en[0]->marker = mrk[0];
    sons[0]->en[1]->bnd = bnd[1];  sons[0]->en[1]->marker = mrk[1];
    sons[0]->vn[1]->bnd = bnd[0];

    sons[1]->en[0]->bnd = bnd[2];  sons[1]->en[0]->marker = mrk[2];
    sons[1]->en[1]->bnd = bnd[3];  sons[1]->en[1]->marker = mrk[3];
    sons[1]->vn[2]->bnd = bnd[1];
  }
  else
  {
    sons[0]->en[0]->bnd = bnd[1];  sons[0]->en[0]->marker = mrk[1];
    sons[0]->en[1]->bnd = bnd[2];  sons[0]->en[1]->marker = mrk[2];
    sons[0]->vn[1]->bnd = bnd[1];

    sons[1]->en[0]->bnd = bnd[3];  sons[1]->en[0]->marker = mrk[3];
    sons[1]->en[1]->bnd = bnd[0];  sons[1]->en[1]->marker = mrk[0];
    sons[1]->vn[2]->bnd = bnd[0];
  }

  //set pointers to parent element for sons
  for(int i = 0; i < 4; i++)
	  if(sons[i] != NULL)
		  sons[i]->parent = e;

  // copy son pointers (could not have been done earlier because of the union)
  memcpy(e->sons, sons, 4 * sizeof(Element*));
}


void Mesh::refine_element_to_triangles(int id)
{
  Element* e = get_element(id);
  if (!e->used) error("Invalid element id number.");
  if (!e->active) error("Attempt to refine element #%d which has been refined already.", e->id);

  if (e->is_triangle())
    return;
  else
    refine_quad_to_triangles(e);

  seq = g_mesh_seq++;
}

void Mesh::load(const char* filename, bool debug)
{
	warn("Deprecated function used. Please update your code to use MeshLoader classes.");

	H2DReader loader;
	loader.load(filename, this);
}


void Mesh::save(const char* filename)
{
	warn("Deprecated function used. Please update your code to use MeshLoader classes.");

	H2DReader loader;
  loader.save(filename, this);
}

//// save_raw, load_raw ////////////////////////////////////////////////////////////////////////////

void Mesh::save_raw(FILE* f)
{
  int i, nn, mm;
  int null = -1;

  assert(sizeof(int) == 4);
  assert(sizeof(double) == 8);

  hermes_fwrite("H2DM\001\000\000\000", 1, 8, f);

  #define output(n, type) \
    hermes_fwrite(&(n), sizeof(type), 1, f)

  output(nbase, int);
  output(ntopvert, int);
  output(nactive, int);

  nn = nodes.get_num_items();
  mm = nodes.get_size();
  output(nn, int);
  output(mm, int);

  // dump all nodes
  Node* n;
  for_all_nodes(n, this)
  {
    output(n->id, int);
    unsigned bits = n->ref | (n->type << 29) | (n->bnd << 30) | (n->used << 31);
    output(bits, unsigned);

    if (n->type == HERMES_TYPE_VERTEX)
    {
      output(n->x, double);
      output(n->y, double);
    }
    else
    {
      output(n->marker, int);
      output(n->elem[0] ? n->elem[0]->id : null, int);
      output(n->elem[1] ? n->elem[1]->id : null, int);
    }

    output(n->p1, int);
    output(n->p2, int);
  }

  nn = elements.get_num_items();
  mm = elements.get_size();
  output(nn, int);
  output(mm, int);

  // dump all elements
  Element* e;
  for (int id = 0; id < get_max_element_id(); id++)
  {
    if ((e = get_element_fast(id))->used || id < nbase)
    {
      output(e->id, int);
      unsigned bits = e->nvert | (e->active << 30) | (e->used << 31);
      output(bits, unsigned);

      if (e->used)
      {
        output(e->marker, int);
        output(e->userdata, int);
        output(e->iro_cache, int);

        for (unsigned i = 0; i < e->nvert; i++)
          output(e->vn[i]->id, int);

        if (e->active)
          for (unsigned i = 0; i < e->nvert; i++)
            output(e->en[i]->id, int);
        else
          for (unsigned i = 0; i < 4; i++)
            output(e->sons[i] ? e->sons[i]->id : null, int);

        if (e->is_curved()) error("Not implemented for curved elements yet.");
      }
    }
  }
  // TODO: curved elements

  #undef output
}


void Mesh::load_raw(FILE* f)
{
  int i, nv, mv, ne, me, id;

  assert(sizeof(int) == 4);
  assert(sizeof(double) == 8);

  // check header
  struct { char magic[4]; int ver; } hdr;
  hermes_fread(&hdr, sizeof(hdr), 1, f);
  if (hdr.magic[0] != 'H' || hdr.magic[1] != '2' || hdr.magic[2] != 'D' || hdr.magic[3] != 'M')
    error("Not a Hermes2D raw mesh file.");
  if (hdr.ver > 1)
    error("Unsupported file version.");

  #define input(n, type) \
    hermes_fread(&(n), sizeof(type), 1, f)

  //printf("Calling Mesh::free() in Mesh::load_raw().\n");
  free();

  input(nbase, int);
  input(ntopvert, int);
  input(nactive, int);

  input(nv, int);
  input(mv, int);
  nodes.force_size(mv);

  // load nodes
  for (i = 0; i < nv; i++)
  {
    input(id, int);
    if (id < 0 || id >= mv) error("Corrupt data.");
    Node* n = &(nodes[id]);
    n->id = id;
    n->used = 1;

    unsigned bits;
    input(bits, unsigned);
    n->ref  =  bits & 0x1fffffff;
    n->type = (bits >> 29) & 0x1;
    n->bnd  = (bits >> 30) & 0x1;

    if (n->type == HERMES_TYPE_VERTEX)
    {
      input(n->x, double);
      input(n->y, double);
    }
    else
    {
      input(n->marker, int);
      n->elem[0] = n->elem[1] = NULL;
      input(n->elem[0], int);
      input(n->elem[1], int);
    }

    input(n->p1, int);
    input(n->p2, int);
  }
  nodes.post_load_scan();

  int hsize = H2D_DEFAULT_HASH_SIZE;
  while (hsize < nv) hsize *= 2;
  HashTable::init(hsize);
  HashTable::rebuild();

  input(ne, int);
  input(me, int);
  elements.force_size(me);

  // load elements
  for (i = 0; i < ne; i++)
  {
    input(id, int);
    if (id < 0 || id >= me) error("Corrupt data.");
    Element* e = &(elements[id]);
    e->id = id;

    unsigned bits;
    input(bits, unsigned);
    e->nvert  =  bits & 0x3fffffff;
    e->active = (bits >> 30) & 0x1;
    e->used   = (bits >> 31) & 0x1;

    if (e->used)
    {
      input(e->marker, int);
      input(e->userdata, int);
      input(e->iro_cache, int);

      // load vertex node ids
      for (unsigned int j = 0; j < e->nvert; j++)
      {
        input(id, int);
        if (id < 0 || id >= mv) error("Corrupt data.");
        e->vn[j] = get_node(id);
      }

      if (e->active)
      {
        // load edge node ids
        for (unsigned int j = 0; j < e->nvert; j++)
        {
          input(id, int);
          if (id < 0 || id >= mv) error("Corrupt data.");
          e->en[j] = get_node(id);
        }
      }
      else
      {
        // load son ids
        for (unsigned int j = 0; j < 4; j++)
        {
          input(id, int);
          if (id < 0)
            e->sons[j] = NULL;
          else if (id < me)
            e->sons[j] = &(elements[id]);
          else
            error("Corrupt data.");
        }
      }
    }
  }
  elements.post_load_scan(nbase);

  // update edge node element pointers
  Node* n;
  for_all_edge_nodes(n, this)
    for (unsigned int j = 0; j < 2; j++)
      if ((int) (long) n->elem[j] == -1)
        n->elem[j] = NULL;
      else
        n->elem[j] = get_element((int) (long) n->elem[j]);

  #undef input
  seq++;
}


Mesh::MarkersConversion::MarkersConversion()
{
  conversion_table_for_element_markers = new std::map<int, std::string>;
  conversion_table_for_boundary_markers = new std::map<int, std::string>;
  conversion_table_for_element_markers_inverse = new std::map<std::string, int>;
  conversion_table_for_boundary_markers_inverse = new std::map<std::string, int>;

  min_boundary_marker_unused = 1;
  min_element_marker_unused = 0;
}

Mesh::MarkersConversion::~MarkersConversion()
{
  delete conversion_table_for_element_markers;
  delete conversion_table_for_boundary_markers;
  delete conversion_table_for_element_markers_inverse;
  delete conversion_table_for_boundary_markers_inverse;
}

void Mesh::MarkersConversion::insert_element_marker(int internal_marker, std::string user_marker)
{
    // First a check that the string value is not already present.
  if(user_marker != "")
    if(conversion_table_for_element_markers_inverse->find(user_marker) != conversion_table_for_element_markers_inverse->end())
      return;
  if(conversion_table_for_element_markers->size() == 0 || conversion_table_for_element_markers->find(internal_marker) == conversion_table_for_element_markers->end()) {
    conversion_table_for_element_markers->insert(std::pair<int, std::string>(internal_marker, user_marker));
    conversion_table_for_element_markers_inverse->insert(std::pair<std::string, int>(user_marker, internal_marker));
    if(user_marker != "")
      this->min_element_marker_unused++;
  }
  return;
}

void Mesh::MarkersConversion::insert_boundary_marker(int internal_marker, std::string user_marker)
{
  // First a check that the string value is not already present.
  if(user_marker != "")
    if(conversion_table_for_boundary_markers_inverse->find(user_marker) != conversion_table_for_boundary_markers_inverse->end())
      return;
  if(conversion_table_for_boundary_markers->size() == 0 || conversion_table_for_boundary_markers->find(internal_marker) == conversion_table_for_boundary_markers->end()) {
    conversion_table_for_boundary_markers->insert(std::pair<int, std::string>(internal_marker, user_marker));
    conversion_table_for_boundary_markers_inverse->insert(std::pair<std::string, int>(user_marker, internal_marker));
    if(user_marker != "")
      this->min_boundary_marker_unused++;
  }
  return;
}

std::string Mesh::MarkersConversion::get_user_element_marker(int internal_marker)
{
  if(conversion_table_for_element_markers->find(internal_marker) == conversion_table_for_element_markers->end())
    error("MarkersConversions class asked for a non existing marker %d", internal_marker);
  return conversion_table_for_element_markers->find(internal_marker)->second;
}
std::string Mesh::MarkersConversion::get_user_boundary_marker(int internal_marker)
{
  if(conversion_table_for_boundary_markers->find(internal_marker) == conversion_table_for_boundary_markers->end())
    error("MarkersConversions class asked for a non existing marker %d", internal_marker);
  return conversion_table_for_boundary_markers->find(internal_marker)->second;
}

int Mesh::MarkersConversion::get_internal_element_marker(std::string user_marker)
{
  if(conversion_table_for_element_markers_inverse->find(user_marker) == conversion_table_for_element_markers_inverse->end())
    error("MarkersConversions class asked for a non existing marker %s", user_marker.c_str());
  return conversion_table_for_element_markers_inverse->find(user_marker)->second;
}
int Mesh::MarkersConversion::get_internal_boundary_marker(std::string user_marker)
{
  if(conversion_table_for_boundary_markers_inverse->find(user_marker) == conversion_table_for_boundary_markers_inverse->end())
    error("MarkersConversions class asked for a non existing marker %s", user_marker.c_str());
  return conversion_table_for_boundary_markers_inverse->find(user_marker)->second;
}

void Mesh::MarkersConversion::check_boundary_marker(int marker)
{
  // This marker is already present.
  if(conversion_table_for_boundary_markers->find(marker) != conversion_table_for_boundary_markers->end())
    // If the present boundary marker is a different one (user-supplied string).
    if(conversion_table_for_boundary_markers->find(marker)->second != "")
    {
      // We have to reassign the markers.
      std::string temp_user_marker = conversion_table_for_boundary_markers->find(marker)->second;
      conversion_table_for_boundary_markers->erase(marker);
      conversion_table_for_boundary_markers_inverse->erase(temp_user_marker);
      insert_boundary_marker(this->min_boundary_marker_unused, temp_user_marker);
    }
  // Now we need to check if the variable min_boundary_marker_unused does not collide with
  // this marker, because we will need to raise the variable if it does as
  // we are about to insert this user-supplied integer marker.
  if(marker == min_boundary_marker_unused)
    min_boundary_marker_unused++;
}

void Mesh::MarkersConversion::check_element_marker(int marker)
{
  // This marker is already present.
  if(conversion_table_for_element_markers->find(marker) != conversion_table_for_element_markers->end())
    // If the present element marker is a different one (user-supplied string).
    if(conversion_table_for_element_markers->find(marker)->second != "")
    {
      // We have to reassign the markers.
      std::string temp_user_marker = conversion_table_for_element_markers->find(marker)->second;
      conversion_table_for_element_markers->erase(marker);
      conversion_table_for_element_markers_inverse->erase(temp_user_marker);
      insert_element_marker(this->min_element_marker_unused, temp_user_marker);
    }
  // Now we need to check if the variable min_element_marker_unused does not collide with
  // this marker, because we will need to raise the variable if it does as
  // we are about to insert this user-supplied integer marker.
  if(marker == min_element_marker_unused)
    min_element_marker_unused++;
}
