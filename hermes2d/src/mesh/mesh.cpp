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

#include "mesh_util.h"
#include "refmap.h"
#include <algorithm>
#include "global.h"
#include "api2d.h"
#include "mesh_reader_h2d.h"
#include "forms.h"

namespace Hermes
{
  namespace Hermes2D
  {
    unsigned g_mesh_seq = 0;
    static const int H2D_DG_INNER_EDGE_INT = -54125631;
    static const std::string H2D_DG_INNER_EDGE = "-54125631";

    Mesh::Mesh() : HashTable(), meshHashGrid(NULL), nbase(0), nactive(0), ntopvert(0), ninitial(0), seq(g_mesh_seq++),
      bounding_box_calculated(0)
    {
    }

    Mesh::~Mesh() 
    {
      free();
    }

    bool Mesh::isOkay() const
    {
      bool okay = true;

      if(this->elements.get_size() < 1)
        okay = false;
      if(this->nodes.get_size() < 1)
        okay = false;
      if(seq < 0)
        okay = false;
      return okay;
    }

    Mesh::ReferenceMeshCreator::ReferenceMeshCreator(MeshSharedPtr coarse_mesh, int refinement) : coarse_mesh(coarse_mesh), refinement(refinement)
    {
    }

    MeshSharedPtr Mesh::ReferenceMeshCreator::create_ref_mesh()
    {
      Mesh* ref_mesh = new Mesh;
      ref_mesh->copy(this->coarse_mesh);
      ref_mesh->refine_all_elements(refinement, false);
      return MeshSharedPtr(ref_mesh);
    }

    void Mesh::initial_single_check()
    {
      RefMap r;
      Element* e;
      Quad2D* quad = &g_quad_2d_std;
      for_all_active_elements(e, this)
      {
        r.set_active_element(e);

        int i, mo = quad->get_max_order(e->get_mode());

        int k = e->is_triangle() ? 2 : 3;

        double const_m[2][2] = 
        { 
          { e->vn[1]->x - e->vn[0]->x,  e->vn[k]->x - e->vn[0]->x },
          { e->vn[1]->y - e->vn[0]->y,  e->vn[k]->y - e->vn[0]->y } 
        };

        double const_jacobian = 0.25 * (const_m[0][0] * const_m[1][1] - const_m[0][1] * const_m[1][0]);
        if(const_jacobian <= 0.0)
          throw Hermes::Exceptions::MeshLoadFailureException("Element #%d is concave or badly oriented in initial_single_check().", e->id);

        // check first the positivity of the jacobian
        double3* pt = quad->get_points(mo, e->get_mode());
        double2x2* m = r.get_inv_ref_map(mo);
        double* jac = r.get_jacobian(mo);
        for (i = 0; i < quad->get_num_points(mo, e->get_mode()); i++)
          if(jac[i] <= 0.0)
            throw Hermes::Exceptions::MeshLoadFailureException("Element #%d is concave or badly oriented in initial_single_check().", e->id);
      }
    }

    void Mesh::initial_multimesh_check(Hermes::vector<MeshSharedPtr > meshes)
    {
    }

    void Mesh::create(int nv, double2* verts, int nt, int3* tris, std::string* tri_markers,
      int nq, int4* quads, std::string* quad_markers, int nm, int2* mark, std::string* boundary_markers)
    {
      free();

      // initialize hash table
      int size = 16;
      while (size < 2*nv) size *= 2;
      init(size);

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
      Element* e;
      for (int i = 0; i < nt; i++)
      {
        e = create_triangle(this->element_markers_conversion.insert_marker(tri_markers[i]), &nodes[tris[i][0]], &nodes[tris[i][1]],
          &nodes[tris[i][2]], NULL);
      }

      // create quads
      for (int i = 0; i < nq; i++)
      {
        e = create_quad(this->element_markers_conversion.insert_marker(quad_markers[i]), &nodes[quads[i][0]], &nodes[quads[i][1]],
          &nodes[quads[i][2]], &nodes[quads[i][3]], NULL);
      }

      // set boundary markers
      for (int i = 0; i < nm; i++)
      {
        Node* en = peek_edge_node(mark[i][0], mark[i][1]);
        if(en == NULL)
          throw Hermes::Exceptions::Exception("Boundary data error (edge does not exist)");

        en->marker = this->boundary_markers_conversion.insert_marker(boundary_markers[i]);

        nodes[mark[i][0]].bnd = 1;
        nodes[mark[i][1]].bnd = 1;
        en->bnd = 1;
      }

      nbase = nactive = ninitial = nt + nq;
      seq = g_mesh_seq++;
    }

    int Mesh::get_num_elements() const
    {
      if(this == NULL) throw Hermes::Exceptions::Exception("this == NULL in Mesh::get_num_elements().");
      if(this->seq < 0)
        return -1;
      else
        return elements.get_num_items();
    }

    /// Returns the number of coarse mesh elements.
    int Mesh::get_num_base_elements() const
    {
      if(this == NULL) throw Hermes::Exceptions::Exception("this == NULL in Mesh::get_num_base_elements().");

      if(this->seq < 0)
        return -1;
      else
        return nbase;
    }

    /// Returns the number of coarse mesh elements.
    int Mesh::get_num_used_base_elements() const
    {
      int toReturn = 0;
      if(this == NULL)
        throw Hermes::Exceptions::Exception("this == NULL in Mesh::get_num_base_elements().");

      if(this->seq < 0)
        return -1;
      else
      {
        Element* e;
        for_all_base_elements(e, this)
          toReturn++;
      }
      return toReturn;
    }

    /// Returns the current number of active elements in the mesh.
    int Mesh::get_num_active_elements() const
    {
      if(this == NULL) 
        throw Hermes::Exceptions::Exception("this == NULL in Mesh::get_num_active_elements().");
      if(this->seq < 0)
        return -1;
      else
        return nactive;
    }

    /// Returns the maximum node id number plus one.
    int Mesh::get_max_element_id() const
    {
      if(this == NULL) 
        throw Hermes::Exceptions::Exception("this == NULL in Mesh::get_max_element_id().");
      if(this->seq < 0)
        return -1;
      else
        return elements.get_size();
    }

    int Mesh::get_num_vertex_nodes() const
    {
      if(this == NULL) 
        throw Hermes::Exceptions::Exception("this == NULL in Mesh::get_num_vertex_nodes().");
      if(this->seq < 0)
        return -1;
      else
      {
        int to_return = 0;
        for (int i = 0; i < this->get_num_nodes(); i++)
        {
          if(get_node(i)->used && !get_node(i)->type)
            to_return++;
        }
        return to_return;
      }
    }

    int Mesh::get_num_edge_nodes() const
    {
      if(this == NULL) 
        throw Hermes::Exceptions::Exception("this == NULL in Mesh::get_num_vertex_nodes().");
      if(this->seq < 0)
        return -1;
      else
      {
        int to_return = 0;
        for (int i = 0; i < this->get_num_nodes(); i++)
        {
          if(get_node(i)->used && get_node(i)->type)
            to_return++;
        }
        return to_return;
      }
    }

    Element* Mesh::get_element(int id) const
    {
      if(id < 0 || id >= elements.get_size())
        throw Hermes::Exceptions::Exception("Invalid element ID %d, current range:[0; %d]", id, elements.get_size());
      return &(elements[id]);
    }

    unsigned Mesh::get_seq() const
    {
      return seq;
    }

    void Mesh::calc_bounding_box()
    {
      // find bounding box of the whole mesh
      bool first = true;
      Node* n;
      for_all_vertex_nodes(n, this)
      {
        if(first)
        {
          this->bottom_left_x = this->top_right_x = n->x;
          this->bottom_left_y = this->top_right_y = n->y;
          first = false;
        }
        else
        {
          if(n->x > this->top_right_x)
            this->top_right_x = n->x;
          if(n->x < this->bottom_left_x)
            this->bottom_left_x = n->x;
          if(n->y > this->top_right_y)
            this->top_right_y = n->y;
          if(n->y < this->bottom_left_y)
            this->bottom_left_y = n->y;
        }
      }
    }

    void Mesh::get_bounding_box(double& bottom_left_x_, double& bottom_left_y_, double& top_right_x_, double& top_right_y_)
    {
      if(!this->bounding_box_calculated)
        this->calc_bounding_box();

      top_right_x_ = this->top_right_x;
      bottom_left_x_ = this->bottom_left_x;
      top_right_y_ = this->top_right_y;
      bottom_left_y_ = this->bottom_left_y;

      bounding_box_calculated = true;
    }

    void Mesh::set_seq(unsigned seq)
    {
      this->seq = seq;
    }

    Element* Mesh::get_element_fast(int id) const
    {
      return &(elements[id]);
    }

    int Mesh::get_edge_sons(Element* e, int edge, int& son1, int& son2) const
    {
      assert(!e->active);

      if(!e->is_triangle())
      {
        if(e->sons[2] == NULL) // horz quad
        {
          if(edge == 0 || edge == 2) { son1 = edge >> 1;   return 1; }
          else if(edge == 1) { son1 = 0; son2 = 1; return 2; }
          else { son1 = 1; son2 = 0; return 2; }
        }
        else if(e->sons[0] == NULL) // vert quad
        {
          if(edge == 1 || edge == 3) { son1 = (edge == 1) ? 3 : 2; return 1; }
          else if(edge == 0) { son1 = 2; son2 = 3; return 2; }
          else { son1 = 3; son2 = 2; return 2; }
        }
      }

      // triangle or 4-son quad
      son1 = edge;
      son2 = e->next_vert(edge);
      return 2;
    }

    Element* Mesh::create_triangle(int marker, Node* v0, Node* v1, Node* v2, CurvMap* cm, int id)
    {
      // create a new element
      Element* e = elements.add();

      if(id != -1)
        e->id = id;

      // initialize the new element
      e->active = 1;
      e->marker = marker;
      e->nvert = 3;
      e->iro_cache = -1;
      e->cm = cm;
      e->parent = NULL;
      e->visited = false;

      // set vertex and edge node pointers
      if(v0 == v1 || v1 == v2 || v2 == v0)
        throw Hermes::Exceptions::MeshLoadFailureException("Some of the vertices of element #%d are identical which is impossible.", e->id);
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

    Element* Mesh::create_quad(int marker, Node* v0, Node* v1, Node* v2, Node* v3,
      CurvMap* cm, int id)
    {
      // create a new element
      Element* e = elements.add();

      if(id != -1)
        e->id = id;

      // initialize the new element
      e->active = 1;
      e->marker = marker;
      e->nvert = 4;
      e->iro_cache = -1;
      e->cm = cm;
      e->parent = NULL;
      e->visited = false;

      // set vertex and edge node pointers
      if(v0 == v1 || v1 == v2 || v2 == v3 || v3 == v0 || v2 == v0 || v3 == v1)
        throw Hermes::Exceptions::MeshLoadFailureException("Some of the vertices of element #%d are identical which is not right.", e->id);
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

    CurvMap* create_son_curv_map(Element* e, int son)
    {
      // if the top three bits of part are nonzero, we would overflow
      // -- make the element non-curvilinear
      if(e->cm->part & 0xe000000000000000ULL)
        return NULL;

      // if the parent element is already almost straight-edged,
      // the son will be even more straight-edged
      if(e->iro_cache == 0)
        return NULL;

      CurvMap* cm = new CurvMap;
      if(e->cm->toplevel == false)
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

    void Mesh::refine_triangle_to_triangles(Element* e, Element** sons_out)
    {
      // remember the markers of the edge nodes
      int bnd[3] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd    };
      int mrk[3] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker };

      // obtain three mid-edge vertex nodes
      Node* x0, *x1, *x2;
      x0 = get_vertex_node(e->vn[0]->id, e->vn[1]->id);
      x1 = get_vertex_node(e->vn[1]->id, e->vn[2]->id);
      x2 = get_vertex_node(e->vn[2]->id, e->vn[0]->id);

      CurvMap* cm[H2D_MAX_NUMBER_EDGES];
      memset(cm, 0, sizeof(cm));

      // adjust mid-edge coordinates if this is a curved element
      if(e->is_curved())
      {
        double2 pt[3] = { { 0.0, -1.0 }, { 0.0, 0.0 }, { -1.0, 0.0 } };
        e->cm->get_mid_edge_points(e, pt, 3);
        x0->x = pt[0][0]; x0->y = pt[0][1];
        x1->x = pt[1][0]; x1->y = pt[1][1];
        x2->x = pt[2][0]; x2->y = pt[2][1];

        // create CurvMaps for sons (pointer to parent element, part)
        for (int i = 0; i < 4; i++)
          cm[i] = create_son_curv_map(e, i);
      }

      // create the four sons
      Element* sons[H2D_MAX_ELEMENT_SONS];
      sons[0] = create_triangle(e->marker, e->vn[0], x0, x2, cm[0]);
      sons[1] = create_triangle(e->marker, x0, e->vn[1], x1, cm[1]);
      sons[2] = create_triangle(e->marker, x2, x1, e->vn[2], cm[2]);
      sons[3] = create_triangle(e->marker, x1, x2, x0, cm[3]);

      // update coefficients of curved reference mapping
      for (int i = 0; i < 4; i++)
        if(sons[i]->is_curved())
          sons[i]->cm->update_refmap_coeffs(sons[i]);

      // deactivate this element and unregister from its nodes
      e->active = 0;
      this->nactive += 3;
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
      {
        if(sons[i] != NULL) sons[i]->parent = e;
      }

      // copy son pointers (could not have been done earlier because of the union)
      memcpy(e->sons, sons, 4 * sizeof(Element*));

      // If sons_out != NULL, copy son pointers there.
      if(sons_out != NULL)
      {
        for(int i = 0; i < 3; i++) sons_out[i] = sons[i];
      }
    }

    Node* get_vertex_node(Node* v1, Node* v2)
    {
      // initialize the new Node
      Node* newnode = new Node();
      newnode->type = HERMES_TYPE_VERTEX;
      newnode->ref = 0;
      newnode->bnd = 0;
      newnode->p1 = -9999;
      newnode->p2 = -9999;
      newnode->x = (v1->x + v2->x) * 0.5;
      newnode->y = (v1->y + v2->y) * 0.5;

      return newnode;
    }

    Node* get_edge_node()
    {
      // initialize the new Node
      Node* newnode = new Node();
      newnode->type = HERMES_TYPE_EDGE;
      newnode->ref = 0;
      newnode->bnd = 0;
      newnode->p1 = -9999;
      newnode->p2 = -9999;
      newnode->marker = 0;
      newnode->elem[0] = newnode->elem[1] = NULL;

      return newnode;
    }

    void Mesh::refine_quad(Element* e, int refinement, Element** sons_out)
    {
      int i, j;
      Element* sons[H2D_MAX_ELEMENT_SONS] = {NULL, NULL, NULL, NULL};

      // remember the markers of the edge nodes
      int bnd[H2D_MAX_NUMBER_EDGES] = { e->en[0]->bnd, e->en[1]->bnd, e->en[2]->bnd, e->en[3]->bnd };
      int mrk[H2D_MAX_NUMBER_EDGES] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker, e->en[3]->marker };

      // deactivate this element and unregister from its nodes
      e->active = false;
      nactive--;
      e->unref_all_nodes(this);

      // now the original edge nodes may no longer exist...
      CurvMap* cm[H2D_MAX_NUMBER_EDGES];
      memset(cm, 0, sizeof(cm));

      // default refinement: one quad to four quads
      if(refinement == 0)
      {
        // obtain four mid-edge vertex nodes and one mid-element vertex node
        Node* x0, *x1, *x2, *x3, *mid;
        x0 = get_vertex_node(e->vn[0]->id, e->vn[1]->id);
        x1 = get_vertex_node(e->vn[1]->id, e->vn[2]->id);
        x2 = get_vertex_node(e->vn[2]->id, e->vn[3]->id);
        x3 = get_vertex_node(e->vn[3]->id, e->vn[0]->id);
        mid = get_vertex_node(x0->id, x2->id);

        // adjust mid-edge coordinates if this is a curved element
        if(e->is_curved())
        {
          double2 pt[5] = { { 0.0, -1.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }, { -1.0, 0.0 }, { 0.0, 0.0 } };
          e->cm->get_mid_edge_points(e, pt, 5);
          x0->x = pt[0][0];  x0->y = pt[0][1];
          x1->x = pt[1][0];  x1->y = pt[1][1];
          x2->x = pt[2][0];  x2->y = pt[2][1];
          x3->x = pt[3][0];  x3->y = pt[3][1];
          mid->x = pt[4][0]; mid->y = pt[4][1];

          // create CurvMaps for sons (pointer to parent element, part)
          for (i = 0; i < H2D_MAX_ELEMENT_SONS; i++)
            cm[i] = create_son_curv_map(e, i);
        }

        // create the four sons
        sons[0] = create_quad(e->marker, e->vn[0], x0, mid, x3, cm[0]);
        sons[1] = create_quad(e->marker, x0, e->vn[1], x1, mid, cm[1]);
        sons[2] = create_quad(e->marker, mid, x1, e->vn[2], x2, cm[2]);
        sons[3] = create_quad(e->marker, x3, mid, x2, e->vn[3], cm[3]);

        // Increase the number of active elements by 4.
        this->nactive += H2D_MAX_ELEMENT_SONS;

        // set correct boundary markers for the new edge nodes
        for (i = 0; i < H2D_MAX_NUMBER_EDGES; i++)
        {
          j = (i > 0) ? i-1 : 3;
          sons[i]->en[j]->bnd = bnd[j];  sons[i]->en[j]->marker = mrk[j];
          sons[i]->en[i]->bnd = bnd[i];  sons[i]->en[i]->marker = mrk[i];
          sons[i]->vn[j]->bnd = bnd[j];
        }
      }

      // refinement '1': one quad to two 'horizontal' quads
      else if(refinement == 1)
      {
        Node* x1, *x3;
        x1 = get_vertex_node(e->vn[1]->id, e->vn[2]->id);
        x3 = get_vertex_node(e->vn[3]->id, e->vn[0]->id);

        // adjust mid-edge coordinates if this is a curved element
        if(e->is_curved())
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

        this->nactive += 2;

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
      else if(refinement == 2)
      {
        Node* x0, *x2;
        x0 = get_vertex_node(e->vn[0]->id, e->vn[1]->id);
        x2 = get_vertex_node(e->vn[2]->id, e->vn[3]->id);

        // adjust mid-edge coordinates if this is a curved element
        if(e->is_curved())
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

        this->nactive += 2;

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
        if(sons[i] != NULL && sons[i]->cm != NULL)
          sons[i]->cm->update_refmap_coeffs(sons[i]);

      // optimization: iro never gets worse
      if(e->iro_cache == 0)
        for (i = 0; i < 4; i++)
          if(sons[i] != NULL)
            sons[i]->iro_cache = 0;

      // set pointers to parent element for sons
      for(int i = 0; i < 4; i++)
        if(sons[i] != NULL) sons[i]->parent = e;

      // copy son pointers (could not have been done earlier because of the union)
      memcpy(e->sons, sons, sizeof(sons));

      // If sons_out != NULL, copy son pointers there.
      if(sons_out != NULL)
      {
        for(int i = 0; i < 4; i++) sons_out[i] = sons[i];
      }
    }

    void Mesh::unrefine_element_internal(Element* e)
    {
      this->refinements.push_back(std::pair<unsigned int, int>(e->id, -1));
      assert(!e->active);
      unsigned int i;
      int s1, s2;

      // obtain markers and bnds from son elements
      int mrk[H2D_MAX_NUMBER_EDGES], bnd[H2D_MAX_NUMBER_EDGES];
      for (unsigned i = 0; i < e->get_nvert(); i++)
      {
        get_edge_sons(e, i, s1, s2);
        assert(e->sons[s1]->active);
        mrk[i] = e->sons[s1]->en[i]->marker;
        bnd[i] = e->sons[s1]->en[i]->bnd;
      }

      // remove all sons
      for (i = 0; i < H2D_MAX_ELEMENT_SONS; i++)
      {
        Element* son = e->sons[i];
        if(son != NULL)
        {
          son->unref_all_nodes(this);
          if(son->cm != NULL) delete son->cm;
          elements.remove(son->id);
          this->nactive--;
        }
      }

      // recreate edge nodes
      for (i = 0; i < e->get_nvert(); i++)
        e->en[i] = this->get_edge_node(e->vn[i]->id, e->vn[e->next_vert(i)]->id);

      e->ref_all_nodes();
      e->active = 1;
      nactive++;

      // restore edge node markers and bnds
      for (i = 0; i < e->get_nvert(); i++)
      {
        e->en[i]->marker = mrk[i];
        e->en[i]->bnd = bnd[i];
      }
    }

    void Mesh::refine_element(Element* e, int refinement)
    {
      this->refinements.push_back(std::pair<unsigned int, int>(e->id, refinement));

      if(e->is_triangle())
      {
        if(refinement == 3)
        {
          refine_triangle_to_quads(e);
        }
        else
        {
          this->refine_triangle_to_triangles(e);
        }
      }
      else refine_quad(e, refinement);

      this->seq = g_mesh_seq++;
    }

    void Mesh::refine_element_id(int id, int refinement)
    {
      if(refinement == -1)
        return;
      Element* e = this->get_element(id);
      if(!e->used) throw Hermes::Exceptions::Exception("Invalid element id number.");
      if(!e->active) throw Hermes::Exceptions::Exception("Attempt to refine element #%d which has been refined already.", e->id);
      this->refine_element(e, refinement);
    }

    void Mesh::refine_all_elements(int refinement, bool mark_as_initial)
    {
      ninitial = this->get_max_element_id();

      if(refinement == -1)
        return;

      elements.set_append_only(true);

      Element* e;

      for_all_active_elements(e, this)
        refine_element(e, refinement);

      elements.set_append_only(false);

      if(mark_as_initial)
        ninitial = this->get_max_element_id();
    }

    static int rtb_marker;
    static bool rtb_aniso;
    static char* rtb_vert;

    void Mesh::refine_by_criterion(int (*criterion)(Element*), int depth, bool mark_as_initial)
    {
      Element* e;
      elements.set_append_only(true);
      for (int r, i = 0; i < depth; i++)
      {
        for_all_active_elements(e, this)
        {
          if((r = criterion(e)) >= 0)
            refine_element_id(e->id, r);
        }
      }
      elements.set_append_only(false);

      if(mark_as_initial)
        ninitial = this->get_max_element_id();
    }

    static int rtv_id;

    static int rtv_criterion(Element* e)
    {
      for (unsigned int i = 0; i < e->get_nvert(); i++)
        if(e->vn[i]->id == rtv_id)
          return 0;
      return -1;
    }

    void Mesh::refine_towards_vertex(int vertex_id, int depth, bool mark_as_initial)
    {
      rtv_id = vertex_id;
      refine_by_criterion(rtv_criterion, depth);
      if(mark_as_initial)
        ninitial = this->get_max_element_id();
    }

    int rtb_criterion(Element* e)
    {
      unsigned int i;
      for (i = 0; i < e->get_nvert(); i++)
      {
        if(e->en[i]->marker == rtb_marker || rtb_vert[e->vn[i]->id])
        {
          break;
        }
      }

      if(i >= e->get_nvert()) return -1;
      // triangle should be split into 3 quads
      //  if(e->is_triangle() && rtb_tria_to_quad) return 3;
      // triangle should be split into 4 triangles or quad should
      // be split into 4 quads
      if(e->is_triangle() || !rtb_aniso) return 0;

      // quads - anisotropic case 1
      if((e->en[0]->marker == rtb_marker && !rtb_vert[e->vn[2]->id] && !rtb_vert[e->vn[3]->id]) ||
        (e->en[2]->marker == rtb_marker && !rtb_vert[e->vn[0]->id] && !rtb_vert[e->vn[1]->id]) ||
        (e->en[0]->marker == rtb_marker && e->en[2]->marker == rtb_marker &&
        e->en[1]->marker != rtb_marker && e->en[3]->marker != rtb_marker)) return 1;

      // quads - anisotropic case 2
      if((e->en[1]->marker == rtb_marker && !rtb_vert[e->vn[3]->id] && !rtb_vert[e->vn[0]->id]) ||
        (e->en[3]->marker == rtb_marker && !rtb_vert[e->vn[1]->id] && !rtb_vert[e->vn[2]->id]) ||
        (e->en[1]->marker == rtb_marker && e->en[3]->marker == rtb_marker &&
        e->en[0]->marker != rtb_marker && e->en[2]->marker != rtb_marker)) return 2;

      return 0;
    }

    void Mesh::refine_towards_boundary(Hermes::vector<std::string> markers, int depth, bool aniso, bool mark_as_initial)
    {
      rtb_aniso = aniso;
      bool refined = true;

      // refinement: refine all elements to quad elements.
      for (int i = 0; i < depth; i++)
      {
        refined = false;
        int size = get_max_node_id() + 1;
        rtb_vert = new char[size];
        memset(rtb_vert, 0, sizeof(char) * size);

        Element* e;
        for_all_active_elements(e, this)
          for (unsigned int j = 0; j < e->get_nvert(); j++)
          {
            bool marker_matched = false;
            for(unsigned int marker_i = 0; marker_i < markers.size(); marker_i++)
              if(e->en[j]->marker == this->boundary_markers_conversion.get_internal_marker(markers[marker_i]).marker)
                marker_matched = true;
            if(marker_matched)
            {
              rtb_vert[e->vn[j]->id] = rtb_vert[e->vn[e->next_vert(j)]->id] = 1;
              refined = true;
            }
          }

          refine_by_criterion(rtb_criterion, 1);
          delete [] rtb_vert;
      }

      if(mark_as_initial)
        ninitial = this->get_max_element_id();
      if(!refined)
        throw Hermes::Exceptions::Exception("None of the markers in Mesh::refine_towards_boundary found in the Mesh.");
    }

    void Mesh::refine_towards_boundary(std::string marker, int depth, bool aniso, bool mark_as_initial)
    {
      if(marker == HERMES_ANY)
        for(std::map<int, std::string>::iterator it = this->boundary_markers_conversion.conversion_table.begin(); it != this->boundary_markers_conversion.conversion_table.end(); ++it)
          refine_towards_boundary(it->second, depth, aniso, mark_as_initial);

      else
      {
        bool refined = true;
        rtb_marker = this->boundary_markers_conversion.get_internal_marker(marker).marker;
        rtb_aniso = aniso;

        // refinement: refine all elements to quad elements.
        for (int i = 0; i < depth; i++)
        {
          refined = false;
          int size = get_max_node_id() + 1;
          rtb_vert = new char[size];
          memset(rtb_vert, 0, sizeof(char) * size);

          Element* e;
          for_all_active_elements(e, this)
            for (unsigned int j = 0; j < e->get_nvert(); j++)
            {
              if(e->en[j]->marker == rtb_marker)
              {
                rtb_vert[e->vn[j]->id] = rtb_vert[e->vn[e->next_vert(j)]->id] = 1;
                refined = true;
              }
            }

            refine_by_criterion(rtb_criterion, 1);
            delete [] rtb_vert;
        }

        if(mark_as_initial)
          ninitial = this->get_max_element_id();
        if(!refined)
          throw Hermes::Exceptions::Exception("None of the markers in Mesh::refine_towards_boundary found in the Mesh.");
      }
    }

    void Mesh::refine_in_area(std::string marker, int depth, bool mark_as_initial)
    {
      Hermes::vector<std::string> markers;
      markers.push_back(marker);
      this->refine_in_areas(markers, depth, mark_as_initial);
    }

    void Mesh::refine_in_areas(Hermes::vector<std::string> markers, int depth, bool mark_as_initial)
    {
      Hermes::vector<int> internal_markers;
      bool any_marker = false;
      for(unsigned int marker_i = 0; marker_i < markers.size(); marker_i++)
      {
        if(markers[marker_i] == HERMES_ANY)
        {
          any_marker = true;
          break;
        }
        int marker = this->element_markers_conversion.get_internal_marker(markers[marker_i]).marker;

        internal_markers.push_back(marker);
      }

      bool refined = true;
      for (int i = 0; i < depth; i++)
      {
        refined = false;
        Element* e;
        if(any_marker)
        {
          for_all_active_elements(e, this)
          {
            this->refine_element(e, 0);
            refined = true;
          }
        }
        else
        {
          for_all_active_elements(e, this)
          {
            for(unsigned int marker_i = 0; marker_i < internal_markers.size(); marker_i++)
              if(e->marker == internal_markers[marker_i])
              {
                this->refine_element(e, 0);
                refined = true;
              }
          }
        }
      }

      if(mark_as_initial)
        ninitial = this->get_max_element_id();
      if(!refined)
        throw Hermes::Exceptions::Exception("None of the markers in Mesh::refine_in_areas found in the Mesh.");
    }

    void Mesh::unrefine_element_id(int id)
    {
      Element* e = get_element(id);
      if(!e->used) throw Hermes::Exceptions::Exception("Invalid element id number.");
      if(e->active) return;

      for (int i = 0; i < 4; i++)
        if(e->sons[i] != NULL)
          unrefine_element_id(e->sons[i]->id);

      unrefine_element_internal(e);
      seq = g_mesh_seq++;
    }

    void Mesh::unrefine_all_elements(bool keep_initial_refinements)
    {
      // find inactive elements with active sons
      Hermes::vector<int> list;
      Element* e;
      for_all_inactive_elements(e, this)
      {
        bool found = true;
        for (unsigned int i = 0; i < 4; i++)
          if(e->sons[i] != NULL &&
            (!e->sons[i]->active || (keep_initial_refinements && e->sons[i]->id < ninitial))
            )
          { found = false; break; }

          if(found) list.push_back(e->id);
      }

      // unrefine the found elements
      for (unsigned int i = 0; i < list.size(); i++)
        unrefine_element_id(list[i]);
    }

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

    double Mesh::vector_length(double a_1, double a_2)
    {
      return sqrt(sqr(a_1) + sqr(a_2));
    }

    bool Mesh::same_line(double p_1, double p_2, double q_1, double q_2, double r_1, double r_2)
    {
      double pq_1 = q_1 - p_1, pq_2 = q_2 - p_2, pr_1 = r_1 - p_1, pr_2 = r_2 - p_2;
      double length_pq = vector_length(pq_1, pq_2);
      double length_pr = vector_length(pr_1, pr_2);
      double sin_angle = (pq_1*pr_2 - pq_2*pr_1)/(length_pq*length_pr);
      if(fabs(sin_angle) < Hermes::Epsilon) return true;
      else return false;
    }

    bool Mesh::is_convex(double a_1, double a_2, double b_1, double b_2)
    {
      if(a_1*b_2 - a_2*b_1 > 0) return true;
      else return false;
    }

    void Mesh::check_triangle(int i, Node *&v0, Node *&v1, Node *&v2)
    {
      // checking that all edges have nonzero length
      double
        length_1 = vector_length(v1->x - v0->x, v1->y - v0->y),
        length_2 = vector_length(v2->x - v1->x, v2->y - v1->y),
        length_3 = vector_length(v0->x - v2->x, v0->y - v2->y);
      if(length_1 < Hermes::epsilon || length_2 < Hermes::epsilon || length_3 < Hermes::epsilon)
        throw Hermes::Exceptions::Exception("Edge of triangular element #%d has length less than Hermes::epsilon.", i);

      // checking that vertices do not lie on the same line
      if(same_line(v0->x, v0->y, v1->x, v1->y, v2->x, v2->y))
        throw Hermes::Exceptions::Exception("Triangular element #%d: all vertices lie on the same line.", i);

      // checking positive orientation. If not positive, swapping vertices
      if(!is_convex(v1->x - v0->x, v1->y - v0->y, v2->x - v0->x, v2->y - v0->y))
      {
        std::swap(v1, v2);
      }
    }

    void Mesh::check_quad(int i, Node *&v0, Node *&v1, Node *&v2, Node *&v3)
    {
      // checking that all edges have nonzero length
      double
        length_1 = vector_length(v1->x - v0->x, v1->y - v0->y),
        length_2 = vector_length(v2->x - v1->x, v2->y - v1->y),
        length_3 = vector_length(v3->x - v2->x, v3->y - v2->y),
        length_4 = vector_length(v0->x - v3->x, v0->y - v3->y);
      if(length_1 < Hermes::epsilon || length_2 < Hermes::epsilon || length_3 < Hermes::epsilon || length_4 < Hermes::epsilon)
        throw Hermes::Exceptions::Exception("Edge of quad element #%d has length less than Hermes::epsilon.", i);

      // checking that both diagonals have nonzero length
      double
        diag_1 = vector_length(v2->x - v0->x, v2->y - v0->y),
        diag_2 = vector_length(v3->x - v1->x, v3->y - v1->y);
      if(diag_1 < Hermes::epsilon || diag_2 < Hermes::epsilon)
        throw Hermes::Exceptions::Exception("Diagonal of quad element #%d has length less than Hermes::epsilon.", i);

      // checking that vertices v0, v1, v2 do not lie on the same line
      if(same_line(v0->x, v0->y, v1->x, v1->y, v2->x, v2->y))
        throw Hermes::Exceptions::Exception("Quad element #%d: vertices v0, v1, v2 lie on the same line.", i);
      // checking that vertices v0, v1, v3 do not lie on the same line
      if(same_line(v0->x, v0->y, v1->x, v1->y, v3->x, v3->y))
        throw Hermes::Exceptions::Exception("Quad element #%d: vertices v0, v1, v3 lie on the same line.", i);
      // checking that vertices v0, v2, v3 do not lie on the same line
      if(same_line(v0->x, v0->y, v2->x, v2->y, v3->x, v3->y))
        throw Hermes::Exceptions::Exception("Quad element #%d: vertices v0, v2, v3 lie on the same line.", i);
      // checking that vertices v1, v2, v3 do not lie on the same line
      if(same_line(v1->x, v1->y, v2->x, v2->y, v3->x, v3->y))
        throw Hermes::Exceptions::Exception("Quad element #%d: vertices v1, v2, v3 lie on the same line.", i);

      // checking that vertex v1 lies on the right of the diagonal v2-v0
      int vertex_1_ok = is_convex(v1->x - v0->x, v1->y - v0->y, v2->x - v0->x, v2->y - v0->y);
      if(!vertex_1_ok) throw Hermes::Exceptions::Exception("Vertex v1 of quad element #%d does not lie on the right of the diagonal v2-v0.", i);
      // checking that vertex v3 lies on the left of the diagonal v2-v0
      int vertex_3_ok = is_convex(v2->x - v0->x, v2->y - v0->y, v3->x - v0->x, v3->y - v0->y);
      if(!vertex_3_ok) throw Hermes::Exceptions::Exception("Vertex v3 of quad element #%d does not lie on the left of the diagonal v2-v0.", i);
      // checking that vertex v2 lies on the right of the diagonal v3-v1
      int vertex_2_ok = is_convex(v2->x - v1->x, v2->y - v1->y, v3->x - v1->x, v3->y - v1->y);
      if(!vertex_2_ok) throw Hermes::Exceptions::Exception("Vertex v2 of quad element #%d does not lie on the right of the diagonal v3-v1.", i);
      // checking that vertex v0 lies on the left of the diagonal v3-v1
      int vertex_0_ok = is_convex(v3->x - v1->x, v3->y - v1->y, v0->x - v1->x, v0->y - v1->y);
      if(!vertex_0_ok) throw Hermes::Exceptions::Exception("Vertex v0 of quad element #%d does not lie on the left of the diagonal v2-v1.", i);
    }

    bool Mesh::rescale(double x_ref, double y_ref)
    {
      // Go through all vertices and rescale coordinates.
      Node* n;
      for_all_vertex_nodes(n, this) {
        n->x /= x_ref;
        n->y /= y_ref;
      }

      // If curvilinear, throw an exception.
      Element* e;
      for_all_elements(e, this)
        if(e->cm != NULL)
        {
          throw CurvedException(e->id);
          return false;
        }
        return true;
    }

    void Mesh::copy(MeshSharedPtr mesh)
    {
      unsigned int i;

      free();

      // copy nodes and elements
      HashTable::copy(mesh.get());
      elements.copy(mesh->elements);

      this->refinements = mesh->refinements;

      Element* e;
      for_all_elements(e, this)
      {
        // update vertex node pointers
        for (i = 0; i < e->get_nvert(); i++)
          e->vn[i] = &nodes[e->vn[i]->id];

        if(e->active)
        {
          // update edge node pointers
          for (i = 0; i < e->get_nvert(); i++)
            e->en[i] = &nodes[e->en[i]->id];
        }
        else
        {
          // update son pointers
          for (i = 0; i < 4; i++)
            if(e->sons[i] != NULL)
              e->sons[i] = &elements[e->sons[i]->id];
        }

        // copy CurvMap, update its parent
        if(e->cm != NULL)
        {
          e->cm = new CurvMap(e->cm);
          if(!e->cm->toplevel)
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
          if(node->elem[i] != NULL)
            node->elem[i] = &elements[node->elem[i]->id];

      nbase = mesh->nbase;
      nactive = mesh->nactive;
      ntopvert = mesh->ntopvert;
      ninitial = mesh->ninitial;
      seq = mesh->seq;
      boundary_markers_conversion = mesh->boundary_markers_conversion;
      element_markers_conversion = mesh->element_markers_conversion;
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

    void Mesh::init(int size)
    {
      HashTable::init(size);
    }

    void Mesh::copy_base(MeshSharedPtr mesh)
    {
      //printf("Calling Mesh::free() in Mesh::copy_base().\n");
      free();
      init();

      // copy top-level vertex nodes
      for (int i = 0; i < mesh->get_max_node_id(); i++)
      {
        Node* node = &(mesh->nodes[i]);
        if(node->ref < TOP_LEVEL_REF) break;
        Node* newnode = nodes.add();
        assert(newnode->id == i && node->type == HERMES_TYPE_VERTEX);
        memcpy(newnode, node, sizeof(Node));
        newnode->ref = TOP_LEVEL_REF;
      }

      // copy base elements
      Element* e;
      for (int id = 0; id < mesh->get_num_base_elements(); id++)
      {
        e = mesh->get_element_fast(id);
        if(!e->used)
        {
          Element* e_temp = elements.add();
          e_temp->used = false;
          continue;
        }
        Element* enew;
        Node *v0 = &nodes[e->vn[0]->id], *v1 = &nodes[e->vn[1]->id], *v2 = &nodes[e->vn[2]->id];
        if(e->is_triangle())
          enew = this->create_triangle(e->marker, v0, v1, v2, NULL);
        else
          enew = this->create_quad(e->marker, v0, v1, v2, &nodes[e->vn[3]->id], NULL);

        // copy edge markers
        for (unsigned int j = 0; j < e->get_nvert(); j++)
        {
          Node* en = get_base_edge_node(e, j);
          enew->en[j]->bnd = en->bnd; // copy bnd data from the active el.
          enew->en[j]->marker = en->marker;
        }

        if(e->is_curved())
          enew->cm = new CurvMap(e->cm);
      }

      this->boundary_markers_conversion = mesh->boundary_markers_conversion;
      this->element_markers_conversion = mesh->element_markers_conversion;

      nbase = nactive = ninitial = mesh->nbase;
      ntopvert = mesh->ntopvert;
      seq = g_mesh_seq++;
    }

    void Mesh::free()
    {
      Element* e;
      for_all_elements(e, this)
      {
        if(e->cm != NULL)
        {
          delete e->cm;
          e->cm = NULL; // fixme!!!
        }
      }
      elements.free();
      HashTable::free();

      if(this->meshHashGrid)
        delete this->meshHashGrid;

      this->boundary_markers_conversion.conversion_table.clear();
      this->boundary_markers_conversion.conversion_table_inverse.clear();
      this->element_markers_conversion.conversion_table.clear();
      this->element_markers_conversion.conversion_table_inverse.clear();
      this->refinements.clear();
      this->seq = -1;
    }

    Element* Mesh::element_on_physical_coordinates(double x, double y)
    {
      // If no hash grid exists, create one.
      if(!this->meshHashGrid)
        this->meshHashGrid = new MeshHashGrid(this);
      else
      {
        // If no the hash grid exists, but the mesh has been refined afterwards, re-create.
        if(this->get_seq() != this->meshHashGrid->get_mesh_seq())
        {
          delete this->meshHashGrid;
          this->meshHashGrid = new MeshHashGrid(this);
        }
      }

      return this->meshHashGrid->getElement(x, y);
    }

    void Mesh::copy_converted(MeshSharedPtr mesh)
    {
      free();
      HashTable::copy(mesh.get());
      this->boundary_markers_conversion = mesh->boundary_markers_conversion;
      this->element_markers_conversion = mesh->element_markers_conversion;

      // clear reference for all nodes
      for(int i = 0; i < nodes.get_size(); i++)
      {
        Node& node = nodes[i];
        if(node.type == HERMES_TYPE_EDGE) { //process only edge nodes
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
        if(e->is_triangle())
        {
          // create a new element
          enew = elements.add();
          enew->active = 1;
          enew->marker = e->marker;
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
        for (unsigned int j = 0; j < e->get_nvert(); j++)
        {
          Node* en = get_base_edge_node(e, j);
          enew->en[j]->bnd = en->bnd;
          enew->en[j]->marker = en->marker;
        }

        if(e->is_curved())
          enew->cm = new CurvMap(e->cm);
      }

      nbase = nactive = ninitial = mesh->nactive;
      ntopvert = mesh->ntopvert = get_num_nodes();
      seq = g_mesh_seq++;
    }

    void Mesh::convert_quads_to_triangles()
    {
      Element* e;

      elements.set_append_only(true);
      for_all_active_elements(e, this)
        refine_element_to_triangles_id(e->id);
      elements.set_append_only(false);

      MeshSharedPtr mesh_tmp_for_convert(new Mesh);
      mesh_tmp_for_convert->copy_converted(MeshSharedPtr(this));

      for (int i = 0; i < mesh_tmp_for_convert->ntopvert; i++)
      {
        if(mesh_tmp_for_convert->nodes[i].type == 1)
        {
          mesh_tmp_for_convert->nodes[i].y = 0.0;
        }
      }
      MeshReaderH2D loader_mesh_tmp_for_convert;
      char* mesh_file_tmp = NULL;
      mesh_file_tmp = tmpnam(NULL);
      loader_mesh_tmp_for_convert.save(mesh_file_tmp, mesh_tmp_for_convert);
      loader_mesh_tmp_for_convert.load(mesh_file_tmp, mesh_tmp_for_convert);
      remove(mesh_file_tmp);
      copy(mesh_tmp_for_convert);
    }

    void Mesh::convert_to_base()
    {
      Element* e;

      elements.set_append_only(true);
      for_all_active_elements(e, this)
        convert_element_to_base_id(e->id);
      elements.set_append_only(false);

      MeshSharedPtr mesh_tmp_for_convert(new Mesh);
      mesh_tmp_for_convert->copy_converted(MeshSharedPtr(this));
      for (int i = 0; i < mesh_tmp_for_convert->ntopvert; i++)
      {
        if(mesh_tmp_for_convert->nodes[i].type == 1)
        {
          mesh_tmp_for_convert->nodes[i].y = 0.0;
        }
      }
      MeshReaderH2D loader_mesh_tmp_for_convert;
      char* mesh_file_tmp = NULL;
      mesh_file_tmp = tmpnam(NULL);
      loader_mesh_tmp_for_convert.save(mesh_file_tmp, mesh_tmp_for_convert);
      loader_mesh_tmp_for_convert.load(mesh_file_tmp, mesh_tmp_for_convert);
      remove(mesh_file_tmp);
      copy(mesh_tmp_for_convert);
    }

    void Mesh::refine_triangle_to_quads(Element* e, Element** sons_out)
    {
      // remember the markers of the edge nodes
      int bnd[3] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd    };
      int mrk[3] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker };

      // obtain three mid-edge and one gravity vertex nodes
      Node* x0 = this->get_vertex_node(e->vn[0]->id, e->vn[1]->id);
      Node* x1 = this->get_vertex_node(e->vn[1]->id, e->vn[2]->id);
      Node* x2 = this->get_vertex_node(e->vn[2]->id, e->vn[0]->id);
      Node* mid = this->get_vertex_node(x0->id, e->vn[1]->id);

      mid->x = (x0->x + x1->x + x2->x)/3;
      mid->y = (x0->y + x1->y + x2->y)/3;

      // check if element e is a internal element.
      bool e_inter = true;
      for (unsigned int n = 0; n < e->get_nvert(); n++)
      {
        if(bnd[n] == 1)
          e_inter = false;
      }

      CurvMap* cm[4];
      memset(cm, 0, sizeof(cm));

      // adjust mid-edge and gravity coordinates if this is a curved element
      if(e->is_curved())
      {
        if(!e_inter)
        {
          double2 pt[4] = { { 0.0, -1.0 }, { 0.0, 0.0 }, { -1.0, 0.0 }, { -0.33333333, -0.33333333 } };
          e->cm->get_mid_edge_points(e, pt, 4);
          x0->x = pt[0][0]; x0->y = pt[0][1];
          x1->x = pt[1][0]; x1->y = pt[1][1];
          x2->x = pt[2][0]; x2->y = pt[2][1];
          mid->x = pt[3][0]; mid->y = pt[3][1];
        }
      }

      // get the boundary edge angle.
      double refinement_angle[3] = {0.0, 0.0, 0.0};
      if(e->is_curved() && (!e_inter))
      {
        // for base element.
        if(e->cm->toplevel == true)
        {
          for (unsigned int n = 0; n < e->get_nvert(); n++)
          {
            if(e->cm->nurbs[n] != NULL)
            {
              //this->info("angle = %f", e->cm->nurbs[n]->angle);
              refinement_angle[n] = e->cm->nurbs[n]->angle;
            }
          }
        }
        else
          // one level refinement.
          if(e->parent->cm->toplevel == true)
          {
            for (unsigned int n = 0; n < e->get_nvert(); n++)
            {
              if(e->parent->cm->nurbs[n] != NULL)
              {
                //this->info("angle = %f", e->parent->cm->nurbs[n]->angle);
                refinement_angle[n] = e->parent->cm->nurbs[n]->angle / 2;
              }
            }
          }
          else
            // two level refinements.
            if(e->parent->parent->cm->toplevel == true)
            {
              for (unsigned int n = 0; n < e->get_nvert(); n++)
              {
                if(e->parent->parent->cm->nurbs[n] != NULL)
                {
                  //this->info("angle = %f", e->parent->parent->cm->nurbs[n]->angle);
                  refinement_angle[n] = e->parent->parent->cm->nurbs[n]->angle / 4;
                }
              }
            }
            else
              // three level refinements.
              if(e->parent->parent->parent->cm->toplevel == true)
              {
                for (unsigned int n = 0; n < e->get_nvert(); n++)
                {
                  if(e->parent->parent->parent->cm->nurbs[n] != NULL)
                  {
                    //this->info("angle = %f", e->parent->parent->parent->cm->nurbs[n]->angle);
                    refinement_angle[n] = e->parent->parent->parent->cm->nurbs[n]->angle / 8;
                  }
                }
              }
              else
                // four level refinements.
                if(e->parent->parent->parent->parent->cm->toplevel == true)
                {
                  for (unsigned int n = 0; n < e->get_nvert(); n++)
                  {
                    if(e->parent->parent->parent->parent->cm->nurbs[n] != NULL)
                    {
                      //this->info("angle = %f", e->parent->parent->parent->parent->cm->nurbs[n]->angle);
                      refinement_angle[n] = e->parent->parent->parent->parent->cm->nurbs[n]->angle / 16;
                    }
                  }
                }
      }

      double angle2 = 0.0;
      int idx = 0;
      // create CurvMaps for sons if this is a curved element
      if((e->is_curved()) && (!e_inter))
      {
        for (idx = 0; idx < 2; idx++)
        {
          if(e->cm->nurbs[idx] != NULL)
          {
            cm[idx] = new CurvMap;
            memset(cm[idx], 0, sizeof(CurvMap));
            cm[idx + 1] = new CurvMap;
            memset(cm[idx + 1], 0, sizeof(CurvMap));
          }
        }

        idx = 0;
        if(e->cm->nurbs[idx] != NULL)
        {
          angle2 = refinement_angle[0] / 2;
          Node* node_temp = this->get_vertex_node(e->vn[idx%3]->id, e->vn[(idx + 1)%3]->id);

          for (int k = 0; k < 2; k++)
          {
            int p1, p2;
            int idx2 = 0;

            if(k == 0)
            {
              p1 = e->vn[(idx)%3]->id;
              p2 = node_temp->id;
              if(idx == 0) idx2 = 0;
              if(idx == 1) idx2 = 1;
              if(idx == 2) continue;
            }
            else
            {
              p1 = node_temp->id;
              p2 = e->vn[(idx + 1)%3]->id;
              idx = (idx + 1)%3;
              if(idx == 0) continue;
              if(idx == 1) idx2 = 0;
              if(idx == 2) idx2 = 0;
            }

            Nurbs* nurbs = new Nurbs;
            bool cricle = true;

            nurbs->arc = cricle;
            nurbs->degree = 2;

            int inner = 1, outer = 0;
            nurbs->np = inner + 2;
            nurbs->pt = new double3[nurbs->np];

            nurbs->pt[0][0] = nodes[p1].x;
            nurbs->pt[0][1] = nodes[p1].y;
            nurbs->pt[0][2] = 1.0;

            nurbs->pt[inner + 1][0] = nodes[p2].x;
            nurbs->pt[inner + 1][1] = nodes[p2].y;
            nurbs->pt[inner + 1][2] = 1.0;

            double angle = angle2;
            double a = (180.0 - angle) / 180.0 * M_PI;
            nurbs->angle = angle;

            // generate one control point
            double x = 1.0 / tan(a * 0.5);
            nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
            nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
            nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

            int i = 0;
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
        if(e->cm->nurbs[idx] != NULL)
        {
          angle2 = refinement_angle[1] / 2;
          Node* node_temp = this->get_vertex_node(e->vn[idx%3]->id, e->vn[(idx + 1)%3]->id);
          for (int k = 0; k < 2; k++)
          {
            int p1, p2;
            int idx2 = 0;
            if(k == 0)
            {
              p1 = e->vn[(idx)%3]->id;
              p2 = node_temp->id;
              if(idx == 0) idx2 = 0;
              if(idx == 1) idx2 = 1;
              if(idx == 2) continue;
            }
            else
            {
              p1 = node_temp->id;
              p2 = e->vn[(idx + 1)%3]->id;
              idx = (idx + 1)%3;
              if(idx == 0) continue;
              if(idx == 1) idx2 = 0;
              if(idx == 2) idx2 = 0;
            }

            Nurbs* nurbs = new Nurbs;
            bool cricle = true;

            nurbs->arc = cricle;
            nurbs->degree = 2;
            int inner = 1, outer = 0;
            nurbs->np = inner + 2;
            nurbs->pt = new double3[nurbs->np];

            nurbs->pt[0][0] = nodes[p1].x;
            nurbs->pt[0][1] = nodes[p1].y;
            nurbs->pt[0][2] = 1.0;

            nurbs->pt[inner + 1][0] = nodes[p2].x;
            nurbs->pt[inner + 1][1] = nodes[p2].y;
            nurbs->pt[inner + 1][2] = 1.0;

            double angle = angle2;
            double a = (180.0 - angle) / 180.0 * M_PI;
            nurbs->angle = angle;

            // generate one control point
            double x = 1.0 / tan(a * 0.5);
            nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
            nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
            nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

            int i = 0;
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
      Element* sons[3];
      sons[0] = this->create_quad(e->marker, e->vn[0], x0, mid, x2, cm[0]);
      sons[1] = this->create_quad(e->marker, x0, e->vn[1], x1, mid, cm[1]);
      sons[2] = this->create_quad(e->marker, x1, e->vn[2], x2, mid, cm[2]);

      // update coefficients of curved reference mapping
      for (int i = 0; i < 3; i++)
        if(sons[i]->is_curved())
          sons[i]->cm->update_refmap_coeffs(sons[i]);

      // deactivate this element and unregister from its nodes
      e->active = 0;
      if(this != NULL)
      {
        this->nactive += 2;
        e->unref_all_nodes(this);
      }
      // now the original edge nodes may no longer exist...

      // set correct boundary status and markers for the new nodes
      sons[0]->en[0]->bnd = bnd[0];  sons[0]->en[0]->marker = mrk[0];
      sons[0]->en[3]->bnd = bnd[2];  sons[0]->en[3]->marker = mrk[2];
      sons[1]->en[0]->bnd = bnd[0];  sons[1]->en[0]->marker = mrk[0];
      sons[1]->en[1]->bnd = bnd[1];  sons[1]->en[1]->marker = mrk[1];
      sons[2]->en[0]->bnd = bnd[1];  sons[2]->en[0]->marker = mrk[1];
      sons[2]->en[1]->bnd = bnd[2];  sons[2]->en[1]->marker = mrk[2];

      //set pointers to parent element for sons
      for(int i = 0; i < 3; i++)
      {
        if(sons[i] != NULL) sons[i]->parent = e;
      }

      // copy son pointers (could not have been done earlier because of the union)
      memcpy(e->sons, sons, 3 * sizeof(Element*));

      // If sons_out != NULL, copy son pointers there.
      if(sons_out != NULL)
      {
        for(int i = 0; i < 3; i++) sons_out[i] = sons[i];
      }
    }

    void Mesh::refine_element_to_quads_id(int id)
    {
      Element* e = get_element(id);
      if(!e->used) throw Hermes::Exceptions::Exception("Invalid element id number.");
      if(!e->active) throw Hermes::Exceptions::Exception("Attempt to refine element #%d which has been refined already.", e->id);

      if(e->is_triangle())
        refine_triangle_to_quads(e);
      else
        refine_quad_to_quads(e);

      seq = g_mesh_seq++;
    }

    void Mesh::refine_quad_to_triangles(Element* e)
    {
      // remember the markers of the edge nodes
      int bnd[H2D_MAX_NUMBER_EDGES] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd,    e->en[3]->bnd };
      int mrk[H2D_MAX_NUMBER_EDGES] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker, e->en[3]->marker };

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

      if((length_x_0_2 + length_y_0_2) > (length_x_1_3 + length_y_1_3))
      {
        bcheck = false;
      }

      double angle2;
      unsigned int idx;
      CurvMap* cm[2];
      memset(cm, 0, sizeof(cm));

      // create CurvMaps for sons if this is a curved element
      if(e->is_curved())
      {
        int i_case2 = 0;
        if(bcheck == true)
        {
          if((e->cm->nurbs[0] != NULL) || (e->cm->nurbs[1] != NULL))
          {
            cm[0] = new CurvMap;
            memset(cm[0], 0, sizeof(CurvMap));
          }
          if((e->cm->nurbs[2] != NULL) || (e->cm->nurbs[3] != NULL))
          {
            cm[1] = new CurvMap;
            memset(cm[1], 0, sizeof(CurvMap));
          }
        }
        else if(bcheck == false)
        {
          if((e->cm->nurbs[1] != NULL) || (e->cm->nurbs[2] != NULL))
          {
            cm[0] = new CurvMap;
            memset(cm[0], 0, sizeof(CurvMap));
          }
          if((e->cm->nurbs[3] != NULL) || (e->cm->nurbs[0] != NULL))
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
            if(e->cm->nurbs[(idx + i_case2)%4] != NULL)
            {
              angle2 = e->cm->nurbs[(idx + i_case2)%4]->angle;

              int p1, p2;

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

              nurbs->pt[inner + 1][0] = nodes[p2].x;
              nurbs->pt[inner + 1][1] = nodes[p2].y;
              nurbs->pt[inner + 1][2] = 1.0;

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
      Element* sons[H2D_MAX_ELEMENT_SONS];
      if(bcheck == true)
      {
        sons[0] = this->create_triangle(e->marker, e->vn[0], e->vn[1], e->vn[2], cm[0]);
        sons[1] = this->create_triangle(e->marker, e->vn[2], e->vn[3], e->vn[0], cm[1]);
        sons[2] = NULL;
        sons[3] = NULL;
      }
      else
      {
        sons[0] = this->create_triangle(e->marker, e->vn[1], e->vn[2], e->vn[3], cm[0]);
        sons[1] = this->create_triangle(e->marker, e->vn[3], e->vn[0], e->vn[1], cm[1]);
        sons[2] = NULL;
        sons[3] = NULL;
      }

      // update coefficients of curved reference mapping
      for (int i = 0; i < 2; i++)
      {
        if(sons[i]->is_curved())
        {
          sons[i]->cm->update_refmap_coeffs(sons[i]);
        }
      }
      nactive += 2;
      // now the original edge nodes may no longer exist...
      // set correct boundary status and markers for the new nodes
      if(bcheck == true)
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
      for(int i = 0; i < H2D_MAX_ELEMENT_SONS; i++)
        if(sons[i] != NULL)
          sons[i]->parent = e;

      // copy son pointers (could not have been done earlier because of the union)
      memcpy(e->sons, sons, H2D_MAX_ELEMENT_SONS * sizeof(Element*));
    }

    void Mesh::refine_element_to_triangles_id(int id)
    {
      Element* e = get_element(id);
      if(!e->used) throw Hermes::Exceptions::Exception("Invalid element id number.");
      if(!e->active) throw Hermes::Exceptions::Exception("Attempt to refine element #%d which has been refined already.", e->id);

      if(e->is_triangle())
        return;
      else
        refine_quad_to_triangles(e);

      seq = g_mesh_seq++;
    }

    void Mesh::convert_element_to_base_id(int id)
    {
      Element* e = get_element(id);
      if(!e->used) throw Hermes::Exceptions::Exception("Invalid element id number.");
      if(!e->active) throw Hermes::Exceptions::Exception("Attempt to refine element #%d which has been refined already.", e->id);

      if(e->is_triangle())
        convert_triangles_to_base(e);
      else
        convert_quads_to_base(e);// FIXME:

      seq = g_mesh_seq++;
    }

    Mesh::MarkersConversion::MarkersConversion() : min_marker_unused(1)
    {
    }

    Mesh::MarkersConversion::StringValid::StringValid()
    {
    }

    Mesh::MarkersConversion::StringValid::StringValid(std::string marker, bool valid) : marker(marker), valid(valid)
    {
    }

    Mesh::MarkersConversion::IntValid::IntValid()
    {
    }

    Mesh::MarkersConversion::IntValid::IntValid(int marker, bool valid) : marker(marker), valid(valid)
    {
    }

    Mesh::ElementMarkersConversion::ElementMarkersConversion()
    {
    }

    Mesh::MarkersConversion::MarkersConversionType Mesh::ElementMarkersConversion::get_type() const
    {
      return HERMES_ELEMENT_MARKERS_CONVERSION;
    }

    Mesh::BoundaryMarkersConversion::BoundaryMarkersConversion()
    {
    }

    Mesh::MarkersConversion::MarkersConversionType Mesh::BoundaryMarkersConversion::get_type() const
    {
      return HERMES_BOUNDARY_MARKERS_CONVERSION;
    }

    const Mesh::ElementMarkersConversion &Mesh::get_element_markers_conversion() const
    {
      return element_markers_conversion;
    }

    const Mesh::BoundaryMarkersConversion &Mesh::get_boundary_markers_conversion() const
    {
      return boundary_markers_conversion;
    }

    Mesh::ElementMarkersConversion &Mesh::get_element_markers_conversion()
    {
      return element_markers_conversion;
    }

    Mesh::BoundaryMarkersConversion &Mesh::get_boundary_markers_conversion()
    {
      return boundary_markers_conversion;
    }

    int Mesh::MarkersConversion::insert_marker(std::string user_marker)
    {
      // First a check that the string value is not already present.
      std::map<std::string, int>::iterator it = conversion_table_inverse.find(user_marker);
      if(it != conversion_table_inverse.end())
        return it->second;
      conversion_table.insert(std::pair<int, std::string>(this->min_marker_unused, user_marker));
      conversion_table_inverse.insert(std::pair<std::string, int>(user_marker, this->min_marker_unused));
      return this->min_marker_unused++;
    }

    int Mesh::MarkersConversion::size() const
    {
      return this->conversion_table.size();
    }

    Mesh::MarkersConversion::StringValid Mesh::MarkersConversion::get_user_marker(int internal_marker) const
    {
      if(internal_marker == H2D_DG_INNER_EDGE_INT)
        return StringValid(H2D_DG_INNER_EDGE, true);

      std::map<int, std::string>::const_iterator marker = conversion_table.find(internal_marker);
      if(marker == conversion_table.end())
        return StringValid("-999", false);
      else
        return StringValid(marker->second, true);
    }

    Mesh::MarkersConversion::IntValid Mesh::MarkersConversion::get_internal_marker(std::string user_marker) const
    {
      if(user_marker == H2D_DG_INNER_EDGE)
        return IntValid(H2D_DG_INNER_EDGE_INT, true);

      std::map<std::string, int>::const_iterator marker = conversion_table_inverse.find(user_marker);
      if(marker == conversion_table_inverse.end())
        return IntValid(-999, false);
      else
        return IntValid(marker->second, true);
    }

    Mesh::CurvedException::CurvedException(int elementId) : elementId(elementId)
    {
      char * msg = new char[150];
      sprintf(msg, "Element id %i is curved, this is not supported in this method.", elementId);
      message = msg;
    }

    Mesh::CurvedException::CurvedException(const CurvedException & e)
    {
      char * msg = new char[strlen(e.what())+1];
      strcpy(msg, e.what());
      message = msg;
      elementId = e.elementId;
    }

    int Mesh::CurvedException::getElementId() const
    {
      return this->elementId;
    }

    void Mesh::convert_triangles_to_base(Element *e)
    {
      // remember the markers of the edge nodes
      int bnd[3] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd    };
      int mrk[3] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker };

      // check if element e is a internal element.
      bool e_inter = true;
      for (unsigned int n = 0; n < e->get_nvert(); n++)
      {
        if(bnd[n] == 1)
          e_inter = false;
      }

      // get the boundary edge angle.
      double refinement_angle[3] = {0.0, 0.0, 0.0};
      if(e->is_curved() && (!e_inter))
      {
        // for base element.
        if(e->cm->toplevel == true)
        {
          for (unsigned int n = 0; n < e->get_nvert(); n++)
          {
            if(e->cm->nurbs[n] != NULL)
            {
              //this->info("angle = %f", e->cm->nurbs[n]->angle);
              refinement_angle[n] = e->cm->nurbs[n]->angle;
            }
          }
        }
        else
          // one level refinement.
          if(e->parent->cm->toplevel == true)
          {
            for (unsigned int n = 0; n < e->get_nvert(); n++)
            {
              if(e->parent->cm->nurbs[n] != NULL)
              {
                //this->info("angle = %f", e->parent->cm->nurbs[n]->angle);
                refinement_angle[n] = e->parent->cm->nurbs[n]->angle / 2;
              }
            }
          }
          else
            // two level refinements.
            if(e->parent->parent->cm->toplevel == true)
            {
              for (unsigned int n = 0; n < e->get_nvert(); n++)
              {
                if(e->parent->parent->cm->nurbs[n] != NULL)
                {
                  //this->info("angle = %f", e->parent->parent->cm->nurbs[n]->angle);
                  refinement_angle[n] = e->parent->parent->cm->nurbs[n]->angle / 4;
                }
              }
            }
            else
              // three level refinements.
              if(e->parent->parent->parent->cm->toplevel == true)
              {
                for (unsigned int n = 0; n < e->get_nvert(); n++)
                {
                  if(e->parent->parent->parent->cm->nurbs[n] != NULL)
                  {
                    //this->info("angle = %f", e->parent->parent->parent->cm->nurbs[n]->angle);
                    refinement_angle[n] = e->parent->parent->parent->cm->nurbs[n]->angle / 8;
                  }
                }
              }
              else
                // four level refinements.
                if(e->parent->parent->parent->parent->cm->toplevel == true)
                {
                  for (unsigned int n = 0; n < e->get_nvert(); n++)
                  {
                    if(e->parent->parent->parent->parent->cm->nurbs[n] != NULL)
                    {
                      //this->info("angle = %f", e->parent->parent->parent->parent->cm->nurbs[n]->angle);
                      refinement_angle[n] = e->parent->parent->parent->parent->cm->nurbs[n]->angle / 16;
                    }
                  }
                }
      }

      // deactivate this element and unregister from its nodes
      e->active = false;
      e->unref_all_nodes(this);

      double angle2 = 0.0;
      int idx = 0;
      CurvMap* cm;
      memset(&cm, 0, sizeof(cm));

      // create CurvMaps for sons if this is a curved element
      if((e->is_curved()) && (!e_inter))
      {
        cm = new CurvMap;
        memset(cm, 0, sizeof(CurvMap));

        for (idx = 0; idx < 3; idx++)
          if((e->cm->nurbs[idx] != NULL) && (bnd[idx] == 1))
          {
            angle2 = refinement_angle[idx];
            int p1, p2;
            p1 = e->en[idx]->p1;
            p2 = e->en[idx]->p2;
            if(p1 > p2) std::swap(p1, p2);

            Nurbs* nurbs = new Nurbs;
            bool cricle = true;

            nurbs->arc = cricle;
            nurbs->degree = 2;

            int inner = 1, outer = 0;
            nurbs->np = inner + 2;
            nurbs->pt = new double3[nurbs->np];

            nurbs->pt[0][0] = nodes[p1].x;
            nurbs->pt[0][1] = nodes[p1].y;
            nurbs->pt[0][2] = 1.0;

            nurbs->pt[inner + 1][0] = nodes[p2].x;
            nurbs->pt[inner + 1][1] = nodes[p2].y;
            nurbs->pt[inner + 1][2] = 1.0;

            double angle = angle2;
            double a = (180.0 - angle) / 180.0 * M_PI;
            nurbs->angle = angle;

            // generate one control point
            double x = 1.0 / tan(a * 0.5);
            nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
            nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
            nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

            int i = 0;
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

            cm->toplevel = 1;
            cm->order = 4;
            cm->nurbs[idx] = nurbs;
            nurbs->ref++;
          }
      }

      // create a new element.
      Element* enew;
      Node *v0 = &nodes[e->vn[0]->id], *v1 = &nodes[e->vn[1]->id], *v2 = &nodes[e->vn[2]->id];
      enew = this->create_triangle(e->marker, v0, v1, v2, cm);

      if(enew->is_curved())
      {
        enew->cm->update_refmap_coeffs(enew);
      }

      // now the original edge nodes may no longer exist...
      // set correct boundary status and markers for the new nodes
      enew->en[0]->bnd = bnd[0];
      enew->en[1]->bnd = bnd[1];
      enew->en[2]->bnd = bnd[2];
      enew->en[0]->marker = mrk[0];
      enew->en[1]->marker = mrk[1];
      enew->en[2]->marker = mrk[2];
      enew->parent = e;
    }

    void Mesh::convert_quads_to_base(Element *e)
    {
      // remember the markers of the edge nodes
      int bnd[H2D_MAX_NUMBER_EDGES] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd,    e->en[3]->bnd    };
      int mrk[H2D_MAX_NUMBER_EDGES] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker, e->en[3]->marker };

      // check if element e is a internal element.
      bool e_inter = true;
      for (unsigned int n = 0; n < e->get_nvert(); n++)
      {
        if(bnd[n] == 1)
          e_inter = false;
      }

      // get the boundary edge angle.
      double refinement_angle[H2D_MAX_NUMBER_EDGES] = {0.0, 0.0, 0.0, 0.0};
      if(e->is_curved() && (!e_inter))
      {
        // for base element.
        if(e->cm->toplevel == true)
        {
          for (unsigned int n = 0; n < e->get_nvert(); n++)
          {
            if((e->cm->nurbs[n] != NULL) && (bnd[n] == 1))
            {
              //this->info("angle = %f", e->cm->nurbs[n]->angle);
              refinement_angle[n] = e->cm->nurbs[n]->angle;
            }
          }
        }
        else
          // one level refinement.
          if(e->parent->cm->toplevel == true)
          {
            for (unsigned int n = 0; n < e->get_nvert(); n++)
            {
              if((e->parent->cm->nurbs[n] != NULL) && (bnd[n] == 1))
              {
                //this->info("angle = %f", e->parent->cm->nurbs[n]->angle);
                refinement_angle[n] = e->parent->cm->nurbs[n]->angle / 2;
              }
            }
          }
          else
            // two level refinements.
            if(e->parent->parent->cm->toplevel == true)
            {
              for (unsigned int n = 0; n < e->get_nvert(); n++)
              {
                if((e->parent->parent->cm->nurbs[n] != NULL) && (bnd[n] == 1))
                {
                  //this->info("angle = %f", e->parent->parent->cm->nurbs[n]->angle);
                  refinement_angle[n] = e->parent->parent->cm->nurbs[n]->angle / 4;
                }
              }
            }
            else
              // three level refinements.
              if(e->parent->parent->parent->cm->toplevel == true)
              {
                for (unsigned int n = 0; n < e->get_nvert(); n++)
                {
                  if((e->parent->parent->parent->cm->nurbs[n] != NULL) && (bnd[n] == 1))
                  {
                    //this->info("angle = %f", e->parent->parent->parent->cm->nurbs[n]->angle);
                    refinement_angle[n] = e->parent->parent->parent->cm->nurbs[n]->angle / 8;
                  }
                }
              }
              else
                // four level refinements.
                if(e->parent->parent->parent->parent->cm->toplevel == true)
                {
                  for (unsigned int n = 0; n < e->get_nvert(); n++)
                  {
                    if((e->parent->parent->parent->parent->cm->nurbs[n] != NULL) && (bnd[n] == 1))
                    {
                      //this->info("angle = %f", e->parent->parent->parent->parent->cm->nurbs[n]->angle);
                      refinement_angle[n] = e->parent->parent->parent->parent->cm->nurbs[n]->angle / 16;
                    }
                  }
                }
      }

      // FIXME:
      if(rtb_aniso)
        for (unsigned int i = 0; i < e->get_nvert(); i++)
          refinement_angle[i] = refinement_angle[i]*2;

      // deactivate this element and unregister from its nodes
      e->active = false;
      e->unref_all_nodes(this);

      double angle2 = 0.0;
      int idx = 0;
      CurvMap* cm;
      memset(&cm, 0, sizeof(cm));

      // create CurvMaps for sons if this is a curved element
      if((e->is_curved()) && (!e_inter))
      {
        bool create_new = false;
        for (unsigned int i = 0; i < e->get_nvert(); i++)
        {
          if(fabs(refinement_angle[i] - 0.0) > 1e-4)
          {
            create_new = true;
          }
        }

        if(create_new)

        {
          cm = new CurvMap;
          memset(cm, 0, sizeof(CurvMap));
        }

        for (idx = 0; idx < 4; idx++)
          if(fabs(refinement_angle[idx] - 0.0) > 1e-4)
          {
            angle2 = refinement_angle[idx];
            int p1, p2;
            p1 = e->en[idx]->p1;
            p2 = e->en[idx]->p2;
            if(p1 > p2) std::swap(p1, p2);

            Nurbs* nurbs = new Nurbs;
            bool cricle = true;

            nurbs->arc = cricle;
            nurbs->degree = 2;

            int inner = 1, outer = 0;
            nurbs->np = inner + 2;
            nurbs->pt = new double3[nurbs->np];

            nurbs->pt[0][0] = nodes[p1].x;
            nurbs->pt[0][1] = nodes[p1].y;
            nurbs->pt[0][2] = 1.0;

            nurbs->pt[inner + 1][0] = nodes[p2].x;
            nurbs->pt[inner + 1][1] = nodes[p2].y;
            nurbs->pt[inner + 1][2] = 1.0;

            double angle = angle2;
            double a = (180.0 - angle) / 180.0 * M_PI;
            nurbs->angle = angle;

            // generate one control point
            double x = 1.0 / tan(a * 0.5);
            nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
            nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
            nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

            int i = 0;
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

            cm->toplevel = 1;
            cm->order = 4;
            cm->nurbs[idx] = nurbs;
            nurbs->ref++;
          }
      }

      // create a new element.
      Element* enew;
      Node *v0 = &nodes[e->vn[0]->id], *v1 = &nodes[e->vn[1]->id], *v2 = &nodes[e->vn[2]->id],  *v3 = &nodes[e->vn[3]->id];
      enew = this->create_quad(e->marker, v0, v1, v2, v3, cm);

      if(enew->is_curved())
      {
        enew->cm->update_refmap_coeffs(enew);
      }

      // now the original edge nodes may no longer exist...
      // set correct boundary status and markers for the new nodes
      enew->en[0]->bnd = bnd[0];
      enew->en[1]->bnd = bnd[1];
      enew->en[2]->bnd = bnd[2];
      enew->en[3]->bnd = bnd[3];
      enew->en[0]->marker = mrk[0];
      enew->en[1]->marker = mrk[1];
      enew->en[2]->marker = mrk[2];
      enew->en[3]->marker = mrk[3];
      enew->parent = e;
    }

    void Mesh::refine_quad_to_quads(Element* e, int refinement)
    {
      // remember the markers of the edge nodes
      int bnd[H2D_MAX_NUMBER_EDGES] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd,    e->en[3]->bnd    };
      int mrk[H2D_MAX_NUMBER_EDGES] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker, e->en[3]->marker };

      // check if element e is a internal element.
      bool e_inter = true;
      for (unsigned int n = 0; n < e->get_nvert(); n++)
      {
        if(bnd[n] == 1)
          e_inter = false;
      }

      // get the boundary edge angle.
      double refinement_angle[H2D_MAX_NUMBER_EDGES] = {0.0, 0.0, 0.0, 0.0};
      if(e->is_curved() && (!e_inter))
      {
        // for base element.
        if(e->cm->toplevel == true)
        {
          for (unsigned int n = 0; n < e->get_nvert(); n++)
          {
            if((e->cm->nurbs[n] != NULL) && (bnd[n] == 1))
            {
              //this->info("angle = %f", e->cm->nurbs[n]->angle);
              refinement_angle[n] = e->cm->nurbs[n]->angle;
            }
          }
        }
        else
          // one level refinement.
          if(e->parent->cm->toplevel == true)
          {
            for (unsigned int n = 0; n < e->get_nvert(); n++)
            {
              if((e->parent->cm->nurbs[n] != NULL) && (bnd[n] == 1))
              {
                //this->info("angle = %f", e->parent->cm->nurbs[n]->angle);
                refinement_angle[n] = e->parent->cm->nurbs[n]->angle / 2;
              }
            }
          }
          else
            // two level refinements.
            if(e->parent->parent->cm->toplevel == true)
            {
              for (unsigned int n = 0; n < e->get_nvert(); n++)
              {
                if((e->parent->parent->cm->nurbs[n] != NULL) && (bnd[n] == 1))
                {
                  //this->info("angle = %f", e->parent->parent->cm->nurbs[n]->angle);
                  refinement_angle[n] = e->parent->parent->cm->nurbs[n]->angle / 4;
                }
              }
            }
            else
              // three level refinements.
              if(e->parent->parent->parent->cm->toplevel == true)
              {
                for (unsigned int n = 0; n < e->get_nvert(); n++)
                {
                  if((e->parent->parent->parent->cm->nurbs[n] != NULL) && (bnd[n] == 1))
                  {
                    //this->info("angle = %f", e->parent->parent->parent->cm->nurbs[n]->angle);
                    refinement_angle[n] = e->parent->parent->parent->cm->nurbs[n]->angle / 8;
                  }
                }
              }
              else
                // four level refinements.
                if(e->parent->parent->parent->parent->cm->toplevel == true)
                {
                  for (unsigned int n = 0; n < e->get_nvert(); n++)
                  {
                    if((e->parent->parent->parent->parent->cm->nurbs[n] != NULL) && (bnd[n] == 1))
                    {
                      //this->info("angle = %f", e->parent->parent->parent->parent->cm->nurbs[n]->angle);
                      refinement_angle[n] = e->parent->parent->parent->parent->cm->nurbs[n]->angle / 16;
                    }
                  }
                }
      }

      // deactivate this element and unregister from its nodes
      e->active = false;
      this->nactive--;
      e->unref_all_nodes(this);

      double angle2 = 0.0;
      int i, j;
      int idx = 0;
      Element* sons[H2D_MAX_ELEMENT_SONS] = {NULL, NULL, NULL, NULL};
      CurvMap* cm[H2D_MAX_ELEMENT_SONS];
      memset(cm, 0, sizeof(cm));

      // default refinement: one quad to four quads
      if(refinement == 0)
      {
        // obtain four mid-edge vertex nodes and one mid-element vetex node
        Node* x0 = get_vertex_node(e->vn[0]->id, e->vn[1]->id);
        Node* x1 = get_vertex_node(e->vn[1]->id, e->vn[2]->id);
        Node* x2 = get_vertex_node(e->vn[2]->id, e->vn[3]->id);
        Node* x3 = get_vertex_node(e->vn[3]->id, e->vn[0]->id);
        Node* mid = get_vertex_node(x0->id, x2->id);

        // adjust mid-edge coordinates if this is a curved element
        if(e->is_curved())
        {
          double2 pt[5] = { { 0.0, -1.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }, { -1.0, 0.0 }, { 0.0, 0.0 } };
          e->cm->get_mid_edge_points(e, pt, 5);
          x0->x = pt[0][0];  x0->y = pt[0][1];
          x1->x = pt[1][0];  x1->y = pt[1][1];
          x2->x = pt[2][0];  x2->y = pt[2][1];
          x3->x = pt[3][0];  x3->y = pt[3][1];
          mid->x = pt[4][0]; mid->y = pt[4][1];
        }

        // create CurvMaps for sons.
        if((e->is_curved()) && (!e_inter))
        {
          //bool create_new = false;
          for (unsigned int i = 0; i < e->get_nvert(); i++)
          {
            if(fabs(refinement_angle[i] - 0.0) > 1e-4)
            {
              cm[i%4] = new CurvMap;
              memset(cm[i%4], 0, sizeof(CurvMap));
              cm[(i + 1)%4] = new CurvMap;
              memset(cm[(i + 1)%4], 0, sizeof(CurvMap));
            }
          }

          for (idx = 0; idx < 4; idx++)
            if(cm[idx] != NULL)
            {
              if((fabs(refinement_angle[idx%4] - 0.0) > 1e-4))
              {
                angle2 = refinement_angle[idx%4] / 2;
                Node* node_temp = this->get_vertex_node(e->vn[idx%4]->id, e->vn[(idx + 1)%4]->id);

                int p1, p2;
                p1 = e->vn[(idx)%4]->id;
                p2 = node_temp->id;

                Nurbs* nurbs = new Nurbs;
                bool cricle = true;

                nurbs->arc = cricle;
                nurbs->degree = 2;

                int inner = 1, outer = 0;
                nurbs->np = inner + 2;
                nurbs->pt = new double3[nurbs->np];

                nurbs->pt[0][0] = nodes[p1].x;
                nurbs->pt[0][1] = nodes[p1].y;
                nurbs->pt[0][2] = 1.0;

                nurbs->pt[inner + 1][0] = nodes[p2].x;
                nurbs->pt[inner + 1][1] = nodes[p2].y;
                nurbs->pt[inner + 1][2] = 1.0;

                double angle = angle2;
                double a = (180.0 - angle) / 180.0 * M_PI;
                nurbs->angle = angle;

                // generate one control point
                double x = 1.0 / tan(a * 0.5);
                nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
                nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
                nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

                int i = 0;
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
                cm[idx]->nurbs[idx%4] = nurbs;
                nurbs->ref++;
              }

              if((fabs(refinement_angle[(idx + 3)%4] - 0.0) > 1e-4))
              {
                angle2 = refinement_angle[(idx + 3)%4]/2;
                Node* node_temp = this->get_vertex_node(e->vn[idx%4]->id, e->vn[(idx + 1)%4]->id);

                int p1, p2;
                p1 = e->vn[(idx)%4]->id;
                p2 = node_temp->id;

                Nurbs* nurbs = new Nurbs;
                bool cricle = true;

                nurbs->arc = cricle;
                nurbs->degree = 2;

                int inner = 1, outer = 0;
                nurbs->np = inner + 2;
                nurbs->pt = new double3[nurbs->np];

                nurbs->pt[0][0] = nodes[p1].x;
                nurbs->pt[0][1] = nodes[p1].y;
                nurbs->pt[0][2] = 1.0;

                nurbs->pt[inner + 1][0] = nodes[p2].x;
                nurbs->pt[inner + 1][1] = nodes[p2].y;
                nurbs->pt[inner + 1][2] = 1.0;

                double angle = angle2;
                double a = (180.0 - angle) / 180.0 * M_PI;
                nurbs->angle = angle;

                // generate one control point
                double x = 1.0 / tan(a * 0.5);
                nurbs->pt[1][0] = 0.5*((nurbs->pt[2][0] + nurbs->pt[0][0]) + (nurbs->pt[2][1] - nurbs->pt[0][1]) * x);
                nurbs->pt[1][1] = 0.5*((nurbs->pt[2][1] + nurbs->pt[0][1]) - (nurbs->pt[2][0] - nurbs->pt[0][0]) * x);
                nurbs->pt[1][2] = cos((M_PI - a) * 0.5);

                int i = 0;
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
                cm[idx]->nurbs[(idx + 3)%4] = nurbs;
                nurbs->ref++;
              }
            }
        }

        // create the four sons
        sons[0] = this->create_quad(e->marker, e->vn[0], x0, mid, x3, cm[0]);
        sons[1] = this->create_quad(e->marker, x0, e->vn[1], x1, mid, cm[1]);
        sons[2] = this->create_quad(e->marker, mid, x1, e->vn[2], x2, cm[2]);
        sons[3] = this->create_quad(e->marker, x3, mid, x2, e->vn[3], cm[3]);

        // Increase the number of active elements by 4.
        this->nactive += 4;

        // set correct boundary markers for the new edge nodes
        for (i = 0; i < 4; i++)
        {
          j = (i > 0) ? i-1 : 3;
          sons[i]->en[j]->bnd = bnd[j];  sons[i]->en[j]->marker = mrk[j];
          sons[i]->en[i]->bnd = bnd[i];  sons[i]->en[i]->marker = mrk[i];
          sons[i]->vn[j]->bnd = bnd[j];
        }
      }
      else assert(0);

      // update coefficients of curved reference mapping
      for (i = 0; i < 4; i++)
        if(sons[i] != NULL && sons[i]->cm != NULL)
          sons[i]->cm->update_refmap_coeffs(sons[i]);

      //set pointers to parent element for sons
      for(int i = 0; i < 4; i++)
        if(sons[i] != NULL)
          sons[i]->parent = e;

      // copy son pointers (could not have been done earlier because of the union)
      memcpy(e->sons, sons, sizeof(sons));
    }

    int Mesh::get_edge_degree(Node* v1, Node* v2)
    {
      int degree = 0;
      Node* v3 = peek_vertex_node(v1->id, v2->id);
      if(v3 != NULL)
      {
        degree = 1 + std::max(get_edge_degree(v1, v3), get_edge_degree(v3, v2));
      }
      return degree;
    }

    void Mesh::regularize_triangle(Element* e)
    {
      int i, k, k1, k2;

      Element* t[3];

      int eo[3] = { get_edge_degree(e->vn[0], e->vn[1]),
        get_edge_degree(e->vn[1], e->vn[2]),
        get_edge_degree(e->vn[2], e->vn[0]) };

      int sum = eo[0] + eo[1] + eo[2];
      if(sum == 3)
      {
        refine_element_id(e->id);
      }
      else if(sum > 0)
      {
        // remember the markers of the edge nodes
        int bnd[3] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd    };
        int mrk[3] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker };

        if(sum == 1)
        {
          Node* v4;
          for(i = 0; i < 3; i++)
            if(eo[i] == 1) k = i;
          k1 = e->next_vert(k);
          k2 = e->prev_vert(k);
          v4 = peek_vertex_node(e->vn[k]->id, e->vn[k1]->id);

          e->active = 0;
          nactive += 1;
          e->unref_all_nodes(this);

          t[0] = this->create_triangle(e->marker, e->vn[k], v4, e->vn[k2], NULL);
          t[1] = this->create_triangle(e->marker, v4, e->vn[k1], e->vn[k2], NULL);

          // set correct boundary status and markers for the new nodes
          t[0]->en[2]->bnd = bnd[k2];
          t[1]->en[1]->bnd = bnd[k1];
          t[0]->en[2]->marker = mrk[k2];
          t[1]->en[1]->marker = mrk[k1];

          e->sons[0] = t[0];
          e->sons[1] = t[1];
          e->sons[2] = NULL;
          e->sons[3] = NULL;
        }
        else if(sum == 2)
        {
          Node *v4, *v5;
          for(i = 0; i < 3; i++)
            if(eo[i] == 0) k = i;
          k1 = e->next_vert(k);
          k2 = e->prev_vert(k);
          v4 = peek_vertex_node(e->vn[k1]->id, e->vn[k2]->id);
          v5 = peek_vertex_node(e->vn[k2]->id, e->vn[k]->id);

          e->active = 0;
          nactive += 2;
          e->unref_all_nodes(this);

          t[0] = this->create_triangle(e->marker, e->vn[k], e->vn[k1], v4,  NULL);
          t[1] = this->create_triangle(e->marker, v4, v5, e->vn[k], NULL);
          t[2] = this->create_triangle(e->marker, v4, e->vn[k2], v5, NULL);

          // set correct boundary status and markers for the new nodes
          t[0]->en[0]->bnd = bnd[k];
          t[0]->en[0]->marker = mrk[k];

          e->sons[0] = t[0];
          e->sons[1] = t[1];
          e->sons[2] = t[2];
          e->sons[3] = NULL;
        }
      }

      // store id of parent
      if(!e->active)
      {
        for (i = 0; i < 4; i++)
          assign_parent(e, i);
      }
    }

    void Mesh::regularize_quad(Element* e)
    {
      int i, k = 0, k1, k2, k3, n = 0, m = 0;
      Node *v4, *v5;
      Element* t[4];

      int eo[4] = { get_edge_degree(e->vn[0], e->vn[1]),
        get_edge_degree(e->vn[1], e->vn[2]),
        get_edge_degree(e->vn[2], e->vn[3]),
        get_edge_degree(e->vn[3], e->vn[0]) };

      int sum = eo[0] + eo[1] + eo[2] + eo[3];
      if(sum == 4)
      {
        refine_element_id(e->id);
      }
      else if(sum > 0)
      {
        // remember the markers of the edge nodes
        int bnd[H2D_MAX_NUMBER_EDGES] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd  ,  e->en[3]->bnd  };
        int mrk[H2D_MAX_NUMBER_EDGES] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker, e->en[3]->marker };

        if(sum == 1)
        {
          for(i = 0; i < 4; i++)
            if(eo[i] == 1) k = i;
          k1 = e->next_vert(k);
          k2 = e->next_vert(k1);
          k3 = e->prev_vert(k);
          v4 = peek_vertex_node(e->vn[k]->id, e->vn[k1]->id);

          e->active = 0;
          nactive += 2;
          e->unref_all_nodes(this);

          t[0] = this->create_triangle(e->marker, e->vn[k], v4, e->vn[k3], NULL);
          t[1] = this->create_triangle(e->marker, v4, e->vn[k1], e->vn[k2], NULL);
          t[2] = this->create_triangle(e->marker, v4, e->vn[k2], e->vn[k3], NULL);

          // set correct boundary status and markers for the new nodes
          t[0]->en[2]->bnd = bnd[k3];
          t[1]->en[1]->bnd = bnd[k1];
          t[2]->en[1]->bnd = bnd[k2];
          t[0]->en[2]->marker = mrk[k3];
          t[1]->en[1]->marker = mrk[k1];
          t[2]->en[1]->marker = mrk[k2];

          e->sons[0] = t[0];
          e->sons[1] = t[1];
          e->sons[2] = t[2];
          e->sons[3] = NULL;
        }
        else if(sum == 2)
        {
          // two hanging nodes opposite to each other
          if(eo[0] == 1 && eo[2] == 1) refine_element_id(e->id, 2);
          else if(eo[1] == 1 && eo[3] == 1) refine_element_id(e->id, 1);
          else // two hanging nodes next to each other
          {
            for(i = 0; i < 4; i++)
              if(eo[i] == 1 && eo[e->next_vert(i)] == 1) k = i;
            k1 = e->next_vert(k);
            k2 = e->next_vert(k1);
            k3 = e->prev_vert(k);
            v4 = peek_vertex_node(e->vn[k]->id, e->vn[k1]->id);
            v5 = peek_vertex_node(e->vn[k1]->id, e->vn[k2]->id);

            e->active = 0;
            nactive += 3;
            e->unref_all_nodes(this);

            t[0] = this->create_triangle(e->marker, e->vn[k1], v5, v4, NULL);
            t[1] = this->create_triangle(e->marker, v5, e->vn[k2], e->vn[k3], NULL);
            t[2] = this->create_triangle(e->marker, v4, v5, e->vn[k3], NULL);
            t[3] = this->create_triangle(e->marker, v4, e->vn[k3], e->vn[k], NULL);

            t[1]->en[1]->bnd = bnd[k2];
            t[3]->en[1]->bnd = bnd[k3];
            t[1]->en[1]->marker = mrk[k2];
            t[3]->en[1]->marker = mrk[k3];

            e->sons[0] = t[0];
            e->sons[1] = t[1];
            e->sons[2] = t[2];
            e->sons[3] = t[3];
          }
        }
        else //sum = 3
        {
          if(eo[0] == 1 && eo[2] == 1)
          {
            refine_element_id(e->id, 2);
            for (i = 0; i < 4; i++)
              assign_parent(e, i);
            n = 2; m = 3;
          }
          else if(eo[1] == 1 && eo[3] == 1)
          {
            refine_element_id(e->id, 1);
            for (i = 0; i < 4; i++)
              assign_parent(e, i);
            n = 0; m = 1;
          }

          regularize_quad(e->sons[n]);
          regularize_quad(e->sons[m]);
        }
      }

      // store id of parent
      if(!e->active)
      {
        for (i = 0; i < 4; i++)
          assign_parent(e, i);
      }
    }

    void Mesh::flatten()
    {
      Node* node;
      for_all_edge_nodes(node, this)
      {
        if(node->elem[0] != NULL) node->elem[0] = (Element*) (node->elem[0]->id + 1);
        if(node->elem[1] != NULL) node->elem[1] = (Element*) (node->elem[1]->id + 1);
      }

      int* idx = new int[elements.get_size() + 1];
      Array<Element> new_elements;
      Element* e;
      for_all_active_elements(e, this)
      {
        Element* ee = new_elements.add();
        int temp = ee->id;
        *ee = *e;
        ee->id = temp;
        idx[e->id] = temp;
        parents[ee->id] = parents[e->id];
      }

      elements.copy(new_elements);
      nbase = nactive = elements.get_num_items();

      for_all_edge_nodes(node, this)
      {
        if(node->elem[0] != NULL) node->elem[0] = &(elements[idx[((int) (long) node->elem[0]) - 1]]);
        if(node->elem[1] != NULL) node->elem[1] = &(elements[idx[((int) (long) node->elem[1]) - 1]]);
      }

      delete [] idx;
    }

    void Mesh::assign_parent(Element* e, int i)
    {
      if(e->sons[i] != NULL)
      {
        if(e->sons[i]->id >= parents_size)
        {
          parents_size = 2 * parents_size;
          parents = (int*) realloc(parents, sizeof(int) * parents_size);
        }

        parents[e->sons[i]->id] = parents[e->id];
      }
    }

    int* Mesh::regularize(int n)
    {
      int j;
      bool ok;
      bool reg = false;
      Element* e;

      if(n < 1)
      {
        n = 1;
        reg = true;
      }

      parents_size = 2*get_max_element_id();
      parents = (int*) malloc(sizeof(int) * parents_size);
      for_all_active_elements(e, this)
        parents[e->id] = e->id;

      do
      {
        ok = true;
        for_all_active_elements(e, this)
        {
          int iso = -1;
          if(e->is_triangle())
          {
            for(unsigned int i = 0; i < e->get_nvert(); i++)
            {
              j = e->next_vert(i);
              if(get_edge_degree(e->vn[i], e->vn[j]) > n)
              { iso = 0; ok = false; break; }
            }
          }
          else
          {
            if(   ((get_edge_degree(e->vn[0], e->vn[1]) > n)  || (get_edge_degree(e->vn[2], e->vn[3]) > n))
              && (get_edge_degree(e->vn[1], e->vn[2]) <= n) && (get_edge_degree(e->vn[3], e->vn[0]) <= n) )
            { iso = 2; ok = false; }
            else if(    (get_edge_degree(e->vn[0], e->vn[1]) <= n)  && (get_edge_degree(e->vn[2], e->vn[3]) <= n)
              && ((get_edge_degree(e->vn[1], e->vn[2]) > n) || (get_edge_degree(e->vn[3], e->vn[0]) > n)) )
            { iso = 1; ok = false; }
            else
            {
              for(unsigned int i = 0; i < e->get_nvert(); i++)
              {
                j = e->next_vert(i);
                if(get_edge_degree(e->vn[i], e->vn[j]) > n)
                { iso = 0; ok = false; break; }
              }
            }
          }

          if(iso >= 0)
          {
            refine_element_id(e->id, iso);
            for (int i = 0; i < 4; i++)
              assign_parent(e, i);
          }
        }
      }
      while (!ok);

      if(reg)
      {
        for_all_active_elements(e, this)
        {
          if(e->is_curved()) throw Hermes::Exceptions::Exception("Regularization of curved elements is not supported.");

          if(e->is_triangle())
            regularize_triangle(e);
          else
            regularize_quad(e);
        }
        flatten();
      }

      return parents;
    }
  }
}
