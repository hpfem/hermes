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


int Mesh::get_edge_degree(Node* v1, Node* v2)
{
  int degree = 0;
  Node* v3 = peek_vertex_node(v1->id, v2->id);
  if (v3 != NULL)
  {
    degree = 1 + std::max(get_edge_degree(v1,v3), get_edge_degree(v3,v2));
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
  if (sum == 3)
  {
    refine_element(e->id);
  }
  else if (sum > 0)
  {
    // remember the markers of the edge nodes
    int bnd[3] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd    };
    int mrk[3] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker };

    if (sum == 1)
    {
      Node* v4;
      for(i = 0; i < 3; i++)
        if (eo[i] == 1) k = i;
      k1 = e->next_vert(k);
      k2 = e->prev_vert(k);
      v4 = peek_vertex_node(e->vn[k]->id, e->vn[k1]->id);

      e->active = 0;
      nactive += 1;
      e->unref_all_nodes(this);

      t[0] = create_triangle(e->marker, e->vn[k], v4, e->vn[k2], NULL);
      t[1] = create_triangle(e->marker, v4, e->vn[k1], e->vn[k2], NULL);

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
    else if (sum == 2)
    {
      Node *v4, *v5;
      for(i = 0; i < 3; i++)
        if (eo[i] == 0) k = i;
      k1 = e->next_vert(k);
      k2 = e->prev_vert(k);
      v4 = peek_vertex_node(e->vn[k1]->id, e->vn[k2]->id);
      v5 = peek_vertex_node(e->vn[k2]->id, e->vn[k]->id);

      e->active = 0;
      nactive += 2;
      e->unref_all_nodes(this);

      t[0] = create_triangle(e->marker, e->vn[k], e->vn[k1], v4,  NULL);
      t[1] = create_triangle(e->marker, v4, v5, e->vn[k], NULL);
      t[2] = create_triangle(e->marker, v4, e->vn[k2], v5, NULL);

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
  if (!e->active)
  {
    for (i = 0; i < 4; i++)
      assign_parent(e, i);
  }
}


void Mesh::regularize_quad(Element* e)
{
  int i, k, k1, k2, k3, n, m;
 Node *v4, *v5;
  Element* t[4];

  int eo[4] = { get_edge_degree(e->vn[0], e->vn[1]),
                get_edge_degree(e->vn[1], e->vn[2]),
                get_edge_degree(e->vn[2], e->vn[3]),
                get_edge_degree(e->vn[3], e->vn[0]) };

  int sum = eo[0] + eo[1] + eo[2] + eo[3];
  if (sum == 4)
  {
    refine_element(e->id);
  }
  else if (sum > 0)
  {
    // remember the markers of the edge nodes
    int bnd[4] = { e->en[0]->bnd,    e->en[1]->bnd,    e->en[2]->bnd  ,  e->en[3]->bnd  };
    int mrk[4] = { e->en[0]->marker, e->en[1]->marker, e->en[2]->marker, e->en[3]->marker };

    if (sum == 1)
    {
      for(i = 0; i < 4; i++)
        if (eo[i] == 1) k = i;
      k1 = e->next_vert(k);
      k2 = e->next_vert(k1);
      k3 = e->prev_vert(k);
      v4 = peek_vertex_node(e->vn[k]->id, e->vn[k1]->id);

      e->active = 0;
      nactive += 2;
      e->unref_all_nodes(this);

      t[0] = create_triangle(e->marker, e->vn[k], v4, e->vn[k3], NULL);
      t[1] = create_triangle(e->marker, v4, e->vn[k1], e->vn[k2], NULL);
      t[2] = create_triangle(e->marker, v4, e->vn[k2], e->vn[k3], NULL);

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
    else if (sum == 2)
    {
      // two hanging nodes opposite to each other
      if (eo[0] == 1 && eo[2] == 1) refine_element(e->id, 2);
      else if (eo[1] == 1 && eo[3] == 1) refine_element(e->id, 1);
      else // two hanging nodes next to each other
      {
        for(i = 0; i < 4; i++)
          if (eo[i] == 1 && eo[e->next_vert(i)] == 1) k = i;
        k1 = e->next_vert(k);
        k2 = e->next_vert(k1);
        k3 = e->prev_vert(k);
        v4 = peek_vertex_node(e->vn[k]->id, e->vn[k1]->id);
        v5 = peek_vertex_node(e->vn[k1]->id, e->vn[k2]->id);

        e->active = 0;
        nactive += 3;
        e->unref_all_nodes(this);

        t[0] = create_triangle(e->marker, e->vn[k1], v5, v4, NULL);
        t[1] = create_triangle(e->marker, v5, e->vn[k2], e->vn[k3], NULL);
        t[2] = create_triangle(e->marker, v4, v5, e->vn[k3], NULL);
        t[3] = create_triangle(e->marker, v4, e->vn[k3], e->vn[k], NULL);

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
       if (eo[0] == 1 && eo[2] == 1)
       {
         refine_element(e->id, 2);
         for (i = 0; i < 4; i++)
           assign_parent(e, i);
         n = 2; m = 3;
       }
       else if (eo[1] == 1 && eo[3] == 1)
       {
         refine_element(e->id, 1);
         for (i = 0; i < 4; i++)
           assign_parent(e, i);
         n = 0; m = 1;
       }

       regularize_quad(e->sons[n]);
       regularize_quad(e->sons[m]);

    }

  }

  // store id of parent
  if (!e->active)
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
    if (node->elem[0] != NULL) node->elem[0] = (Element*) (node->elem[0]->id + 1);
    if (node->elem[1] != NULL) node->elem[1] = (Element*) (node->elem[1]->id + 1);
  }

  AUTOLA_OR(int, idx, elements.get_size()+1);
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
    if (node->elem[0] != NULL) node->elem[0] = &(elements[idx[((int) (long) node->elem[0]) - 1]]);
    if (node->elem[1] != NULL) node->elem[1] = &(elements[idx[((int) (long) node->elem[1]) - 1]]);
  }
}


void Mesh::assign_parent(Element* e, int i)
{
  if (e->sons[i] != NULL)
  {
    if (e->sons[i]->id >= parents_size)
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
  int iso = 0;
  Element* e;

  if (n < 1)
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
      if (e->is_triangle())
      {
        for(unsigned int i = 0; i < e->nvert; i++)
        {
          j = e->next_vert(i);
          if (get_edge_degree(e->vn[i], e->vn[j]) > n)
            { iso = 0; ok = false; break; }
        }
      }
      else
      {
        if (   ((get_edge_degree(e->vn[0], e->vn[1]) > n)  || (get_edge_degree(e->vn[2], e->vn[3]) > n))
            && (get_edge_degree(e->vn[1], e->vn[2]) <= n) && (get_edge_degree(e->vn[3], e->vn[0]) <= n) )
          { iso = 2; ok = false; }
        else if (    (get_edge_degree(e->vn[0], e->vn[1]) <= n)  && (get_edge_degree(e->vn[2], e->vn[3]) <= n)
                  && ((get_edge_degree(e->vn[1], e->vn[2]) > n) || (get_edge_degree(e->vn[3], e->vn[0]) > n)) )
          { iso = 1; ok = false; }
        else
        {
          for(unsigned int i = 0; i < e->nvert; i++)
          {
            j = e->next_vert(i);
            if (get_edge_degree(e->vn[i], e->vn[j]) > n)
              { iso = 0; ok = false; break; }
          }
        }
      }

      if (iso >= 0)
      {
        refine_element(e->id, iso);
        for (int i = 0; i < 4; i++)
          assign_parent(e, i);
      }
    }
  }
  while (!ok);


  if (reg)
  {
    for_all_active_elements(e,this)
    {
      if (e->is_curved()) error("Regularization of curved elements is not supported.");

      if (e->is_triangle())
        regularize_triangle(e);
      else
        regularize_quad(e);
    }
    flatten();
  }

  return parents;

}
