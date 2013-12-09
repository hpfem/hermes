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

#include "global.h"
#include "mesh.h"
#include "hash.h"

namespace Hermes
{
  namespace Hermes2D
  {
    HashTable::HashTable()
    {
      v_table = NULL; e_table = NULL;
    }

    HashTable::~HashTable()
    {
      free();
    }

    void HashTable::init(int size)
    {
      v_table = e_table = NULL;

      mask = size-1;
      if(size & mask) throw Hermes::Exceptions::Exception("Parameter 'size' must be a power of two.");

      // allocate and initialize the hash tables
      v_table = new Node*[size];
      e_table = new Node*[size];

      memset(v_table, 0, size * sizeof(Node*));
      memset(e_table, 0, size * sizeof(Node*));
    }

    void HashTable::copy_list(Node** ptr, Node* node)
    {
      while (node != NULL)
      {
        *ptr = &nodes[node->id];
        ptr = &((*ptr)->next_hash);
        node = node->next_hash;
      }
      *ptr = NULL;
    }

    Node* HashTable::get_node(int id) const
    {
      return &(nodes[id]);
    }

    /// Returns the total number of nodes stored.
    int HashTable::get_num_nodes() const
    {
      return nodes.get_num_items();
    }

    /// Returns the maximum node id number plus one.
    int HashTable::get_max_node_id() const
    {
      return nodes.get_size();
    }

    void HashTable::copy(const HashTable* ht)
    {
      free();
      nodes.copy(ht->nodes);
      mask = ht->mask;

      v_table = new Node*[mask + 1];
      e_table = new Node*[mask + 1];
      for (int i = 0; i <= mask; i++)
      {
        copy_list(v_table + i, ht->v_table[i]);
        copy_list(e_table + i, ht->e_table[i]);
      }
    }

    void HashTable::rebuild()
    {
      memset(v_table, 0, (mask + 1) * sizeof(Node*));
      memset(e_table, 0, (mask + 1) * sizeof(Node*));

      Node* node;
      for_all_nodes(node, this)
      {
        int p1 = node->p1, p2 = node->p2;
        if(p1 > p2) std::swap(p1, p2);
        int idx = hash(p1, p2);

        if(node->type == HERMES_TYPE_VERTEX)
        {
          node->next_hash = v_table[idx];
          v_table[idx] = node;
        }
        else
        {
          node->next_hash = e_table[idx];
          e_table[idx] = node;
        }
      }
    }

    void HashTable::free()
    {
      nodes.free();
      if(v_table != NULL)
      {
        delete [] v_table;
        v_table = NULL;
      }
      if(e_table != NULL)
      {
        delete [] e_table;
        e_table = NULL;
      }
    }

    inline Node* HashTable::search_list(Node* node, int p1, int p2) const
    {
      while (node != NULL)
      {
        if(node->p1 == p1 && node->p2 == p2)
          return node;
        node = node->next_hash;
      }
      return NULL;
    }

    Node* HashTable::get_vertex_node(int p1, int p2)
    {
      // search for the node in the vertex hashtable
      if(p1 > p2) std::swap(p1, p2);
      int i = hash(p1, p2);
      Node* node = search_list(v_table[i], p1, p2);
      if(node != NULL) 
        return node;

      // not found - create a new one
      Node* newnode = nodes.add();

      // initialize the new Node
      newnode->type = HERMES_TYPE_VERTEX;
      newnode->ref = 0;
      newnode->bnd = 0;
      newnode->p1 = p1;
      newnode->p2 = p2;
      assert(nodes[p1].type == HERMES_TYPE_VERTEX && nodes[p2].type == HERMES_TYPE_VERTEX);
      newnode->x = (nodes[p1].x + nodes[p2].x) * 0.5;
      newnode->y = (nodes[p1].y + nodes[p2].y) * 0.5;

      // insert into hashtable
      newnode->next_hash = v_table[i];
      v_table[i] = newnode;

      return newnode;
    }

    Node* HashTable::get_edge_node(int p1, int p2)
    {
      // search for the node in the edge hashtable
      if(p1 > p2) std::swap(p1, p2);
      int i = hash(p1, p2);
      Node* node = search_list(e_table[i], p1, p2);
      if(node != NULL) return node;

      // not found - create a new one
      Node* newnode = nodes.add();

      // initialize the new node
      newnode->type = HERMES_TYPE_EDGE;
      newnode->ref = 0;
      newnode->bnd = 0;
      newnode->p1 = p1;
      newnode->p2 = p2;
      newnode->marker = 0;
      newnode->elem[0] = newnode->elem[1] = NULL;

      // insert into hashtable
      newnode->next_hash = e_table[i];
      e_table[i] = newnode;

      return newnode;
    }

    Node* HashTable::peek_vertex_node(int p1, int p2) const
    {
      if(p1 > p2) std::swap(p1, p2);
      return search_list(v_table[hash(p1, p2)], p1, p2);
    }

    Node* HashTable::peek_edge_node(int p1, int p2) const
    {
      if(p1 > p2) std::swap(p1, p2);
      return search_list(e_table[hash(p1, p2)], p1, p2);
    }

    void HashTable::remove_vertex_node(int id)
    {
      // remove the node from the hash table
      int i = hash(nodes[id].p1, nodes[id].p2);
      Node** ptr = v_table + i;
      Node* node = *ptr;
      while (node != NULL)
      {
        if(node->id == id)
        {
          *ptr = node->next_hash;
          break;
        }
        ptr = &node->next_hash;
        node = *ptr;
      }

      // remove node from the array
      nodes.remove(id);
    }

    void HashTable::remove_edge_node(int id)
    {
      // remove the node from the hash table
      int i = hash(nodes[id].p1, nodes[id].p2);
      Node** ptr = e_table + i;
      Node* node = *ptr;
      while (node != NULL)
      {
        if(node->id == id)
        {
          *ptr = node->next_hash;
          break;
        }
        ptr = &node->next_hash;
        node = *ptr;
      }

      // remove node from the array
      nodes.remove(id);
    }
  }
}