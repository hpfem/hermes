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

#include "sub_element_map.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename T>
    SubElementMap<T>::SubElementMap() : root(NULL)
    {
    }

    template<typename T>
    SubElementMap<T>::~SubElementMap()
    {
      this->clear();
    }

    template<typename T>
    void SubElementMap<T>::clear()
    {
      this->root.clear_subtree();
      memset(this->root.children, 0, 8 * sizeof(Node*));
    }

    template<typename T>
    T* SubElementMap<T>::get(uint64_t sub_idx, bool& to_add)
    {
      if(sub_idx == 0)
        return root.data;
      Node* node = &root;
      while (sub_idx > 0)
      {
        int current_son = (sub_idx - 1) & 7;
        Node* child_node = node->children[current_son];
        if(child_node == NULL)
          if(to_add)
            node->children[current_son] = new Node(NULL);
          else
            return NULL;
        sub_idx = (sub_idx - 1) >> 3;
        node = child_node;
      }
      if(to_add)
        to_add = false;
      return node->data;
    }

    template<typename T>
    void SubElementMap<T>::add(uint64_t sub_idx, T* data)
    {
      Node* node = &root;
      while (sub_idx > 0)
      {
        int current_son = (sub_idx - 1) & 7;
        Node* child_node = node->children[current_son];
        if(child_node == NULL)
        {
          node->children[current_son] = new Node(NULL);
        }
        sub_idx = (sub_idx - 1) >> 3;
        node = child_node;
      }
      node->data = data;
    }

    template<typename T>
    SubElementMap<T>::Node::Node(T* data) : data(data)
    {
      memset(this->children, 0, 8 * sizeof(Node*));
    }

    template<typename T>
    void SubElementMap<T>::Node::clear_subtree()
    {
      for(int i = 0; i < 8; i++)
      {
        if(this->children[i] != NULL)
        {
          this->children[i]->clear_subtree;
          delete this->children[i];
        }
      }
    }
  }
}