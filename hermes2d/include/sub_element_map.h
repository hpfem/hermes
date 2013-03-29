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
#ifndef __H2D_SUB_ELEMENT_MAP_H
#define __H2D_SUB_ELEMENT_MAP_H

#include "exceptions.h"
#include "function/transformable.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Sub-Element maps - used for replacing std::map<uint64_t, T>.
    template<typename T>
    class SubElementMap
    {
      template<typename Scalar> friend class Function;
      friend class RefMap;
      typedef void (*ProcessingFunction)(T* data);
    public:
      SubElementMap()
      {
      }
      ~SubElementMap()
      {
        this->clear();
      }
      void clear()
      {
        this->root.clear_subtree();
        memset(this->root.children, 0, 8 * sizeof(Node*));
        this->root.data = NULL;
      }

      void run_for_all(ProcessingFunction f)
      {
        run_for_node(&this->root, f);
        this->root.data = NULL;
      }

    private:
      class Node
      {
      public:
        Node() : data(NULL)
        {
          init();
        }
        Node(T* data) : data(data)
        {
          init();
        }
        T* data;
        Node* children[8];
        void clear_subtree()
        {
          for(int i = 0; i < 8; i++)
          {
            if(this->children[i] != NULL)
            {
              this->children[i]->clear_subtree();
              delete this->children[i];
            }
          }
        }
      private:
        void init()
        {
          memset(this->children, 0, 8 * sizeof(Node*));
        }
      };

    public:
      void run_for_node(Node* node, ProcessingFunction f)
      {
        for(int i = 0; i < 8; i++)
        {
          if(node->children[i] != NULL)
            run_for_node(node->children[i], f);
        }
        if(node->data != NULL)
          f(node->data);
      }

      Node* get(uint64_t sub_idx, bool& to_add)
      {
        if(sub_idx == 0)
        {
          if(root.data != NULL)
            to_add = false;
          return &root;
        }
        Node* node = &root;
        bool added = false;
        int sons[Transformable::H2D_MAX_TRN_LEVEL];
        int sons_count = 0;
        while (sub_idx > 0)
        {
          sons[sons_count++] = (sub_idx - 1) & 7;
          sub_idx = (sub_idx - 1) >> 3;
        }
        
        for(int i = sons_count - 1; i >= 0 ; i--)
        {
          Node* child_node = node->children[sons[i]];
          if(child_node == NULL)
            if(to_add)
            {
              added = true;
              child_node = node->children[sons[i]] = new Node(NULL);
            }
            else
              return NULL;
          node = child_node;
        }
        if(!added && node->data != NULL)
          to_add = false;
        return node;
      }
    
      Node root;
    };
  }
}

#endif
