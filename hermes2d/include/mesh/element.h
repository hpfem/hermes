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

#ifndef __H2D_ELEMENT_H
#define __H2D_ELEMENT_H

#include "../global.h"
#include "curved.h"
#include "hash.h"

namespace Hermes
{
  namespace Hermes2D
  {
    class HashTable;

    enum ///< node types
    {
      HERMES_TYPE_VERTEX = 0,
      HERMES_TYPE_EDGE = 1
    };

    class Element;

    /// \brief Stores one node of a mesh.
    ///
    /// There are are two variants of this structure, depending on the value of
    /// the member 'type':
    /// <ol> <li> HERMES_TYPE_VERTEX -- vertex node. Has physical coordinates x, y.
    ///&nbsp;   <li> HERMES_TYPE_EDGE   -- edge node. Only stores edge marker and two element pointers.
    /// </ol>
    ///
    struct HERMES_API Node
    {
      /// node id number
      int id;
      /// the number of elements using the node
      unsigned ref : 29;
      /// 0 = vertex node; 1 = edge node
      unsigned type : 1;
      /// 1 = boundary node; 0 = inner node
      unsigned bnd : 1;
      /// array item usage flag
      unsigned used : 1;

      union
      {
        struct
        {
          /// vertex node coordinates
          double x, y;
        };
        struct
        {
          /// edge marker
          int marker;
          /// elements sharing the edge node
          Element* elem[2];
        };
      };

      /// parent id numbers
      int p1, p2;
      /// next node in hash synonym list
      Node* next_hash;

      /// Returns true if the (vertex) node is constrained.
      bool is_constrained_vertex() const;

      void ref_element(Element* e = nullptr);
      void unref_element(HashTable* ht, Element* e = nullptr);
    };

    /// \brief Stores one element of a mesh.
    ///
    /// The element can be a triangle or a quad (nvert == 3 or nvert = 4), active or inactive.
    ///
    /// Vertex/node index number
    ///&nbsp;    [2]
    ///&nbsp;(3)-------(2)
    ///&nbsp; |         |
    ///[3]|  quad.  |[1]
    ///&nbsp; |         |
    ///&nbsp;(0)-------(1)
    ///&nbsp;    [0]
    /// Active elements are actual existing elements in the mesh, which take part in the
    /// computation. Inactive elements are those which have once been active, but were refined,
    /// ie., replaced by several other (smaller) son elements. The purpose of the union is the
    /// following. Active elements store pointers to their vertex and edge nodes. Inactive
    /// elements store pointers to thier son elements and to vertex nodes.
    ///
    /// If an element has curved edges, the member 'cm' points to an associated CurvMap structure,
    /// otherwise it is nullptr.
    ///
    class HERMES_API Element
    {
    public:
      Element();
      /// element id number
      int id;
      /// 0 = active, no sons; 1 = inactive (refined), has sons
      bool active;
      /// array item usage flag
      bool used;
      /// pointer to the parent element for the current son
      Element* parent;
      /// true if the element has been visited during assembling
      bool visited;

      /// Calculates the area of the element.
      /// \param[in] precise_for_curvature If curved elements should be evaluated exactly. \
                                          /// This takes much longer.
      void calc_area(bool precise_for_curvature = false);

      /// Calculates the diameter.
      void calc_diameter();

      /// Returns the center of gravity.
      void get_center(double& x, double& y);

      /// vertex node pointers
      Node* vn[H2D_MAX_NUMBER_VERTICES];
      union
      {
        /// edge node pointers
        Node* en[H2D_MAX_NUMBER_EDGES];
        /// son elements (up to four)
        Element* sons[H2D_MAX_ELEMENT_SONS];
      };

      /// element marker
      int marker;

      // returns the edge orientation. This works for the unconstrained edges.
      bool get_edge_orientation(int ie) const;

      ElementMode2D  get_mode() const {
        return (nvert == 3) ? HERMES_MODE_TRIANGLE : HERMES_MODE_QUAD;
      }

      inline bool is_triangle() const {
        return nvert == 3;
      }
      inline bool is_quad() const {
        return nvert == 4;
      }
      inline bool is_curved() const {
        return cm != nullptr;
      }
      inline unsigned char get_nvert() const {
        return this->nvert;
      }

      inline bool is_parallelogram() const
      {
        if (this->nvert == 3)
          return false;
        else if (this->id == -1)
          return true;

        const double eps = 1e-14;
        return fabs(this->vn[2]->x - (this->vn[1]->x + this->vn[3]->x - this->vn[0]->x)) < eps &&
          fabs(this->vn[2]->y - (this->vn[1]->y + this->vn[3]->y - this->vn[0]->y)) < eps;
      }

      inline bool has_const_ref_map() const
      {
        return (this->nvert == 3 || is_parallelogram()) && (!this->cm);
      }

      bool hsplit() const;
      bool vsplit() const;
      bool bsplit() const;

      /// curved mapping, nullptr if not curvilinear
      CurvMap* cm;
      /// Serves for saving the once calculated area of this element.
      double area;

      bool center_set;
      double x_center, y_center;

      /// Serves for saving the once calculated diameter of this element.
      double diameter;

      /// Increase in integration order, see RefMap::calc_inv_ref_order()
      unsigned short iro_cache;

      /// Helper functions to obtain the index of the next or previous vertex/edge
      inline unsigned char next_vert(unsigned char i) const
      {
        return ((i + 1) % nvert);
      }

      inline unsigned char prev_vert(unsigned char i) const
      {
        return ((i + nvert - 1) % nvert);
      }

      /// Returns a pointer to the neighboring element across the edge 'ie', or
      /// nullptr if it does not exist or is across an irregular edge.
      Element* get_neighbor(int ie) const;

      /// Internal.
      void ref_all_nodes();
      /// Internal.
      void unref_all_nodes(HashTable* ht);

      /// number of vertices (3 or 4)
      unsigned char nvert;
    };

    static Node* get_edge_node();
    static Node* get_vertex_node(Node* v1, Node* v2);

    const int TOP_LEVEL_REF = 123456;
  }
}
#endif
