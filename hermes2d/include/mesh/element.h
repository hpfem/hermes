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
      int id;          ///< node id number
      unsigned ref:29; ///< the number of elements using the node
      unsigned type:1; ///< 0 = vertex node; 1 = edge node
      unsigned bnd:1;  ///< 1 = boundary node; 0 = inner node
      unsigned used:1; ///< array item usage flag

      union
      {
        struct
        {
          double x, y; ///< vertex node coordinates
        };
        struct
        {
          int marker;       ///< edge marker
          Element* elem[2]; ///< elements sharing the edge node
        };
      };

      int p1, p2; ///< parent id numbers
      Node* next_hash; ///< next node in hash synonym list

      /// Returns true if the (vertex) node is constrained.
      bool is_constrained_vertex() const;

      void ref_element(Element* e = NULL);
      void unref_element(HashTable* ht, Element* e = NULL);
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
    /// otherwise it is NULL.
    ///
    class HERMES_API Element
    {
    public:
      Element();
      int id;              ///< element id number
      bool active;   ///< 0 = active, no sons; 1 = inactive (refined), has sons
      bool used;     ///< array item usage flag
      Element* parent;     ///< pointer to the parent element for the current son
      bool visited;        ///< true if the element has been visited during assembling

      /// Calculates the area of the element.
      /// \param[in] precise_for_curvature If curved elements should be evaluated exactly. \
      /// This takes much longer.
      double get_area(bool precise_for_curvature = false);

      /// Returns the length of the longest edge for triangles, and the
      /// length of the longer diagonal for quads. Ignores element curvature.
      double get_diameter();

      /// Returns the center of gravity.
      void get_center(double& x, double& y);

      Node* vn[H2D_MAX_NUMBER_VERTICES];   ///< vertex node pointers
      union
      {
        Node* en[H2D_MAX_NUMBER_EDGES];      ///< edge node pointers
        Element* sons[H2D_MAX_ELEMENT_SONS]; ///< son elements (up to four)
      };

      int marker;        ///< element marker

      // returns the edge orientation. This works for the unconstrained edges.
      int get_edge_orientation(int ie) const;
      ElementMode2D  get_mode() const;

      bool is_triangle() const;
      bool is_quad() const;
      bool is_curved() const;
      int get_nvert() const;

      bool hsplit() const;
      bool vsplit() const;
      bool bsplit() const;

      CurvMap* cm; ///< curved mapping, NULL if not curvilinear
      /// Serves for saving the once calculated area of this element.
      bool areaCalculated;
      /// Serves for saving the once calculated area of this element.
      double area;

      bool center_set;
      double x_center, y_center;

      /// Serves for saving the once calculated diameter of this element.
      bool diameterCalculated;
      /// Serves for saving the once calculated diameter of this element.
      double diameter;

      /// Increase in integration order, see RefMap::calc_inv_ref_order()
      int iro_cache;

      /// Helper functions to obtain the index of the next or previous vertex/edge
      int next_vert(int i) const;
      int prev_vert(int i) const;

      /// Returns a pointer to the neighboring element across the edge 'ie', or
      /// NULL if it does not exist or is across an irregular edge.
      Element* get_neighbor(int ie) const;

      /// Internal.
      void ref_all_nodes();
      /// Internal.
      void unref_all_nodes(HashTable* ht);

      unsigned nvert:30; ///< number of vertices (3 or 4)
    };

    static Node* get_edge_node();
    static Node* get_vertex_node(Node* v1, Node* v2);

    const int TOP_LEVEL_REF = 123456;
  }
}
#endif
