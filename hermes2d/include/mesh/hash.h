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

#ifndef __H2D_HASH_H
#define __H2D_HASH_H

#include "hermes_common.h"

namespace Hermes
{
  namespace Hermes2D
  {
    struct Node;

    namespace Views
    {
      class MeshView;
      class ScalarView;
      class Linearizer;
      class Vectorizer;
      class Orderizer;
    };

    /// \brief Stores and searches node tables.
    ///
    /// HashTable is a base class for Mesh. It serves as a container for all nodes
    /// of a mesh. Moreover, it has node searching functions based on hash tables.
    ///
    class HERMES_API HashTable : public Hermes::Mixins::Loggable
    {
    public:
      /// Retrieves a node by its id number.
      Node* get_node(int id) const;

      /// Returns the total number of nodes stored.
      int get_num_nodes() const;

      /// Returns the maximum node id number plus one.
      int get_max_node_id() const;

    protected:
      HashTable();
      ~HashTable();

      /// Returns a vertex node with parent id's p1 and p2 if it exists, NULL otherwise.
      Node* peek_vertex_node(int p1, int p2) const;

      /// Returns an edge node with parent id's p1 and p2 if it exists, NULL otherwise.
      Node* peek_edge_node(int p1, int p2) const;

      /// Central function: obtains a vertex node pointer given the id
      /// numbers of its parents. If the vertex node does not exist, it is
      /// created first.
      Node* get_vertex_node(int p1, int p2);

      /// Central function: obtains an edge node pointer given the id
      /// numbers of its parents. If the edge node does not exist, it is
      /// created first.
      Node* get_edge_node(int p1, int p2);

      static const int H2D_DEFAULT_HASH_SIZE = 0x8000; // 32K entries

      Array<Node> nodes; ///< Array storing all nodes

      /// Initializes the hash table.
      /// \param size[in] Hash table size; must be a power of two.
      void init(int size = H2D_DEFAULT_HASH_SIZE);

      /// Copies another hash table contents
      void copy(const HashTable* ht);

      /// Reconstructs the hashtable, after, e.g., the nodes have been loaded from a file.
      void rebuild();

      /// Frees all memory used by the instance.
      void free();

      /// Removes a vertex node with parent id's p1 and p2.
      void remove_vertex_node(int id);

      /// Removes an edge node with parent id's p1 and p2.
      void remove_edge_node(int id);

      // Internal members
    private:

      Node** v_table; ///< Vertex node hash table
      Node** e_table; ///< Edge node hash table

      int mask;

      inline int hash(int p1, int p2) const { return (984120265*p1 + 125965121*p2) & mask; }

      /// Searches a list of hash synonyms given the first list item.
      /// Returns the node matching the parent ids p1 and p2.
      Node* search_list(Node* node, int p1, int p2) const;

      /// Creates a copy of a hash synonym list.
      void copy_list(Node** ptr, Node* node);

      friend struct Node;
      friend class MeshReaderH2D;
      template<typename Scalar> friend class NeighborSearch;
      template<typename Scalar> friend class Space;
      template<typename Scalar> friend class H1Space;
      template<typename Scalar> friend class L2Space;
      template<typename Scalar> friend class HcurlSpace;
      template<typename Scalar> friend class HdivSpace;
      friend class Views::ScalarView;
      friend class Views::Linearizer;
    };
  }
}
#endif