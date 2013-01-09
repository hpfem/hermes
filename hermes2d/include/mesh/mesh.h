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

#ifndef __H2D_MESH_H
#define __H2D_MESH_H

#include "../global.h"
#include "curved.h"
#include "hash.h"
#include "../mixins2d.h"

namespace Hermes
{
  namespace Hermes2D
  {
    enum ///< node types
    {
      HERMES_TYPE_VERTEX = 0,
      HERMES_TYPE_EDGE = 1
    };

    class Element;
    class HashTable;

    template<typename Scalar> class Space;
    template<typename Scalar> class KellyTypeAdapt;
    struct MItem;
    struct Rect;
    extern unsigned g_mesh_seq;

    namespace RefinementSelectors
    {
      template<typename Scalar> class Selector;
      template<typename Scalar> class POnlySelector;
      template<typename Scalar> class HOnlySelector;
      template<typename Scalar> class OptimumSelector;
      template<typename Scalar> class ProjBasedSelector;
      template<typename Scalar> class H1ProjBasedSelector;
      template<typename Scalar> class L2ProjBasedSelector;
      template<typename Scalar> class HcurlProjBasedSelector;
    }

    namespace Views
    {
      class MeshView;
    }

    /// \brief Stores one node of a mesh.
    ///
    /// There are are two variants of this structure, depending on the value of
    /// the member 'type':
    /// <ol> <li> HERMES_TYPE_VERTEX -- vertex node. Has physical coordinates x, y.
    ///      <li> HERMES_TYPE_EDGE   -- edge node. Only stores edge marker and two element pointers.
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
    ///       [2]
    ///   (3)-------(2)
    ///    |         |
    ///[3]|  quad.  |[1]
    ///    |         |
    ///   (0)-------(1)
    ///       [0]
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
      
      /// Calculates the area of the element. For curved elements, this is only
      /// an approximation: the curvature is not accounted for.
      double get_area();

      /// Returns the length of the longest edge for triangles, and the
      /// length of the longer diagonal for quads. Ignores element curvature.
      double get_diameter();

      Node* vn[4];   ///< vertex node pointers
      union
      {
        Node* en[4];      ///< edge node pointers
        Element* sons[4]; ///< son elements (up to four)
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

    protected:
      CurvMap* cm; ///< curved mapping, NULL if not curvilinear
      /// Serves for saving the once calculated area of this element.
      bool areaCalculated;
      /// Serves for saving the once calculated area of this element.
      double area;

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
    private:
      unsigned nvert:30; ///< number of vertices (3 or 4)

      friend class Mesh;
      friend class MeshReader;
      friend class MeshReaderH2D;
      friend class MeshReaderH1DXML;
      friend class MeshReaderH2DXML;
      friend class PrecalcShapeset;
      template<typename Scalar> friend class Space;
      template<typename Scalar> friend class Adapt;
      template<typename Scalar> friend class H1Space;
      template<typename Scalar> friend class HcurlSpace;
      template<typename Scalar> friend class HdivSpace;
      template<typename Scalar> friend class L2Space;
      template<typename Scalar> friend class KellyTypeAdapt;
      template<typename Scalar> friend class DiscreteProblem;
      template<typename Scalar> friend class Solution;
      template<typename Scalar> friend class NeighborSearch;
      template<typename Scalar> friend class Filter;
      template<typename Scalar> friend class MeshFunction;
      template<typename Scalar> friend class Global;
      template<typename Scalar> friend class RefinementSelectors::Selector;
      template<typename Scalar> friend class RefinementSelectors::POnlySelector;
      template<typename Scalar> friend class RefinementSelectors::HOnlySelector;
      template<typename Scalar> friend class RefinementSelectors::OptimumSelector;
      template<typename Scalar> friend class RefinementSelectors::ProjBasedSelector;
      template<typename Scalar> friend class RefinementSelectors::H1ProjBasedSelector;
      template<typename Scalar> friend class RefinementSelectors::L2ProjBasedSelector;
      template<typename Scalar> friend class RefinementSelectors::HcurlProjBasedSelector;
      friend class Views::ScalarView;
      friend class Views::Linearizer;
      friend class RefMap;
      friend class Traverse;
      friend class Transformable;
      friend class CurvMap;
      friend class Views::Orderizer;
      friend class Views::Vectorizer;
      friend bool is_twin_nurbs(Element* e, int i);
      friend int rtb_criterion(Element* e);
      friend CurvMap* create_son_curv_map(Element* e, int son);
    };

    /// \brief Represents a finite element mesh.
    ///
    class HERMES_API Mesh : public HashTable, public Hermes::Hermes2D::Mixins::StateQueryable
    {
    public:
      Mesh();
      virtual ~Mesh();

      /// Initializes the mesh.
      /// \param size[in] Hash table size; must be a power of two.
      void init(int size = H2D_DEFAULT_HASH_SIZE);

      /// State querying helpers.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "Mesh"; }

      /// Rescales the mesh.
      bool rescale(double x_ref, double y_ref);

      /// Creates a copy of another mesh.
      void copy(const Mesh* mesh);

      /// Copies the coarsest elements of another mesh.
      void copy_base(Mesh* mesh);

      /// Copies the active elements of a converted mesh.
      void copy_converted(Mesh* mesh);

      /// Creates a mesh from given vertex, triangle, quad, and marker arrays
      void create(int nv, double2* verts, int nt, int3* tris, std::string* tri_markers,
                  int nq, int4* quads, std::string* quad_markers, int nm, int2* mark, std::string* boundary_markers);

      /// Frees all data associated with the mesh.
      void free();

      /// Retrieves an element by its id number.
      Element* get_element(int id) const;

      /// Returns the total number of elements stored.
      int get_num_elements() const;

      /// Returns the number of coarse mesh elements.
      int get_num_base_elements() const;

      /// Returns the current number of active elements in the mesh.
      int get_num_active_elements() const;

      /// Returns the maximum node id number plus one.
      int get_max_element_id() const;

      /// Refines an element.
      /// \param id[in] Element id number.
      /// \param refinement[in] Ignored for triangles. If the element
      /// is a quad, 0 means refine in both directions, 1 means refine
      /// horizontally (with respect to the reference domain), 2 means
      /// refine vertically.
      void refine_element_id(int id, int refinement = 0);

      /// Refines all elements.
      /// \param refinement[in] Same meaning as in refine_element_id().
      void refine_all_elements(int refinement = 0, bool mark_as_initial = false);

      /// Selects elements to refine according to a given criterion and
      /// performs 'depth' levels of refinements. The criterion function
      /// receives a pointer to an element to be considered.
      /// It must return -1 if the element is not to be refined, 0 if it
      /// should be refined uniformly, 1 if it is a quad and should be split
      /// horizontally, 2 if it is a quad and should be split vertically,
      /// and 3 if it is a triangle and should be split into three quads.
      void refine_by_criterion(int (*criterion)(Element* e), int depth = 1, bool mark_as_initial = false);

      /// Performs repeated refinements of elements containing the given vertex.
      /// A mesh graded towards the vertex is created.
      void refine_towards_vertex(int vertex_id, int depth = 1, bool mark_as_initial = false);

      /// Performs repeated refinements of elements touching a part of the
      /// boundary marked by 'marker'. Elements touching both by an edge or
      /// by a vertex are refined. 'aniso' allows or disables anisotropic
      /// splits of quads.
      void refine_towards_boundary(Hermes::vector<std::string> markers, int depth = 1, bool aniso = true, bool mark_as_initial = false);
      void refine_towards_boundary(std::string marker, int depth = 1, bool aniso = true, bool mark_as_initial = false);

      /// Refines all element sharing the marker passed.
      void refine_in_area(std::string marker, int depth = 1, bool mark_as_initial = false);
      /// Refines all element sharing the markers passed.
      void refine_in_areas(Hermes::vector<std::string> markers, int depth = 1, bool mark_as_initial = false);

      /// Regularizes the mesh by refining elements with hanging nodes of
      /// degree more than 'n'. As a result, n-irregular mesh is obtained.
      /// If n = 0, completely regular mesh is created. In this case, however,
      /// due to incompatible refinements, the element refinement hierarchy
      /// is removed and all elements become top-level elements. Also, total
      /// regularization does not work on curved elements.
      /// Returns an array of new element parents which can be passed to
      /// Space::distribute_orders(). The array must be deallocated with ::free().
      int* regularize(int n);

      /// Recursively removes all son elements of the given element and
      /// makes it active.
      void unrefine_element_id(int id);

      /// Unrefines all elements with immediate active sons. In effect, this
      /// shaves off one layer of refinements from the mesh. If done immediately
      /// after refine_all_elements(), this function reverts the mesh to its
      /// original state. However, it is not exactly an inverse to
      /// refine_all_elements().
      void unrefine_all_elements(bool keep_initial_refinements = true);

      /// For internal use.
      Element* get_element_fast(int id) const;

      /// For internal use.
      unsigned get_seq() const;

      /// For internal use.
      void set_seq(unsigned seq);

      /// Class for creating reference mesh.
      class HERMES_API ReferenceMeshCreator
      {
      public:
        /// Constructor.
        /// \param[in] coarse_mesh The coarse (original) mesh.
        /// \param refinement[in] Ignored for triangles. If the element
        /// is a quad, 0 means refine in both directions, 1 means refine
        /// horizontally (with respect to the reference domain), 2 means
        /// refine vertically.
        ReferenceMeshCreator(Mesh* coarse_mesh, int refinement = 0);

        /// Method that does the creation.
        /// THIS IS THE METHOD TO OVERLOAD FOR CUSTOM CREATING OF A REFERENCE MESH.
        virtual Mesh* create_ref_mesh();

      private:
        /// Storage.
        Mesh* coarse_mesh;
        int refinement;
      };

    private:
      /// For internal use.
      int get_edge_sons(Element* e, int edge, int& son1, int& son2) const;

      /// Refines all quad elements to triangles.
      /// It refines a quadrilateral element into two triangles.
      /// Note: this function creates a base mesh.
      void convert_quads_to_triangles();

      /// Convert all active elements to a base mesh.
      void convert_to_base();

      void refine_element_to_quads_id(int id);

      void refine_triangle_to_quads(Mesh* mesh, Element* e, Element** elems_out = NULL);

      void refine_element_to_triangles_id(int id);

      void refine_quad_to_triangles(Element* e);

      /// Refines one quad element into four quad elements.
      /// The difference between refine_quad_to_quads() and refine_quad()
      /// is that all the internal edges of the former's son elements are
      /// straight edges.
      void refine_quad_to_quads(Element* e, int refinement = 0);

      void convert_element_to_base_id(int id);
      void convert_triangles_to_base(Element* e);
      void convert_quads_to_base(Element* e);

      Array<Element> elements;
      int nactive;
      unsigned seq;

      int nbase, ntopvert;
      int ninitial;

      void unrefine_element_internal(Element* e);

      /// Returns a NURBS curve with reversed control points and inverted knot vector.
      /// Used for curved edges inside a mesh, where two mirror Nurbs have to be created
      /// for the adjacent elements
      ///
      Nurbs* reverse_nurbs(Nurbs* nurbs);
      Node*  get_base_edge_node(Element* base, int edge);

      int* parents;
      int parents_size;

      int  get_edge_degree(Node* v1, Node* v2);
      void assign_parent(Element* e, int i);
      void regularize_triangle(Element* e);
      void regularize_quad(Element* e);
      void flatten();

      class HERMES_API MarkersConversion
      {
      public:
        MarkersConversion();

        /// Info about the maximum marker used so far, used in determining
        /// of the internal marker for a user-supplied std::string identification for
        /// the purpose of disambiguity.
        int min_marker_unused;

        /// Function inserting a marker into conversion_table_for_element_markers.
        /// This function controls if this user_marker x internal_marker is already
        /// present, and if not, it inserts the std::pair.
        void insert_marker(int internal_marker, std::string user_marker);

        /// Struct for return type of get_user_marker().
        struct StringValid
        {
          StringValid();
          StringValid(std::string marker, bool valid);
          std::string marker;
          bool valid;
        };

        /// Struct for return type of get_internal_marker().
        struct IntValid
        {
          IntValid();
          IntValid(int marker, bool valid);
          int marker;
          bool valid;
        };

        /// Lookup functions.
        /// Find a user marker for this internal marker.
        StringValid get_user_marker(int internal_marker) const;

        /// Find an internal marker for this user_marker.
        IntValid get_internal_marker(std::string user_marker) const;

        enum MarkersConversionType {
          HERMES_ELEMENT_MARKERS_CONVERSION = 0,
          HERMES_BOUNDARY_MARKERS_CONVERSION = 1
        };

        virtual MarkersConversionType get_type() const = 0;

      protected:
        /// Conversion tables between the std::string markers the user sets and
        /// the markers used internally as members of Elements, Nodes.
        std::map<int, std::string> conversion_table;

        /// Inverse tables, so that it is possible to search using either
        /// the internal representation, or the user std::string value.
        std::map<std::string, int> conversion_table_inverse;

        friend class Space<double>;
        friend class Space<std::complex<double> >;
        friend class Mesh;
      };

      /// \brief Curved element exception.
      /// Exception occurs when there is a curved element where we only process not curved.
      class HERMES_API CurvedException : public Hermes::Exceptions::Exception
      {
      public:
        /// Constructor
        /// \param[in] element_id Id of the element that is curved.
        CurvedException(int elementId);
        CurvedException(const CurvedException & e);

        int getElementId() const;
      private:
        int elementId;
      };

      class ElementMarkersConversion : public MarkersConversion
      {
      public:
        ElementMarkersConversion();
        virtual MarkersConversionType get_type() const;
      };

      class BoundaryMarkersConversion : public MarkersConversion
      {
      public:
        BoundaryMarkersConversion();
        virtual MarkersConversionType get_type() const;
      };

      ElementMarkersConversion element_markers_conversion;
      BoundaryMarkersConversion boundary_markers_conversion;

      friend class MeshReaderH2D;
      friend class MeshReaderH2DXML;
      friend class MeshReaderH1DXML;
      friend class MeshReaderExodusII;
      friend class DiscreteProblem<double>;
      friend class DiscreteProblem<std::complex<double> >;
      friend class WeakForm<double>;
      friend class WeakForm<std::complex<double> >;
      template<typename Scalar> friend class Adapt;
      friend class KellyTypeAdapt<double>;
      template<typename Scalar> friend class Global;
      friend class KellyTypeAdapt<std::complex<double> >;
      template<typename Scalar> friend class Solution;
      template<typename Scalar> friend class Filter;
      template<typename Scalar> friend class MeshFunction;
      friend class RefMap;
      friend class Traverse;
      friend class Transformable;
      friend class Curved;
      friend class Views::MeshView;
      template<typename Scalar> friend class RefinementSelectors::Selector;
      template<typename Scalar> friend class RefinementSelectors::POnlySelector;
      template<typename Scalar> friend class RefinementSelectors::HOnlySelector;
      template<typename Scalar> friend class RefinementSelectors::OptimumSelector;
      template<typename Scalar> friend class RefinementSelectors::ProjBasedSelector;
      template<typename Scalar> friend class RefinementSelectors::H1ProjBasedSelector;
      template<typename Scalar> friend class RefinementSelectors::L2ProjBasedSelector;
      template<typename Scalar> friend class RefinementSelectors::HcurlProjBasedSelector;
      friend class PrecalcShapeset;
      template<typename Scalar> friend class Space;
      template<typename Scalar> friend class H1Space;
      template<typename Scalar> friend class HcurlSpace;
      template<typename Scalar> friend class HdivSpace;
      template<typename Scalar> friend class L2Space;
      friend class Views::ScalarView;
      friend class Views::Orderizer;
    public:
      ElementMarkersConversion &get_element_markers_conversion();
      BoundaryMarkersConversion &get_boundary_markers_conversion();

      /*  node and son numbering on a triangle:

        -Triangle to triangles refinement

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

        -Triangle to quads refinement

                        vn[2]                                     vn[2]

                          *                                        *
                          / \                                      / \
                        /   \                                    /   \
                        /     \                                  /     \
                      /       \                          vn[3] * son[2]* vn[1]
                      /         \                       vn[3] *  \     /  * vn[2]
            en[2]   *           *   en[1]                   / \  \   /  / \
                    /             \                         /   \ vn[0] /   \
                  /               \                       /     \  *  /     \
                  /                 \                     /       \   /       \
                /         *         \                   /   vn[2] * * vn[3]   \
                /                     \                 /          | |          \
              /                       \               /  son[0]   | |  son[1]   \
              /                         \             /            | |            \
            *-------------*-------------*           *-------------* *-------------*
                                                  vn[0]      vn[1]   vn[0]        vn[1]
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

      Element* create_quad(int marker, Node* v0, Node* v1, Node* v2, Node* v3, CurvMap* cm, int id = -1);
      Element* create_triangle(int marker, Node* v0, Node* v1, Node* v2, CurvMap* cm, int id = -1);
      void refine_element(Element* e, int refinement);

      /// Vector for storing refinements in order to be able to save/load meshes with identical element IDs.
      /// Refinement "-1" stands for unrefinement.
      Hermes::vector<std::pair<unsigned int, int> > refinements;

      /// Refines a quad element into four quads, or two quads (horizontally or
      /// vertically. If mesh != NULL, the new elements are incorporated into
      /// the mesh. The option mesh == NULL is used to perform adaptive numerical
      /// quadrature. If sons_out != NULL, pointers to the new elements will be
      /// saved there.
      void refine_quad(Element* e, int refinement, Element** sons_out = NULL);
      void refine_triangle_to_triangles(Element* e, Element** sons = NULL);

      /// Computing vector length.
      static double vector_length(double a_1, double a_2);

      /// Checking whether the points p, q, r lie on the same line.
      static bool same_line(double p_1, double p_2, double q_1, double q_2, double r_1, double r_2);

      /// Checking whether the angle of vectors 'a' and 'b' is between zero and Pi.
      static bool is_convex(double a_1, double a_2, double b_1, double b_2);
      static void check_triangle(int i, Node *&v0, Node *&v1, Node *&v2);
      static void check_quad(int i, Node *&v0, Node *&v1, Node *&v2, Node *&v3);
    };

    static Node* get_edge_node();
    static Node* get_vertex_node(Node* v1, Node* v2);

    /// Helper macros for easy iteration through all elements, nodes etc. in a Mesh.
    #define for_all_elements(e, mesh) \
            for (int _id = 0, _max = (mesh)->get_max_element_id(); _id < _max; _id++) \
              if(((e) = (mesh)->get_element_fast(_id))->used)

    #define for_all_base_elements(e, mesh) \
            for (int _id = 0; _id < (mesh)->get_num_base_elements(); _id++) \
              if(((e) = (mesh)->get_element_fast(_id))->used)

    #define for_all_base_elements_incl_inactive(e, mesh) \
            for (int _id = 0; _id < (mesh)->get_num_base_elements(); _id++) \
              if(((e) = (mesh)->get_element_fast(_id))->used || !((e) = (mesh)->get_element_fast(_id))->used)

    #define for_all_active_elements(e, mesh) \
            for (int _id = 0, _max = (mesh)->get_max_element_id(); _id < _max; _id++) \
              if(((e) = (mesh)->get_element_fast(_id))->used) \
                if((e)->active)

    #define for_all_inactive_elements(e, mesh) \
            for (int _id = 0, _max = (mesh)->get_max_element_id(); _id < _max; _id++) \
              if(((e) = (mesh)->get_element_fast(_id))->used) \
                if(!(e)->active)

    #define for_all_nodes(n, mesh) \
            for (int _id = 0, _max = (mesh)->get_max_node_id(); _id < _max; _id++) \
              if(((n) = (mesh)->get_node(_id))->used)

    #define for_all_vertex_nodes(n, mesh) \
            for (int _id = 0, _max = (mesh)->get_max_node_id(); _id < _max; _id++) \
              if(((n) = (mesh)->get_node(_id))->used) \
                if(!(n)->type)

    #define for_all_edge_nodes(n, mesh) \
            for (int _id = 0, _max = (mesh)->get_max_node_id(); _id < _max; _id++) \
              if(((n) = (mesh)->get_node(_id))->used) \
                if((n)->type)

    const int TOP_LEVEL_REF = 123456;
  }
}
#endif
