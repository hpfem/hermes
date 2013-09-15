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

#include "element.h"
#include "mesh_util.h"
#include "hash.h"
#include "../mixins2d.h"

#ifdef _WINDOWS
typedef std::shared_ptr<Hermes::Hermes2D::Mesh> MeshSharedPtr;
#else
typedef std::tr1::shared_ptr<Hermes::Hermes2D::Mesh> MeshSharedPtr;
#endif

namespace Hermes
{
  namespace Hermes2D
  {
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
    /// \brief Represents a finite element mesh.
    /// Typical usage:
    /// MeshSharedPtr mesh;
    /// Hermes::Hermes2D::MeshReaderH2DXML mloader;
    /// try
    /// {
    ///&nbsp;mloader.load("mesh.xml", &mesh);
    /// }
    /// catch(Exceptions::MeshLoadFailureException& e)
    /// {
    ///&nbsp;e.print_msg();
    ///&nbsp;return -1;
    /// }
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
      void copy(MeshSharedPtr mesh);

      /// Copies the coarsest elements of another mesh.
      void copy_base(MeshSharedPtr mesh);

      /// Copies the active elements of a converted mesh.
      void copy_converted(MeshSharedPtr mesh);

      /// Creates a mesh from given vertex, triangle, quad, and marker arrays
      void create(int nv, double2* verts, int nt, int3* tris, std::string* tri_markers,
        int nq, int4* quads, std::string* quad_markers, int nm, int2* mark, std::string* boundary_markers);

#pragma region MeshHashGrid
      /// Returns the element pointer located at physical coordinates x, y.
      /// \param[in] x Physical x-coordinate.
      /// \param[in] y Physical y-coordinate.
      /// \param[in] x_reference Optional parameter, in which the x-coordinate of x in the reference domain will be returned.
      /// \param[in] y_reference Optional parameter, in which the y-coordinate of y in the reference domain will be returned.
      Element* element_on_physical_coordinates(double x, double y);

      MeshHashGrid* meshHashGrid;
#pragma endregion

#pragma region MarkerArea
      double get_marker_area(int marker);

      std::map<int, MarkerArea*> marker_areas;
#pragma endregion

#pragma region getters
      /// Retrieves an element by its id number.
      Element* get_element(int id) const;

      /// Returns the total number of elements stored.
      int get_num_elements() const;

      /// Returns the number of base mesh elements.
      int get_num_base_elements() const;

      /// Returns the number of used base mesh elements.
      int get_num_used_base_elements() const;

      /// Returns the current number of active elements in the mesh.
      int get_num_active_elements() const;

      /// Returns the maximum node id number plus one.
      int get_max_element_id() const;

      /// Returns the number of vertex nodes.
      int get_num_vertex_nodes() const;

      /// Returns the number of edge nodes.
      int get_num_edge_nodes() const;
      
      /// Get the mesh bounding box.
      /// \param [out] bottom_left_x Bottom left corner - x coordinate.
      /// \param [out] bottom_left_y Bottom left corner - y coordinate.
      /// \param [out] top_right_x Top right corner - x coordinate.
      /// \param [out] top_right_y Top right corner - y coordinate.
      void get_bounding_box(double& bottom_left_x, double& bottom_left_y, double& top_right_x, double& top_right_y);

      /// For internal use.
      Element* get_element_fast(int id) const;

      /// For internal use.
      unsigned get_seq() const;

      static const std::string eggShellInnerMarker;
      static const std::string eggShell1Marker;
      static const std::string eggShell0Marker;
      
      /// Return the "Egg-shell".
      /// Finds all the elements that neighbor an area with a marker marker.
      /// \param[in] mesh The source mesh
      /// \param[in/out] elements The array where the elements will be returned.
      /// \param[in/out] n_elements Size of the array.
      /// \param[in] marker The marker
      /// \param[in] n_element_guess(optional) Approximate number of elements that will be in this method. Used as an allocation hint. -1 for not-known.
      static MeshSharedPtr get_egg_shell(MeshSharedPtr mesh, std::string marker, unsigned int levels, int n_element_guess = -1);
#pragma endregion

#pragma region refinements
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
      void refine_in_area(std::string marker, int depth = 1, int refinement = 0, bool mark_as_initial = false);
      /// Refines all element sharing the markers passed.
      void refine_in_areas(Hermes::vector<std::string> markers, int depth = 1, int refinement = 0, bool mark_as_initial = false);

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
#pragma endregion

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
        ReferenceMeshCreator(MeshSharedPtr coarse_mesh, int refinement = 0);

        /// Method that does the creation.
        /// THIS IS THE METHOD TO OVERLOAD FOR CUSTOM CREATING OF A REFERENCE MESH.
        virtual MeshSharedPtr create_ref_mesh();

      private:
        /// Storage.
        MeshSharedPtr coarse_mesh;
        int refinement;
      };

      class HERMES_API MarkersConversion
      {
      public:
        MarkersConversion();

        /// Function inserting a marker into conversion_table_for_element_markers.
        /// This function controls if this user_marker x internal_marker is already
        /// present, and if not, it inserts the std::pair.
        /// \return The internal marker where user_marker was inserted.
        int insert_marker(std::string user_marker);

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

        int size() const;

      protected:

        /// Info about the maximum marker used so far, used in determining
        /// of the internal marker for a user-supplied std::string identification for
        /// the purpose of disambiguity.
        int min_marker_unused;

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

      /// Frees all data associated with the mesh.
      void free();

      /// For internal use.
      void set_seq(unsigned seq);

    private:

      /// Internal.
      /// Return the "Egg-shell" internal structures.
      /// Finds all the elements that neighbor an area with a marker marker.
      /// \param[in] mesh The target mesh
      /// \param[in/out] elements The array where the elements will be returned.
      /// \param[in/out] n_elements Size of the array.
      /// \param[in] marker The marker
      /// \param[in] n_element_guess(optional) Approximate number of elements that will be in this method. Used as an allocation hint. -1 for not-known.
      static void get_egg_shell_structures(MeshSharedPtr target_mesh, Element**& elements, int& n_elements, std::string marker, unsigned int levels, int n_element_guess = -1);
      
      
      /// Internal.
      /// Return the "Egg-shell" mesh.
      /// Finds all the elements that neighbor an area with a marker marker.
      /// \param[in/out] mesh The target mesh
      /// \param[in] elements The array from get_egg_shell_structures.
      /// \param[in] n_elements Size of the array from get_egg_shell_structures.
      static void make_egg_shell_mesh(MeshSharedPtr target_mesh, Element** elements, int n_elements);

      /// For internal use.
      void initial_single_check();
      static void initial_multimesh_check(Hermes::vector<MeshSharedPtr > meshes);

      /// For internal use.
      int get_edge_sons(Element* e, int edge, int& son1, int& son2) const;

      /// Refines all quad elements to triangles.
      /// It refines a quadrilateral element into two triangles.
      /// Note: this function creates a base mesh.
      void convert_quads_to_triangles();

      /// Convert all active elements to a base mesh.
      void convert_to_base();

      void refine_element_to_quads_id(int id);

      void refine_triangle_to_quads(Element* e, Element** elems_out = NULL);

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

      /// Bounding box.
      double bottom_left_x, bottom_left_y, top_right_x, top_right_y;
      /// Bounding box calculated.
      bool bounding_box_calculated;
      /// Bounding box calculation.
      void calc_bounding_box();

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

      friend class MeshHashGrid;
      friend class MeshReaderH2D;
      friend class MeshReaderH2DBSON;
      friend class MeshReaderH2DXML;
      friend class MeshReaderH1DXML;
      friend class MeshReaderExodusII;
      friend class DiscreteProblem<double>;
      friend class DiscreteProblem<std::complex<double> >;
      template<typename Scalar> friend class DiscreteProblemDGAssembler;
      friend class WeakForm<double>;
      friend class WeakForm<std::complex<double> >;
      template<typename Scalar> friend class Adapt;
      friend class KellyTypeAdapt<double>;
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
      const ElementMarkersConversion &get_element_markers_conversion() const;
      const BoundaryMarkersConversion &get_boundary_markers_conversion() const;
      ElementMarkersConversion &get_element_markers_conversion();
      BoundaryMarkersConversion &get_boundary_markers_conversion();

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
  }
}
#endif
