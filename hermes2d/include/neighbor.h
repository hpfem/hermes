#ifndef __H2D_NEIGHBOR_H
#define __H2D_NEIGHBOR_H

#include "mesh/mesh.h"
#include "quadrature/quad.h"
#include "function/solution.h"
#include "forms.h"
#include "mesh/refmap.h"
#include "asmlist.h"
#include "space/space.h"
namespace Hermes
{
  namespace Hermes2D
  {
    /*** Class NeighborSearch. ***/

    /*!\class NeighborSearch neighbor.h "src/neighbor.h"
    * \brief This class characterizes a neighborhood of a given edge in terms of adjacent elements
    *        and provides methods for getting limit values of discontinuous functions from both sides of the edge.
    *
    * Each instance of the class is connected to a mesh and its active element. The active element becomes the
    * <em> central element </em> of the neighborhood and all adjacent elements the <em> neighbors </em>.
    *
    * In order to search for the neighboring elements, one selects a particular edge of the central element and calls
    * the method \c set_active_edge. This enumerates the neighbors and fills in the array \c transformations, which
    * will be needed later for getting function values at matching points from both sides of the selected (active) edge.
    *
    * The actual procedure depends on the relative size of the central element with respect to the neighbor element(s)
    * across the active edge:
    *    - If active edge is shared by two elements of same size, then \c Element::get_neighbor is used to retrieve
    *      pointer to the neighboring element. Moreover, no transformations are needed to obtain values of any given
    *      function from either side of active edge.
    *    - If the neighbor element is bigger than the central element, then we "go up" on the central element, until
    *      we find its parent that has the same size as the neighbor. We keep track of the visited intermediate parents
    *      and after the final one has been found, we use them in reverse order to fill in the \c transformations array.
    *      These transformations will have to be pushed to a function on the neighboring (bigger) element in order to
    *      obtain values at points matching those from the central element's side.
    *    - If the neighbor element is smaller than the central element, it means that it is one of several neighbors
    *      across the active edge. Hence, we "go down" in the central element in order to find a (virtual) sub-element
    *      matching the currently processed neighbor and store the corresponding transformations in the neighbor's field
    *      of the array \c transformations. This way, we obtain for each neighbor a set of transformations which have to
      *      be applied on the central element to transform integration points to its correct part matching the neighbor.
    *
    * There may be either one or more neighbors (see below). In the second case, each of the neighbors shares with the
    * central element a unique segment of the active edge. The next step is therefore a loop through all neighbors,
    * setting the appropriate active segment through \c set_active_segment at the beginning of every iteration. This also
    * handles the "go-down" transformations of the attached central element's pss and creates a pss on the
    * neighbor element for use in evaluating discontinuous bilinear forms (see below).
    *
    * If only external functions are supposed to be discontinuous across the active edge (e.g. in the case of linear
    * forms), the \c init_ext_fn methods may now be used to create pointers to new \c DiscontinuousFunc objects - first to
    * one used for evaluating integration order and then, after computing the quadrature order for given form and setting
    * it via \c set_quad_order, to that with the actual values at the matching integration points from both sides of the
    * active segment.
    *
    * If also the test functions need to be considered discontinuous (e.g when assembling DG bilinear forms), the local
    * precalculated shapesets on the central element (attached from the assembling procedure) and on the current neighbor
    * element are extended by zero to the whole neighborhood of the two elements (see \c create_extended_shapeset).
    * A so-called <em> extended shapeset </em> thus created is stored in the variable \c supported_shapes and may be queried
    * for the count of all contained extended shapes and their global DOF numbers. Assembling is done over all shape
    * functions from the extended shapes. The user gets a pointer to the current extended shape function by calling
    * \c supported_shapes.get_extended_shape_fn and may use it then to obtain \c DiscontinuousFunc objects with the shape
    * function's order/values on both sides of active edge as in the previous case of external functions.
    *
    */

    template<typename Scalar>
    class HERMES_API NeighborSearch
    {
    public:

      /// Constructor.
      ///
      /// \param[in]  el    Central element of the neighborhood (current active element in the assembling procedure).
      /// \param[in]  mesh  Mesh on which we search for the neighbors.
      ///
      NeighborSearch(Element* el, const Mesh* mesh);
      NeighborSearch(const NeighborSearch& ns);

      /// Destructor.
      ~NeighborSearch();

      /*** Methods for changing active state for further calculations. ***/

      /// Set active edge and compute all information about the neighbors.
      ///
      /// In particular, it fills the \c neighbors and \c neighbor_edges vectors and the \c transformations array used
      /// for transforming either the central or the neighboring elements to the same size. It also sets \c neighborhood_type,
      /// according to whether we can find a vertex in the middle of the active edge - if we can, then we go down and the
      /// corresponding transformations will be filled for the bigger central element. This is performed by the method
      /// \c find_act_elem_down. Otherwise, we go up and will later push the transformations from \c transfomations to
      /// functions on the bigger neighboring element. Way up is Realized by the method \c find_act_elem_up.
      ///
      /// \param[in] edge Local (element dependent) number of the edge.
      /// \param[in] ignore_visited_segments   If true, it indicates that an edge-based discontinuous Galerkin formulation
      ///                 is assumed. This means that once an edge is visited during assembling, integrals on that edge from
      ///                 both its sides are assembled at once and the edge should be ignored when the neighbor is assembled.
      ///                 A sum of edge-integrals over edges represents interelement coupling across the edge in the math.
      ///                 formulation.
      ///
      ///                 If false, an element-based DG formulation is assumed, where each edge of every element is assembled
      ///                 (a sum of edge-integrals over elements appears in the mathematical formulation). Note that this may
      ///                 be also use in the edge-based form., e.g. to specify linear forms.
      ///
      void set_active_edge(int edge);

      /// Enhancement of set_active_edge for multimesh assembling.
      bool set_active_edge_multimesh(const int& edge);

      /// Extract transformations in the correct direction from the provided sub_idx.
      Hermes::vector<unsigned int> get_transforms(uint64_t sub_idx) const;

      /// Gives an info if edge is an intra- or inter- element edge.
      bool is_inter_edge(const int& edge, const Hermes::vector<unsigned int>& transformations) const;

      /// Update according to the subelement mapping of the central element.
      void update_according_to_sub_idx(const Hermes::vector<unsigned int>& transformations);

      /// Special function for handling subelement transformations in the case of more than one neighboring active elements.
      void handle_sub_idx_way_down(const Hermes::vector<unsigned int>& transformations);

      /// Give the info if the two transformations are correct, w.r.t. the edge.
      /// Simply compares a to b in case of triangles, does more work in case of quads.
      bool compatible_transformations(unsigned int a, unsigned int b, int edge) const;

      /// Clear the initial_sub_idxs from the central element transformations of NeighborSearches with multiple neighbors.
      /// Does nothing in the opposite case.
      void clear_initial_sub_idx();

      /// In case we determine a neighbor is not correct due to a subelement mapping, we delete it.
      void delete_neighbor(unsigned int position);

      /// Fill arrays with function values / derivatives at both sides of the active segment of active edge.
      /// Assumes that integration order has been set by \c set_quad_order.
      ///
      /// \param[in] fu MeshFunction whose values are requested.
      /// \return Pointer to a discontinuous function allowing to access the values from each side of the active edge.
      ///
      DiscontinuousFunc<Scalar>* init_ext_fn(MeshFunction<Scalar>* fu);

      class ExtendedShapeset;
      ExtendedShapeset *supported_shapes; ///< Object allowing to set/get a particular shape function from the extended
      ///< shapeset and retrieve global assembly information for it.

      /// Form the extended shapeset.
      ///
      /// Initializes the class attribute \c supported_shapes for the first time and updates it whenever the active segment
      /// or edge change. The extended shapeset consists of shape functions from both the central element and current
      /// neighbor element, extended by zero to the union of these elements.
      ///
      /// \param[in]  space Space from which the local shapesets are drawn.
      /// \param[in]  al    Assembly list for the central element.
      /// \return     Number of shape functions in the extended shapeset (sum of central and neighbor elems' local counts).
      ///
      ExtendedShapeset* create_extended_asmlist(const Space<Scalar>* space, AsmList<Scalar>* al);
      ExtendedShapeset* create_extended_asmlist_multicomponent(const Space<Scalar>* space, AsmList<Scalar>* al);

      /*** Methods for working with quadrature on the active edge. ***/

      /// Sets the quadrature order to be used for obtaining integration points and weights in both neighbors.
      ///
      /// \param[in] order  The integration order.
      ///
      void set_quad_order(int order);

      /// Get the integration pseudo-order for the active edge (assumes the true order has been set by \c set_quad_order).
      ///
      /// \param[in] on_neighbor  If true, order is returned for the neighbor el. (using its local active edge number).
      /// \return    The edge "pseudo-order" (determined by the true order and local edge number).
      ///
      int get_quad_eo(bool on_neighbor = false) const;

      /*** Methods for retrieving additional information about the neighborhood. ***/

      /// Return the number of elements (neighbors) adjacent to the central element across the common (active) edge.
      ///
      /// \return The number of neighbors.
      ///
      int get_num_neighbors() const;

      /// Retrieve all neighbor elements.
      ///
      /// \return pointer to the vector of neighboring elements.
      ///
      const Hermes::vector<Element*>* get_neighbors() const;

      /// Frees the memory occupied by the extended shapeset.
      void clear_supported_shapes();

      /// Function that sets the variable ignore_errors. See the variable description.
      void set_ignore_errors(bool value);

      /// This class represents the extended shapeset, consisting of shape functions from both the central element and
      /// current neighbor element, extended by zero to the union of these elements.
      class ExtendedShapeset
      {
      public:
        /// Constructor.
        ///
        /// \param[in]  neighborhood  Neighborhood on which the extended shapeset is defined.
        /// \param[in]  central_al    Assembly list for the currently assembled edge on the central element.
        /// \param[in]  space         Space from which the neighbor's assembly list will be obtained.
        ///
        ExtendedShapeset(NeighborSearch<Scalar>* neighborhood, AsmList<Scalar>* central_al, const Space<Scalar>*space);

        ExtendedShapeset(const ExtendedShapeset & other);

        /// Destructor.
        ~ExtendedShapeset();

        void free_central_al();

        /// Create assembly list for the extended shapeset by joining central and neighbor element's assembly lists.
        void combine_assembly_lists();

        /// Update the extended shapeset when active segment or active edge is changed (i.e. there will be a new neighbor).
        ///
        /// \param[in]  neighborhood  Neighborhood on which the extended shapeset is defined.
        /// \param[in]  space         Space from which the neighbor's assembly list will be obtained.
        ///
        void update(NeighborSearch* neighborhood, const Space<Scalar>* space);

      public:
        int cnt;  ///< Number of shape functions in the extended shapeset.
        int *dof; ///< Array of global DOF numbers of shape functions in the extended shapeset.

        bool has_support_on_neighbor(unsigned int index) const;

        AsmList<Scalar>* central_al;                    ///< Assembly list for the currently assembled edge on the central elem.
        AsmList<Scalar>* neighbor_al;                   ///< Assembly list for the currently assembled edge on the neighbor elem.

        friend class NeighborSearch; // Only a NeighborSearch is allowed to create an ExtendedShapeset.
      };

      /*** Neighborhood information. ***/
      /// Structure containing all the needed information about the active edge from the neighbor's side.
      class NeighborEdgeInfo
      {
      public:
        NeighborEdgeInfo() : local_num_of_edge(-1), orientation(-1) {};

        int local_num_of_edge;  ///< Local number of the edge on neighbor element.
        int orientation;        ///< Relative orientation of the neighbor edge with respect to the active edge
                                ///< (0 - same orientation, 1 - reverse orientation).
      };

      /// When creating sparse structure of a matrix using this class, we want to ignore errors
      /// and do nothing instead when set_active_edge() function is called for a non-boundary edge.
      bool ignore_errors;

      /// Returns the current active segment.
      int get_active_segment() const;

      /// Sets the active segment, neighbor element, and neighbor edge accordingly.
      void set_active_segment(unsigned int index);

      /// Returns the current neighbor element according to the current active segment.
      Element* get_neighb_el() const;

      /// Returns the current active neighbor edge according to the current active segment.
      NeighborEdgeInfo get_neighbor_edge() const;

      /// Returns the number(depth) of the current central transformation according to the current active segment.
      unsigned int get_central_n_trans(unsigned int index) const;

      /// Returns the current central transformations according to the current active segment.
      unsigned int get_central_transformations(unsigned int index_1, unsigned int index_2) const;

      /// Returns the number(depth) of the current neighbor transformation according to the current active segment.
      unsigned int get_neighbor_n_trans(unsigned int index) const;

      /// Returns the current neighbor transformations according to the current active segment.
      unsigned int get_neighbor_transformations(unsigned int index_1, unsigned int index_2) const;

      /// Transformations of an element to one of its neighbors.
      struct Transformations
      {
      public:
        static const int max_level = Transformable::H2D_MAX_TRN_LEVEL; ///< Number of allowed transformations (or equiv. number of neighbors
                                                                       ///< in a go-down neighborhood) - see Transformable::push_transform.

        unsigned int transf[max_level];   ///< Array holding the transformations at subsequent levels.
        unsigned int num_levels;          ///< Number of transformation levels actually used in \c transf.

        Transformations();
        Transformations(const Transformations* t);
        Transformations(const Hermes::vector<unsigned int>& t);

        void copy_from(const Hermes::vector<unsigned int>& t);

        void copy_from(const Transformations* t);

        void copy_to(Hermes::vector<unsigned int>* t);

        void reset();

        void strip_initial_transformations(unsigned int number_of_stripped);

        void apply_on(Transformable* tr) const;

        void apply_on(const Hermes::vector<Transformable*>& tr) const;

        template<typename T> friend class NeighborSearch;
        template<typename T> friend class KellyTypeAdapt;
        template<typename T> friend class Adapt;
        template<typename T> friend class Func;
        template<typename T> friend class DiscontinuousFunc;
        template<typename T> friend class DiscreteProblem;
        template<typename T> friend class DiscreteProblemLinear;
      };

    private:

      const Mesh* mesh;

      /*** Transformations. ***/

      LightArray< Transformations* > central_transformations;     ///< Array of transformations of the central element to each neighbor
                                                                  ///< (in a go-down neighborhood; stored as on \c Transformation structure
                                                                  ///< for each neighbor).
      LightArray< Transformations* > neighbor_transformations;    ///< Array of transformations of the neighbor to the central element (go-up).

      uint64_t original_central_el_transform;                  ///< Sub-element transformation of any function that comes from the
                                                               ///< assembly, before transforms from \c transformations are pushed to it.

      /*** Significant objects of the neighborhood. ***/
      Element* central_el;          ///< Central (currently assembled) element.
      Element* neighb_el;           ///< Currently selected neighbor element (on the other side of active segment).

      int active_edge;               ///< Local number of the currently assembled edge, w.r.t. the central element.
      NeighborEdgeInfo neighbor_edge;///< Assembled edge, w.r.t. the element on the other side.
      int active_segment;            ///< Part of the active edge shared by central and neighbor elements.

      Hermes::vector<NeighborEdgeInfo> neighbor_edges;   ///< Active edge information from each neighbor.
      Hermes::vector<Element*> neighbors;                ///< Vector with pointers to the neighbor elements.
      unsigned int n_neighbors;                          ///< Number of neighbors (>1 for a go-down neighborhood, 1 otherwise).

      /// Possible neighborhood types, according to which way we went on the neighbor element in order to get to the
      /// other side of the neighbor. The way is characterized by transformations needed to be pushed either on the
      /// central or the neighbor element.
      enum NeighborhoodType
      {
        H2D_DG_NOT_INITIALIZED = -1,
        H2D_DG_NO_TRANSF = 0, ///< Active edge has same-sized active elements on either side. Central element already
        ///< matches the neighbor, no transformation is needed.
        H2D_DG_GO_DOWN = 1,   ///< Central element is bigger than the neighbor element, we need to go down to its part with
        ///< the same size as the neighbor. Transformations will be pushed on the central element.
        H2D_DG_GO_UP = 2      ///< Neighbor element is bigger than the central element, we need to go up to the matching
        ///< parent. Corresponding transformations will be pushed in reverse order to the neighbor.
      };
      NeighborhoodType neighborhood_type;

      /// Find the neighbor element to a smaller central element.
      ///
      /// Central element is neccessarily a descendant of one or more inactive elements in this case. We go up
      /// through these parents and check their edge with the same local number as the active edge. For each
      /// inactive intermediate parent, this edge will not be used on the mesh (\c peek_edge_node return NULL).
      /// Once a used edge is found, its actual owner is the active neighbor element, but it shares it with the
      /// parent of the central element we were looking for. Transformation of the central element to this parent
      /// is determined from the sequence of middle vertices of the intermediate parent edges. However, we actually use
      /// the inverse transformation on the neighbor element.
      ///
      /// \param[in] elem             Pointer to a parent of the element from previous step.
      /// \param[in] orig_vertex_id   Array of oriented vertices of the active edge.
      /// \param[in] par_mid_vertices Array of vertices between those in \c orig_vertex_id visited on the way up.
      /// \param[in] n_parents        Number of intermediate parents visited on the way up.
      ///
      void find_act_elem_up( Element* elem, int* orig_vertex_id, Node** par_mid_vertices, int n_parents);

      /// Find all neighbors to a bigger central element.
      ///
      /// This is a recursive bisection of the active edge, until its segment is found that is also used as another
      /// edge on the mesh - by one of the active neighbor elements. The sequence of visited middle vertices is used
      /// to define the transformation path for the central element through those of its (virtual) sub-elements that
      /// lead to one completely adjacent to the found neighbor. Then we go back in the recursion tree and continue
      /// in another branch down to the next neighbor, until all of them are found.
      ///
      /// \param[in] vertex             Pointer to a middle vertex of the edge from previous step.
      /// \param[in] bounding_verts_id  Array of id's of vertices bounding the edge from previous step. They are passed
      ///                               explicitly in order to keep information about the orientation of the edge.
      /// \param[in] sons               Array that identifies virtual sons of the central element that must be visited in
      ///                               order to get to one matching the actual neighbor. This defines the row of the
      ///                               \c transformations array corresponding to that neighbor.
      /// \param[in] n_sons             Number of sons that lead to the current neighbor's counterpart.
      ///
      void find_act_elem_down( Node* vertex, int* bounding_verts_id, int* sons, unsigned int n_sons);

      /// Determine relative orientation of the neighbor edge w.r.t. the active edge.
      ///
      /// When this method is called from \c find_act_elem_up, \c bounding_vert1 and \c bounding_vert2 represent vertices
      /// bounding the edge of the central element's first inactive parent which completely matches the neighbor. Argument
      /// \c segment is always zero in this case.
      ///
      /// When this method is called from \c find_act_elem_down, \c bounding_vert1 and \c bounding_vert2 represent vertices
      /// which bound the edge of one of the inactive parents of the active neighbor element, with a <em>middle vertex</em>
      /// in between. Which of these vertices is the startpoint and which the endpoint of the neighbor's active edge is
      /// determined by the argument \c segment: neighbor edge spans from \c bounding_vert1 to the <em>middle vertex</em> if
      /// <tt>segment == 0</tt>, or from the <em>middle vertex</em> to \c bounding_vert2 if <tt>segment == 1</tt>.
      ///
      /// \param[in] bounding_vert1 ID of one endpoint of the edge (see above).
      /// \param[in] bounding_vert2 ID of the other endpoint of the edge (see above).
      ///
      /// \return 1 if the orientation of the neighbor's edge is reversed w.r.t. the central el.'s edge,
      ///         0 otherwise.
      ///
      int neighbor_edge_orientation(int bounding_vert1, int bounding_vert2, int segment) const;

      /// Cleaning of internal structures before a new edge is set as active.
      void reset_neighb_info();

      /*** Quadrature on the active edge. ***/
      Quad2D* quad;

      int central_quad_order;  ///< Quadrature data of the active edge with respect to the central element.
      int neighb_quad_order;   ///< Quadrature data of the active edge with respect to the element on the other side.

      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class Adapt;
      template<typename T> friend class Func;
      template<typename T> friend class DiscontinuousFunc;
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
    };
  }
}
#endif