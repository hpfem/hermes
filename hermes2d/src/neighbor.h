#ifndef NEIGHBOR_H_
#define NEIGHBOR_H_

#include "../hermes_common/common.h"
#include "mesh/mesh.h"
#include "quadrature/quad.h"
#include "function/solution.h"
#include "weakform/forms.h"
#include "mesh/refmap.h"

/*** Global sanity checks for the NeighborSearch class ***/

#define ensure_active_edge(obj) \
  assert_msg( obj->active_edge >= 0 && obj->active_edge < obj->central_el->nvert,\
              "Wrong active edge or active edge not set." )

#define ensure_active_segment(obj) \
  assert_msg( obj->active_segment >= 0 && obj->active_segment < obj->n_neighbors,\
              "Active segment of the active edge has not been set or exceeds the number of neighbors" )

#define ensure_central_pss_rm(obj) \
  assert_msg( obj->central_pss != NULL && obj->central_pss->get_active_element() == obj->central_el &&\
              obj->central_rm != NULL && obj->central_rm->get_active_element() == obj->central_el,\
              "Precalculated shapeset and refmap have not been attached or have a wrong active element." )

#define ensure_central_rm(obj) \
  assert_msg( obj->central_rm != NULL && obj->central_rm->get_active_element() == obj->central_el,\
              "Reference mapping has not been attached or has a wrong active element." )

//TODO: Add a test for overshooting the maximum allowed edge order.
#define ensure_set_quad_order(obj) \
  assert_msg( obj.eo > 0 && obj.np > 0 && obj.pt != NULL, \
              "Quadrature order must be set before calculating geometry and function values." )

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
 * If values of a test function will be needed, a precalculated shapeset that will receive the transformations
 * must also be attached via \c attach_pss_and_rm (so that a call to \c detach_pss_and_rm may be performed afterwards
 * to return the pss to its original state).
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

class HERMES_API NeighborSearch
{
public:

  /// Constructor.
  ///
  /// \param[in]  el    Central element of the neighborhood (current active element in the assembling procedure).
  /// \param[in]  mesh  Mesh on which we search for the neighbors.
  ///
  NeighborSearch(Element* el, Mesh* mesh);

/*** Methods for changing active state for further calculations. ***/

  /// Set active edge and compute all information about the neighbors.
  ///
  /// In particular, it fills the \c neighbors and \c neighbor_edges vectors and the \c transformations array used
  /// for transforming either the central or the neighboring elements to the same size. It also sets \c neighborhood_type,
  /// according to whether we can find a vertex in the middle of the active edge - if we can, then we go down and the
  /// corresponding transformations will be filled for the bigger central element. This is performed by the method
  /// \c find_act_elem_down. Otherwise, we go up and will later push the transformations from \c transfomations to
  /// functions on the bigger neighboring element. Way up is realized by the method \c find_act_elem_up.
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
  void set_active_edge(int edge, bool ignore_visited_segments = true);

  /// Set the part of active edge shared by the central element and a given neighbor.
  ///
  /// \param[in] neighbor   Number of the neighbor as enumerated in \c set_active_edge (i.e. its index in \c neighbors).
  /// \param[in] with_neighbor_pss  If true, also creates and/or sets the transformation of a PrecalcShapeset on the
  ///                               given neighbor element; this is then used for forming the extended shapeset.
  /// \return false if NeighborSearch is set to ignore already visited segments and the given one has already been visited
  /// \return true otherwise.
  ///
  bool set_active_segment(int neighbor, bool with_neighbor_pss = true);

/*** Methods for computing values of external functions. ***/

  /// Fill arrays with function values / derivatives at both sides of the active segment of active edge.
  /// Assumes that integration order has been set by \c set_quad_order.
  ///
  /// \param[in] fu MeshFunction whose values are requested.
  /// \return Pointer to a discontinuous function allowing to access the values from each side of the active edge.
  ///
  DiscontinuousFunc<scalar>* init_ext_fn(MeshFunction* fu);

  /// Initialize the polynomial orders of the given function at both sides of the active segment of active edge.
  ///
  /// \param[in] fu MeshFunction whose order is requested.
  /// \return Pointer to a discontinuous function allowing to access the order from each side of the active edge.
  ///
  DiscontinuousFunc<Ord>* init_ext_fn_ord(Solution* fu);

  /// Initialize the polynomial orders of the given function at both sides of the active segment of active edge.
  ///
  /// The order can be computed more precisely for Solutions than for general MeshFunctions, since we know the individual
  /// orders of the solution on both the central and neighbor elements. Note that we could even obtain the appropriate
  /// axial orders for the active edge, if we could obtain the pointer to the Solution's approximation space (this is
  /// currently impossible).
  ///
  /// \param[in] fu Solution whose order is requested.
  /// \return Pointer to a discontinuous function allowing to access the order from each side of the active edge.
  ///
  DiscontinuousFunc<Ord>* init_ext_fn_ord(MeshFunction* fu);

/*** Methods for working with shape functions. ***/

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
  int create_extended_shapeset(Space* space, AsmList* al);

  /// Assign a precalculated shapeset and associated reference mapping to the central element.
  ///
  /// Apart from setting the corresponding pointers, it also stores the original transformation of the pss and refmap.
  /// This is required for a <em>way down</em>-type neighborhood (big central element with several neighbors), where
  /// transforms from the \c transfomations array need to be pushed to the active element's pss and refmap, in order to
  /// later allow resetting the pss and refmap to their original active element.
  ///
  /// \param[in]  pss Pointer to a PrecalcShapeset object (typically from the assembling procedure).
  /// \param[in]  rm  Pointer to a RefMap object associated with the pss.
  ///
  void attach_pss_and_rm(PrecalcShapeset* pss, RefMap* rm);

  /// Assign a reference mapping to the central element.
  ///
  /// If geometric data about the neighborhood is needed (see function \c init_geometry), a reference mapping with
  /// the correctly pushed <em>way down</em> transformations is required. If the extended shapesets will not be used,
  /// this method may be used for this purpose instead of \c attach_pss_and_rm.
  ///
  /// \param[in]  rm  Pointer to a RefMap object associated with the pss.
  ///
  void attach_rm(RefMap* rm);

  /// Restore the transformation set for central element's pss and refmap to that before their attachment to the
  /// NeighborSearch.
  void detach_pss_and_rm();

  /// Restore the transformation set for central element's refmap to that before its attachment to the NeighborSearch.
  void detach_rm();

/*** Methods for working with quadrature on the active edge. ***/

  /// Sets the quadrature order to be used for obtaining integration points and weights in both neighbors.
  ///
  /// \param[in] order  The integration order.
  ///
  void set_quad_order(int order);

  /// Get quadrature points on the active edge (assumes that integration order has been set by \c set_quad_order).
  ///
  /// \param[in] on_neighbor  If true, points are returned for the neighbor el. (using its local active edge number).
  /// \return    Pointer to the array of quadrature points.
  ///
  double3* get_quad_pt(bool on_neighbor = false);

  /// Get the integration pseudo-order for the active edge (assumes the true order has been set by \c set_quad_order).
  ///
  /// \param[in] on_neighbor  If true, order is returned for the neighbor el. (using its local active edge number).
  /// \return    The edge "pseudo-order" (determined by the true order and local edge number).
  ///
  int get_quad_eo(bool on_neighbor = false);


  /// Get the number of integration points on the active edge (assumes that integration order has been set by
  /// \c set_quad_order).
  ///
  /// \return  The number of integration points corresponding to the set order. Note that it is the same for both
  ///          neighboring elements.
  ///
  int get_quad_np();

  /// Get the precalculated shapeset on the central element.
  PrecalcShapeset * get_pss() {return this->central_pss;}

  /// Get the refmap on the central element.
  RefMap * get_rm() {return this->central_rm;}


/*** Methods for getting geometry and integration data. ***/

  /// Initialize the geometry data for the active segment.
  ///
  /// Calculates physical points, normals, etc. at the integration points, using cached values whenever possible.
  /// This function assumes that integration order has been set by \c set_quad_order.
  ///
  /// \param[in,out]  ext_cache_e   Geometry cache from the assembling procedure. A key to this cache is the sole
  ///                               edge integration pseudo-order (which contains information about the true order
  ///                               and local edge number).
  ///
  ///                               If no transformations are to be pushed on the central element (not a
  ///                               <em>way-down</em> neighborhood), the returned value is either fetched from this
  ///                               cache or a new entry is created and the cache is updated.
  ///
  ///                               When central element requires transformation to a son matching the current
  ///                               neighbor, an internal cache is used which is queried by both edge pseudo-order
  ///                               and active segment (which uniquely determines the transformation).
  ///
  ///                               If ext_cache_e == NULL, no cache is used whatsoever.
  ///
  /// \param[in] ep Active edge data required by the \c init_geom_surf function.
  /// \return Pointer to a structure holding the geometry data as well as diameter, id and marker of elements on both
  ///         sides of the edge.
  ///
  Geom<double>* init_geometry(Geom< double >** ext_cache_e, SurfPos* ep);

  /// Initialize the products of Jacobian and quadrature weights.
  ///
  /// This function uses both internal cache and cache from the assembling procedure, in the same way as
  /// \c init_geometry.
  /// This function assumes that integration order has been set by \c set_quad_order.
  ///
  /// \param[in,out]  ext_cache_jwt   Cache from the assembling procedure (may be set to NULL to bypass the caching).
  /// \return         The array of Jacobian * weights.
  ///
  double* init_jwt(double** ext_cache_jwt);

/*** Methods for retrieving additional information about the neighborhood. ***/

  /// Return the number of elements (neighbors) adjacent to the central element across the common (active) edge.
  ///
  /// \return The number of neighbors.
  ///
  int get_num_neighbors() { return n_neighbors; }

  /// Retrieve all neighbor elements.
  ///
  /// \return pointer to the vector of neighboring elements.
  ///
  std::vector<Element*>* get_neighbors() { return &neighbors; }

  /// Retrieve the central element in the neighborhood defined by current active segment.
  ///
  /// \return pointer to the central element.
  ///
  Element* get_current_central_element()
  {
    ensure_active_segment(this);
    return central_el;
  }

  /// Retrieve the neighbor element in the neighborhood defined by current active segment.
  ///
  /// \return pointer to the neighbor of the central element.
  ///
  Element* get_current_neighbor_element()
  {
    ensure_active_segment(this);
    return neighb_el;
  }

  /// Get transfomations that must be pushed to a Transormable either on the central or neighbor element in order for
  /// its values from both sides to match.
  ///
  /// \param[in] segment Part of active edge shared by the central element and selected neighbor.
  /// \return the array of transformation sub-indices.
  ///
  int* get_transformations(int segment) {
    ensure_active_edge(this);
    assert(segment >= 0 && segment < n_neighbors);
    return transformations[segment];
  }

  /// Return local number of the selected active edge segment relatively to the neighboring element.
  ///
  /// \param[in] segment Part of active edge shared by the central element and selected neighbor.
  /// \return Local number of the active edge in the neighboring element.
  ///
  int get_neighb_edge_number(int segment);

  /// Return orientation of the selected active edge segment with respect to both adjacent elements.
  ///
  /// \param[in] active_segment  Part of active edge shared by the central element and selected neighbor.
  /// \return 0 if orientation of the common edge is the same on both adjacent elements,
  /// \return 1 otherwise.
  ///
  int get_neighb_edge_orientation(int segment);

/*** Methods for cleaning up. ***/

  /// Frees the memory occupied by the extended shapeset.
  void clear_supported_shapes() {
    if (supported_shapes != NULL) delete supported_shapes; supported_shapes = NULL;
  }

  /// Frees the memory occupied by the neighbor element's precalc. shapeset and reference mapping.
  void clear_neighbor_pss() {
    if (neighb_pss != NULL) delete neighb_pss;
    if (neighb_rm != NULL) delete neighb_rm;
  }

  /// Frees the memory occupied by the internal geometric and jac*wt caches.
  void clear_caches() {
    for (std::map<Key, Geom<double>*, Compare>::iterator it = cache_e.begin(); it != cache_e.end(); it++)
    {
      (it->second)->free();
      delete it->second;
    }
    cache_e.clear();

    for (std::map<Key, double*, Compare>::iterator it = cache_jwt.begin(); it != cache_jwt.end(); it++)
      delete [] it->second;
    cache_jwt.clear();
  }

  /// Destructor.
  ~NeighborSearch();

  /// This variable has the meaning how many neighbors have been used for a single edge so far,
  /// and it is used for the allocation of the arrays NeighborSearch::transformations and NeighborSearch::n_trans.
  static int max_neighbors;

  /// Function that sets the variable ignore_errors. See the variable description.
  void set_ignore_errors(bool value) {this->ignore_errors = value;};

  /// Functionality for caching of NeighborSearch instances.
  struct MainKey
  {
    int element_id, isurf;
    MainKey(int element_id, int isurf) : element_id(element_id), isurf(isurf) {};
  };

  struct MainCompare
  {
    bool operator()(MainKey a, MainKey b) const
    {
      if (a.element_id < b.element_id)
        return true;
      else if (a.element_id > b.element_id)
        return false;
      else
        return (a.isurf < b.isurf);
    }
  };

  /// Two caches of NeighborSearch class instances.
  static std::map<MainKey, NeighborSearch*, MainCompare> main_cache_m;
  static std::map<MainKey, NeighborSearch*, MainCompare> main_cache_n;

  /// Method to empty the above caches.
  static void empty_main_caches();

private:

  Mesh* mesh;

/*** Transformations. ***/

  static const int max_n_trans = 20;              ///< Number of allowed transformations (or equiv. number of neighbors
                                                  ///< in a go-down neighborhood) - see Transformable::push_transform.
                                                  ///< TODO: Revise this for multimesh.
  std::vector<int *> transformations;             ///< Vector of transformations of the central element to each neighbor
                                                  ///< (in a go-down neighborhood; stored row-wise for each neighbor)
                                                  ///< or of the neighbor to the central element (go-up).
  std::vector<int> n_trans;                       ///< Number of transforms stored in each row of \c transformations.
  unsigned int original_central_el_transform;              ///< Sub-element transformation of any function that comes from the
                                                  ///< assembly, before transforms from \c transformations are pushed
                                                  ///< to it.


/*** Significant objects of the neighborhood. ***/

  Element* central_el;          ///< Central (currently assembled) element.
  Element* neighb_el;           ///< Currently selected neighbor element (on the other side of active segment).
  RefMap* central_rm;           ///< Reference mapping of the central element.
  RefMap* neighb_rm;            ///< Reference mapping of the neighbor element.
  PrecalcShapeset *central_pss; ///< Precalculated shapeset on the central element.
  PrecalcShapeset *neighb_pss;  ///< Precalculated shapeset on the neighbor element.


/*** Neighborhood information. ***/

  int active_edge;              ///< Local number of the currently assembled edge, w.r.t. the central element.
  int neighbor_edge;            ///< Local number of the currently assembled edge, w.r.t. the element on the other side.
  int active_segment;           ///< Part of the active edge shared by central and neighbor elements.
  bool ignore_visited_segments; ///< True if each edge should be assembled only from one side.

  /// Structure containing all the needed information about the active edge from the neighbor's side.
  struct NeighborEdgeInfo
  {
    int local_num_of_edge;  ///< Local number of the edge on neighbor element.
    int orientation;        ///< Relative orientation of the neighbor edge with respect to the active edge
                            ///< (0 - same orientation, 1 - reverse orientation).
    NeighborEdgeInfo() : local_num_of_edge(-1), orientation(-1) {};
  };

  std::vector<NeighborEdgeInfo> neighbor_edges;   ///< Active edge information from each neighbor.
  std::vector<Element*> neighbors;                ///< Vector with pointers to the neighbor elements.
  int n_neighbors;                                ///< Number of neighbors (>1 for a go-down neighborhood, 1 otherwise).

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
  /// inactive intermediate parent,this edge will not be used on the mesh (\c peek_edge_node return NULL).
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
  void find_act_elem_down( Node* vertex, int* bounding_verts_id, int* sons, int n_sons);


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
  int neighbor_edge_orientation(int bounding_vert1, int bounding_vert2, int segment);

  /// Cleaning of internal structures before a new edge is set as active.
  void reset_neighb_info();


/*** Quadrature on the active edge. ***/

  Quad2D* quad;

  /// Structure specifying the quadrature on the active segment.
  struct QuadInfo
  {
    int eo;       ///< Edge pseudo-order.
    int np;       ///< Number of integration points on the active edge (segment).
    double3* pt;  ///< Array of physical points and weights on the active segment.

    QuadInfo() : eo(0), np(0), pt(NULL) {};

    void init(Quad2D* quad, int eo) {
      this->eo = eo;
      this->np = quad->get_num_points(eo);
      this->pt = quad->get_points(eo);
    }
  };
  QuadInfo central_quad;  ///< Quadrature data of the active edge with respect to the central element.
  QuadInfo neighb_quad;   ///< Quadrature data of the active edge with respect to the element on the other side.


/*** Caching data for the poosibly transformed central element ***/

  struct Key
  {
    int order, active_segment;
    Key(int order, int active_segment) : order(order), active_segment(active_segment) {};
  };

  struct Compare
  {
    bool operator()(Key a, Key b) const
    {
      if (a.order < b.order) return true;
      else if (a.order > b.order) return false;
      else return (a.active_segment < b.active_segment);
    }
  };

  std::map<Key, Geom<double>*, Compare> cache_e;
  std::map<Key, double*, Compare> cache_jwt;

/*** Geometric calculations. ***/

  double* calculate_jwt(int edge_order);

public:
  /// This class represents the extended shapeset, consisting of shape functions from both the central element and
  /// current neighbor element, extended by zero to the union of these elements.
  class ExtendedShapeset
  {
    public:
      /// This class represents a concrete shape function from the extended shapeset.
      class ExtendedShapeFunction
      {
        private:
          NeighborSearch *neibhood;   ///< Neighborhood on which the extended shape function is defined.
          PrecalcShapeset *active_pss;///< Pointer to the precalc. shapeset from which the active shape function
                                      ///< is drawn. Depending on the index of the active shape fn. within the
                                      ///< extended shapeset, it may be the local shapeset precalculated either on the
                                      ///< central element or on the neighbor element (see \c activate).
          RefMap *active_rm;          ///< Reference mapping associated with the \c active_pss.

          int order;                  ///< Polynomial order of the active shape function on the active edge
                                      ///< (for quads vertical or horizontal as appropriate).
          bool support_on_neighbor;   ///< True if the active shape fn. is nonzero on the neighbor element (and thus
                                      ///< zero on the central element), false otherwise.
          bool reverse_neighbor_side; ///< True if the orientation of the neighbor edge w.r.t. the active edge is
                                      ///< reversed.

          ///
          void activate(int index, AsmList* central_al, AsmList* neighb_al);

          /// Constructor.
          /// \param[in] neighborhood Neighborhood on which the extended shape function is defined.
          ///
          ExtendedShapeFunction(NeighborSearch* neighborhood) : neibhood(neighborhood) {};

        public:

          /*** Assembly list information for the active shape ***/

          int idx;      ///< shape function index
          int dof;      ///< basis function number
          scalar coef;  ///< coefficient

          /*** Get methods ***/

          RefMap* get_activated_refmap() { return active_rm; }
          PrecalcShapeset* get_activated_pss() { return active_pss; }
          int get_fn_order() { return order; }

          /// Get function values, derivatives, etc. of the extended shape function in quadrature points along the
          /// active segment.
          ///
          /// \param[in,out]  ext_cache_fn  Reference to the cache of shape functions obtained from Discrete/FeProblem.
          /// \return         Pointer to a \c DiscontinuousFunc object which may be queried for values on either side
          ///                 of the active segment.
          ///
          DiscontinuousFunc<double>* get_fn(std::map< PrecalcShapeset::Key, Func< double >*, PrecalcShapeset::Compare >& ext_cache_fn);

          /// Get \c DiscontinuousFunc representation of the active shape function's polynomial order.
          DiscontinuousFunc<Ord>* get_fn_ord() {
            int inc = (active_pss->get_num_components() == 2) ? 1 : 0;
            return extend_by_zero( init_fn_ord(this->order + inc) );
          }

          /// Extend by zero the active shape function to the other element.
          ///
          /// It uses the \c support_on_neighbor and \c reverse_neighbor_side of the activated shape function to
          /// correctly set the zero and the non-zero part of the resulting discontinuous shape function.
          ///
          /// \param[in]  fu  \c Func representation of the active shape function as obtained from
          ///                 \c DiscreteProblem::get_fn
          /// \return     Pointer to a \c DiscontinuousFunc object which may be queried for values on either side
          ///             of the discontinuity.
          ///
          DiscontinuousFunc<double>* extend_by_zero(Func<double>* fu) {
            return new DiscontinuousFunc<double>(fu, support_on_neighbor, reverse_neighbor_side);
          }

          /// Extend by zero the \c Func representation of the active shape's polynomial order to the other element.
          DiscontinuousFunc<Ord>* extend_by_zero(Func<Ord>* fu) {
            return new DiscontinuousFunc<Ord>(fu, support_on_neighbor);
          }

          // Only an ExtendedShapeset is allowed to create an ExtendedShapeFunction.
          friend class NeighborSearch::ExtendedShapeset;
      };

      /// Set active element, push neccessary transforms and set the active shape for \c active_pss.
      ///
      /// If the given index is lower than the number of shape functions in central element's assembly list for active
      /// edge, the corresponding function from central element's precalculated shapeset is selected and go-down
      /// transformations (if any) are pushed to it. The resulting extended shape function will be equal to this
      /// function on the central element and will be zero on the neighbor.
      ///
      /// If the index is greater than the number of shape functions in central element's assembly list and lower than
      /// that in neighbor element's assembly list, active shape function is drawn from the neighbor element's pss.
      /// The resulting extended shape function will attain zero values on the central element and be equal to the
      /// active shape function on the neighbor.
      ///
      /// \param[in]  index index of selected extended shape function
      ///                   (ranging from zero to <tt>central_al->cnt + neighbor_al->cnt</tt>)
      /// \return     Pointer to the active extended shape function.
      ///
      ExtendedShapeFunction* get_extended_shape_fn(int index) {
        active_shape->activate(index, central_al, neighbor_al);
        return active_shape;
      }

    private:
      ExtendedShapeFunction *active_shape;    ///< Extended shape function with activated \c active_pss.
      AsmList* central_al;                    ///< Assembly list for the currently assembled edge on the central elem.
      AsmList* neighbor_al;                   ///< Assembly list for the currently assembled edge on the neighbor elem.

      /// Create assembly list for the extended shapeset by joining central and neighbor element's assembly lists.
      void combine_assembly_lists();

      /// Update the extended shapeset when active segment or active edge is changed (i.e. there will be a new neighbor).
      ///
      /// \param[in]  neighborhood  Neighborhood on which the extended shapeset is defined.
      /// \param[in]  space         Space from which the neighbor's assembly list will be obtained.
      ///
      void update(NeighborSearch* neighborhood, Space* space) {
        delete [] this->dof;
        space->get_boundary_assembly_list(neighborhood->neighb_el, neighborhood->neighbor_edge, neighbor_al);
        combine_assembly_lists();
      }

      /// Constructor.
      ///
      /// \param[in]  neighborhood  Neighborhood on which the extended shapeset is defined.
      /// \param[in]  central_al    Assembly list for the currently assembled edge on the central element.
      /// \param[in]  space         Space from which the neighbor's assembly list will be obtained.
      ///
      ExtendedShapeset(NeighborSearch* neighborhood, AsmList* central_al, Space *space);

      /// Destructor.
      ~ExtendedShapeset() {
        delete [] dof; delete active_shape; delete neighbor_al;
      }

    public:
      int cnt;  ///< Number of shape functions in the extended shapeset.
      int *dof; ///< Array of global DOF numbers of shape functions in the extended shapeset.

      friend class NeighborSearch; // Only a NeighborSearch is allowed to create an ExtendedShapeset.
  };

  /// When creating sparse structure of a matrix using this class, we want to ignore errors
  /// and do nothing instead when set_active_edge() funciton is called for a non-boundary edge.
  bool ignore_errors;
};

typedef NeighborSearch::ExtendedShapeset::ExtendedShapeFunction* ExtendedShapeFnPtr;


#endif /* NEIGHBOR_H_ */
