#ifndef NEIGHBOR_H_
#define NEIGHBOR_H_

#include "common.h"
#include "mesh.h"
#include "quad.h"
#include "solution.h"
#include "forms.h"
#include "refmap.h"

// Global sanity checks for the NeighborSearch class

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

/*!\class NeighborSearch neighbor.h "src/neighbor.h"
 * \brief This class serves to find and offer all information from neighbors at a given(active) element and its concrete edge.
 *
 * How works the finding neighbors
 * We will call the active  element as a central element.
 * If we have irregular mesh, there can be three options in relation of central element and neighbor element at common edge.
 * First, the neighbor is same "size" as central, so the edge has active elements on both sides. This option is tested by function get_neighbor().
 * Second, the neighbor is "bigger" then central. Then we have to go "way up" and use method finding_act_elem_up().
 * Third, the neighbor is "smaller", we have more neighbors against the edge. This solves "way down" by method finding_act_elem_down().
 * The choice between way up or way down is made by testing if we can find vertex in the middle of the common edge. If we can
 * then we go way down.

 * Also at every way we fill function values, derivatives and etc. of central and neighbor elements threw method set_fn_values(). Last step is
 * possible change of order of neighbor's function values to correspond function values of central element at same points.

 * We also need transform solution either on neighbor or central element to get points at correct part of the edge. We use method "push_transform"
 * and use only range [0-3]. These types of transformation are common for triangles and quads and choosing right transformation can
 * be derived from local numbers of edges.

 * For numbering and ordering of edges, vertices and sons of an element look into mesh.cpp
 */

class H2D_API NeighborSearch
{
public:

	///\brief Common constructor. 
	NeighborSearch(Element* el, Mesh* mesh, bool ignore_if_visited = true);
  
  // Methods for changing active state for further calculations.
  
	///\brief Set active edge and compute all information about the neighbors.
	/// \param[in] edge Local (element dependent) number of the edge.
	void set_active_edge(int edge);
  
  ///\brief Set the part of active edge shared by the central element and a given neighbor.
  /// \param[in] neighbor   Number of the neighbor as enumerated in set_active_edge (i.e. its index in the vector neighbors).
  /// \param[in] with_neighbor_pss  If true, also creates and/or sets the transformation of a PrecalcShapeset on the given 
  ///                               neighbor element; this is then used for forming the extended shapeset.
  /// \return false if NeighborSearch is set to ignore already visited segments and the one given has already been visited,
  /// \return true otherwise.
  bool set_active_segment(int neighbor, bool with_neighbor_pss = true);
  
  // Methods for computing values of external functions.

  ///\brief Fill arrays with function values / derivatives at both sides of the active segment of active edge.
  /// \param[in] fu MeshFunction whose values are requested.
  /// \return Pointer to a discontinuous function allowing to access the values from each side of the active edge.
  DiscontinuousFunc<scalar>* init_ext_fn(MeshFunction* fu);
  
  /// Initialize the polynomial orders of the given function at both sides of the active segment of active edge.
  /// \param[in] fu MeshFunction whose order is requested.
  /// \return Pointer to a discontinuous function allowing to access the order from each side of the active edge.
  DiscontinuousFunc<Ord>* init_ext_fn_ord(Solution* fu);
  
  ///\brief Initialize the polynomial orders of the given function at both sides of the active segment of active edge.
  /// The order can be computed more precisely for Solutions than for general MeshFunctions, since we know the individual 
  /// orders of the solution on both the central and neighbor elements. Note that we could even obtain the appropriate 
  /// axial orders for the active edge, if we could obtain the pointer to the Solution's approximation space (this is
  /// currently impossible).
  /// \param[in] fu Solution whose order is requested.
  /// \return Pointer to a discontinuous function allowing to access the order from each side of the active edge.
  DiscontinuousFunc<Ord>* init_ext_fn_ord(MeshFunction* fu);
  
  // Methods for working with shape functions.

  class ExtendedShapeset;
  ExtendedShapeset *supported_shapes; 
  
  int extend_attached_shapeset(Space* space, AsmList* al);
  void attach_pss(PrecalcShapeset* pss, RefMap* rm);
  void detach_pss();
  
  // Methods for working with quadrature on the active edge.
  
  void set_quad_order(int order)
  {
    ensure_active_segment(this);
    quad->set_mode(central_el->get_mode());
    central_quad.init(quad, quad->get_edge_points(active_edge, order));
    quad->set_mode(neighb_el->get_mode());
    neighb_quad.init(quad, quad->get_edge_points(neighbor_edge, order));
  }
  
  double3* get_quad_pt(bool on_neighbor = false) { return on_neighbor ? neighb_quad.pt : central_quad.pt; }
  int get_quad_eo(bool on_neighbor = false) { return on_neighbor ? neighb_quad.eo : central_quad.eo; }
  int get_quad_np() { return central_quad.np; }

  // Methods for working with quadrature on the active edge.

  Geom<double>* init_geometry(Geom< double >** ext_cache_e, EdgePos* ep);
  double* init_jwt(double** ext_cache_jwt);
  
  // Methods for retrieving additional information about the neighborhood.

  ///\brief Return the number of elements (neighbors) adjacent to the active (central) element across a common (active) edge.
  /// \return The number of neighbors.
  int get_number_of_neighbs() { return n_neighbors; }

  /// Return pointer to the vector of neighboring elements.
  std::vector<Element*>* get_neighbors() { return &neighbors; }

  ///\brief Get transfomations that must be pushed to a Transormable either on the central or neighbor element in order for its
  ///\brief values from both sides to match.
  /// \param[in] segment Part of active edge shared by the central element and selected neighbor.
  /// \return the array of transformation sub-indices.
  int* get_transformations(int segment) { 
    ensure_active_edge(this);
    assert(segment >= 0 && segment < n_neighbors);
    return transformations[segment]; 
  }

  ///\brief Return local number of the selected active edge segment relatively to the neighboring element.
  /// \param[in] segment Part of active edge shared by the central element and selected neighbor.
  /// \return Local number of the active edge in the neighboring element.
  int get_neighb_edge_number(int segment);

  ///\brief Return orientation of the selected active edge segment with respect to both adjacent elements.
  /// \param[in] active_segment  Part of active edge shared by the central element and selected neighbor.
  /// \return 0 if orientation of the common edge is the same on both adjacent elements,
  /// \return 1 otherwise.
  int get_neighb_edge_orientation(int segment);
  
  // Methods for cleaning up.
  
  void clear_supported_shapes() {  
    if (supported_shapes != NULL) delete supported_shapes; supported_shapes = NULL; 
  }
  void clear_neighbor_pss() {
    if (neighb_pss != NULL) delete neighb_pss;
    if (neighb_rm != NULL) delete neighb_rm;
  }
  void clear_geometry_cache() {
    for (std::map<Key, Geom<double>*, Compare>::iterator it = cache_e.begin(); it != cache_e.end(); it++)
      (it->second)->free();
    cache_e.clear();
  }
  
  ~NeighborSearch();
  
  
private:  
  static const int max_n_trans = 20;    ///< Number of allowed transformations, see "push_transform" in transform.h.

  bool ignore_if_visited;
	int n_neighbors; ///< Number of neighbors.
	
	Mesh* mesh;
 
  Quad2D* quad;
  struct QuadInfo
  { 
    int eo; 
    int np; 
    double3* pt;
    
    void init(Quad2D* quad, int eo) {
      this->eo = eo;
      this->np = quad->get_num_points(eo);
      this->pt = quad->get_points(eo);
    }
  };  
  QuadInfo central_quad, neighb_quad;
  
  RefMap *central_rm, *neighb_rm;
  PrecalcShapeset *central_pss, *neighb_pss;
  
	Element* central_el; ///< Central element.
	Element* neighb_el;  ///< Actual neighbor element we are working with,
     
	int transformations[max_n_trans][max_n_trans];	///< Table of transformations for all neighbors.
	int n_trans[max_n_trans]; ///< Number of transformations for every neighbor.
	int active_edge;          ///< Edge where we are searching for neighbors.
	int neighbor_edge;        ///< Edge of the working neighbor corresponding to active_edge.
  int active_segment;

	int way_flag; ///< This flag holds which way was used on the active edge. So is equal to one of the members of Trans_flag.

  int original_central_el_transform;

	
  

  /// This serves for distinguish which way was used for finding neighbors and according the way how are obtain values in method
  /// set_fn_values().
  ///
  enum Trans_flag{
    H2D_NO_TRANSF = 0,  ///< Don't use any transformation, the edge has on both sides active element.
    H2D_WAY_DOWN = 1,   ///< Transformation of central element, against the edge neighbor has some sons.
    H2D_WAY_UP = 2      ///< Transformation of neighbor element, central element is son.
  };
  
	/** \brief Method "way up" for finding neighbor element, from smaller to larger.
   * 
   * We use recurrence in this way.
   * If the neighbor is "bigger" then this means central element is descendant of some inactive elements. We go threw this parents and
   * stop when against an edge, which has same local number as the original edge, we have active element.
   * Important is that all sons have same orientation as parent, so local number of the edge is same.
   *
   * \param[in] elem The pointer to parent element of the element from previous step.
   * \param[in] orig_vertex_id Array containing oriented vertices of the active edge.
   * \param[in] road_vertices Array of vertices which we used in finding the active neighbor.
   * \param[in] n_road_vertices Number of vertices in array road_vertices.
  **/
	void finding_act_elem_up( Element* elem, int* orig_vertex_id, Node** road_vertices, int n_road_vertices);

	/** \brief Method "Way down" for finding neighbor elements, from larger to smaller.
	 *
   *  Again we use recurrence in this way. In every step we take middle vertex of the edge (starting with active edge). This vertex split the edge
   *  on two parts. On every part (an edge) we test if the new edge is active. If not, the middle vertex is found and the method is called
   *  again with this new vertex on this part.
   *
	 * \param[in] vertex The pointer to a vertex which was in the middle of the edge we were working with in previous step.
	 * \param[in] par_vertex_id Array containing id of vertices which are "parents" of the middle vertex (first parameter). They are passed separated because we need to conserve orientation of the vertices.
	 * \param[in] road Array which serves for storing codes of transformations (code is equal to number of son on which we will transform solution).
	 * \param[in] n_road Number of valid members in array road.
	 * \param[in] use_edge Just local name for active edge.
	 * \param[in] n_vert Number of vertices of the central element.
	**/
	void finding_act_elem_down( Node* vertex, int* par_vertex_id, int* road, int n_road);


  /// Structure containing all needed information about neighbor's edge. Initial values of both members are invalid.
	struct NeighborEdgeInfo
  {
		NeighborEdgeInfo(){
			local_num_of_edge = -1;
			orientation = -1;
		}

		/// Local number of the edge on neighbor element.
		int local_num_of_edge;

    /// Relative orientation of the neighbor edge. If equal to 0 then the neighbor edge has same orientation as the active edge,
    /// otherwise is equal to 1.
		int orientation;
	};
  
  /// Vector containing all neighbor edges information corresponding to active edge.
  std::vector<NeighborEdgeInfo> neighbor_edges;
  
  ///  Vector containing pointers to  all neighbors. 
  std::vector<Element*> neighbors; 
  

	/** Find the orientation of the neighbor edge in relation with the active edge. All input paramaters depend on in which way we
   * we are. For all ways parent1 and parent2 are ids of vertices oriented in same direction as are vertices of central element.
   * Parameter active_segment is 0 or 1 and actualy is used only in way down. In other ways is set to 0.
   * For way up parameters parent1 and parent2 are vertices which define the edge of inactive parent of central element. The edge has on other
   * side active neighbor.
   * For way down parent1 and parent2 are vertices of the edge, which is part of active edge. Both vertices belong to direct parent element
   * of active neighbor element. This means that one of them for sure belongs to active neighbor. The second vertex serves for finding the other(middle) vertex
   * which define neighbor edge corresponding to active edge.
   * Parameter active_segment is 0 if the neighbor edge of active neighbor has vertices parent1 and middle vertex and equal to 1 if the vertices which define
   * the neighbor edge are middle vertex and parent2.
   * \param[out] edge_info The relative orientation is copied into the struct NeighborEdgeInfo.
  **/
	int direction_neighbor_edge(int parent1, int parent2, int segment);
  
  /// Cleaning of internal structures before a new edge is set as active.
  void reset_neighb_info();
  
  
  // Key for caching transformed function values on elements
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
  
  // Caching transformed geometric for element
  std::map<Key, Geom<double>*, Compare> cache_e;
  
public:
  class ExtendedShapeset
  {
    public:
      class ExtendedShapeFunction
      {
        private:
          NeighborSearch *neibhood;
          PrecalcShapeset *active_pss;
          RefMap *active_rm;
          
          int order;
          bool support_on_neighbor;
          bool reverse_neighbor_side;
          
          void activate(int index, AsmList* central_al, AsmList* neighb_al);
          
          ExtendedShapeFunction(NeighborSearch* neighborhood) : neibhood(neighborhood) {};
          
        public:
          int idx;      ///< shape function index
          int dof;      ///< basis function number
          scalar coef;  ///< coefficient
          
          RefMap* get_activated_refmap() { return active_rm; }
          PrecalcShapeset* get_activated_pss() { return active_pss; }
          int get_fn_order() { return order; }            
          
          Func<double>* get_active_func(int eo) {
            return init_fn(active_pss, active_rm, eo);
          }
          Func<Ord>* get_active_func_ord() {
            int inc = (active_pss->get_num_components() == 2) ? 1 : 0;
            return init_fn_ord(this->order + inc);
          }
          DiscontinuousFunc<double>* make_discontinuous(Func<double>* fu) {
            return new DiscontinuousFunc<double>(fu, support_on_neighbor, reverse_neighbor_side);
          }
          DiscontinuousFunc<Ord>* make_discontinuous(Func<Ord>* fu) {
            return new DiscontinuousFunc<Ord>(fu, support_on_neighbor);
          }
          int get_quad_eo() { return neibhood->get_quad_eo(support_on_neighbor); }
          
        friend class NeighborSearch::ExtendedShapeset; // Only an ExtendedShapeset is allowed to create an ExtendedShapeFunction.
      };
      
      ExtendedShapeFunction* get_extended_shape_fn(int index) { 
        active_shape->activate(index, central_al, neighbor_al); 
        return active_shape;
      }
      
    private:
      ExtendedShapeFunction *active_shape;
      AsmList* central_al;
      AsmList* neighbor_al;
      
      void combine_assembly_lists() {
        assert(central_al != NULL && neighbor_al != NULL);
        
        cnt = central_al->cnt + neighbor_al->cnt;
        dof = new int [cnt];
        memcpy(dof, central_al->dof, sizeof(int)*central_al->cnt);
        memcpy(dof + central_al->cnt, neighbor_al->dof, sizeof(int)*neighbor_al->cnt);
      }
      
      void update(NeighborSearch* neighborhood, Space* space) {
        delete [] this->dof;
        space->get_edge_assembly_list(neighborhood->neighb_el, neighborhood->neighbor_edge, neighbor_al);
        combine_assembly_lists();
      }
      
      ExtendedShapeset(NeighborSearch* neighborhood, AsmList* central_al, Space *space) : central_al(central_al)
      {
        neighbor_al = new AsmList();
        space->get_edge_assembly_list(neighborhood->neighb_el, neighborhood->neighbor_edge, neighbor_al);
        combine_assembly_lists();
        active_shape = new ExtendedShapeFunction(neighborhood);
      }
      
      ~ExtendedShapeset() {  delete [] dof; delete active_shape; delete neighbor_al; }
      
    public:
      int cnt;
      int *dof;
      
    friend class NeighborSearch; // Only a NeighborSearch is allowed to create an ExtendedShapeset.
  };
};

typedef NeighborSearch::ExtendedShapeset::ExtendedShapeFunction* ExtendedShapeFnPtr;


#endif /* NEIGHBOR_H_ */
