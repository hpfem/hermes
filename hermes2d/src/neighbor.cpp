#include "neighbor.h"

const char* NeighborSearch::ERR_ACTIVE_SEGMENT_NOT_SET = "Active segment of the active edge has not been set.";

void NeighborSearch::ExtendedShapeset::ExtendedShapeFunction::activate(int index, AsmList* central_al, AsmList* neighb_al)
{
  assert_msg(neibhood->active_segment != -1, ERR_ACTIVE_SEGMENT_NOT_SET);
  assert_msg(index >= 0, "Wrong shape function index.");
  assert_msg(neibhood->neighb_pss != NULL && neibhood->central_pss != NULL, 
             "Cannot activate extended shape function. PrecalcShapesets not set.");

  if (index >= central_al->cnt) 
  {
    active_pss = neibhood->neighb_pss;
    active_rm = neibhood->neighb_rm;
    
    int idx_loc = index - central_al->cnt;
    idx = neighb_al->idx[idx_loc];
    dof = neighb_al->dof[idx_loc];
    coef = neighb_al->coef[idx_loc];
   
    support_on_neighbor = true;
  } 
  else 
  {
    active_pss = neibhood->central_pss;
    active_rm = neibhood->central_rm;
    
    idx = central_al->idx[index];
    dof = central_al->dof[index];
    coef = central_al->coef[index];
    
    support_on_neighbor = false;
  }
  
  active_pss->set_active_shape(idx);
  reverse_neighbor_side = (neibhood->neighbor_edges[neibhood->active_segment].orientation == 1);
  order = active_pss->get_edge_fn_order(neibhood->active_edge);
}

NeighborSearch::NeighborSearch(Element* el, Space* space, bool ignore_if_visited) : 
  central_el(el), central_pss(NULL), central_rm(NULL),
  neighb_el(NULL), neighb_pss(NULL), neighb_rm(NULL),
  space(space), mesh(space->get_mesh()), 
  supported_shapes(NULL), quad(&g_quad_2d_std), ignore_if_visited(ignore_if_visited)
{
  assert_msg(central_el != NULL, "You must set an active element before using RefMap in NeighborSearch.");  
  neighbors.reserve(max_n_trans);
  neighbor_edges.reserve(max_n_trans);
}; 
 

NeighborSearch::~NeighborSearch()
{
	neighbor_edges.clear();
	neighbors.clear();
  clear_geometry_cache();
  clear_supported_shapes();
  clear_neighbor_pss();
  detach_pss();
};

void NeighborSearch::attach_pss(PrecalcShapeset *pss, RefMap *rm)
{
  central_pss = pss;

  /* NOTE: Workaround for a possible segfault caused by pushing transforms to the slave pss fv and a bug in Judy usage.
  central_pss = new PrecalcShapeset(pss->get_shapeset());
  central_pss->set_quad_2d(pss->get_quad_2d());
  central_pss->set_active_element(central_el);
  central_pss->set_transform(pss->get_transform());
  */
  central_rm = rm;
  
  original_central_el_transform = pss->get_transform();
  assert_msg( original_central_el_transform == rm->get_transform(), 
              "Cannot use a RefMap that is transformed differently than the supplied PrecalcShapeset." );
}

void NeighborSearch::detach_pss()
{
  if (central_pss == NULL) return;
  
  if (way_flag == H2D_WAY_DOWN) {
    if (central_pss->get_transform() != original_central_el_transform) {
      central_pss->set_transform(original_central_el_transform);
      central_rm->force_transform(central_pss->get_transform(), central_pss->get_ctm());
    }
  }
  
  // NOTE: Workaround for a possible segfault caused by pushing transforms to the slave pss fv and a bug in Judy usage.
  // delete central_pss; central_pss = NULL;
}

void NeighborSearch::set_active_edge(int edge)
{
	// Erase all data from previous edge or element.
	reset_neighb_info();
  
	active_edge = edge;

	//debug_log("central element: %d", central_el->id);
	if (central_el->en[active_edge]->bnd == 0)
	{
		// Finding active element.
		neighb_el = central_el->get_neighbor(active_edge);

		//First case : The neighboring element is of the same size as the central one.
		if (neighb_el != NULL)
		{
			//debug_log("active neighbor el: %d", neighb_el->id);
			// Finding the correct edge.
			for (int j = 0; j < neighb_el->nvert; j++)
				if (central_el->en[active_edge] == neighb_el->en[j])
				{
					neighbor_edge = j;
					break;
				}

			// Get orientation of neighbor edge.
			int p1 = central_el->vn[active_edge]->id;
			int p2 = central_el->vn[(active_edge + 1) % central_el->nvert]->id;

			NeighborEdgeInfo local_edge_info;
			local_edge_info.local_num_of_edge = neighbor_edge;
			local_edge_info.orientation = direction_neighbor_edge(p1, p2, 0);

			// Add the local_edge_info into the vector.
			neighbor_edges.push_back(local_edge_info);

			// There is only one neighbor in this case.
			n_neighbors = 1;

			// Add neighbor to the neighbors vector.
			neighbors.push_back(neighb_el);
			
			// No need for transformation, since the neighboring element is of the same size.
			way_flag = H2D_NO_TRANSF;
		} 
		else
		{
			Node* vertex = mesh->peek_vertex_node(central_el->en[active_edge]->p1,	central_el->en[active_edge]->p2);
			int orig_vertex_id[2];
			orig_vertex_id[0] = central_el->vn[active_edge]->id;
			orig_vertex_id[1]	= central_el->vn[(active_edge + 1) % central_el->nvert]->id;
			if (vertex == NULL)
			{
				// way up
				Element* parent = central_el->parent;

				// array containing vertices we went through
				Node** road_vertices = new Node*[max_n_trans]; 
				// number of used vertices
				int n_road_vertices = 0; 

				for (int j = 0; j < max_n_trans; j++)
					road_vertices[j] = NULL;

				finding_act_elem_up(parent, orig_vertex_id, road_vertices, n_road_vertices);

				delete[] road_vertices;

				way_flag = H2D_WAY_UP;
			} 
			else
			{
				// way down
				int road[max_n_trans]; // array for temporal transformation
				int n_road = 0; // number of used transformations

				finding_act_elem_down( vertex, orig_vertex_id, road, n_road);

				way_flag = H2D_WAY_DOWN;
        
				//debug_log("number of neighbors on the way down: %d ", n_neighbors);
			}
		}
	}
	else
		error("The given edge isn't inner");
};



//way up
/*! \brief Function for finding "bigger" neighbor.
 *
 * We use recurrence in this way.
 * If the neighbor is "bigger" then this means central element is descendant of some inactive elements. We go threw this parents and
 * stop when against an edge, which has same local number as the original edge, we have active element.
 * Important is that all sons have same orientation as parent, so local number of the edge is same.
 */
void NeighborSearch::finding_act_elem_up( Element* elem, int* orig_vertex_id, Node** road_vertices, int n_road_vertices)
{
	Node* edge = NULL;
	Node* vertex = NULL;
	int p1, p2; // id of parents of the edge

	// Order parents in direction of parent element (needed for transformation of solution).
	p1 = elem->vn[active_edge]->id;
	p2 = elem->vn[(active_edge + 1) % elem->nvert]->id;

	int id_of_par_orient_1 = p1;
	int id_of_par_orient_2 = p2;

	// Find if between parents p1 and p2 is used edge (is used by neighbor element).
	edge = mesh->peek_edge_node(p1, p2);

	// When we are on parent, we take middle vertex on the edge and add it to road_vertices. This is for consequent transformation of solution
	// on neighbor element.
	vertex = mesh->peek_vertex_node(p1, p2);
	if(vertex != NULL)
	{
		if (n_road_vertices == 0)
			road_vertices[n_road_vertices++] = vertex;
		else
			if (n_road_vertices == max_n_trans - 1)
				error("n_road_vertices == max_n_trans - 1 in NeighborSearch::finding_act_elem_up");
			else
				if(road_vertices[n_road_vertices - 1]->id != vertex->id)
					road_vertices[n_road_vertices++] = vertex;
	}
	// If we have not found the parent neighboring an active element.
	if ((edge == NULL) || (central_el->en[active_edge]->id == edge->id))
		finding_act_elem_up(elem->parent, orig_vertex_id, road_vertices, n_road_vertices);
	else
		for (int i = 0; i < 2; i++)
		{
			// this condition test if on one of sides is some element and if the element is active, because it may happen that
			// something is found even thought it's not an active element.
			if ((edge->elem[i] != NULL) && (edge->elem[i]->active == 1))
			{
				// Getting to the correct edge.
        //debug_log("way up neighbor: %d", edge->elem[i]->id);
				neighb_el = edge->elem[i];

				// Finding the correct edge.
				neighbor_edge = -1;
				for(int j = 0; j < neighb_el->nvert; j++)
					if(neighb_el->en[j] == edge)
					{
						neighbor_edge = j;
						break;
					}
				if(neighbor_edge == -1) error("edge wasn't found");
			
				Node* n = NULL;
				// Filling the array containing number of transformation needed for this neighbor.
				// There is always only one neighbor (since the central element is the smaller of the two neighboring ones).
				assert(n_neighbors == 0);
				n_trans[n_neighbors] = n_road_vertices;
	
				//debugging
        /*
				for(int k = 0 ; k < n_road_vertices; k++)
					debug_log("vertices on the way: %d", road_vertices[k]->id);
				debug_log("\n");
				*/
        
				// Go between the inactive parent neighboring an active element and the central 
				// element and save correct transformation along the way down.
				for(int j = n_road_vertices - 1; j > 0; j-- )
				{
						n = mesh->peek_vertex_node(road_vertices[j]->id, p1);
						if(n == NULL)
						{
							transformations[n_neighbors][n_road_vertices - j - 1] = neighbor_edge;
							p1 = road_vertices[j]->id;
						}
						else
						{
								if(n->id == road_vertices[j-1]->id)
								{
									transformations[n_neighbors][n_road_vertices - j - 1] = (neighbor_edge + 1) % neighb_el->nvert;
									p2 = road_vertices[j]->id;
								}
								else
								{
									transformations[n_neighbors][n_road_vertices - j - 1] = neighbor_edge;
									p1 = road_vertices[j]->id;
								}
						}
				}

				// Final transformation on central element.
				int test = 0;
				if (orig_vertex_id[0] == road_vertices[0]->id)
					test = 1;

				if(test == 1){
					transformations[n_neighbors][n_road_vertices - 1] = neighbor_edge;
				}
				else{
					transformations[n_neighbors][n_road_vertices - 1] = (neighbor_edge + 1) % neighb_el->nvert;
				}

				NeighborEdgeInfo local_edge_info;
				local_edge_info.local_num_of_edge = neighbor_edge;
				local_edge_info.orientation = direction_neighbor_edge(id_of_par_orient_1, id_of_par_orient_2, 0);

				// Add the local_edge_info into the vector.
				neighbor_edges.push_back(local_edge_info);

				// There is only one neighbor.
				n_neighbors = 1;

				// Add neighbor to the neighbors vector.
				neighbors.push_back(neighb_el);
			}
		}
};

/*! \brief On active edge we have more neighbors. Gives us information from all neighbors.
 *
 *	Again we use recurrence in this way. In every step we take middle vertex of the edge (starting with active edge). This vertex split the edge
 *	on two parts. On every part (an edge) we test if the new edge is active. If not, the middle vertex is found and the method is called
 *	again with this new vertex on this part.
 */
//way down
void NeighborSearch::finding_act_elem_down( Node* vertex, int* par_vertex_id, int* road, int n_road)
{
	int son = vertex->id;
	int parents[2];
	parents[0] = par_vertex_id[0];
	parents[1] = par_vertex_id[1];

	for (int i = 0; i < 2; i++)
	{
		road[n_road] = (active_edge + i) % central_el->nvert;

		Node* edge = mesh->peek_edge_node(son, parents[i]);

		// Test if edge is inactive. Means there is no active element on each side.
		if (edge == NULL)
		{
			Node * n = mesh->peek_vertex_node(son, parents[i]);
			if(n == NULL)
				error("wasn't able to find middle vertex");
			else
			{
				if(i == 0) 
					par_vertex_id[1] = son;
				else 
					par_vertex_id[0] = son;
				
				finding_act_elem_down( n, par_vertex_id, road, n_road+1);

				par_vertex_id[0] = parents[0];
				par_vertex_id[1] = parents[1];
			}
		} 
		else
		{
			// Test if on one of sides is active element.
			for (int j = 0; j < 2; j++)
			{
				if (edge->elem[j] != NULL)
					if (edge->elem[j]->active == 1)
					{
							neighb_el = mesh->get_element(edge->elem[j]->id);

						  //debug_log("way down neighbor: %d", edge->elem[j]->id);

							// Getting to correct edge.
							neighbor_edge = -1;
							for(int k = 0; k < neighb_el->nvert; k++)
								if(neighb_el->en[k] == edge)
								{
									neighbor_edge = k;
									break;
								}
							if(neighbor_edge == -1) error("edge wasn't found");

							// Filling transformation.
							for(int k = 0; k <= n_road; k++) transformations[n_neighbors][k] = road[k];

							// + 1 is because how to n_road is computed it's one less then number of transformations.
							n_trans[n_neighbors] = n_road + 1;

							NeighborEdgeInfo local_edge_info;
							local_edge_info.local_num_of_edge = neighbor_edge;
							local_edge_info.orientation = direction_neighbor_edge(parents[0], parents[1], i);

							// Add the local_edge_info into the vector.
							neighbor_edges.push_back(local_edge_info);

							// Raise number of neighbors.
							n_neighbors++;

							// Add neighbors to vectors.
							neighbors.push_back(neighb_el);
					}
			}
		}
	}
};

bool NeighborSearch::set_active_segment(int neighbor, bool with_neighbor_pss)
{
  active_segment = neighbor;
  neighb_el = neighbors[active_segment];
  neighbor_edge = neighbor_edges[active_segment].local_num_of_edge;
  
  if (neighb_el->visited && ignore_if_visited)
    return false;
 
  if(central_pss->get_transform() != original_central_el_transform)
    central_pss->set_transform(original_central_el_transform);
  
  if (way_flag == H2D_WAY_DOWN) {
    for(int i = 0; i < n_trans[active_segment]; i++)
      central_pss->push_transform(transformations[active_segment][i]);
    central_rm->force_transform(central_pss->get_transform(), central_pss->get_ctm()); 
  }
  
  if (with_neighbor_pss) {
    if (neighb_pss == NULL) {
      neighb_pss = new PrecalcShapeset(central_pss->get_shapeset());
      neighb_pss->set_quad_2d(quad);
      neighb_pss->set_active_element(neighb_el);
      
      neighb_rm = new RefMap();
      neighb_rm->set_quad_2d(quad);
      neighb_rm->set_active_element(neighb_el);
    }
    
    neighb_pss->reset_transform();
    
    if (way_flag == H2D_WAY_UP) {
      for(int i = 0; i < n_trans[active_segment]; i++)        
        neighb_pss->push_transform(transformations[active_segment][i]);
      neighb_rm->force_transform(neighb_pss->get_transform(), neighb_pss->get_ctm());
    }
  }
  /*
  if (neighb_el->visited)
    debug_log("%d | %d skipped.\n", central_el->id, neighb_el->id);
  */
  return true;
}

int NeighborSearch::extend_attached_shapeset(AsmList* al)
{
  assert_msg(active_segment != -1, ERR_ACTIVE_SEGMENT_NOT_SET);
    
  if (supported_shapes == NULL) 
    supported_shapes = new ExtendedShapeset(this, al);
  else
    supported_shapes->update(this);
  
  return supported_shapes->cnt;
}

// Init geometry and jacobian*weights.
Geom<double>* NeighborSearch::init_geometry(Geom<double>** ext_cache_e, EdgePos *ep)
{
  int eo = get_quad_eo();
  
  if (n_neighbors == 1) 
  {
    if (ext_cache_e[eo] == NULL)
      ext_cache_e[eo] = init_geom_surf(central_rm, ep, eo);
    return ext_cache_e[eo];
  } 
  
  Key key(eo, active_segment);
  if (cache_e[key] == NULL)
    cache_e[key] = init_geom_surf(central_rm, ep, eo);
  return cache_e[key];
}

double* NeighborSearch::init_jwt(double** ext_cache_jwt)
{
  int eo = get_quad_eo();
  
  if (ext_cache_jwt[eo] == NULL) {
    int np = get_quad_np();
    double3* pt = get_quad_pt();
    double3* tan = central_rm->get_tangent(active_edge, eo);
    
    ext_cache_jwt[eo] = new double[np];
    for(int i = 0; i < np; i++)
      ext_cache_jwt[eo][i] = pt[i][2] * tan[i][2];
  }
  
  return ext_cache_jwt[eo];
}


/*! \brief Fill function values / derivatives for the central element and one of its neighbors.
 *
 */
//NOTE: This is not very efficient, since it sets the active elements and possibly pushes transformations
// for each mesh function in each cycle of the innermost assembly loop. This is neccessary because only in
// this innermost cycle (in function LinSystem::eval_form), we know the quadrature order (dependent on 
// the actual basis and test function), which is needed for creating the Func<scalar> objects via init_fn.
// The reason for storing the central and neighbor values of any given function in these objects is that otherwise
// we would have to have one independent copy of the function for each of the neighboring elements. However, it 
// could unify the way PrecalcShapesets and MeshFunctions are treated in NeighborSearch and maybe these additional 
// deep memory copying, performed only after setting the active edge part (before the nested loops over basis and 
// test functions), would be actually more efficient than this. This would require implementing the copy operation for Filters.
//WARNING: Contrary to LinSystem::init_ext_fn, the order here is the real edge order (like 'order'), not the pseudo-order (like 'eo').
DiscontinuousFunc<scalar>* NeighborSearch::init_ext_fn(MeshFunction* fu)
{
  assert_msg(active_segment != -1, ERR_ACTIVE_SEGMENT_NOT_SET);
      
  fu->set_active_element(neighb_el);
  
  if (way_flag == H2D_WAY_UP) 
    // Transform neighbor element on appropriate part.
    for(int i = 0; i < n_trans[active_segment]; i++)        
      fu->push_transform(transformations[active_segment][i]);
      
  Func<scalar>* fn_neighbor = init_fn(fu, fu->get_refmap(), get_quad_eo(true));  

  fu->set_active_element(central_el);
  
  if (way_flag == H2D_WAY_DOWN)
    // Transform central element on appropriate part.
    for(int i = 0; i < n_trans[active_segment]; i++)
      fu->push_transform(transformations[active_segment][i]);
  
  Func<scalar>* fn_central = init_fn(fu, fu->get_refmap(), get_quad_eo(false));

  // Reset transformations.
  fu->reset_transform();

  return new DiscontinuousFunc<scalar>(fn_central, fn_neighbor, (neighbor_edges[active_segment].orientation == 1));    
}

DiscontinuousFunc<Ord>* NeighborSearch::init_ext_fn_ord(MeshFunction* fu)
{
  assert_msg(active_segment != -1, ERR_ACTIVE_SEGMENT_NOT_SET);
  Func<Ord>* fo = init_fn_ord(fu->get_edge_fn_order(active_edge));
  return new DiscontinuousFunc<Ord>(fo, fo);
}

DiscontinuousFunc<Ord>* NeighborSearch::init_ext_fn_ord(Solution* fu)
{   
  assert_msg(active_segment != -1, ERR_ACTIVE_SEGMENT_NOT_SET);
  int inc = (space->get_shapeset()->get_num_components() == 2) ? 1 : 0;
  int central_order = fu->get_edge_fn_order(active_edge, space, central_el) + inc;
  int neighbor_order = fu->get_edge_fn_order(neighbor_edge, space, neighb_el) + inc;
  return new DiscontinuousFunc<Ord>(init_fn_ord(central_order), init_fn_ord(neighbor_order));
}

// Fills the info about orientation into edge_info.
// Here segment can be only 0 or 1, 1 is meaningful only for way down.
// Returns 1 if the orientation is reversed, 0 otherwise.
int NeighborSearch::direction_neighbor_edge(int parent1, int parent2, int segment)
{
  assert_msg(segment != -1, ERR_ACTIVE_SEGMENT_NOT_SET);
  
	if (segment == 0)
	{
		if (neighb_el->vn[neighbor_edge]->id != parent1)
			return 1;
	}
	else
		if (neighb_el->vn[neighbor_edge]->id == parent2)
			return 1;
	return 0;
}

void NeighborSearch::reset_neighb_info()
{
  // Reset information about active state.
  active_segment = -1;
	active_edge = -1;
	neighb_el = NULL;
	neighbor_edge = -1;

  // Clear vectors with neighbor elements and their edge info for the active edge.
  neighbor_edges.clear();
  neighbors.clear();
  n_neighbors = 0;
  
  // Reset transformations.
	for(int i = 0; i < max_n_trans; i++)
	{
		n_trans[i] = 0;
		for(int j = 0; j < max_n_trans; j++)
			transformations[i][j] = -1;
	}
  way_flag = -1;
}


int NeighborSearch::get_neighb_edge_number(int active_segment)
{
	if(active_segment >= neighbor_edges.size())
		error("given number is bigger than actual number of neighbors ");
	else
		return neighbor_edges[active_segment].local_num_of_edge;
};

int NeighborSearch::get_neighb_edge_orientation(int active_segment)
{
	if(active_segment >= neighbor_edges.size())
		error("given number is bigger than actual number of neighbors ");
	else
		return neighbor_edges[active_segment].orientation;
};
