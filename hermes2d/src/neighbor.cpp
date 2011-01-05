#include "neighbor.h"


int NeighborSearch::max_neighbors = 1;

std::map<NeighborSearch::MainKey, NeighborSearch*, NeighborSearch::MainCompare> NeighborSearch::main_cache_m = 
  *new std::map<NeighborSearch::MainKey, NeighborSearch*, NeighborSearch::MainCompare>(); 

std::map<NeighborSearch::MainKey, NeighborSearch*, NeighborSearch::MainCompare> NeighborSearch::main_cache_n = 
  *new std::map<NeighborSearch::MainKey, NeighborSearch*, NeighborSearch::MainCompare>();

void NeighborSearch::empty_main_caches()
{
  std::map<NeighborSearch::MainKey, NeighborSearch*, NeighborSearch::MainCompare>::iterator it;
  for(it = NeighborSearch::main_cache_m.begin(); it != NeighborSearch::main_cache_m.end(); it++) {
    delete (*it).second->central_pss;
    (*it).second->central_pss = NULL;
    delete (*it).second->central_rm;
    (*it).second->central_rm = NULL;
    delete (*it).second;
  }
  // This maybe is not needed.
  NeighborSearch::main_cache_m.clear();

  for(it = NeighborSearch::main_cache_n.begin(); it != NeighborSearch::main_cache_n.end(); it++) {
    delete (*it).second->central_pss;
    (*it).second->central_pss = NULL;
    delete (*it).second->central_rm;
    (*it).second->central_rm = NULL;
    delete (*it).second;
  }
  // This maybe is not needed.
  NeighborSearch::main_cache_n.clear();
}


NeighborSearch::NeighborSearch(Element* el, Mesh* mesh) :
  supported_shapes(NULL),
  mesh(mesh),
  central_el(el),
  neighb_el(NULL),
  central_rm(NULL),
  neighb_rm(NULL),
  central_pss(NULL),
  neighb_pss(NULL),
  ignore_visited_segments(true),
  quad(&g_quad_2d_std)
{
  transformations.reserve(NeighborSearch::max_neighbors * 2);
  for(int i = 0; i < NeighborSearch::max_neighbors * 2; i++)
    transformations.push_back(new int[NeighborSearch::max_n_trans]);

  n_trans.reserve(NeighborSearch::max_neighbors*2);
  for(int i = 0; i < NeighborSearch::max_neighbors * 2; i++)
    n_trans.push_back(0);

  assert_msg(central_el != NULL && central_el->active == 1, 
             "You must pass an active element to the NeighborSearch constructor.");  
  neighbors.reserve(NeighborSearch::max_neighbors * 2);
  neighbor_edges.reserve(NeighborSearch::max_neighbors * 2);

  ignore_errors = false;
} 
 

NeighborSearch::~NeighborSearch()
{
	neighbor_edges.clear();
	neighbors.clear();
  clear_caches();
  clear_supported_shapes();
  clear_neighbor_pss();
  detach_pss();
  for(unsigned int i = 0; i < transformations.size(); i++)
    delete [] transformations.at(i);
  transformations.clear();
  n_trans.clear();
}

void NeighborSearch::reset_neighb_info()
{  
  // Reset information about the neighborhood's active state.
  active_segment = -1;
  active_edge = -1;
  neighb_el = NULL;
  neighbor_edge = -1;  
  
  // Clear vectors with neighbor elements and their edge info for the active edge.
  neighbor_edges.clear();
  neighbors.clear();
  n_neighbors = 0;
  
  // Reset transformations.
  for(int i = 0; i < NeighborSearch::max_neighbors; i++)
  {
    n_trans[i] = 0;
    for(int j = 0; j < max_n_trans; j++)
      transformations[i][j] = -1;
  }
  neighborhood_type = H2D_DG_NOT_INITIALIZED;
  
  // Clear geometry and jac*wt caches.
  clear_caches();
}

void NeighborSearch::set_active_edge(int edge, bool ignore_visited)
{
	// Erase all data from previous edge or element.
	reset_neighb_info();
  
	active_edge = edge;
  ignore_visited_segments = ignore_visited;

	//debug_log("central element: %d", central_el->id);
	if (central_el->en[active_edge]->bnd == 0)
	{
		neighb_el = central_el->get_neighbor(active_edge);

		// First case : The neighboring element is of the same size as the central one.
		if (neighb_el != NULL)
		{
			//debug_log("active neighbor el: %d", neighb_el->id);
      
      // Get local number of the edge used by the neighbor.
			for (int j = 0; j < neighb_el->nvert; j++)
				if (central_el->en[active_edge] == neighb_el->en[j])
				{
					neighbor_edge = j;
					break;
				}
			
			NeighborEdgeInfo local_edge_info;
			local_edge_info.local_num_of_edge = neighbor_edge;
      
      // Query the orientation of the neighbor edge relative to the central el.
      int p1 = central_el->vn[active_edge]->id;
      int p2 = central_el->vn[(active_edge + 1) % central_el->nvert]->id;
			local_edge_info.orientation = neighbor_edge_orientation(p1, p2, 0);
      
			neighbor_edges.push_back(local_edge_info);

			// There is only one neighbor in this case.
			n_neighbors = 1;
			neighbors.push_back(neighb_el);
			
			// No need for transformation, since the neighboring element is of the same size.
			neighborhood_type = H2D_DG_NO_TRANSF;
		} 
		else
		{
      // Peek the vertex in the middle of the active edge (if there is none, vertex will be NULL).
			Node* vertex = mesh->peek_vertex_node(central_el->en[active_edge]->p1,	central_el->en[active_edge]->p2);
      
      // Endpoints of the active edge.
			int orig_vertex_id[2];
			orig_vertex_id[0] = central_el->vn[active_edge]->id;
			orig_vertex_id[1]	= central_el->vn[(active_edge + 1) % central_el->nvert]->id;
      
			if (vertex == NULL)
			{
        neighborhood_type = H2D_DG_GO_UP;
        
				Element* parent = central_el->parent;

				// Array of middle-point vertices of the intermediate parent edges that we climb up to the correct parent element.
				Node** par_mid_vertices = new Node*[max_n_trans]; 
				// Number of visited intermediate parents.
				int n_parents = 0; 

				for (int j = 0; j < max_n_trans; j++)
					par_mid_vertices[j] = NULL;

				find_act_elem_up(parent, orig_vertex_id, par_mid_vertices, n_parents);

				delete[] par_mid_vertices;
			} 
			else
			{
        neighborhood_type = H2D_DG_GO_DOWN;
        
				int sons[max_n_trans]; // array of virtual sons of the central el. visited on the way down to the neighbor
				int n_sons = 0; // number of used transformations

        // Start the search by going down to the first son.
				find_act_elem_down( vertex, orig_vertex_id, sons, n_sons+1);

				//debug_log("number of neighbors on the way down: %d ", n_neighbors);
			}
		}
	}
	else
    if(!ignore_errors)
		  error("The given edge isn't inner");
  
  // Reset the neighbor information to force the user to set the correct active segment.
  neighb_el = NULL;
  neighbor_edge = -1;
}


void NeighborSearch::find_act_elem_up( Element* elem, int* orig_vertex_id, Node** par_mid_vertices, int n_parents)
{
	Node* edge = NULL;
	Node* vertex = NULL;
	
  // IDs of vertices bounding the current intermediate parent edge.
	int p1 = elem->vn[active_edge]->id;
	int p2 = elem->vn[(active_edge + 1) % elem->nvert]->id;

	int id_of_par_orient_1 = p1;
	int id_of_par_orient_2 = p2;

	// Find if p1 and p2 bound a used edge (used by the neighbor element).
	edge = mesh->peek_edge_node(p1, p2);

	// Add the vertex in the middle of the parent edge to the array of intermediate parent vertices. This is for  
  // consequent transformation of functions on neighbor element.
	vertex = mesh->peek_vertex_node(p1, p2);
	if(vertex != NULL)
	{
		if (n_parents == 0)
			par_mid_vertices[n_parents++] = vertex;
		else
			if (n_parents == max_n_trans - 1)
				error("Maximum number of intermediate parents exceeded in NeighborSearch::finding_act_elem_up");
			else
				if(par_mid_vertices[n_parents - 1]->id != vertex->id)
					par_mid_vertices[n_parents++] = vertex;
	}
	
	if ((edge == NULL) || (central_el->en[active_edge]->id == edge->id))
  {
    // We have not yet found the parent of the central element completely adjacent to the neighbor.
		find_act_elem_up(elem->parent, orig_vertex_id, par_mid_vertices, n_parents);
  }
	else
  {
		for (int i = 0; i < 2; i++)
		{
			// Get a pointer to the active neighbor element.
			if ((edge->elem[i] != NULL) && (edge->elem[i]->active == 1))
			{
				neighb_el = edge->elem[i];  //debug_log("way up neighbor: %d", neighb_el->id);
        
				// Get local number of the edge used by the neighbor.
				neighbor_edge = -1;
				for(int j = 0; j < neighb_el->nvert; j++)
					if(neighb_el->en[j] == edge)
					{
						neighbor_edge = j;
						break;
					}
				if(neighbor_edge == -1) error("Neighbor edge wasn't found");
			
				Node* n = NULL;
        
				// Add to the array of transformations one that transforms central el. to its parent completely
        // adjacent to the single big neighbor.
				assert(n_neighbors == 0);
				n_trans[n_neighbors] = n_parents;
        if(n_neighbors > NeighborSearch::max_neighbors)
          NeighborSearch::max_neighbors = n_neighbors;
	
        /*
				for(int k = 0 ; k < n_parents; k++)
					debug_log("vertices on the way: %d", parents[k]->id);
				debug_log("\n");
				*/
        
				// Go back through the intermediate inactive parents down to the central element and stack corresponding
        // transformations into the array 'transformations'. 
				for(int j = n_parents - 1; j > 0; j-- )
				{
          n = mesh->peek_vertex_node(par_mid_vertices[j]->id, p1);
          if(n == NULL)
          {
            transformations[n_neighbors][n_parents - j - 1] = neighbor_edge;
            p1 = par_mid_vertices[j]->id;
          }
          else
          {
            if(n->id == par_mid_vertices[j-1]->id)
            {
              transformations[n_neighbors][n_parents - j - 1] = (neighbor_edge + 1) % neighb_el->nvert;
              p2 = par_mid_vertices[j]->id;
            }
            else
            {
              transformations[n_neighbors][n_parents - j - 1] = neighbor_edge;
              p1 = par_mid_vertices[j]->id;
            }
          }
				}

				// Final transformation to the central element itself.
        if (orig_vertex_id[0] == par_mid_vertices[0]->id)
					transformations[n_neighbors][n_parents - 1] = neighbor_edge;
				else
					transformations[n_neighbors][n_parents - 1] = (neighbor_edge + 1) % neighb_el->nvert;

        
				NeighborEdgeInfo local_edge_info;
				local_edge_info.local_num_of_edge = neighbor_edge;
        // Query the orientation of the neighbor edge relative to the central el.
				local_edge_info.orientation = neighbor_edge_orientation(id_of_par_orient_1, id_of_par_orient_2, 0);
				
        neighbor_edges.push_back(local_edge_info);

				// There is only one neighbor,...
				n_neighbors = 1;

				// ...add it to the vector of neighbors.
				neighbors.push_back(neighb_el);
			}
		}
  }
}


void NeighborSearch::find_act_elem_down( Node* vertex, int* bounding_verts_id, int* sons, int n_sons)
{
	int mid_vert = vertex->id; // ID of vertex in between vertices from par_vertex_id.
	int bnd_verts[2];
	bnd_verts[0] = bounding_verts_id[0];
	bnd_verts[1] = bounding_verts_id[1];

	for (int i = 0; i < 2; i++)
	{
		sons[n_sons-1] = (active_edge + i) % central_el->nvert;

    // Try to get a pointer to the edge between the middle vertex and one of the vertices bounding the previously
    // tested segment.
		Node* edge = mesh->peek_edge_node(mid_vert, bnd_verts[i]);

    if (edge == NULL) // The edge is not used, i.e. there is no active element on either side.
		{
      // Get the middle vertex of this edge and try again on the segments into which this vertex splits the edge.
			Node * n = mesh->peek_vertex_node(mid_vert, bnd_verts[i]);
			if(n == NULL)
				error("wasn't able to find middle vertex");
			else
			{
        // Make sure the next visited segment has the same orientation as the original central element's active edge.
				if(i == 0) 
					bounding_verts_id[1] = mid_vert;
				else 
					bounding_verts_id[0] = mid_vert;
				
				find_act_elem_down( n, bounding_verts_id, sons, n_sons+1);

				bounding_verts_id[0] = bnd_verts[0];
				bounding_verts_id[1] = bnd_verts[1];
			}
		} 
    else  // We have found a used edge, the active neighbor we are looking for is on one of its sides.
		{
			for (int j = 0; j < 2; j++)
			{
        if ((edge->elem[j] != NULL) && (edge->elem[j]->active == 1))
        {
          neighb_el = mesh->get_element(edge->elem[j]->id);  //debug_log("way down neighbor: %d", edge->elem[j]->id);

          // Get local number of the edge used by the neighbor.
          neighbor_edge = -1;
          for(int k = 0; k < neighb_el->nvert; k++)
            if(neighb_el->en[k] == edge)
            {
              neighbor_edge = k;
              break;
            }
          if(neighbor_edge == -1) error("Neighbor edge wasn't found");

          // Construct the transformation path to the current neighbor.
          for(int k = 0; k < n_sons; k++) 
            transformations[n_neighbors][k] = sons[k];

          n_trans[n_neighbors] = n_sons;

          
          NeighborEdgeInfo local_edge_info;
          local_edge_info.local_num_of_edge = neighbor_edge;
          // Query the orientation of the neighbor edge relative to the central el.
          local_edge_info.orientation = neighbor_edge_orientation(bnd_verts[0], bnd_verts[1], i);

          neighbor_edges.push_back(local_edge_info);

          // Append the new neighbor.
          n_neighbors++;
          if(n_neighbors > NeighborSearch::max_neighbors)
            NeighborSearch::max_neighbors = n_neighbors;
          neighbors.push_back(neighb_el);
        }
			}
		}
	}
}

int NeighborSearch::neighbor_edge_orientation(int bounding_vert1, int bounding_vert2, int segment)
{
  if (segment == 0)
  {
    // neighbor edge goes from parent1 to middle vertex
    if (neighb_el->vn[neighbor_edge]->id != bounding_vert1)
      return 1; // orientation reversed
  }
  else
  {
    // neighbor edge goes from middle vertex to parent2
    if (neighb_el->vn[neighbor_edge]->id == bounding_vert2)
      return 1; // orientation reversed
  }
  return 0;
}


bool NeighborSearch::set_active_segment(int neighbor, bool with_neighbor_pss)
{
  ensure_active_edge(this);
  
  active_segment = neighbor; 
  ensure_active_segment(this);
  
  neighb_el = neighbors[active_segment];
  neighbor_edge = neighbor_edges[active_segment].local_num_of_edge;
  
  // TODO: If two spaces share one mesh, can this cause troubles (if the parameter
  // ignore_visited_segments is forgotten??
  if (neighb_el->visited && ignore_visited_segments)
    return false;
 
  if(central_pss != NULL) {
    ensure_central_pss_rm(this);
    
    // Reset the transformation of the pss on central element.
    if(central_pss->get_transform() != original_central_el_transform)
      central_pss->set_transform(original_central_el_transform);
    
    // Push the central element's transformation to the central pss and refmap.
    if (neighborhood_type == H2D_DG_GO_DOWN)
      for(int i = 0; i < n_trans[active_segment]; i++)
        central_pss->push_transform(transformations[active_segment][i]);
    
    central_rm->force_transform(central_pss->get_transform(), central_pss->get_ctm()); 
  }  
  
  if (with_neighbor_pss) // If the extended shapeset is needed...
  { 
    if (neighb_pss == NULL) // Create the neighbor objects for the first time.
    {
      neighb_pss = new PrecalcShapeset(central_pss->get_shapeset());
      neighb_pss->set_quad_2d(quad);
      neighb_rm = new RefMap();
      neighb_rm->set_quad_2d(quad);
    }
    
    // Set active element for the neighbor objects or update it in the case of a go-down neighborhood.
    neighb_pss->set_active_element(neighb_el);
    neighb_rm->set_active_element(neighb_el);
    // Reset the transformation on the neighbor (this has an effect only when the active edge is changed
    // and the previous one defined a go-up neighborhood).
    neighb_pss->reset_transform();
    neighb_rm->reset_transform();
    neighb_pss->set_active_shape(central_pss->get_active_shape());
    
    // Push the neighbor element's transformations in the case of a go-up neighborhood.
    if (neighborhood_type == H2D_DG_GO_UP) {
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

void NeighborSearch::attach_pss(PrecalcShapeset *pss, RefMap *rm)
{
  central_pss = pss;
  central_rm = rm;
  
  original_central_el_transform = pss->get_transform();
  assert_msg( original_central_el_transform == rm->get_transform(), 
              "Cannot use a RefMap that is transformed differently than the supplied PrecalcShapeset." );
}

void NeighborSearch::detach_pss()
{
  if (central_pss == NULL) return;
  
  if (neighborhood_type == H2D_DG_GO_DOWN) {
    if (central_pss->get_transform() != original_central_el_transform) {
      central_pss->set_transform(original_central_el_transform);
      central_rm->force_transform(central_pss->get_transform(), central_pss->get_ctm());
    }
  }
}

int NeighborSearch::create_extended_shapeset(Space *space, AsmList* al)
{
  ensure_central_pss_rm(this);
  ensure_active_segment(this);
    
  if (supported_shapes == NULL) 
    supported_shapes = new ExtendedShapeset(this, al, space);
  else
    supported_shapes->update(this, space);
  
  return supported_shapes->cnt;
}

void NeighborSearch::set_quad_order(int order)
{
  ensure_active_segment(this);
  
  quad->set_mode(central_el->get_mode());
  central_quad.init(quad, quad->get_edge_points(active_edge, order));
  quad->set_mode(neighb_el->get_mode());
  neighb_quad.init(quad, quad->get_edge_points(neighbor_edge, order));
}

double3* NeighborSearch::get_quad_pt(bool on_neighbor) 
{ 
  if (on_neighbor) 
  { 
    ensure_set_quad_order(neighb_quad);
    return neighb_quad.pt;
  } 
  else 
  {
    ensure_set_quad_order(central_quad);
    return central_quad.pt; 
  }
}

int NeighborSearch::get_quad_eo(bool on_neighbor) 
{ 
  if (on_neighbor) 
  { 
    ensure_set_quad_order(neighb_quad);
    return neighb_quad.eo;
  } 
  else 
  {
    ensure_set_quad_order(central_quad);
    return central_quad.eo; 
  }
}

int NeighborSearch::get_quad_np() 
{ 
  ensure_set_quad_order(central_quad);
  return central_quad.np; 
}


Geom<double>* NeighborSearch::init_geometry(Geom<double>** ext_cache_e, SurfPos *ep)
{
  ensure_central_pss_rm(this);
  ensure_active_segment(this);
  
  int eo = get_quad_eo();
  
  if (n_neighbors == 1) // go-up or no-transf neighborhood
  {
    // Do the same as if assembling standard (non-DG) surface forms.
    if (ext_cache_e[eo] == NULL)
      ext_cache_e[eo] = new InterfaceGeom<double> (init_geom_surf(central_rm, ep, eo), 
                                                   neighb_el->marker, neighb_el->id, neighb_el->get_diameter());
    return ext_cache_e[eo];
  } 
  else // go-down neighborhood
  {
    // Also take into account the transformations of the central element.
    Key key(eo, active_segment);
    if (cache_e[key] == NULL)
      cache_e[key] = new InterfaceGeom<double> (init_geom_surf(central_rm, ep, eo), 
                                                neighb_el->marker, neighb_el->id, neighb_el->get_diameter());
    return cache_e[key];
  }    
}

double* NeighborSearch::init_jwt(double** ext_cache_jwt)
{
  ensure_central_pss_rm(this);
  ensure_active_segment(this);
  
  int eo = get_quad_eo();
  
  if (n_neighbors == 1) // go-up or no-transf neighborhood
  {
    // Do the same as if assembling standard (non-DG) surface forms.
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
  else  // go-down neighborhood
  {
    // Also take into account the transformations of the central element.
    Key key(eo, active_segment);
    if (cache_jwt[key] == NULL) {
      int np = get_quad_np();
      double3* pt = get_quad_pt();
      double3* tan = central_rm->get_tangent(active_edge, eo);
      
      cache_jwt[key] = new double[np];
      for(int i = 0; i < np; i++)
        cache_jwt[key][i] = pt[i][2] * tan[i][2];
    }
    return cache_jwt[key];
  }
}

DiscontinuousFunc<Ord>* NeighborSearch::init_ext_fn_ord(MeshFunction* fu)
{
  ensure_active_segment(this);
  Func<Ord>* fo1 = init_fn_ord(fu->get_edge_fn_order(active_edge));
  Func<Ord>* fo2 = init_fn_ord(fu->get_edge_fn_order(active_edge));
  return new DiscontinuousFunc<Ord>(fo1, fo2);
}

DiscontinuousFunc<Ord>* NeighborSearch::init_ext_fn_ord(Solution* fu)
{   
  ensure_active_segment(this);
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  int central_order = fu->get_edge_fn_order(active_edge) + inc;
  int neighbor_order = fu->get_edge_fn_order(neighbor_edge) + inc;
  return new DiscontinuousFunc<Ord>(init_fn_ord(central_order), init_fn_ord(neighbor_order));
}

DiscontinuousFunc<scalar>* NeighborSearch::init_ext_fn(MeshFunction* fu)
{
  ensure_active_segment(this);
  ensure_set_quad_order(central_quad);
  ensure_set_quad_order(neighb_quad);
  
  if(fu->get_active_element() != central_el) {
    // Just in case - the input function should have the active element set to the central element from get_next_state
    // called in the assembling procedure.  
    fu->set_active_element(central_el);
    fu->set_transform(original_central_el_transform);
  }
        
  if (neighborhood_type == H2D_DG_GO_DOWN)
    // Shrink the bigger central element to match the smaller neighbor (shrinks the active edge to the active segment).
    for(int i = 0; i < n_trans[active_segment]; i++)
      fu->push_transform(transformations[active_segment][i]);
  
  Func<scalar>* fn_central = init_fn(fu, fu->get_refmap(), get_quad_eo(false));

  // Change the active element of the function. Note that this also resets the transformations on the function.
  fu->set_active_element(neighb_el);  
  
  if (neighborhood_type == H2D_DG_GO_UP) 
    // Shrink the bigger neighbor to match the smaller central element.
    for(int i = 0; i < n_trans[active_segment]; i++)        
      fu->push_transform(transformations[active_segment][i]);
    
  Func<scalar>* fn_neighbor = init_fn(fu, fu->get_refmap(), get_quad_eo(true));  
  
  // Restore the original function.
  fu->set_active_element(central_el);
  fu->set_transform(original_central_el_transform);

  return new DiscontinuousFunc<scalar>(fn_central, fn_neighbor, (neighbor_edges[active_segment].orientation == 1));
  
  //NOTE: This function is not very efficient, since it sets the active elements and possibly pushes transformations
  // for each mesh function in each cycle of the innermost assembly loop. This is neccessary because only in
  // this innermost cycle (in function DiscreteProblem::eval_form), we know the quadrature order (dependent on 
  // the actual basis and test function), which is needed for creating the Func<scalar> objects via init_fn.
  // The reason for storing the central and neighbor values of any given function in these objects is that otherwise
  // we would have to have one independent copy of the function for each of the neighboring elements. However, it 
  // could unify the way PrecalcShapesets and MeshFunctions are treated in NeighborSearch and maybe these additional 
  // deep memory copying, performed only after setting the active edge part (before the nested loops over basis and 
  // test functions), would be actually more efficient than this. This would require implementing copy for Filters.
}

int NeighborSearch::get_neighb_edge_number(int segment)
{
    ensure_active_edge(this);
	if( (unsigned) segment >= neighbor_edges.size())
		error("given number is bigger than actual number of neighbors ");
	else
		return neighbor_edges[segment].local_num_of_edge;
    return 0;
}

int NeighborSearch::get_neighb_edge_orientation(int segment)
{
    ensure_active_edge(this);
	if( (unsigned) segment >= neighbor_edges.size())
		error("given number is bigger than actual number of neighbors ");
	else
		return neighbor_edges[segment].orientation;
    return 0;
}



/*** _____________________________________________ EXTENDED SHAPESET _____________________________________________ ***/


NeighborSearch::ExtendedShapeset::ExtendedShapeset(NeighborSearch* neighborhood, AsmList* central_al, Space* space) :
  central_al(central_al)
{
  neighbor_al = new AsmList();
  space->get_boundary_assembly_list(neighborhood->neighb_el, neighborhood->neighbor_edge, neighbor_al);
  combine_assembly_lists();
  active_shape = new ExtendedShapeFunction(neighborhood);
}

void NeighborSearch::ExtendedShapeset::combine_assembly_lists() 
{
  assert(central_al != NULL && neighbor_al != NULL);        
  cnt = central_al->cnt + neighbor_al->cnt;
  dof = new int [cnt];
  memcpy(dof, central_al->dof, sizeof(int)*central_al->cnt);
  memcpy(dof + central_al->cnt, neighbor_al->dof, sizeof(int)*neighbor_al->cnt);
}

void NeighborSearch::ExtendedShapeset::ExtendedShapeFunction::activate(int index, AsmList* central_al, AsmList* neighb_al)
{
  ensure_active_segment(neibhood);
  ensure_central_pss_rm(neibhood);
  assert_msg(neibhood->neighb_pss != NULL, "Cannot activate extended shape function."
                                           "PrecalcShapeset for neighbor has not been set."  );
  assert_msg(index >= 0, "Wrong shape function index.");
  
  if (index >= central_al->cnt) 
  {
    // Active shape is nonzero on the neighbor element
    support_on_neighbor = true;
    
    active_pss = neibhood->neighb_pss;
    active_rm = neibhood->neighb_rm;
    
    // AsmList entries for the active shape, taken from neighbor.
    int idx_loc = index - central_al->cnt;
    idx = neighb_al->idx[idx_loc];
    dof = neighb_al->dof[idx_loc];
    coef = neighb_al->coef[idx_loc];
  } 
  else
  {
    // Active shape is nonzero on the central element
    support_on_neighbor = false;
    
    active_pss = neibhood->central_pss;
    active_rm = neibhood->central_rm;
    
    // AsmList entries for the active shape, taken from central.
    idx = central_al->idx[index];
    dof = central_al->dof[index];
    coef = central_al->coef[index];
  }
  
  active_pss->set_active_shape(idx);
  reverse_neighbor_side = (neibhood->neighbor_edges[neibhood->active_segment].orientation == 1);
  order = active_pss->get_edge_fn_order(neibhood->active_edge);
}

DiscontinuousFunc<double>* 
NeighborSearch::ExtendedShapeset::ExtendedShapeFunction::get_fn(  std::map< PrecalcShapeset::Key,
                                                                            Func< double >*, 
                                                                            PrecalcShapeset::Compare >& ext_cache_fn  )
{
  int eo = neibhood->get_quad_eo(support_on_neighbor);
  PrecalcShapeset::Key key( 256 - active_pss->get_active_shape(),
                            eo,
                            active_pss->get_transform(),
                            active_pss->get_shapeset()->get_id()  );
  
  if (ext_cache_fn[key] == NULL)
    ext_cache_fn[key] = init_fn(active_pss, active_rm, eo);
  
  return extend_by_zero( ext_cache_fn[key] );
}
