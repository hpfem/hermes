#include "hermes2d.h"

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
  quad(&g_quad_2d_std)
{
  central_transformations.reserve(NeighborSearch::max_neighbors * 2);
  for(int i = 0; i < NeighborSearch::max_neighbors * 2; i++)
    central_transformations.push_back(new int[NeighborSearch::max_n_trans]);

  central_n_trans.reserve(NeighborSearch::max_neighbors*2);
  for(int i = 0; i < NeighborSearch::max_neighbors * 2; i++)
    central_n_trans.push_back(0);

  neighbor_transformations.reserve(NeighborSearch::max_neighbors * 2);
  for(int i = 0; i < NeighborSearch::max_neighbors * 2; i++)
    neighbor_transformations.push_back(new int[NeighborSearch::max_n_trans]);

  neighbor_n_trans.reserve(NeighborSearch::max_neighbors*2);
  for(int i = 0; i < NeighborSearch::max_neighbors * 2; i++)
    neighbor_n_trans.push_back(0);

  assert_msg(central_el != NULL && central_el->active == 1,
             "You must pass an active element to the NeighborSearch constructor.");
  neighbors.reserve(NeighborSearch::max_neighbors * 2);
  neighbor_edges.reserve(NeighborSearch::max_neighbors * 2);

  ignore_errors = false;
}

NeighborSearch::NeighborSearch(const NeighborSearch& ns) :
  supported_shapes(NULL),
  mesh(ns.mesh),
  central_el(ns.central_el),
  neighb_el(NULL),
  central_rm(NULL),
  neighb_rm(NULL),
  central_pss(NULL),
  neighb_pss(NULL),
  quad(&g_quad_2d_std),
  ignore_errors(ns.ignore_errors),
  n_neighbors(ns.n_neighbors),
  neighborhood_type(ns.neighborhood_type),
  original_central_el_transform(ns.original_central_el_transform),
  active_edge(ns.active_edge),
  neighbor_edge(ns.neighbor_edge),
  active_segment(ns.active_segment)
{
  central_transformations.reserve(NeighborSearch::max_neighbors * 2);
  central_n_trans.reserve(NeighborSearch::max_neighbors*2);
  neighbor_transformations.reserve(NeighborSearch::max_neighbors * 2);
  neighbor_n_trans.reserve(NeighborSearch::max_neighbors*2);
  neighbors.reserve(NeighborSearch::max_neighbors * 2);
  neighbor_edges.reserve(NeighborSearch::max_neighbors * 2);

  for(int i = 0; i < ns.central_transformations.size(); i++) {
  this->central_transformations.push_back(new int[ns.central_n_trans[i]]);
    for(unsigned int j = 0; j < ns.central_n_trans[i]; j++)
      this->central_transformations[i][j] = ns.central_transformations[i][j];
  }

  for(int i = 0; i < ns.central_n_trans.size(); i++)
    this->central_n_trans.push_back(ns.central_n_trans[i]);

  for(int i = 0; i < ns.neighbor_transformations.size(); i++) {
    this->neighbor_transformations.push_back(new int[ns.neighbor_n_trans[i]]);
    for(unsigned int j = 0; j < ns.neighbor_n_trans[i]; j++)
      this->neighbor_transformations[i][j] = ns.neighbor_transformations[i][j];
  }

  for(int i = 0; i < ns.neighbor_n_trans.size(); i++)
    this->neighbor_n_trans.push_back(ns.neighbor_n_trans[i]);

  assert_msg(central_el != NULL && central_el->active == 1,
             "You must pass an active element to the NeighborSearch constructor.");

  for(int i = 0; i < ns.neighbors.size(); i++)
    this->neighbors.push_back(ns.neighbors[i]);
  for(int i = 0; i < ns.neighbor_edges.size(); i++)
    this->neighbor_edges.push_back(ns.neighbor_edges[i]);


}

NeighborSearch::~NeighborSearch()
{
	neighbor_edges.clear();
	neighbors.clear();
  clear_caches();
  clear_supported_shapes();
  clear_neighbor_pss();
  detach_pss_and_rm();
  for(unsigned int i = 0; i < central_transformations.size(); i++)
    delete [] central_transformations.at(i);
  central_transformations.clear();
  central_n_trans.clear();
  for(unsigned int i = 0; i < neighbor_transformations.size(); i++)
    delete [] neighbor_transformations.at(i);
  neighbor_transformations.clear();
  neighbor_n_trans.clear();
}

void NeighborSearch::reset_neighb_info()
{
  // Reset information about the neighborhood's active state.
  active_segment = -1;
  active_edge = -1;
  neighb_el = NULL;
  neighbor_edge.local_num_of_edge = -1;

  // Clear vectors with neighbor elements and their edge info for the active edge.
  neighbor_edges.clear();
  neighbors.clear();
  n_neighbors = 0;

  // Reset transformations.
  for(int i = 0; i < NeighborSearch::max_neighbors; i++) {
    central_n_trans[i] = 0;
    for(int j = 0; j < max_n_trans; j++)
      central_transformations[i][j] = -1;
  }
  for(int i = 0; i < NeighborSearch::max_neighbors; i++) {
    neighbor_n_trans[i] = 0;
    for(int j = 0; j < max_n_trans; j++)
      neighbor_transformations[i][j] = -1;
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

	//debug_log("central element: %d", central_el->id);
	if (central_el->en[active_edge]->bnd == 0)
	{
		neighb_el = central_el->get_neighbor(active_edge);

		// First case : The neighboring element is of the same size as the central one.
		if (neighb_el != NULL)
		{
			//debug_log("active neighbor el: %d", neighb_el->id);

      // Get local number of the edge used by the neighbor.
			for (unsigned int j = 0; j < neighb_el->nvert; j++)
				if (central_el->en[active_edge] == neighb_el->en[j])
				{
          neighbor_edge.local_num_of_edge = j;
					break;
				}

			NeighborEdgeInfo local_edge_info;
      local_edge_info.local_num_of_edge = neighbor_edge.local_num_of_edge;

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
}

void NeighborSearch::set_active_edge_multimesh(const int& edge)
{
  Hermes::vector<unsigned int> transformations = get_transforms(original_central_el_transform);
  // Inter-element edge.
  if(is_inter_edge(edge, transformations)) {
    set_active_edge(edge);
    update_according_to_sub_idx(transformations);
  }
  // Intra-element edge.
  else {
    neighb_el = central_el;
    
    for(unsigned int i = 0; i < transformations.size(); i++)
      neighbor_transformations[0][i] = transformations[i];
    neighbor_n_trans[0] = transformations.size();

    neighbor_edge.local_num_of_edge = active_edge = edge;
    NeighborEdgeInfo local_edge_info;
    local_edge_info.local_num_of_edge = neighbor_edge.local_num_of_edge;
	  //! The "opposite" view of the same edge has the same orientation.
	  local_edge_info.orientation = 0;
	  neighbor_edges.push_back(local_edge_info);

	  n_neighbors = 1;
	  neighbors.push_back(neighb_el);
	  neighborhood_type = H2D_DG_NO_TRANSF;
  }
  return;
}

Hermes::vector<unsigned int> NeighborSearch::get_transforms(uint64_t sub_idx)
{
  Hermes::vector<unsigned int> transformations_backwards;
  int i = 0;
  while(sub_idx >> 3 > 0)
    transformations_backwards[i++] = sub_idx % 8;
  
  Hermes::vector<unsigned int> transformations;
  for(unsigned int i = 0; i < transformations_backwards.size(); i++)
    transformations[i] = transformations_backwards[transformations_backwards.size() - i];

  return transformations;
}

bool NeighborSearch::is_inter_edge(const int& edge, const Hermes::vector<unsigned int>& transformations)
{
  // No subelements => of course this edge is an inter-element one.
	if(transformations.size() == 0)
		return true;

	// Triangles.
  for(unsigned int i = 0; i < transformations.size(); i++)
	  if(central_el->get_mode() == HERMES_MODE_TRIANGLE) {
      if ((edge == 0 && (transformations[i] == 3 || transformations[i] == 4)) ||
          (edge == 1 && (transformations[i] == 1 || transformations[i] == 4)) ||
          (edge == 2 && (transformations[i] == 2 || transformations[i] == 4)))
        return false;
    }
	  // Quads.
	  else {
      if ((edge == 0 && (transformations[i] == 3 || transformations[i] == 4 || transformations[i] == 6)) ||
          (edge == 1 && (transformations[i] == 1 || transformations[i] == 4 || transformations[i] == 7)) ||
          (edge == 2 && (transformations[i] == 1 || transformations[i] == 2 || transformations[i] == 5)) ||
          (edge == 3 && (transformations[i] == 2 || transformations[i] == 3 || transformations[i] == 0)))
				return false;
    }
	return true;
}

void NeighborSearch::update_according_to_sub_idx(const Hermes::vector<unsigned int>& transformations)
{
  if(neighborhood_type == H2D_DG_NO_TRANSF || neighborhood_type == H2D_DG_GO_UP) {
    for(unsigned int i = 0; i < transformations.size(); i++)
      // Triangles.
      if(central_el->get_mode() == HERMES_MODE_TRIANGLE)
        if ((active_edge == 0 && transformations[i] == 1) ||
            (active_edge == 1 && transformations[i] == 2) ||
            (active_edge == 2 && transformations[i] == 3))
          neighbor_transformations[0][neighbor_n_trans[0]++] = (!neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 3);
				else
          neighbor_transformations[0][neighbor_n_trans[0]++] = (neighbor_edges[0].orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 3);
	    // Quads.
	    else
        if ((active_edge == 0 && (transformations[i] == 1 || transformations[i] == 7)) ||
            (active_edge == 1 && (transformations[i] == 2 || transformations[i] == 5)) ||
            (active_edge == 2 && (transformations[i] == 3 || transformations[i] == 0)) ||
            (active_edge == 3 && (transformations[i] == 4 || transformations[i] == 6)))
          neighbor_transformations[0][neighbor_n_trans[0]++] = (!neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 4);
        else
          neighbor_transformations[0][neighbor_n_trans[0]++] = (neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 4);
  }
  else handle_sub_idx_way_down(transformations);
}

void NeighborSearch::handle_sub_idx_way_down(const Hermes::vector<unsigned int>& transformations)
{
  // We basically identify the neighbors that are not compliant with the current sub-element mapping on the central element.
  for(int level = 0; level < transformations.size(); level++) {
    // In case of bigger (i.e. ~ way up) neighbor, there will be one neighbor left, so this saves time.
    // Commented out for debugging, as the check should pass.
    //if(n_neighbors == 1)
    //  break;
    for(int i = 0; i < n_neighbors; i++) {
      // If the found neighbor is not a neighbor of this subelement.
      if(!compatible_transformations(central_transformations[i][level], transformations[level], active_edge))
        delete_neighbor(i);
      else
        // We want to use the transformations from assembling, because set_active_edge only uses bsplit.
        central_transformations[i][level] = transformations[level];
        // If we are already on a bigger (i.e. ~ way up) neighbor.
        if(central_n_trans[i] == level + 1) {
          for(unsigned int i = level + 1; i < transformations.size(); i++)
            // Triangles.
            if(central_el->get_mode() == HERMES_MODE_TRIANGLE)
              if ((active_edge == 0 && transformations[i] == 1) ||
                  (active_edge == 1 && transformations[i] == 2) ||
                  (active_edge == 2 && transformations[i] == 3))
                neighbor_transformations[0][neighbor_n_trans[0]++] = (!neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 3);
				      else
                neighbor_transformations[0][neighbor_n_trans[0]++] = (neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 3);
	          // Quads.
	          else
              if ((active_edge == 0 && (transformations[i] == 1 || transformations[i] == 7)) ||
                  (active_edge == 1 && (transformations[i] == 2 || transformations[i] == 5)) ||
                  (active_edge == 2 && (transformations[i] == 3 || transformations[i] == 0)) ||
                  (active_edge == 3 && (transformations[i] == 4 || transformations[i] == 6)))
                neighbor_transformations[0][neighbor_n_trans[0]++] = (!neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 4);
              else
                neighbor_transformations[0][neighbor_n_trans[0]++] = (neighbor_edge.orientation ? neighbor_edge.local_num_of_edge : (neighbor_edge.local_num_of_edge + 1) % 4);
        }
      }
  }
}

bool NeighborSearch::compatible_transformations(unsigned int a, unsigned int b, int edge)
{
  if(a == b)
    return true;
  if(edge == 0)
    if ((a == 1 && b == 7) ||
        (a == 2 && b == 0))
      return true;
    else
      return false;
  if(edge == 1)
    if ((a == 2 && b == 5) ||
        (a == 3 && b == 6))
      return true;
    else
      return false;
  if(edge == 2)
    if ((a == 3 && b == 0) ||
        (a == 4 && b == 7))
      return true;
    else
      return false;
  if(edge == 3)
    if ((a == 4 && b == 6) ||
        (a == 1 && b == 5))
      return true;
    else
      return false;
  return false;
}

void NeighborSearch::clear_initial_sub_idx()
{
  if(!neighborhood_type == H2D_DG_GO_DOWN)
    return;
  // Obtain the transformations sequence.
  Hermes::vector<unsigned int> transformations = get_transforms(original_central_el_transform);
  // Test for active element.
  if(transformations.empty())
    return;
  for(unsigned int i = 0; i < n_neighbors; i++) {
    // Find the index where the additional subelement mapping (on top of the initial one from assembling) starts.
    int j = 0;
    // Note that we do not have to test if central_transformations is empty or how long it is, because it has to be
    // longer than transformations (and that is tested).
    // Also the function compatible_transformations() does not have to be used, as now the array central_transformations
    // has been adjusted so that it contains the array transformations.
    while(central_transformations[i][j] == transformations[j])
      if(++j > transformations.size() - 1)
        break;
    // Create a new array of transformations.
    int* shifted_trfs = new int[NeighborSearch::max_n_trans];
    // Move the old one to the new one.
    for(unsigned int k = j; k < central_n_trans[i]; k++)
      shifted_trfs[k -j] = central_transformations[i][k];
    // Point to the new one and delete the old one.
    delete central_transformations[i];
    central_transformations[i] = shifted_trfs;
  }
}

void NeighborSearch::delete_neighbor(unsigned int position)
{
  delete [] central_transformations[position];
  central_transformations.erase (central_transformations.begin()+position);
  central_n_trans.erase (central_n_trans.begin()+position);
  neighbor_edges.erase (neighbor_edges.begin()+position);
  neighbors.erase (neighbors.begin()+position);
  n_neighbors--;
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
        neighbor_edge.local_num_of_edge = -1;
				for(unsigned int j = 0; j < neighb_el->nvert; j++)
					if(neighb_el->en[j] == edge)
					{
            neighbor_edge.local_num_of_edge = j;
						break;
					}
        if(neighbor_edge.local_num_of_edge == -1) error("Neighbor edge wasn't found");

				Node* n = NULL;

				// Add to the array of neighbor_transformations one that transforms central el. to its parent completely
        // adjacent to the single big neighbor.
				assert(n_neighbors == 0);
				neighbor_n_trans[n_neighbors] = n_parents;
        if(n_neighbors > NeighborSearch::max_neighbors)
          NeighborSearch::max_neighbors = n_neighbors;

        /*
				for(int k = 0 ; k < n_parents; k++)
					debug_log("vertices on the way: %d", parents[k]->id);
				debug_log("\n");
				*/

				// Go back through the intermediate inactive parents down to the central element and stack corresponding
        // neighbor_transformations into the array 'neighbor_transformations'.
				for(int j = n_parents - 1; j > 0; j-- )
				{
          n = mesh->peek_vertex_node(par_mid_vertices[j]->id, p1);
          if(n == NULL)
          {
            neighbor_transformations[n_neighbors][n_parents - j - 1] = neighbor_edge.local_num_of_edge;
            p1 = par_mid_vertices[j]->id;
          }
          else
          {
            if(n->id == par_mid_vertices[j-1]->id)
            {
              neighbor_transformations[n_neighbors][n_parents - j - 1] = (neighbor_edge.local_num_of_edge + 1) % neighb_el->nvert;
              p2 = par_mid_vertices[j]->id;
            }
            else
            {
              neighbor_transformations[n_neighbors][n_parents - j - 1] = neighbor_edge.local_num_of_edge;
              p1 = par_mid_vertices[j]->id;
            }
          }
				}

				// Final transformation to the central element itself.
        if (orig_vertex_id[0] == par_mid_vertices[0]->id)
          neighbor_transformations[n_neighbors][n_parents - 1] = neighbor_edge.local_num_of_edge;
				else
          neighbor_transformations[n_neighbors][n_parents - 1] = (neighbor_edge.local_num_of_edge + 1) % neighb_el->nvert;


				NeighborEdgeInfo local_edge_info;
        local_edge_info.local_num_of_edge = neighbor_edge.local_num_of_edge;
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
          neighbor_edge.local_num_of_edge = -1;
          for(unsigned int k = 0; k < neighb_el->nvert; k++)
            if(neighb_el->en[k] == edge)
            {
              neighbor_edge.local_num_of_edge = k;
              break;
            }
          if(neighbor_edge.local_num_of_edge == -1) error("Neighbor edge wasn't found");

          // Construct the transformation path to the current neighbor.
          for(int k = 0; k < n_sons; k++)
            central_transformations[n_neighbors][k] = sons[k];

          central_n_trans[n_neighbors] = n_sons;


          NeighborEdgeInfo local_edge_info;
          local_edge_info.local_num_of_edge = neighbor_edge.local_num_of_edge;
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
    if (neighb_el->vn[neighbor_edge.local_num_of_edge]->id != bounding_vert1)
      return 1; // orientation reversed
  }
  else
  {
    // neighbor edge goes from middle vertex to parent2
    if (neighb_el->vn[neighbor_edge.local_num_of_edge]->id == bounding_vert2)
      return 1; // orientation reversed
  }
  return 0;
}


bool NeighborSearch::set_active_segment(int neighbor, bool with_neighbor_pss)
{
  _F_
  ensure_active_edge(this);
  active_segment = neighbor;
  ensure_active_segment(this);

  neighb_el = neighbors[active_segment];
  neighbor_edge = neighbor_edges[active_segment];

  if(central_pss != NULL)
  {
    ensure_central_pss_rm(this);

    // Reset the transformation of the pss on central element.
    if(central_pss->get_transform() != original_central_el_transform)
      central_pss->set_transform(original_central_el_transform);

    // Push the central element's transformation to the central pss and refmap.
    if (neighborhood_type == H2D_DG_GO_DOWN)
      for(int i = 0; i < central_n_trans[active_segment]; i++)
        central_pss->push_transform(central_transformations[active_segment][i]);

    central_rm->force_transform(central_pss->get_transform(), central_pss->get_ctm());
  }
  else if (central_rm != NULL)
  {
    // Push the central element's transformation to the central refmap.
    if (neighborhood_type == H2D_DG_GO_DOWN)
      for(int i = 0; i < central_n_trans[active_segment]; i++)
        central_rm->push_transform(central_transformations[active_segment][i]);
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
      for(int i = 0; i < neighbor_n_trans[active_segment]; i++)
        neighb_pss->push_transform(neighbor_transformations[active_segment][i]);
      neighb_rm->force_transform(neighb_pss->get_transform(), neighb_pss->get_ctm());
    }
  }
  /*
  if (neighb_el->visited)
    debug_log("%d | %d skipped.\n", central_el->id, neighb_el->id);
  */
  return true;
}

void NeighborSearch::attach_pss_and_rm(PrecalcShapeset *pss, RefMap *rm)
{
  _F_
  central_pss = pss;
  central_rm = rm;

  original_central_el_transform = pss->get_transform();
  assert_msg( original_central_el_transform == rm->get_transform(),
              "Cannot use a RefMap that is transformed differently than the supplied PrecalcShapeset." );
}

void NeighborSearch::attach_rm(RefMap *rm)
{
  central_rm = rm;
  original_central_el_transform = rm->get_transform();
}

void NeighborSearch::detach_pss_and_rm()
{
  _F_
  if (central_pss == NULL)
  {
    if (central_rm != NULL)
      detach_rm();
    return;
  }

  if (neighborhood_type == H2D_DG_GO_DOWN) {
    if (central_pss->get_transform() != original_central_el_transform) {
      central_pss->set_transform(original_central_el_transform);
      central_rm->force_transform(central_pss->get_transform(), central_pss->get_ctm());
    }
  }
}

void NeighborSearch::detach_rm()
{
  if (central_rm != NULL && neighborhood_type == H2D_DG_GO_DOWN && central_rm->get_transform() != original_central_el_transform)
    central_rm->set_transform(original_central_el_transform);
}

NeighborSearch::ExtendedShapeset* NeighborSearch::create_extended_asmlist(Space *space, AsmList* al)
{
  _F_
  if (supported_shapes == NULL)
    supported_shapes = new ExtendedShapeset(this, al, space);
  else
    supported_shapes->update(this, space);

  return supported_shapes;
}

void NeighborSearch::set_quad_order(int order)
{
  _F_
  quad->set_mode(neighbors[active_segment]->get_mode());
  neighb_quad.init(quad, quad->get_edge_points(neighbor_edge.local_num_of_edge, order));
  quad->set_mode(central_el->get_mode());
  central_quad.init(quad, quad->get_edge_points(active_edge, order));
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
  int eo = get_quad_eo();

  // Do not use the caches at all.
  if (ext_cache_e == NULL)
    return new InterfaceGeom<double> (init_geom_surf(central_rm, ep, eo),
                                      neighb_el->marker, neighb_el->id, neighb_el->get_diameter());

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
  int eo = get_quad_eo();

  // Do not use the cache at all.
  if (ext_cache_jwt == NULL)
    return calculate_jwt(eo);

  if (n_neighbors == 1) // go-up or no-transf neighborhood
  {
    // Do the same as if assembling standard (non-DG) surface forms.
    if (ext_cache_jwt[eo] == NULL)
      ext_cache_jwt[eo] = calculate_jwt(eo);
    return ext_cache_jwt[eo];
  }
  else  // go-down neighborhood
  {
    // Also take into account the transformations of the central element.
    Key key(eo, active_segment);
    if (cache_jwt[key] == NULL)
      cache_jwt[key] = calculate_jwt(eo);
    return cache_jwt[key];
  }
}

double* NeighborSearch::calculate_jwt(int edge_order)
{
  int np = get_quad_np();
  double3* pt = get_quad_pt();
  double3* tan = central_rm->get_tangent(active_edge, edge_order);

  double *jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * tan[i][2];

  return jwt;
}

DiscontinuousFunc<Ord>* NeighborSearch::init_ext_fn_ord(MeshFunction* fu)
{
  _F_
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  int central_order = fu->get_edge_fn_order(active_edge) + inc;
  int neighbor_order = fu->get_edge_fn_order(neighbor_edge.local_num_of_edge) + inc;
  return new DiscontinuousFunc<Ord>(init_fn_ord(central_order), init_fn_ord(neighbor_order));
}

DiscontinuousFunc<Ord>* NeighborSearch::init_ext_fn_ord(Solution* fu)
{
  _F_
  int inc = (fu->get_num_components() == 2) ? 1 : 0;
  int central_order = fu->get_edge_fn_order(active_edge) + inc;
  int neighbor_order = fu->get_edge_fn_order(neighbor_edge.local_num_of_edge) + inc;
  return new DiscontinuousFunc<Ord>(init_fn_ord(central_order), init_fn_ord(neighbor_order));
}

DiscontinuousFunc<scalar>* NeighborSearch::init_ext_fn(MeshFunction* fu)
{
  _F_
  Func<scalar>* fn_central = init_fn(fu, get_quad_eo(false));

  uint64_t original_transform = fu->get_transform();

  // Change the active element of the function. Note that this also resets the transformations on the function.
  fu->set_active_element(neighbors[active_segment]);
  
  for(int i = 0; i < neighbor_n_trans[active_segment]; i++)
    fu->push_transform(neighbor_transformations[active_segment][i]);

  Func<scalar>* fn_neighbor = init_fn(fu, get_quad_eo(true));

  // Restore the original function.
  fu->set_active_element(central_el);
  fu->set_transform(original_transform);

  return new DiscontinuousFunc<scalar>(fn_central, fn_neighbor, (neighbor_edge.orientation == 1));

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
  _F_
  neighbor_al = new AsmList();
  space->get_boundary_assembly_list(neighborhood->neighb_el, neighborhood->neighbor_edge.local_num_of_edge, neighbor_al);
  combine_assembly_lists();
}

void NeighborSearch::ExtendedShapeset::combine_assembly_lists()
{
  assert(central_al != NULL && neighbor_al != NULL);
  cnt = central_al->cnt + neighbor_al->cnt;
  dof = new int [cnt];
  memcpy(dof, central_al->dof, sizeof(int)*central_al->cnt);
  memcpy(dof + central_al->cnt, neighbor_al->dof, sizeof(int)*neighbor_al->cnt);
}
