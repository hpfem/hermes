// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _MESH_H_
#define _MESH_H_

class Space;
class MeshLoader;


#include "h3d_common.h"

// refinement type
#define H3D_REFT_HEX_NONE							0x0000
#define H3D_REFT_HEX_X								0x0001
#define H3D_REFT_HEX_Y								0x0002
#define H3D_REFT_HEX_Z								0x0003
#define H3D_H3D_REFT_HEX_XY							0x0004
#define H3D_H3D_REFT_HEX_XZ							0x0005
#define H3D_H3D_REFT_HEX_YZ							0x0006
#define H3D_H3D_H3D_REFT_HEX_XYZ						0x0007

// refinements on facets
#define H3D_REFT_FACE_NONE							0x0000
#define H3D_REFT_QUAD_HORZ							0x0001
#define H3D_REFT_QUAD_VERT							0x0002
#define H3D_REFT_QUAD_BOTH							0x0003

// splits
#define H3D_SPLIT_NONE								0x0000
#define H3D_SPLIT_HEX_X								0x0001
#define H3D_SPLIT_HEX_Y								0x0002
#define H3D_SPLIT_HEX_Z								0x0004
#define H3D_H3D_SPLIT_HEX_XY							H3D_SPLIT_HEX_X | H3D_SPLIT_HEX_Y
#define H3D_H3D_SPLIT_HEX_XZ							H3D_SPLIT_HEX_X | H3D_SPLIT_HEX_Z
#define H3D_H3D_SPLIT_HEX_YZ							H3D_SPLIT_HEX_Y | H3D_SPLIT_HEX_Z
#define H3D_H3D_H3D_SPLIT_HEX_XYZ						H3D_SPLIT_HEX_X | H3D_SPLIT_HEX_Y | H3D_SPLIT_HEX_Z

/// Represents a vertex in 3D
///
///
class HERMES_API Vertex {
public:
	static const int NUM_COORDS = 3;

	Vertex();
	Vertex(double _x, double _y, double _z);
	Vertex(const Vertex &o);			// copy-constructor
	virtual ~Vertex();

	virtual Vertex *copy();

	// for debugging
	void dump();

	// data
	double x, y, z;						// x-, y-, z-coordinates
};

/// Represents an edge in 3D
///
///
class HERMES_API Edge {
public:
	// geometry
	static const int NUM_VERTICES = 2;

	unsigned bnd:1;							// 1 - edge on the boundary, 0 - edge not on the boundary
	unsigned ref:31;							// number of active elements referring to this edge

	// for debugging
	void dump();

  Edge(const Edge &o) { this->bnd = o.bnd; this->ref = o.ref; }

  Edge & operator=(const Edge &o) 
  { 
    this->bnd = o.bnd;
    this->ref = o.ref; 
    return (*this);
  }

  struct Key
  {
    unsigned int * vtcs;
    unsigned int size;
    Key()
    {
      vtcs = NULL;
      size = 0;
    }
   Key(unsigned int vtcs_ [], unsigned int size_)
    {
      this->size = size_;
      if(size > 0)
        this->vtcs = new unsigned int [size];
      for(unsigned int i = 0; i < size; i++) {
        unsigned int temp_place = i;
        for(unsigned int j = i + 1; j < size; j++)
          if(vtcs_[j] < vtcs_[temp_place])
            temp_place = j;
        this->vtcs[i] = vtcs_[temp_place];
        vtcs_[temp_place] = vtcs_[i];
      }
    };
    Key(const Key &b)
    {
      size = b.size;
      if(size > 0)
        this->vtcs = new unsigned int [size];
      for(unsigned int i = 0; i < size; i++)
        vtcs[i] = b.vtcs[i];
    };
    Key & operator =(const Key &b)
    {
      if(size > 0)
        delete [] vtcs;
      size = b.size;
      if(size > 0)
        this->vtcs = new unsigned int [size];
      for(unsigned int i = 0; i < size; i++)
        vtcs[i] = b.vtcs[i];
      return *this;
    };
    ~Key()
    {
      if(size > 0)
        delete [] vtcs;
    };
    bool operator <(const Key & other) const
    {
      if(this->size < other.size)
        return true;
      else if(this->size > other.size)
        return false;
      else
        for(unsigned int i = 0; i < this->size; i++)
          if(this->vtcs[i] < other.vtcs[i])
            return true;
          else if(this->vtcs[i] > other.vtcs[i])
            return false;

      return false;
    };
    bool operator ==(const Key & other) const
    {
      if(this->size < other.size)
        return false;
      else if(this->size > other.size)
        return false;
      else
        for(unsigned int i = 0; i < this->size; i++)
          if(this->vtcs[i] < other.vtcs[i])
            return false;
          else if(this->vtcs[i] > other.vtcs[i])
            return false;

      return true;
    };

    bool operator !=(const Key & other) const
    {
      return (!((*this)==other));
    };
  };
  static Key invalid_key;
  Edge();
	
  bool is_active() { return ref > 0; }
};


/// Represents a triangle in 3D
///
///
class HERMES_API Tri {
public:
	static const int NUM_VERTICES = 3;
	static const int NUM_EDGES = 3;
};


/// Represents a quadrilateral in 3D
///
///
class HERMES_API Quad {
public:
	static const int NUM_VERTICES = 4;
	static const int NUM_EDGES = 4;
};

#ifndef INVALID_IDX
#define INVALID_IDX					((unsigned int) -1)
#endif


/// Stores information about neighboring elements
///
///
class HERMES_API Facet {
public:
	static const int MAX_SONS = 4;		/// Maximum number of sons on a facet

	enum Type {
		INNER = 0,		// inner facet (left is element, right is element)
		OUTER			// outer facet (left is element, right is boundary)
	};

	Facet(ElementMode2D mode);
	Facet(const Facet &o);
	virtual ~Facet();

	virtual Facet *copy();
	virtual Facet *copy_base();

	void set_left_info(unsigned int elem_id, int face_num = -1) {
		this->left = elem_id;
		this->left_face_num = face_num;
		this->lactive = elem_id != INVALID_IDX;
	}

	void set_right_info(unsigned int elem_id, int face_num = -1) {
		this->right = elem_id;
		this->right_face_num = face_num;
		this->ractive = elem_id != INVALID_IDX;
	}

	/// @return TRUE if we have a constraint on the face
	/// @param idx[in] - element ID
	/// @param iface[in] - local face number on the element 'idx'
	bool ced(unsigned int idx, int iface);

	// for debugging
	virtual void dump();

public:
	Type type;					/// type of the facet
	ElementMode2D mode;				/// mode of the facet (TRI, QUAD)
	int left;				/// ID of an element on the "left"
	int right;				/// ID of an element or a boundary on the "right" (depending on the type of the facet)

	signed left_face_num:4;		/// local facet number with respect to left element facet numbering (-1 = not set).
	signed right_face_num:4;	/// local facet number with respect to right element facet numbering (-1 = not set). Only for INNER facets.
	unsigned lactive:1;			/// information for the left is active; 1 - active; 0 - inactive
	unsigned ractive:1;			/// information for the right is active; 1 - active; 0 - inactive
	unsigned ref_mask:2;		/// how is the facet divided (0 - not divived, 1 - horz, 2 - vert, 3 - both)

  struct Key
  {
    unsigned int * vtcs;
    unsigned int size;
    Key()
    {
      vtcs = NULL;
      size = 0;
    }
    Key(unsigned int vtcs_ [], unsigned int size_)
    {
      this->size = size_;
      if(size > 0)
        this->vtcs = new unsigned int [size];
      for(unsigned int i = 0; i < size; i++) {
        unsigned int temp_place = i;
        for(unsigned int j = i + 1; j < size; j++)
          if(vtcs_[j] < vtcs_[temp_place])
            temp_place = j;
        this->vtcs[i] = vtcs_[temp_place];
        vtcs_[temp_place] = vtcs_[i];
      }
    };
    ~Key()
    {
      if(size > 0)
        delete [] vtcs;
    };
    Key(const Key &b)
    {
      size = b.size;
      if(size > 0)
        this->vtcs = new unsigned int [size];
      for(unsigned int i = 0; i < size; i++)
        vtcs[i] = b.vtcs[i];
    };
    Key & operator =(const Key &b)
    {
      if(size > 0)
        delete [] vtcs;
      size = b.size;
      if(size > 0)
        this->vtcs = new unsigned int [size];
      for(unsigned int i = 0; i < size; i++)
        vtcs[i] = b.vtcs[i];
      return *this;
    };
    bool operator <(const Key & other) const
    {
      if(this->size < other.size)
        return true;
      else if(this->size > other.size)
        return false;
      else
        for(unsigned int i = 0; i < this->size; i++)
          if(this->vtcs[i] < other.vtcs[i])
            return true;
          else if(this->vtcs[i] > other.vtcs[i])
            return false;

      return false;
    };
    bool operator ==(const Key & other) const
    {
      if(this->size < other.size)
        return false;
      else if(this->size > other.size)
        return false;
      else
        for(unsigned int i = 0; i < this->size; i++)
          if(this->vtcs[i] < other.vtcs[i])
            return false;
          else if(this->vtcs[i] > other.vtcs[i])
            return false;
      return true;
    };

    bool operator !=(const Key & other) const
    {
      return (!(*this == other));
    };

  };
  static Key invalid_key;
  
	Key parent;		/// Key of the parent facet
	Key sons[MAX_SONS];		/// Key of child facets (interpretation depend on ref_mask)
};


/// Base class for elements (abstract).
///
///
class Element {
public:
	Element();
	Element(const Element &o);
	virtual ~Element();

	virtual ElementMode3D get_mode() const = 0;
	virtual int get_num_vertices() const = 0;
	virtual int get_num_edges() const = 0;
	virtual int get_num_faces() const = 0;
	virtual int get_num_surf() const = 0;

	virtual unsigned int get_vertex(int vertex_num) const = 0;
	virtual void get_vertices(unsigned int *vtcs) const = 0;

	// these 2 should not be overloaded
	virtual int get_edge_vertices(int edge_num, unsigned int *vtcs) const = 0;
	virtual const int *get_edge_vertices(int edge_num) const = 0;

	virtual int get_edge_orientation(int edge_num) const = 0;

	virtual ElementMode2D get_face_mode(int face_num) const = 0;
	virtual int get_num_face_vertices(int face_num) const = 0;
	virtual int get_num_face_edges(int face_num) const = 0;

	// these 2 should not be overloaded
	virtual int get_face_vertices(int face_num, unsigned int *vtcs) const = 0;
	virtual const int *get_face_vertices(int face_num) const = 0;

	virtual const int *get_face_edges(int face_num) const = 0;
	virtual int get_face_orientation(int face_num) const = 0;

	virtual Element *copy() = 0;
	virtual Element *copy_base() = 0;

	// for debugging
	virtual void dump();

	//
	virtual void ref_all_nodes() = 0;
	virtual void unref_all_nodes() = 0;

	//
	virtual unsigned int get_son(int son_idx) { return INVALID_IDX; }
	virtual int get_num_sons() { return -1; }

public:
	unsigned int id;							// id of an element
	int marker;

	unsigned active:1;					// 0 = not active (refined, has sons); 1 = active (no sons)
	unsigned used:1;					// 1 is used, otherwise unused

	int reft;							// refinement

protected:
	int iro_cache;						// inverse refmap order - cached

	friend class RefMap;
};


/// Represents hexahedron in 3D
///
///
class HERMES_API Hex : public Element {
public:
	// geometry
	static const int NUM_VERTICES = 8;
	static const int NUM_FACES = 6;
	static const int NUM_EDGES = 12;

	static const int NUM_SONS = 8;

	Hex();
	Hex(unsigned int v[]);
	Hex(unsigned int v1, unsigned int v2, unsigned int v3, unsigned int v4, unsigned int v5, unsigned int v6, unsigned int v7, unsigned int v8);
	Hex(const Hex &o);
	virtual ~Hex();

	//
	virtual ElementMode3D get_mode() const { return HERMES_MODE_HEX; }
	virtual int get_num_vertices() const { return NUM_VERTICES; }
	virtual int get_num_edges() const { return NUM_EDGES; }
	virtual int get_num_faces() const { return NUM_FACES; }
	virtual int get_num_surf() const { return NUM_FACES; }

	virtual unsigned int get_vertex(int vertex_num) const { return vtcs[vertex_num]; }
	virtual void get_vertices(unsigned int *vtcs) const { memcpy(vtcs, this->vtcs, sizeof(this->vtcs)); }

	// FIXME
	virtual int get_edge_vertices(int edge_num, unsigned int *vtcs) const;
	virtual const int *get_edge_vertices(int edge_num) const;

	virtual int get_edge_orientation(int edge_num) const;

	virtual ElementMode2D get_face_mode(int face_num) const { return HERMES_MODE_QUAD; }
	virtual int get_num_face_vertices(int face_num) const { return Quad::NUM_VERTICES; }
	virtual int get_num_face_edges(int face_num) const { return Quad::NUM_EDGES; }

	// FIXME
	virtual int get_face_vertices(int face_num, unsigned int *vtcs) const;
	virtual const int *get_face_vertices(int face_num) const;

	virtual const int *get_face_edges(int face_num) const;
	virtual int get_face_orientation(int face_num) const;

	virtual Element *copy();
	virtual Element *copy_base();

	virtual void ref_all_nodes();
	virtual void unref_all_nodes();

	// adaptivity
	virtual unsigned int get_son(int son_idx) { assert(son_idx >= 0 && son_idx < NUM_SONS); return sons[son_idx]; }
	virtual int get_num_sons() { return NUM_SONS; }

	// for debugging
	virtual void dump();

protected:
	unsigned int vtcs[NUM_VERTICES];					// array of vertex indices that build up the hexahedron
	int sons[NUM_SONS];						// indices of son elements

	friend class Mesh;
};


/// Represents tetrahedron in 3D
///
///
class HERMES_API Tetra : public Element {
public:
	static const int NUM_VERTICES = 4;
	static const int NUM_FACES = 4;
	static const int NUM_EDGES = 6;

	Tetra();
	Tetra(unsigned int v[]);
	Tetra(unsigned int v1, unsigned int v2, unsigned int v3, unsigned int v4);
	Tetra(const Tetra &o);
	virtual ~Tetra();

	//
	virtual ElementMode3D get_mode() const { return HERMES_MODE_TET; }
	virtual int get_num_vertices() const { return NUM_VERTICES; }
	virtual int get_num_edges() const { return NUM_EDGES; }
	virtual int get_num_faces() const { return NUM_FACES; }
	virtual int get_num_surf() const { return NUM_FACES; }

	virtual unsigned int get_vertex(int vertex_num) const { return vtcs[vertex_num]; }
	virtual void get_vertices(unsigned int *vtcs) const { memcpy(vtcs, this->vtcs, sizeof(this->vtcs)); }

	// FIXME
	virtual int get_edge_vertices(int edge_num, unsigned int *vtcs) const;
	virtual const int *get_edge_vertices(int edge_num) const;

	virtual int get_edge_orientation(int edge_num) const;

	virtual ElementMode2D get_face_mode(int face_num) const { return HERMES_MODE_TRIANGLE; }
	virtual int get_num_face_vertices(int face_num) const  { return Tri::NUM_VERTICES; }
	virtual int get_num_face_edges(int face_num) const { return Tri::NUM_EDGES; }
	// FIXME
	virtual int get_face_vertices(int face_num, unsigned int *vtcs) const;
	virtual const int *get_face_vertices(int face_num) const;

	virtual const int *get_face_edges(int face_num) const;
	virtual int get_face_orientation(int face_num) const;

	virtual Element *copy();
	virtual Element *copy_base();

	virtual void ref_all_nodes();
	virtual void unref_all_nodes();

	// for debugging
	virtual void dump();

protected:
	unsigned int vtcs[NUM_VERTICES];					// array of vertex indices that build up the tetrahedron
};


/// Represents prism in 3D
///
/// currrently not used
class HERMES_API Prism : public Element {
public:
	// geometry
	static const int NUM_VERTICES = 6;
	static const int NUM_FACES = 5;
	static const int NUM_EDGES = 9;
	// adaptivity

	Prism();
	Prism(unsigned int v[]);
	Prism(unsigned int v1, unsigned int v2, unsigned int v3, unsigned int v4, unsigned int v5, unsigned int v6);
	Prism(const Prism &o);
	virtual ~Prism();

	//
	virtual ElementMode3D get_mode() const { return HERMES_MODE_PRISM; }
	virtual int get_num_vertices() const { return NUM_VERTICES; }
	virtual int get_num_edges() const { return NUM_EDGES; }
	virtual int get_num_faces() const { return NUM_FACES; }
	virtual int get_num_surf() const { return NUM_FACES; }

	virtual unsigned int get_vertex(int vertex_num) const { return vtcs[vertex_num]; }
	virtual void get_vertices(unsigned int *vtcs) const { memcpy(vtcs, this->vtcs, sizeof(this->vtcs)); }

	// FIXME
	virtual int get_edge_vertices(int edge_num, unsigned int *vtcs) const;
	virtual const int *get_edge_vertices(int edge_num) const;

	virtual int get_edge_orientation(int edge_num) const;

	virtual ElementMode2D get_face_mode(int face_num) const;
	virtual int get_num_face_vertices(int face_num) const;
	virtual int get_num_face_edges(int face_num) const;
	// FIXME
	virtual int get_face_vertices(int face_num, unsigned int *vtcs) const;
	virtual const int *get_face_vertices(int face_num) const;

	virtual const int *get_face_edges(int face_num) const;
	virtual int get_face_orientation(int face_num) const;

	virtual Element *copy();
	virtual Element *copy_base();

	virtual void ref_all_nodes();
	virtual void unref_all_nodes();

	// for debugging
	virtual void dump();

protected:
	unsigned int vtcs[NUM_VERTICES];					// array of vertex indices that build up the hexahedron
};


/// Base class for boundaries of all types
///
///
class HERMES_API Boundary {
public:
	Boundary(int marker);
	Boundary(const Boundary &o);
	virtual ~Boundary();

	virtual Boundary *copy() = 0;

	// for debugging
	virtual void dump();

	unsigned int id;										// boundary id
	int marker;
};


/// Triangular boundary
///
///
class HERMES_API BoundaryTri : public Boundary {
public:
	static const int NUM_VERTICES = 3;
	static const int NUM_EDGES = 3;

	BoundaryTri(int marker);
	BoundaryTri(const BoundaryTri &o);
	virtual ~BoundaryTri();

	virtual Boundary *copy();
};


/// Quadrilateral boundary
///
class HERMES_API BoundaryQuad : public Boundary {
public:
	static const int NUM_VERTICES = 4;
	static const int NUM_EDGES = 4;

	BoundaryQuad(int marker);
	BoundaryQuad(const BoundaryQuad &o);
	virtual ~BoundaryQuad();

	virtual Boundary *copy();
};


/// Represents the geometry of a mesh
///
///
class HERMES_API Mesh {
//	Mesh(const Mesh &o);
public:
  Mesh();
  virtual ~Mesh();
  /// Frees all data associated with the mesh.
  void free();

  /// Creates a copy of another mesh.
  void copy(const Mesh &mesh);
  /// Copies the coarsest elements of another mesh.
  void copy_base(const Mesh &mesh);

  /// Returns the total number of elements stored.
  unsigned int get_num_elements() const { return elements.size(); }
  /// Returns the number of coarse mesh elements.
  unsigned int get_num_base_elements() const { return nbase; }
  /// Returns the current number of active elements in the mesh.
  unsigned int get_num_active_elements() const { return nactive; }
  /// Returns the maximum node id number plus one.
  unsigned int get_max_element_id() const 
  { 
    unsigned int temp_max = 0;
    for(std::map<unsigned int, Element *>::const_iterator it = elements.begin(); it != elements.end(); it++)
      if(it->first > temp_max)
        temp_max = it->first;
      return temp_max;
  }

	/// Checks wether it is possible to refine an element.
	/// @return true if it posible to apply the refinement, otherwise false
	/// @param[in] eid Element id number.
	/// @param[in] refinement The refinement that is going to be applied
	bool can_refine_element(unsigned int eid, int reft) const;

	/// Refines an element.
	/// @param[in] id Element id number.
	/// @param[in] refinement Refinement to apply
	bool refine_element(unsigned int id, int refinement);

	/// Refines all elements.
	/// @param refinement [in] Same meaning as in refine_element().
	void refine_all_elements(int refinement);

	/// Keeps refining all elements of the mesh until the given criterion
	/// is met for all elements of the mesh. The criterion function
	/// receives a pointer to an element to be considered.
	/// It must return -1 if the element is not to be refined, 0 if it
	/// should be refined uniformly, 1 if it is a quad and should be split
	/// horizontally or 2 if it is a quad and should be split vertically.
	/// Exactly 'depth' levels of refinements are performed.
	void refine_by_criterion(int (*criterion)(Element* e), int depth);

	/// Performs repeated refinements of elements touching a part of the
	/// boundary marked by 'marker'.
	void refine_towards_boundary(int marker, int depth);

	/// Recursively removes all son elements of the given element and makes it active.
	void unrefine_element(unsigned int id);

	/// Unrefines all elements with immediate active sons. In effect, this
	/// shaves off one layer of refinements from the mesh. If done immediately
	/// after refine_all_elements(), this function reverts the mesh to its
	/// original state. However, it is not exactly an inverse to
	/// refine_all_elements().
	void unrefine_all_elements();

	/// Regularize mesh (only 1-irregularity rule implemented)
	void regularize();

	Facet::Key get_facet_id(Element *e, int face_num) const;

	/// Get ID of the edge
	/// @param e Element
	/// @param edge_num Local number of an edge on element e
	/// @return ID of the edge on the element
  Edge::Key get_edge_id(Element *e, int edge_num) const;
	/// Get ID of the edge between a and b
	/// @param[in] a Index (ID) of the first vertex
	/// @param[in] b Index (ID) of the second vertex
	/// @return ID of the edge, INVALID_IDX is edge does not exist
	Edge::Key get_edge_id(unsigned int a, unsigned int b) const;
	/// Get state of an edge
	/// @param eidx ID of an edge
	/// @return true is edge is active, otherwise false

	void dump();

	/// Gets an index of a midpoint between a and b.
	/// @param[in] a index of the first vertex
	/// @param[in] b index of the second vertex
	/// @return index of the midpoint
	unsigned int peek_midpoint(unsigned int a, unsigned int b) const;

	/// Adds an vertex
	unsigned int add_vertex(double x, double y, double z);
	/// Adds an element
	/// @param[in] e An element to add
	Tetra *add_tetra(unsigned int vtcs[]);
	Hex *add_hex(unsigned int vtcs[]);
	Prism *add_prism(unsigned int vtcs[]);

	/// Get the facing facet
	/// @param[in] fid ID of the facet
	/// @param[in] elem_id ID of the element
	Facet::Key get_facing_facet(Facet::Key fid, unsigned int elem_id);

	Boundary *add_tri_boundary(unsigned int vtcs[], int marker);
	Boundary *add_quad_boundary(unsigned int vtcs[], int marker);

	void ugh();

  /// Create faces (to be used after the mesh is read from file).
  void create_faces();

  // data
  std::map<unsigned int, Vertex *>   vertices;
	std::map<Edge::Key, Edge *>      edges;
	std::map<unsigned int, Element *>  elements;
	std::map<unsigned int, Boundary *> boundaries;
  std::map<Facet::Key, Facet *> facets;

protected:

	
  unsigned int nbase;							/// number of base elements
	unsigned int nactive;						/// number of active elements

	Tetra *create_tetra(unsigned int vtcs[]);
	Hex *create_hex(unsigned int vtcs[]);
	Prism *create_prism(unsigned int vtcs[]);

	bool can_refine_hex(Hex *elem, int refinement) const;

	/// Apply a refinement to a hex
	/// @param[in] elem An element to refine
	/// @param[in] refinement Type of refinement
	bool refine_hex(Hex *elem, int refinement);
	bool refine_hex_2(Hex *parent, int refinement);
	bool refine_hex_4(Hex *parent, int refinement);
	bool refine_hex_8(Hex *parent, int refinement);

	/// Refine quad facet
	///
	/// @param[in] parent Parent element
	/// @param[in] iface Local number of a face whose facet will be refined
	/// @param[in] face_refinement How to refine the face (see REFT_QUAD_XXX for possible values)
	/// @param[in] eid ID of son element
	bool refine_quad_facet(Hex *parent, int iface, unsigned int face_refinement, unsigned int eid);
	bool refine_quad_facet(Hex *parent, int iface, unsigned int face_refinement, unsigned int eid0, 
                               unsigned int eid1);
	bool refine_quad_facet(Hex *parent, int iface, unsigned int face_refinement, unsigned int eid0, 
                               unsigned int eid1, unsigned int eid2, unsigned int eid3);

	/// @return true if we can refine the face, otherwise false
	/// @param[in] facet - facet to refine
	/// @param[in] reft - the refinement that is going to be applied
	bool is_compatible_quad_refinement(Facet *facet, int reft) const;

	/// Add a quadrilateral-shaped facet
	///
	/// @param[in] left_elem ID of the element on the right
	/// @param[in] left_iface Local face number on the element on the left
	/// @param[in] right_elem ID of the element on the right
	/// @param[in] right_iface Local face number on the element on the right
	/// @return Pointer to the newly created facet
	Facet *add_quad_facet(Facet::Type type, unsigned int left_elem, int left_iface, unsigned int right_elem, int right_iface);

  struct MidPointKey {
    unsigned int a;
    unsigned int b;
    MidPointKey() {
      this->a = 0;
      this->b = 0;
    };
    MidPointKey(unsigned int a_, unsigned int b_) {
      if(a_ < b_) {
        this->a = a_;
        this->b = b_;
      }
      else {
        this->a = b_;
        this->b = a_;
      }
    };
    bool operator<(const MidPointKey & other) const {
      if(this->a < other.a)
        return true;
      else if(this->a > other.a)
        return false;
      else
        if(this->b < other.b)
          return true;
        else
          return false;
    };
  };

	// midpoints
	std::map<MidPointKey, unsigned int> midpoints;

	/// Adds a midpoint as a vertex
	/// @param[in] a index of the first vertex
	/// @param[in] b index of the second vertex
	/// @return index of the newly created vertex
	unsigned int create_midpoint(unsigned int a, unsigned int b);
	/// Gets an index of a midpoint between a and b. If midpoint does not exists, it is created
	/// @param[in] a index of the first vertex
	/// @param[in] b index of the second vertex
	/// @return index of the midpoint
	/// FIXME: better name
	unsigned int get_midpoint(unsigned int a, unsigned int b);
	/// Sets the midpoints index for a midpoint
	/// @param[in] a index of the first vertex
	/// @param[in] b index of the second vertex
	void set_midpoint(unsigned int a, unsigned int b, unsigned int idx);

	/// referencing edges
	void ref_edges(Element *e);
	void unref_edges(Element *e);

	/// Check element orientations
	void check_elem_oris();

	// Seq (internal use only)
	int seq;
	/// For internal use.
	int get_seq() const { return seq; }
	/// For internal use.
	void set_seq(unsigned seq) { this->seq = seq; }

	friend class Space;
	friend class WeakForm;
};


#endif
