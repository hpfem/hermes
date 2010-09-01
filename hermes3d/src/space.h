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

#ifndef __SPACE_H
#define __SPACE_H

#include "mesh.h"
#include "shapeset.h"
#include "asmlist.h"
#include "quad.h"
#include "order.h"

#include <common/bitarray.h>

/// @defgroup spaces Spaces
///
/// TODO: description

// for iterating through nodes in space classes

///< Iterates over all vertex node indices.
///< \param idx Vertex node hash table index.
#define FOR_ALL_VERTEX_NODES(idx) \
		for (Word_t (idx) = vn_data.first(); (idx) != INVALID_IDX; (idx) = vn_data.next((idx)))

///< Iterates over all edge node indices.
///< \param idx Edge node hash table index.
#define FOR_ALL_EDGE_NODES(idx) \
		for (Word_t (idx) = en_data.first(); (idx) != INVALID_IDX; (idx) = en_data.next((idx)))

///< Iterates over all face node indices.
///< \param idx Face node hash table index.
#define FOR_ALL_FACE_NODES(idx) \
		for (Word_t (idx) = fn_data.first(); (idx) != INVALID_IDX; (idx) = fn_data.next((idx)))

///< Iterates over all element node indices.
///< \param idx Element node hash table index.
#define FOR_ALL_ELEMENT_NODES(idx) \
		for (Word_t (idx) = elm_data.first(); (idx) != INVALID_IDX; (idx) = elm_data.next((idx)))


/// Possible return values for bc_type_callback():
enum BCType 
{
  BC_ESSENTIAL, /// Essential (Dirichlet) BC.
  BC_NATURAL,   /// Natural (Neumann, Newton) BC.
  BC_NONE       /// Hermes will not attempt to evaluate any boundary
                /// integrals on this part of the boundary.
};


#define H3D_MARKER_UNDEFINED				-1

#define H3D_DOF_UNASSIGNED					-2
#define H3D_DOF_NOT_ANALYZED				-3

/// Base class for all spaces
///
/// The Space class represents a finite element space over a domain defined
/// by 'mesh', consisting of basis functions constructed using 'shapeset'.
///
/// @ingroup spaces
class Space {
public:
	Space(Mesh *mesh, Shapeset *shapeset);
	virtual ~Space();

	virtual Space *dup(Mesh *mesh) const = 0;

	ESpaceType get_type() { return type; }

	void set_bc_types(BCType (*bc_type_callback)(int marker));
	void set_essential_bc_values(scalar (*bc_value_callback_by_coord)(int ess_bdy_marker, double x, double y, double z));
	// TODO: different callback: void (*bc_vec_value_callback_by_coord)(int marker, double x, double y, double z, scalar3 &result)
	void set_essential_bc_values(scalar3 &(*bc_vec_value_callback_by_coord)(int ess_bdy_marker, double x, double y, double z));

	void set_element_order(Word_t eid, order3_t order);
	order3_t get_element_order(Word_t eid) const;
	void set_uniform_order(order3_t order);
	void copy_orders(const Space &space, int inc = 0);

	virtual void enforce_minimum_rule();
	virtual int assign_dofs(int first_dof = 0, int stride = 1);
	int get_dof_count() const { return (next_dof - first_dof) / stride; }
	int get_max_dof() const { return next_dof - stride; }

	Shapeset *get_shapeset() const { return shapeset; }
	Mesh *get_mesh() const { return mesh; }

   	virtual void get_element_assembly_list(Element *e, AsmList *al) = 0;
	virtual void get_boundary_assembly_list(Element *e, int face, AsmList *al) = 0;

	void dump();

	/// Returns true if the space is ready for computation, false otherwise.
	bool is_up_to_date() const { return was_assigned && mesh_seq == mesh->get_seq(); }

protected:
	Mesh *mesh;
	Shapeset *shapeset;
	ESpaceType type;

	int first_dof, next_dof, first_bubble;
	int stride;

	int seq, mesh_seq;
	bool was_assigned;

	// CED
	struct BaseVertexComponent {
		int dof;
		scalar coef;
	};

	struct BaseEdgeComponent {
		Word_t edge_id;							/// ID of the constraining edge
		int ori;								/// the orientation of the constraining edge
		Part part;								/// part of the edge that is constrained
		scalar coef;
	};

	struct BaseFaceComponent {
		Word_t face_id;							/// ID of a constraining face
		unsigned ori:3;							/// the orientation of constraining face
		unsigned dir:1;							/// the orientation of ???
		unsigned iface:4;						/// local number of constraining face
		Part part;								/// part of the face that is constrained
		scalar coef;

		BaseFaceComponent() {
		}

		~BaseFaceComponent() {
		}
	};

	struct NodeData {
		int marker;
		BCType bc_type;

		NodeData() {
			marker = H3D_MARKER_UNDEFINED;
			bc_type = BC_NONE;
		}
	};

	struct VertexData : public NodeData {
		unsigned ced:1;								/// 0 = normal vertex, 1 = constrained vertex
		union {
			/// normal node
			struct {
				int dof;
				int n;								/// number of DOFs
			};
			/// CED node
			struct {
				int ncomponents;
				BaseVertexComponent *baselist;
			};
		};

		scalar bc_proj;								/// projection of boundary condition

		VertexData() {
			bc_proj = 0.0;
		}

		~VertexData() {
			if (ced) ::free(baselist);
		}

		void dump(int id);
	};

	struct EdgeData : public NodeData {
		unsigned ced:1;								/// 1 = is constrained
		union {
			/// normal node
			struct {
				order1_t order;   					/// polynomial order
				int dof;
				int n;								/// number of DOFs
			};
			/// CED
			struct {
				BaseEdgeComponent *edge_baselist;
				int edge_ncomponents;
				BaseFaceComponent *face_baselist;
				int face_ncomponents;
			};
		};

		scalar *bc_proj;

		EdgeData() {
			bc_proj = NULL;
		}

		virtual ~EdgeData() {
			delete [] bc_proj;
			if (ced) {
				::free(edge_baselist);
				::free(face_baselist);
			}
		}

		void dump(int id);
	};

	struct FaceData : public NodeData  {
		unsigned ced:1;								/// 1 = is constrained
		order2_t order;   							/// polynomial order
		union {
			struct {								/// normal node
				int dof;
				int n;								/// number of DOFs
			};
			struct {								/// CED node
				Word_t facet_id;					/// ID of a facing facet
				int ori;							/// orientation of facing facet
				Part part;
			};
		};

		scalar *bc_proj;

		FaceData() {
			bc_proj = NULL;
		}

		virtual ~FaceData() {
			delete [] bc_proj;
		}

		void dump(int id);
	};

	struct ElementData {
		order3_t order;								/// Polynomial degree associated to the element node (interior).
		int dof;									/// The number of the first degree of freedom belonging to the node.
		int n;										/// Total number of degrees of freedom belonging to the node.

		ElementData() {
			order = -1;
			dof = H3D_DOF_NOT_ANALYZED;
			n = -1;
		}

		void dump(int id);
	};

	ArrayPtr<VertexData> vn_data;					/// Vertex node hash table
	ArrayPtr<EdgeData> en_data;						/// Edge node hash table
	ArrayPtr<FaceData> fn_data;						/// Face node hash table
	ArrayPtr<ElementData> elm_data;					/// Element node hash table

	void set_order_recurrent(Word_t eid, order3_t order);

	virtual int get_vertex_ndofs() = 0;
	virtual int get_edge_ndofs(order1_t order) = 0;
	virtual int get_face_ndofs(order2_t order) = 0;
	virtual int get_element_ndofs(order3_t order) = 0;

	virtual void assign_vertex_dofs(Word_t vid);
	virtual void assign_edge_dofs(Word_t eid);
	virtual void assign_face_dofs(Word_t fid);
	virtual void assign_bubble_dofs(Word_t eid);

	virtual void assign_dofs_internal() = 0;

	virtual void get_vertex_assembly_list(Element *e, int ivertex, AsmList *al);
	virtual void get_edge_assembly_list(Element *e, int iedge, AsmList *al);
	virtual void get_face_assembly_list(Element *e, int iface, AsmList *al);
	virtual void get_bubble_assembly_list(Element *e, AsmList *al);

	virtual void calc_boundary_projections();
	virtual void calc_vertex_boundary_projection(Element *elem, int ivertex) = 0;
	virtual void calc_edge_boundary_projection(Element *elem, int iedge) = 0;
	virtual void calc_face_boundary_projection(Element *elem, int iface) = 0;

	void set_bc_info(NodeData *node, BCType bc, int marker);
	void set_bc_information();
	void copy_callbacks(const Space *space);

	void init_data_tables();
	void free_data_tables();

	// CED
	struct FaceInfo {
		Word_t elem_id;
		int face;

		unsigned type:1;				// 1 - quad, 0 - triangle
		union {
			struct {					// quad part
				int h_part;				// horizontal part
				int v_part;				// vertical part
				double h_lo, h_hi;		// limits in horizontal direction
				double v_lo, v_hi;		// limits in vertical direction
			};
			struct {					// triangle part
				// TODO: implement me
			};
		};

		FaceInfo(EMode2D mode, Word_t elem_id, int face) {
			this->type = mode == MODE_QUAD;
			this->elem_id = elem_id;
			this->face = face;
			this->h_part = this->v_part = 0;
			this->h_lo = this->v_lo = -1.0;
			this->h_hi = this->v_hi =  1.0;
		}
	};

	void update_constraints();

	// find constraints
	void find_constraints();
	void fc_base(Word_t eid, int iface);
	/// @param[in] eid - ID of the element
	/// @param[in] iface - local number of the face on the element eid
	void fc_face(Word_t eid, int iface, bool ced);
	/// @param[in] fid - ID of the facet
	void fc_face_left(Word_t fid);
	/// @param[in] fid - ID of the facet
	void fc_face_right(Word_t fid);
	/// @param[in] idx - ID of the element
	void fc_element(Word_t idx);
	BitArray face_ced;

	// update constraints
	void uc_element(Word_t idx);
	void uc_face(Word_t eid, int iface);
	void uc_dep(Word_t eid);
	BitArray uc_deps;

	Array<FaceInfo *> fi_data;

	VertexData *create_vertex_node_data(Word_t vid, bool ced);
	EdgeData *create_edge_node_data(Word_t eid, bool ced);
	FaceData *create_face_node_data(Word_t fid, bool ced);

	void output_component(BaseVertexComponent *&current, BaseVertexComponent *&last, BaseVertexComponent *min, bool add);
	void output_component(BaseEdgeComponent *&current, BaseEdgeComponent *&last, BaseEdgeComponent *min, bool add);
	void output_component(BaseFaceComponent *&current, BaseFaceComponent *&last, BaseFaceComponent *min, bool add);
	void output_component_over(BaseFaceComponent *&current, BaseFaceComponent *min, BaseFaceComponent *m);

	BaseVertexComponent *merge_baselist(BaseVertexComponent *l1, int n1, BaseVertexComponent *l2, int n2, int &ncomponents, bool add);
	BaseEdgeComponent *merge_baselist(BaseEdgeComponent *l1, int n1, BaseEdgeComponent *l2, int n2, int &ncomponents, bool add);
	BaseFaceComponent *merge_baselist(BaseFaceComponent *l1, int n1, BaseFaceComponent *l2, int n2, int &ncomponents, Word_t fid, bool add);

	// all these work for hexahedra
	void calc_vertex_vertex_ced(Word_t vtx1, Word_t vtx2);
	void calc_vertex_edge_ced(Word_t vtx, Word_t edge_id, int ori, int part);
	void calc_vertex_face_ced(Word_t vtx, Word_t fid, int ori, int iface, int hpart, int vpart);
	void calc_edge_edge_ced(Word_t seid, Word_t eid, int ori, int epart, int part);
	void calc_edge_face_ced(Word_t mid_eid, Word_t eid[], Word_t fid, int ori, int iface, int part_ori, int fpart, int epart);
	void calc_face_face_ced(Word_t sfid, Word_t fid, int ori, int hpart, int vpart);

	void calc_mid_vertex_vertex_ced(Word_t mid, Word_t vtx1, Word_t vtx2, Word_t vtx3, Word_t vtx4);
	void calc_mid_vertex_edge_ced(Word_t vtx, Word_t fmp, Word_t eid, int ori, int part);
	void calc_mid_edge_edge_ced(Word_t meid, Word_t eid[], int ori[], int epart, int part);

public:
	BCType (*bc_type_callback)(int);

	// value callbacks for dirichlet
	scalar (*bc_value_callback_by_coord)(int ess_bdy_marker, double x, double y, double z);
	scalar3 &(*bc_vec_value_callback_by_coord)(int ess_bdy_marker, double x, double y, double z);

protected:
	/// Internal. Used by LinProblem to detect changes in the space.
	int get_seq() const { return seq; }

	friend class DiscreteProblem;
	friend class LinearProblem;
};


#endif
