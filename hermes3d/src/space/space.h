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
#include "../../../hermes_common/vector.h"
#include "shapeset/shapeset.h"
#include "asmlist.h"
#include "quad.h"
#include "order.h"

#include "bctypes.h"

/// @defgroup spaces Spaces
///
/// TODO: description
#define H3D_MARKER_UNDEFINED				-1

#define H3D_DOF_UNASSIGNED					-2
#define H3D_DOF_NOT_ANALYZED				-3

/// Base class for all spaces
///
/// The Space class represents a finite element space over a domain defined
/// by 'mesh', consisting of basis functions constructed using 'shapeset'.
///
/// @ingroup spaces
class HERMES_API Space {
public:
  Space(Mesh *mesh, Shapeset *shapeset, BCType (*bc_type_callback)(int), 
        scalar (*bc_value_callback_by_coord)(int, double, double, double), Ord3 p_init = Ord3(1,1,1));
	virtual ~Space();

  virtual Space *dup(Mesh *mesh) const = 0;

  ESpaceType get_type() { return type; }

  void set_bc_types(BCType (*bc_type_callback)(int marker));
  void set_bc_types_init(BCType (*bc_type_callback)(int marker));
  void set_essential_bc_values(scalar (*bc_value_callback_by_coord)(int ess_bdy_marker, double x, double y, double z));
  // TODO: different callback: void (*bc_vec_value_callback_by_coord)(int marker, double x, double y, double z, scalar3 &result)
  void set_essential_bc_values(scalar3 &(*bc_vec_value_callback_by_coord)(int ess_bdy_marker, double x, double y, double z));

  /// Sets the shapeset.
  virtual void set_shapeset(Shapeset* shapeset) = 0;
  void set_element_order(unsigned int eid, Ord3 order);
  Ord3 get_element_order(unsigned int eid) const;
  /// Sets the same polynomial order for all elements in the mesh. Intended for 
  /// the user and thus assign_dofs() is called at the end of this function.
  void set_uniform_order(Ord3 order, int marker = HERMES_ANY_INT);
  /// Sets the same polynomial order for all elements in the mesh. Does not 
  /// call assign_dofs(). For internal use.
  void set_uniform_order_internal(Ord3 order, int marker = HERMES_ANY_INT);

  void copy_orders(const Space &space, int inc = 0);

  virtual void enforce_minimum_rule();
  virtual int assign_dofs(int first_dof = 0, int stride = 1);

  /// \brief Returns the number of basis functions contained in the space.
  int get_num_dofs() { return ndof; }

  int get_dof_count() const { return (next_dof - first_dof) / stride; }
  /// \brief Returns the DOF number of the last basis function.
  int get_max_dof() const { return next_dof - stride; }

  Shapeset *get_shapeset() const { return shapeset; }
  Mesh *get_mesh() const { return mesh; }

  virtual void get_element_assembly_list(Element *e, AsmList *al) = 0;
  virtual void get_boundary_assembly_list(Element *e, int face, AsmList *al) = 0;

  void dump();

  /// Returns true if the space is ready for computation, false otherwise.
  bool is_up_to_date() const { return was_assigned && mesh_seq == mesh->get_seq(); }

  /// Number of degrees of freedom (dimension of the space)
  int ndof;

  /// FE mesh
  Mesh *mesh;

  /// Create globally refined space.
  static Hermes::vector<Space *>* construct_refined_spaces(Hermes::vector<Space *> coarse, int order_increase = 1);
  static Space* construct_refined_space(Space* coarse, int order_increase = 1);

  static int get_num_dofs(Hermes::vector<Space *> spaces);
  static int get_num_dofs(Space* space);
  static int assign_dofs(Hermes::vector<Space*> spaces);

protected:
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
    Edge::Key edge_id;							/// ID of the constraining edge
    int ori;								/// the orientation of the constraining edge
    Part part;								/// part of the edge that is constrained
    scalar coef;
  };

  struct BaseFaceComponent {
    Facet::Key face_id;							/// ID of a constraining face
    unsigned ori:3;							/// the orientation of constraining face
    unsigned dir:1;							/// the orientation of ???
    unsigned iface:4;						/// local number of constraining face
    Part part;								/// part of the face that is constrained
    scalar coef;

    BaseFaceComponent() {}

    ~BaseFaceComponent() {}
  };

  struct NodeData {
    int marker;
    BCType bc_type;

    NodeData() {
      marker = H3D_MARKER_UNDEFINED;
      bc_type = H3D_BC_NONE;
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
	      Ord1 order;   					/// polynomial order
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

    void dump(Edge::Key id);
  };

  struct FaceData : public NodeData  {

    unsigned ced:1;								/// 1 = is constrained
    Ord2 order;   							/// polynomial order
    int dof;
    int n;								/// number of DOFs
    Facet::Key facet_id;					/// ID of a facing facet
    int ori;							/// orientation of facing facet
    Part part;

    scalar *bc_proj;

    FaceData() {
      bc_proj = NULL;
    }

    virtual ~FaceData() {
      delete [] bc_proj;
    }

    void dump(Facet::Key id);
  };

  struct ElementData {
    Ord3 order;	      /// Polynomial degree associated with the element node (interior).
    int dof;	      /// The number of the first degree of freedom belonging to the node.
    int n;		      /// Total number of degrees of freedom belonging to the node.
    int marker;           /// Material marker.  

    ElementData() {
      order = -1;
      dof = H3D_DOF_NOT_ANALYZED;
      n = -1;
      marker = -1;
    }

    void dump(int id);
  };

  std::map<unsigned int, VertexData *> vn_data;		/// Vertex node hash table
  std::map<Edge::Key, EdgeData *> en_data;		/// Edge node hash table
  std::map<Facet::Key, FaceData *> fn_data;		/// Face node hash table
  std::map<unsigned int, ElementData *> elm_data;		/// Element node hash table

  void set_order_recurrent(unsigned int eid, Ord3 order);

  virtual int get_vertex_ndofs() = 0;
  virtual int get_edge_ndofs(Ord1 order) = 0;
  virtual int get_face_ndofs(Ord2 order) = 0;
  virtual int get_element_ndofs(Ord3 order) = 0;

  virtual void assign_vertex_dofs(unsigned int vid);
  virtual void assign_edge_dofs(Edge::Key eid);
  virtual void assign_face_dofs(Facet::Key fid);
  virtual void assign_bubble_dofs(unsigned int eid);

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
    unsigned int elem_id;
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

    FaceInfo(ElementMode2D mode, unsigned int elem_id, int face) {
      this->type = mode == HERMES_MODE_QUAD;
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
  void fc_base(unsigned int eid, int iface);
  /// @param[in] eid - ID of the element
  /// @param[in] iface - local number of the face on the element eid
  void fc_face(unsigned int eid, int iface, bool ced);
  /// @param[in] fid - ID of the facet
  void fc_face_left(Facet::Key fid);
  /// @param[in] fid - ID of the facet
  void fc_face_right(Facet::Key fid);
  /// @param[in] idx - ID of the element
  void fc_element(unsigned int idx);
  std::map<Facet::Key, bool> face_ced;

  // update constraints
  void uc_element(unsigned int idx);
  void uc_face(unsigned int eid, int iface);
  void uc_dep(unsigned int eid);
  std::map<unsigned int, bool> uc_deps;

  std::map<Facet::Key, FaceInfo *> fi_data;

  VertexData *create_vertex_node_data(unsigned int vid, bool ced);
  EdgeData *create_edge_node_data(Edge::Key eid, bool ced);
  FaceData *create_face_node_data(Facet::Key fid, bool ced);

  void output_component(BaseVertexComponent *&current, BaseVertexComponent *&last, BaseVertexComponent *min, bool add);
  void output_component(BaseEdgeComponent *&current, BaseEdgeComponent *&last, BaseEdgeComponent *min, bool add);
  void output_component(BaseFaceComponent *&current, BaseFaceComponent *&last, BaseFaceComponent *min, bool add);
  void output_component_over(BaseFaceComponent *&current, BaseFaceComponent *min, BaseFaceComponent *m);

  BaseVertexComponent *merge_baselist(BaseVertexComponent *l1, int n1, BaseVertexComponent *l2, int n2, int &ncomponents, bool add);
  BaseEdgeComponent *merge_baselist(BaseEdgeComponent *l1, int n1, BaseEdgeComponent *l2, int n2, int &ncomponents, bool add);
  BaseFaceComponent *merge_baselist(BaseFaceComponent *l1, int n1, BaseFaceComponent *l2, int n2, int &ncomponents, Facet::Key fid, bool add);

  // all these work for hexahedra
  void calc_vertex_vertex_ced(unsigned int vtx1, unsigned int vtx2);
  void calc_vertex_edge_ced(unsigned int vtx, Edge::Key edge_id, int ori, int part);
  void calc_vertex_face_ced(unsigned int vtx, Facet::Key fid, int ori, int iface, int hpart, int vpart);
  void calc_edge_edge_ced(Edge::Key seid, Edge::Key eid, int ori, int epart, int part);
  void calc_edge_face_ced(Edge::Key mid_eid, Edge::Key eid[], Facet::Key fid, int ori, int iface, int part_ori, int fpart, int epart);
  void calc_face_face_ced(Facet::Key sfid, Facet::Key fid, int ori, int hpart, int vpart);

  void calc_mid_vertex_vertex_ced(unsigned int mid, unsigned int vtx1, unsigned int vtx2, unsigned int vtx3, unsigned int vtx4);
  void calc_mid_vertex_edge_ced(unsigned int vtx, unsigned int fmp, Edge::Key eid, int ori, int part);
  void calc_mid_edge_edge_ced(Edge::Key meid, Edge::Key eid[], int ori[], int epart, int part);

public:
  BCType (*bc_type_callback)(int);

  // value callbacks for dirichlet
  scalar (*bc_value_callback_by_coord)(int ess_bdy_marker, double x, double y, double z);
  scalar3 &(*bc_vec_value_callback_by_coord)(int ess_bdy_marker, double x, double y, double z);

protected:
  /// Internal. Used by DiscreteProblem to detect changes in the space.
  int get_seq() const { return seq; }

  friend class DiscreteProblem;
};

#endif
