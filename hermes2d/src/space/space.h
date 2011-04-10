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

#ifndef __H2D_SPACE_H
#define __H2D_SPACE_H

#include "../mesh/mesh.h"
#include "../mesh/traverse.h"
#include "../shapeset/shapeset.h"
#include "../asmlist.h"
#include "../shapeset/precalc.h"
#include "../quadrature/quad_all.h"
#include "../boundaryconditions/essential_bcs.h"

/// \brief Represents a finite element space over a domain.
///
/// The Space class represents a finite element space over a domain defined by 'mesh', spanned
/// by basis functions constructed using 'shapeset'. It serves as a base class for H1Space,
/// HcurlSpace, HhivSpace and L2Space, since most of the functionality is common for all these spaces.
///
/// There are four main functions the Space class provides:
/// <ol>
///    <li> It handles the Dirichlet (essential) boundary conditions. The user provides a pointer
///         to an instance of the EssentialBCs class that determines
///         which markers represent the Dirichlet part of the boundary, and which markers represent
///         the Neumann and Newton parts. It is also possible to use the value BC_NONE
///         which supresses all BC processing on such part of the boundary.
///
///         The default BC type is BC_NATURAL for the whole boundary. The default BC value is zero for all
///         markers.
///
///    <li> It stores element polynomial degrees, or 'orders'. All active elements need to have
///         an order set for the Space to be valid. Individual orders can be set by calling
///         set_element_order(). You can also set the same order for all elements using
///         set_uniform_order(). Quadrilateral elements can have different orders in the vertical
///         and horizontal directions. It is therefore necessary to form the order using the macro
///         H2D_MAKE_QUAD_ORDER() when calling the aforementioned functions.
///
///    <li> It builds and enumerates the basis functions. After all element orders have been set,
///         you must call the function assign_dofs(). This function assigns the DOF (degree-of-
///         freedom) numbers to basis functions, starting with 'first_dof' (optional parameter).
///         It also determines constraining relationships in the mesh due to hanging nodes and
///         builds constrained basis functions. The total number of basis functions can then
///         be obtained by calling  (). Standard basis functions are assigned positive
///         numbers from 'first_dof' to ('first_dof' + (get_num_dofs() - 1) * 'stride'). All
///         shape functions belonging to the Dirichlet lift are assigned DOF number of -1. This
///         way the Dirichlet lift becomes a (virtual) basis function. This simplifies assembling.
///
///    <li> Finally, and most importantly, the Space is able to deliver a list of shape functions
///         existing on each element. Such a list is called an "assembly list" and is represented
///         by the class AsmList. The assembly list contains the triplets (idx, dof, coef).
///         'idx' is the shape function index, as understood by the Shapeset. 'dof' is the number
///         of the basis function, whose part the shape function forms. 'coef' is a real constant,
///         with which the shape function must be multiplied in order to fit into the basis
///         function. This is typically 1, but can be less than that for constrained functions.
///         Constrained vertex functions can belong to more than one basis functions. This
///         results in more triplets with the same 'idx'. However, the assembling procedure or
///         the Solution class do not have to worry about that. For instance, the Solution class
///         simply obtains the values of all shape functions contained in the list, multiplies them
///         by 'coef', and forms a linear combination of them by taking the solution vector values
///         at the 'dof' positions. This way the solution to the PDE is obtained.
/// </ol>
///
/// Space is an abstract class and cannot be instatiated. Use one of the specializations H1Space,
/// HcurlSpace or L2Space instead.
///
/// The handling of irregular meshes is desribed in H1Space and HcurlSpace.
///
class Ord2;

class HERMES_API Space
{
public:
  Space(Mesh* mesh, Shapeset* shapeset, EssentialBCs* essential_bcs, Ord2 p_init);

  virtual ~Space();
  virtual void free();

  /// Sets the boundary condition.
  void set_essential_bcs(EssentialBCs* essential_bcs);
  /// Sets element polynomial order. Can be called by the user. Should not be called
  /// for many elements at once, since assign_dofs() is called at the end of this function.
  virtual void set_element_order(int id, int order);
  /// Sets polynomial order to all elements.
  virtual void set_element_orders(int* elem_orders);
  /// Sets element polynomial order. This version does not call assign_dofs() and is
  /// intended primarily for internal use.
  virtual void set_element_order_internal(int id, int order);
  /// Returns element polynomial order.
  int get_element_order(int id) const;
  /// Sets the same polynomial order for all elements in the mesh. Intended for
  /// the user and thus assign_dofs() is called at the end of this function.
  void set_uniform_order(int order, std::string marker = HERMES_ANY);
  /// Sets the same polynomial order for all elements in the mesh. Does not
  /// call assign_dofs(). For internal use.
  void set_uniform_order_internal(Ord2 order, int marker);
  /// Sets the order automatically assigned to all newly created elements.
  /// (The order of these is normally undefined and has to be set explicitly.)
  void set_default_order(int tri_order, int quad_order = -1);
  /// Set the element order relative to the current order.
  /// The parameter min_order prevents decreasing polynomial order below this threshold.
  void adjust_element_order(int order_change, int min_order);
  /// Version for quads.
  void adjust_element_order(int horizontal_order_change, int vertical_order_change, unsigned int horizontal_min_order, unsigned int vertical_min_order);
  
  /// Sets the shapeset.
  virtual void set_shapeset(Shapeset* shapeset) = 0;
  /// Copies element orders from another space. 'inc' is an optional order
  /// increase. If the source space has a coarser mesh, the orders are distributed
  /// recursively. This is useful for reference solution spaces.
  void copy_orders(const Space* space, int inc = 0);
  /// Internal. Obtains the order of an edge, according to the minimum rule.
  virtual int get_edge_order(Element* e, int edge);

  /// \brief Builds basis functions and assigns DOF numbers to them.
  /// \details This functions must be called \b after assigning element orders, and \b before
  /// using the space in a computation, otherwise an error will occur.
  /// \param first_dof [in] The DOF number of the first basis function.
  /// \param stride [in] The difference between the DOF numbers of successive basis functions.
  /// \return The number of basis functions contained in the space.
  virtual int assign_dofs(int first_dof = 0, int stride = 1);

  /// \brief Returns the number of basis functions contained in the space.
  int get_num_dofs() { return ndof; }
  /// \brief Returns the DOF number of the last basis function.
  int get_max_dof() const { return next_dof - stride; }

  Shapeset* get_shapeset() const { return shapeset; }
  Mesh* get_mesh() const { return mesh; }
  void set_mesh(Mesh* mesh);

  /// Creates a copy of the space, increases order of all elements by
  /// "order_increase".
  virtual Space* dup(Mesh* mesh, int order_increase = 0) const = 0;

  /// Returns true if the space is ready for computation, false otherwise.
  bool is_up_to_date() const {
      return was_assigned && mesh_seq == (int) mesh->get_seq();
  }

  /// Sets polynomial orders to elements created by Mesh::regularize() using "parents".
  void distribute_orders(Mesh* mesh, int* parents);

public:

  /// Number of degrees of freedom (dimension of the space).
  int ndof;

  /// Obtains an boundary conditions
  inline EssentialBCs* get_essential_bcs() { return essential_bcs; }

  /// Obtains an assembly list for the given element.
  virtual void get_element_assembly_list(Element* e, AsmList* al);

  /// Obtains an edge assembly list (contains shape functions that are nonzero on the specified edge).
  void get_boundary_assembly_list(Element* e, int surf_num, AsmList* al);

  /// Updates essential BC values. Typically used for time-dependent
  /// essnetial boundary conditions.
  void update_essential_bc_values();

  /// \brief Returns the number of basis functions contained in the spaces.
  static int get_num_dofs(Hermes::vector<Space *> spaces);
  static int get_num_dofs(Space* space);

  /// \brief Assings the degrees of freedom to all Spaces in the Hermes::vector.
  static int assign_dofs(Hermes::vector<Space*> spaces);

protected:
  static const int H2D_UNASSIGNED_DOF = -2; ///< DOF which was not assigned yet.
  static const int H2D_CONSTRAINED_DOF = -1; ///< DOF which is constrained.

  Shapeset* shapeset;
  bool own_shapeset;  ///< true if default shapeset is created in the constructor, false if shapeset is supplied by user.

  // boundary conditions
  EssentialBCs* essential_bcs;

  /// FE mesh
  Mesh* mesh;

  int default_tri_order, default_quad_order;
  int first_dof, next_dof;
  int stride;
  int seq, mesh_seq;
  bool was_assigned;

  struct BaseComponent
  {
    int dof;
    scalar coef;
  };

  union NodeData
  {
    struct // regular node
    {
      int dof;
      union {
        scalar* edge_bc_proj;
        scalar* vertex_bc_coef;
      };
      int n; ///< Number of dofs. Temporarily used during assignment
             ///< of DOFs to indicate nodes which were not processed yet.
    };
    struct // constrained vertex node
    {
      BaseComponent* baselist;
      int ncomponents;
    };
    struct // constrained edge node
    {
      Node* base;
      int part;
    };
    NodeData() : dof(0), edge_bc_proj(NULL) {}
  };

  struct ElementData
  {
    int order;
    int bdof, n;
  };

public:
  NodeData* ndata;    ///< node data table
  int nsize, ndata_allocated; ///< number of items in ndata, allocated space
  ElementData* edata; ///< element data table
  int esize;

protected: //debugging support
  virtual int get_edge_order_internal(Node* en);

  /// \brief Updates internal node and element tables.
  /// \details Since meshes only contain geometric information, the Space class keeps two
  /// tables with FEM-related information. The first one, 'ndata', contains DOF numbers
  /// and other things for each node. The second table, 'edata', holds element orders
  /// and bubble DOF numbers. Both tables are directly indexed by the node and element
  /// IDs. The function resize_tables() is called to check whether the tables are large
  /// enough to contain all node and element id's, and to reallocate them if not.
  virtual void resize_tables();

  void H2D_CHECK_ORDER(int order);
  void copy_orders_recurrent(Element* e, int order);

  virtual void reset_dof_assignment(); ///< Resets assignment of DOF to an unassigned state.
  virtual void assign_vertex_dofs() = 0;
  virtual void assign_edge_dofs() = 0;
  virtual void assign_bubble_dofs() = 0;

  virtual void get_vertex_assembly_list(Element* e, int iv, AsmList* al) = 0;
  virtual void get_boundary_assembly_list_internal(Element* e, int surf_num, AsmList* al) = 0;
  virtual void get_bubble_assembly_list(Element* e, AsmList* al);

  double** proj_mat;
  double*  chol_p;

  void copy_callbacks(const Space* space);
  void precalculate_projection_matrix(int nv, double**& mat, double*& p);
  virtual scalar* get_bc_projection(SurfPos* surf_pos, int order) = 0;
  void update_edge_bc(Element* e, SurfPos* surf_pos);

  /// Called by Space to update constraining relationships between shape functions due
  /// to hanging nodes in the mesh. As this is space-specific, this function is reimplemented
  /// in H1Space and HcurlSpace.
  virtual void update_constraints() {}

  /// Auxiliary function the descendants may implement to perform additional tasks after
  /// the DOFs have been assigned.
  virtual void post_assign() {}

  std::vector<void*> extra_data;
  void free_extra_data();

  void propagate_zero_orders(Element* e);

public:
  /// Internal. Used by DiscreteProblem to detect changes in the space.
  int get_seq() const { return seq; }
  int set_seq(int seq_) { seq = seq_; return seq;}

  /// Internal. Return type of this space (H1 = HERMES_H1_SPACE, Hcurl = HERMES_HCURL_SPACE,
  /// Hdiv = HERMES_HDIV_SPACE, L2 = HERMES_L2_SPACE)
  virtual ESpaceType get_type() const = 0;

  /// Create globally refined space.
  static Hermes::vector<Space *>* construct_refined_spaces(Hermes::vector<Space *> coarse, int order_increase = 1);
  static Space* construct_refined_space(Space* coarse, int order_increase = 1);

  // updating time-dependent essential (Dirichlet) boundary conditions
  static void update_essential_bc_values(Hermes::vector<Space*> spaces, double time);  // multiple spaces
  static void update_essential_bc_values(Space *s, double time);    // one space

};


class Ord2
{
  public:
    Ord2(int order_h, int order_v) : order_h(order_h), order_v(order_v) {};
    Ord2(int order) : order_h(order), order_v(order) {};
    int order_h;
    int order_v;
};


#endif
