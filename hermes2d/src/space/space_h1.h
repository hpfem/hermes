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

#ifndef __H2D_SPACE_H1_H
#define __H2D_SPACE_H1_H

#include "space.h"


/// H1Space represents a space of continuous scalar functions over a domain (mesh).
///
///
///
class HERMES_API H1Space : public Space
{
public:
  H1Space(Mesh* mesh, EssentialBCs* boundary_conditions, int p_init = 1,
          Shapeset* shapeset = NULL);

  H1Space(Mesh* mesh, int p_init = 1,
          Shapeset* shapeset = NULL);

  // Common code for the constructors.
  void init(Shapeset* shapeset, Ord2 p_init);

  virtual ~H1Space();

  virtual void set_shapeset(Shapeset* shapeset);

  /// Removes the degree of freedom from a vertex node with the given id (i.e., its number
  /// in the mesh file) and makes it part of the Dirichlet lift with the given value.
  /// This is a special-purpose function which normally should not be needed.
  /// It is intended for fixing the solution of a system which would otherwise be singular
  /// and for some reason a standard Dirichlet condition (with non-zero measure on the
  /// boundary) is not suitable.
  void fix_vertex(int id, scalar value = 0.0);

  virtual Space* dup(Mesh* mesh, int order_increase = 0) const;

  virtual ESpaceType get_type() const { return HERMES_H1_SPACE; }

protected:

  virtual void assign_vertex_dofs();
  virtual void assign_edge_dofs() {};
  virtual void assign_bubble_dofs() {};

  virtual void get_vertex_assembly_list(Element* e, int iv, AsmList* al);
  virtual void get_boundary_assembly_list_internal(Element* e, int ie, AsmList* al);

  static double** h1_proj_mat;
  static double*  h1_chol_p;
  static int      h1_proj_ref;

  virtual scalar* get_bc_projection(SurfPos* surf_pos, int order);

  struct EdgeInfo
  {
    Node* node;
    int part;
    int ori;
    double lo, hi;
  };

  inline void output_component(BaseComponent*& current, BaseComponent*& last, BaseComponent* min,
                               Node*& edge, BaseComponent*& edge_dofs);

  BaseComponent* merge_baselists(BaseComponent* l1, int n1, BaseComponent* l2, int n2,
                                 Node* edge, BaseComponent*& edge_dofs, int& ncomponents);

  void update_constrained_nodes(Element* e, EdgeInfo* ei0, EdgeInfo* ei1, EdgeInfo* ei2, EdgeInfo* ei3);
  virtual void update_constraints();

  struct FixedVertex
  {
    int id;
    scalar value;
  };

  std::vector<FixedVertex> fixed_vertices;

  inline bool is_fixed_vertex(int id) const;
  virtual void post_assign();

  //void dump_baselist(NodeData& nd);

};


#endif
