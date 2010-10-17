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

#ifndef __H2D_SPACE_HDIV_H
#define __H2D_SPACE_HDIV_H

#include "space.h"


/// HdivSpace represents a space of vector functions with continuous normal
/// components over a domain (mesh).
///
///
///
class H2D_API HdivSpace : public Space
{
public:
  // For backward compatibility.
  HdivSpace(Mesh* mesh, BCType (*bc_type_callback)(int), 
                 scalar (*bc_value_callback_by_coord)(int, double, double), int p_init, 
                 Shapeset* shapeset = NULL);
  HdivSpace(Mesh* mesh = NULL, BCType (*bc_type_callback)(int) = NULL, 
                 scalar (*bc_value_callback_by_coord)(int, double, double) = NULL, Ord2 p_init = Ord2(1,1), 
                 Shapeset* shapeset = NULL);
  
  virtual ~HdivSpace();

  virtual Space* dup(Mesh* mesh) const;

  virtual void set_shapeset(Shapeset* shapeset);

  virtual int get_type() const { return 2; }

protected:

  virtual void assign_vertex_dofs() {}
  virtual void assign_edge_dofs();
  virtual void assign_bubble_dofs();

  virtual void get_vertex_assembly_list(Element* e, int iv, AsmList* al) {}
  virtual void get_boundary_assembly_list_internal(Element* e, int surf_num, AsmList* al);
  virtual void get_bubble_assembly_list(Element* e, AsmList* al);

  static double** hdiv_proj_mat;
  static double*  hdiv_chol_p;
  static int      hdiv_proj_ref;

  virtual scalar* get_bc_projection(SurfPos* surf_pos, int order);

  struct EdgeInfo
  {
    Node* node;
    int part;
    int ori;
    double lo, hi;
  };

  void update_constrained_nodes(Element* e, EdgeInfo* ei0, EdgeInfo* ei1, EdgeInfo* ei2, EdgeInfo* ei3);
  virtual void update_constraints();

};



#endif
