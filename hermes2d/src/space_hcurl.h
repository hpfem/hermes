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

#ifndef __H2D_SPACE_HCURL_H
#define __H2D_SPACE_HCURL_H

#include "space.h"


/// HcurlSpace represents a space of vector functions with continuous tangent
/// components over a domain (mesh).
///
///
///
class H2D_API HcurlSpace : public Space
{
public:

  HcurlSpace(Mesh* mesh = NULL, BCType (*bc_type_callback)(int) = NULL, 
                 scalar (*bc_value_callback_by_coord)(int, double, double) = NULL, int p_init = 1, 
                 Shapeset* shapeset = NULL);
  virtual ~HcurlSpace();

  virtual Space* dup(Mesh* mesh) const;

  virtual int get_type() const { return 1; }

  /// Sets element polynomial order and calls assign_dofs(). Intended for the user.
  virtual void set_element_order(int id, int order);

  /// Sets element polynomial order without calling assign_dofs(). For internal use.
  virtual void set_element_order_internal(int id, int order);


protected:

  virtual void assign_vertex_dofs() {}
  virtual void assign_edge_dofs();
  virtual void assign_bubble_dofs();

  virtual void get_vertex_assembly_list(Element* e, int iv, AsmList* al) {}
  virtual void get_edge_assembly_list_internal(Element* e, int ie, AsmList* al);

  static double** hcurl_proj_mat;
  static double*  hcurl_chol_p;
  static int      hcurl_proj_ref;

  virtual scalar* get_bc_projection(EdgePos* ep, int order);

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
