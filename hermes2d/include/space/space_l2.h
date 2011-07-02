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

#ifndef __H2D_SPACE_L2
#define __H2D_SPACE_L2

#include "space.h"
namespace Hermes
{
  namespace Hermes2D
  {
    /// L2Space represents a space of Scalar functions with discontinuities along
    /// mesh edges.
    ///
    ///
    template<typename Scalar>
    class HERMES_API L2Space : public Space<Scalar>
    {
    public:
      L2Space(Mesh* mesh, int p_init = 0,
        Shapeset* shapeset = NULL);

      /// Common code for the constructors.
      void init(Shapeset* shapeset, Ord2 p_init);

      virtual ~L2Space();

      virtual Space<Scalar>* dup(Mesh* mesh, int order_increase = 0) const;

      virtual int get_edge_order(Element* e, int edge) {
        // There are no continuity constraints on shape functions in L2.
        return Global<Scalar>::make_edge_order( e->get_mode(), edge, this->edata[e->id].order );
      }

      virtual void set_shapeset(Shapeset* shapeset);

      virtual SpaceType get_type() const { return HERMES_L2_SPACE; }

      virtual void get_element_assembly_list(Element* e, AsmList<Scalar>* al);
      
      // FIXME: This function should probably not be used at all.
      virtual Scalar* get_bc_projection(SurfPos* surf_pos, int order);

    protected:

      struct L2Data
      {
        int vdof[4];
        int edof[4];
      };

      L2Data* ldata;
      int lsize;

      virtual void resize_tables();

      virtual void assign_vertex_dofs() {}
      virtual void assign_edge_dofs() {}
      virtual void assign_bubble_dofs();

      virtual void get_vertex_assembly_list(Element* e, int iv, AsmList<Scalar>* al);
      virtual void get_boundary_assembly_list_internal(Element* e, int surf_num, AsmList<Scalar>* al);
      virtual void get_bubble_assembly_list(Element* e, AsmList<Scalar>* al);
    };
  }
}
#endif
