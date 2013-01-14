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
		/// @ingroup spaces
		/// L2Space represents a space of Scalar functions with discontinuities along
		/// mesh edges.
		template<typename Scalar>
		class HERMES_API L2Space : public Space<Scalar>
		{
		public:
			L2Space();
			L2Space(const Mesh* mesh, int p_init = 0,
				Shapeset* shapeset = NULL);

			virtual ~L2Space();

			virtual Scalar* get_bc_projection(SurfPos* surf_pos, int order, EssentialBoundaryCondition<Scalar> *bc);

			virtual void get_element_assembly_list(Element* e, AsmList<Scalar>* al, unsigned int first_dof = 0) const;

			/// Copy from Space instance 'space'
			virtual void copy(const Space<Scalar>* space, Mesh* new_mesh);
		protected:
			virtual int get_edge_order(Element* e, int edge) const {
				return H2D_MAKE_EDGE_ORDER(e->get_mode(), edge, this->edata[e->id].order);
			}

			virtual void set_shapeset(Shapeset* shapeset);

			virtual SpaceType get_type() const { return HERMES_L2_SPACE; }

			/// Common code for the constructors.
			void init(Shapeset* shapeset, int p_init);

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

			virtual void get_vertex_assembly_list(Element* e, int iv, AsmList<Scalar>* al) const;
			virtual void get_boundary_assembly_list_internal(Element* e, int surf_num, AsmList<Scalar>* al) const;
			virtual void get_bubble_assembly_list(Element* e, AsmList<Scalar>* al) const;
			template<typename T> friend class Space<T>::ReferenceSpaceCreator;
			friend class Space<Scalar>;
		};
	}
}
#endif
