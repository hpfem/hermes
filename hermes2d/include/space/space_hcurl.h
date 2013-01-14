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
namespace Hermes
{
	namespace Hermes2D
	{
		/// @ingroup spaces
		/// HcurlSpace represents a space of vector functions with continuous tangent
		/// components over a domain (mesh).
		template<typename Scalar>
		class HERMES_API HcurlSpace : public Space<Scalar>
		{
		public:
			HcurlSpace();
			HcurlSpace(const Mesh* mesh, EssentialBCs<Scalar>* boundary_conditions, int p_init = 1,
				Shapeset* shapeset = NULL);

			HcurlSpace(const Mesh* mesh, int p_init = 1,
				Shapeset* shapeset = NULL);

			virtual ~HcurlSpace();

			virtual void set_shapeset(Shapeset* shapeset);

			virtual Scalar* get_bc_projection(SurfPos* surf_pos, int order, EssentialBoundaryCondition<Scalar> *bc);

			/// Copy from Space instance 'space'
			virtual void copy(const Space<Scalar>* space, Mesh* new_mesh);
		protected:

			virtual SpaceType get_type() const { return HERMES_HCURL_SPACE; }

			/// Common code for the constructors.
			void init(Shapeset* shapeset, int p_init);

			virtual void assign_vertex_dofs() {}
			virtual void assign_edge_dofs();
			virtual void assign_bubble_dofs();

			virtual void get_vertex_assembly_list(Element* e, int iv, AsmList<Scalar>* al) const {}
			virtual void get_boundary_assembly_list_internal(Element* e, int surf_num, AsmList<Scalar>* al) const;

			struct EdgeInfo
			{
				Node* node;
				int part;
				int ori;
				double lo, hi;
			};

			void update_constrained_nodes(Element* e, EdgeInfo* ei0, EdgeInfo* ei1, EdgeInfo* ei2, EdgeInfo* ei3);
			virtual void update_constraints();
			friend class Space<Scalar>;
		};
	}
}
#endif
