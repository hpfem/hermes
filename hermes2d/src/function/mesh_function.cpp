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

#include "mesh_function.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    MeshFunction<Scalar>::MeshFunction()
      : Function<Scalar>()
    {
      refmap = new RefMap;
      mesh = NULL;
      this->element = NULL;
    }

    template<typename Scalar>
    MeshFunction<Scalar>::MeshFunction(Mesh *mesh) :
    Function<Scalar>()
    {
      this->mesh = mesh;
      this->refmap = new RefMap;
      // FIXME - this was in H3D: MEM_CHECK(this->refmap);
      this->element = NULL;		// this comes with Transformable
    }

    template<typename Scalar>
    MeshFunction<Scalar>::~MeshFunction()
    {
      delete refmap;
      if(this->overflow_nodes != NULL)
      {
        for(unsigned int i = 0; i < this->overflow_nodes->get_size(); i++)
          if(this->overflow_nodes->present(i))
            ::free(this->overflow_nodes->get(i));
        delete this->overflow_nodes;
      }
    }

    template<typename Scalar>
    void MeshFunction<Scalar>::init()
    {
    }

    template<typename Scalar>
    void MeshFunction<Scalar>::reinit()
    {
      this->free();
      init();
    }

    template<typename Scalar>
    int MeshFunction<Scalar>::get_edge_fn_order(int edge)
    {
      return Function<Scalar>::get_edge_fn_order(edge);
    }

    template<typename Scalar>
    Mesh* MeshFunction<Scalar>::get_mesh() const
    {
      return mesh;
    }

    template<typename Scalar>
    RefMap* MeshFunction<Scalar>::get_refmap(bool update)
    {
      if(update)
        this->update_refmap();
      return refmap;
    }

    template<typename Scalar>
    void MeshFunction<Scalar>::set_refmap(RefMap* refmap_to_set)
    {
      this->refmap = refmap_to_set;
    }

    template<typename Scalar>
    void MeshFunction<Scalar>::set_quad_2d(Quad2D* quad_2d)
    {
      _F_
      if (quad_2d==NULL) throw Exceptions::NullException(1);
      Function<Scalar>::set_quad_2d(quad_2d);
      refmap->set_quad_2d(quad_2d);
    }


    template<typename Scalar>
    void MeshFunction<Scalar>::set_active_element(Element* e)
    {
      Transformable::set_active_element(e);
      mode = e->get_mode();
      refmap->set_active_element(e);
      Function<Scalar>::reset_transform();
    }

    template<typename Scalar>
    void MeshFunction<Scalar>::handle_overflow_idx()
    {
      if(this->overflow_nodes != NULL) {
        for(unsigned int i = 0; i < this->overflow_nodes->get_size(); i++)
          if(this->overflow_nodes->present(i))
            ::free(this->overflow_nodes->get(i));
        delete this->overflow_nodes;
      }
      this->nodes = new LightArray<typename Function<Scalar>::Node *>;
      this->overflow_nodes = this->nodes;
    }

    template<typename Scalar>
    void MeshFunction<Scalar>::push_transform(int son)
    {
      Transformable::push_transform(son);
      Function<Scalar>::update_nodes_ptr();
    }

    template<typename Scalar>
    void MeshFunction<Scalar>::pop_transform()
    {
      Transformable::pop_transform();
      Function<Scalar>::update_nodes_ptr();
    }

    template<typename Scalar>
    void MeshFunction<Scalar>::force_transform(MeshFunction<Scalar>* mf)
    {
      Function<Scalar>::force_transform(mf->get_transform(), mf->get_ctm());
    }

    template<typename Scalar>
    void MeshFunction<Scalar>::update_refmap()
    {
      refmap->force_transform(this->sub_idx, this->ctm);
    }

    template<typename Scalar>
    void MeshFunction<Scalar>::force_transform(uint64_t sub_idx, Trf* ctm)
    {
      this->sub_idx = sub_idx;
      this->ctm = ctm;
    }

    template class HERMES_API MeshFunction<double>;
    template class HERMES_API MeshFunction<std::complex<double> >;
  }
}
