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

#include "function.h"
#include "../mesh/element.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    int Function<Scalar>::idx2mask[H2D_NUM_FUNCTION_VALUES][2] =
    {
      { H2D_FN_VAL_0, H2D_FN_VAL_1 }, { H2D_FN_DX_0, H2D_FN_DX_1 }, { H2D_FN_DY_0, H2D_FN_DY_1 },
#ifdef H2D_USE_SECOND_DERIVATIVES
      { H2D_FN_DXX_0, H2D_FN_DXX_1 }, { H2D_FN_DYY_0, H2D_FN_DYY_1 }, { H2D_FN_DXY_0, H2D_FN_DXY_1 }
#endif
    };

    template<typename Scalar>
    Function<Scalar>::Function()
      : Transformable()
    {
      order = -1;
      memset(quads, 0, H2D_MAX_QUADRATURES*sizeof(Quad2D*));
      this->values_valid = false;
    }

    template<typename Scalar>
    Function<Scalar>::~Function()
    {
    }

    template<typename Scalar>
    int Function<Scalar>::get_fn_order() const
    {
      return order;
    }

    template<typename Scalar>
    int Function<Scalar>::get_edge_fn_order(int edge) const
    {
      return order;
    }

    template<typename Scalar>
    void Function<Scalar>::set_quad_order(unsigned int order, int mask)
    {
      precalculate(order, mask);
      this->order = order;
    }

    template<typename Scalar>
    void Function<Scalar>::set_active_element(Element* e)
    {
      Transformable::set_active_element(e);
      this->invalidate_values();
    }

    template<typename Scalar>
    const Scalar* Function<Scalar>::get_values(int a, int b) const
    {
      return values[a][b];
    }

    template<typename Scalar>
    void Function<Scalar>::set_quad_2d(Quad2D* quad_2d)
    {
      int i;
      this->invalidate_values();

      // check to see if we already have the quadrature
      for (i = 0; i < H2D_MAX_QUADRATURES; i++)
      {
        if (quads[i] == quad_2d)
        {
          cur_quad = i;
          return;
        }
      }

      // if not, add the quadrature to a free slot
      for (i = 0; i < H2D_MAX_QUADRATURES; i++)
      {
        if (quads[i] == nullptr)
        {
          quads[i] = quad_2d;
          cur_quad = i;
          return;
        }
      }

      throw Hermes::Exceptions::Exception("too many quadratures.");
    }

    template<typename Scalar>
    void Function<Scalar>::precalculate(int order, int mask)
    {
      this->values_valid = true;
    }

    template<typename Scalar>
    void Function<Scalar>::push_transform(int son)
    {
      Transformable::push_transform(son);
      this->invalidate_values();
    }

    template<typename Scalar>
    void Function<Scalar>::pop_transform()
    {
      Transformable::pop_transform();
      this->invalidate_values();
    }

    template<typename Scalar>
    Quad2D* Function<Scalar>::get_quad_2d() const
    {
      return quads[cur_quad];
    }

    template<typename Scalar>
    void Function<Scalar>::invalidate_values()
    {
      this->values_valid = false;
    }

    template<typename Scalar>
    void Function<Scalar>::force_transform(uint64_t sub_idx, Trf* ctm)
    {
      Transformable::force_transform(sub_idx, ctm);
      this->invalidate_values();
    }

    template<typename Scalar>
    void Function<Scalar>::reset_transform()
    {
      Transformable::reset_transform();
      this->invalidate_values();
    }

    template<typename Scalar>
    void Function<Scalar>::set_transform(uint64_t idx)
    {
      Transformable::set_transform(idx);
      this->invalidate_values();
    }

    template<typename Scalar>
    int Function<Scalar>::get_num_components() const
    {
      return num_components;
    }

    template<typename Scalar>
    const Scalar* Function<Scalar>::get_fn_values(int component) const
    {
      assert(this->values_valid);
      return &values[component][0][0];
    }

    template<typename Scalar>
    const Scalar* Function<Scalar>::get_dx_values(int component) const
    {
      assert(this->values_valid);
      return &values[component][1][0];
    }

    template<typename Scalar>
    const Scalar* Function<Scalar>::get_dy_values(int component) const
    {
      assert(this->values_valid);
      return &values[component][2][0];
    }

#ifdef H2D_USE_SECOND_DERIVATIVES
    template<typename Scalar>
    const Scalar* Function<Scalar>::get_dxx_values(int component) const
    {
      assert(this->values_valid);
      return &values[component][3][0];
    }

    template<typename Scalar>
    const Scalar* Function<Scalar>::get_dyy_values(int component) const
    {
      assert(this->values_valid);
      return &values[component][4][0];
    }

    template<typename Scalar>
    const Scalar* Function<Scalar>::get_dxy_values(int component) const
    {
      assert(this->values_valid);
      return &values[component][5][0];
    }
#endif

    template<typename Scalar>
    Scalar* Function<Scalar>::deep_copy_array(int component, int item) const
    {
      assert(this->values_valid);
      int np = this->quads[this->cur_quad]->get_num_points(this->order, this->element->get_mode());
      Scalar* toReturn = malloc_with_check<Scalar>(np);
      memcpy(toReturn, this->values[component][item], sizeof(Scalar)* np);
      return toReturn;
    }

    template class HERMES_API Function<double>;
    template class HERMES_API Function<std::complex<double> >;
  }
}
