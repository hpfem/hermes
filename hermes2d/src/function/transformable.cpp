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

#include "transformable.h"
#include "mesh.h"
namespace Hermes
{
  namespace Hermes2D
  {
#define H2D_IDENTIFY_TRF { 1.0,  1.0 }, { 0.0, 0.0 } ///< Identity transformation.

    Trf tri_trf[H2D_TRF_NUM] =
    {
      { { 0.5,  0.5 }, { -0.5, -0.5 } }, // son 0
      { { 0.5,  0.5 }, {  0.5, -0.5 } }, // son 1
      { { 0.5,  0.5 }, { -0.5,  0.5 } }, // son 2
      { {-0.5, -0.5 }, { -0.5, -0.5 } }, // son 3
      { H2D_IDENTIFY_TRF },              // identity
      { H2D_IDENTIFY_TRF },              // identity
      { H2D_IDENTIFY_TRF },              // identity
      { H2D_IDENTIFY_TRF },              // identity
      { H2D_IDENTIFY_TRF }               // identity
    };

    Trf quad_trf[H2D_TRF_NUM] =
    {
      { { 0.5, 0.5 }, { -0.5, -0.5 } }, // son 0
      { { 0.5, 0.5 }, {  0.5, -0.5 } }, // son 1
      { { 0.5, 0.5 }, {  0.5,  0.5 } }, // son 2
      { { 0.5, 0.5 }, { -0.5,  0.5 } }, // son 3
      { { 1.0, 0.5 }, {  0.0, -0.5 } }, // horz son 0
      { { 1.0, 0.5 }, {  0.0,  0.5 } }, // horz son 1
      { { 0.5, 1.0 }, { -0.5,  0.0 } }, // vert son 0
      { { 0.5, 1.0 }, {  0.5,  0.0 } }, // vert son 1
      { H2D_IDENTIFY_TRF }              // identity
    };

    Transformable::Transformable()
    {
      memset(stack, 0, sizeof(stack));
      reset_transform();
      element = NULL;
    }

    Transformable::~Transformable() {}

    Element* Transformable::get_active_element() const { return element; }

    void Transformable::pop_transform()
    {
      assert(top > 0);
      ctm = stack + (--top);
      sub_idx = (sub_idx - 1) >> 3;
    }

    uint64_t Transformable::get_transform() const { return sub_idx; }

    void Transformable::reset_transform()
    {
      stack[0].m[0] = stack[0].m[1] = 1.0;
      stack[0].t[0] = stack[0].t[1] = 0.0;
      ctm = stack;
      sub_idx = top = 0;
    }

    void Transformable::set_active_element(Element* e)
    {
      if(e==NULL) throw Exceptions::NullException(1);
      element = e;
      this->reset_transform();
    }

    void Transformable::set_transform(uint64_t idx)
    {
      int son[25];
      int i = 0;
      while (idx > 0)
      {
        son[i++] = (idx - 1) & 7;
        idx = (idx - 1) >> 3;
      }
      reset_transform();
      for (int k = i-1; k >= 0; k--)
        push_transform(son[k]);
    }

    void Transformable::push_transform(int son)
    {
      assert(element != NULL);
      if(top >= H2D_MAX_TRN_LEVEL)
        throw Hermes::Exceptions::Exception("Too deep transform.");

      Trf* mat = stack + (++top);
      Trf* tr = (element->is_triangle() ? tri_trf + son : quad_trf + son);

      mat->m[0] = ctm->m[0] * tr->m[0];
      mat->m[1] = ctm->m[1] * tr->m[1];
      mat->t[0] = ctm->m[0] * tr->t[0] + ctm->t[0];
      mat->t[1] = ctm->m[1] * tr->t[1] + ctm->t[1];

      ctm = mat;
      sub_idx = (sub_idx << 3) + son + 1; // see traverse.cpp if this changes
    }

    void Transformable::push_transforms(std::set<Transformable *>& transformables, int son)
    {
      for(std::set<Transformable *>::iterator it = transformables.begin(); it != transformables.end(); ++it)
        if(*it != NULL)
          (*it)->push_transform(son);
    }
    void Transformable::pop_transforms(std::set<Transformable *>& transformables)
    {
      for(std::set<Transformable *>::iterator it = transformables.begin(); it != transformables.end(); ++it)
        if(*it != NULL)
          (*it)->pop_transform();
    }
  }
}