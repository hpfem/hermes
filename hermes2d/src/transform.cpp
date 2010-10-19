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

#include "h2d_common.h"
#include "transform.h"

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
