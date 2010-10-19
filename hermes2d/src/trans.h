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

#ifndef __H2D_TRANS_H
#define __H2D_TRANS_H

#include "h2d_common.h"
#include "mesh.h"
#include "shapeset/shapeset.h"
#include "shapeset/shapeset_h1_all.h"
#include "quad.h"
#include "quad_all.h"

extern HERMES_API double2 *transform(Element *e);
extern HERMES_API void element_polygonal_boundary(Element *e, double2 **tp, int *n);

#endif
