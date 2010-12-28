// This file is part of Hermes3D.
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D. If not, see <http://www.gnu.org/licenses/>.

#ifndef __H3D_COMMON_H_
#define __H3D_COMMON_H_

#include "../../hermes_common/common.h"

// H3D-specific error codes.
#define H3D_ERR_FACE_INDEX_OUT_OF_RANGE         "Face index out of range."
#define H3D_ERR_EDGE_INDEX_OUT_OF_RANGE         "Edge index out of range."
#define H3D_ERR_TETRA_NOT_COMPILED              "hermes3d was not built with tetra elements."
#define H3D_ERR_HEX_NOT_COMPILED                "hermes3d was not built with hex elements."
#define H3D_ERR_PRISM_NOT_COMPILED              "hermes3d was not built with prism elements."

// Maximum polynomial order of elements.
#define H3D_MAX_ELEMENT_ORDER							10

#endif
