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

#ifndef __H2D_NORM_H
#define __H2D_NORM_H

#include "solution.h"
#include "refmap.h"

// main functions used for error calculation in Hermes
extern H2D_API double calc_rel_error(MeshFunction* sln1, MeshFunction* sln2, int norm_type); // Note: coarse mesh sln has to be first, then 
                                                                                             // ref_sln.
extern H2D_API double calc_abs_error(MeshFunction* sln1, MeshFunction* sln2, int norm_type);
extern H2D_API double calc_norm(MeshFunction* sln, int norm_type);

// helper functions
extern H2D_API double calc_abs_error(double (*fn)(MeshFunction*, MeshFunction*, RefMap*, RefMap*), MeshFunction* sln1, MeshFunction* sln2);
extern H2D_API double calc_norm(double (*fn)(MeshFunction*, RefMap*), MeshFunction* sln);

extern H2D_API double error_fn_l2(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv);
extern H2D_API double norm_fn_l2(MeshFunction* sln, RefMap* ru);
//extern H2D_API double l2_error(MeshFunction* sln1, MeshFunction* sln2);
//extern H2D_API double l2_norm(MeshFunction* sln);

extern H2D_API double error_fn_h1(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv);
extern H2D_API double norm_fn_h1(MeshFunction* sln, RefMap* ru);
//extern H2D_API double h1_error(MeshFunction* sln1, MeshFunction* sln2);
//extern H2D_API double h1_norm(MeshFunction* sln);

extern H2D_API double error_fn_hc(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv);
extern H2D_API double norm_fn_hc(MeshFunction* sln, RefMap* ru);
//extern H2D_API double hcurl_error(MeshFunction* sln1, MeshFunction* sln2);
//extern H2D_API double hcurl_norm(MeshFunction* sln);

extern H2D_API double error_fn_hcl2(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv);
extern H2D_API double norm_fn_hcl2(MeshFunction* sln, RefMap* ru);
extern H2D_API double hcurl_l2error(MeshFunction* sln1, MeshFunction* sln2);
extern H2D_API double hcurl_l2norm(MeshFunction* sln);

extern H2D_API double error_fn_hdiv(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv);
extern H2D_API double norm_fn_hdiv(MeshFunction* sln, RefMap* ru);

#endif
