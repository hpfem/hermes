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

// DEPRECATED
// main functions used for error calculation in Hermes
extern HERMES_API double calc_rel_error(MeshFunction* sln1, MeshFunction* sln2, int norm_type); // Note: coarse mesh sln has to be first, then 
// DEPRECATED                                                                                             // ref_sln.
extern HERMES_API double calc_abs_error(MeshFunction* sln1, MeshFunction* sln2, int norm_type);
// DEPRECATED
extern HERMES_API double calc_norm(MeshFunction* sln, int norm_type);

// DEPRECATED
// function calculating errors between solutions in right and left vectors, returning all necessary parameters
// returns correct parameters only if the return value is true
// coarse mesh sln has to be first, then ref_sln
HERMES_API bool calc_errors(Hermes::Tuple<Solution* > left, Hermes::Tuple<Solution *> right, Hermes::Tuple<double> & err_abs, Hermes::Tuple<double> & norm_vals, 
      double & err_abs_total, double & norm_total, double & err_rel_total, Hermes::Tuple<ProjNormType> norms = Hermes::Tuple<ProjNormType>());

// helper functions
// DEPRECATED
extern HERMES_API double calc_abs_error(double (*fn)(MeshFunction*, MeshFunction*, RefMap*, RefMap*), MeshFunction* sln1, MeshFunction* sln2);
// DEPRECATED
extern HERMES_API double calc_norm(double (*fn)(MeshFunction*, RefMap*), MeshFunction* sln);

extern HERMES_API double error_fn_l2(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv);
extern HERMES_API double norm_fn_l2(MeshFunction* sln, RefMap* ru);
//extern HERMES_API double l2_error(MeshFunction* sln1, MeshFunction* sln2);
//extern HERMES_API double l2_norm(MeshFunction* sln);

extern HERMES_API double error_fn_h1(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv);
extern HERMES_API double norm_fn_h1(MeshFunction* sln, RefMap* ru);
//extern HERMES_API double h1_error(MeshFunction* sln1, MeshFunction* sln2);
//extern HERMES_API double h1_norm(MeshFunction* sln);

extern HERMES_API double error_fn_hc(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv);
extern HERMES_API double norm_fn_hc(MeshFunction* sln, RefMap* ru);
//extern HERMES_API double hcurl_error(MeshFunction* sln1, MeshFunction* sln2);
//extern HERMES_API double hcurl_norm(MeshFunction* sln);

extern HERMES_API double error_fn_hcl2(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv);
extern HERMES_API double norm_fn_hcl2(MeshFunction* sln, RefMap* ru);
extern HERMES_API double hcurl_l2error(MeshFunction* sln1, MeshFunction* sln2);
extern HERMES_API double hcurl_l2norm(MeshFunction* sln);

extern HERMES_API double error_fn_hdiv(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv);
extern HERMES_API double norm_fn_hdiv(MeshFunction* sln, RefMap* ru);

#endif
