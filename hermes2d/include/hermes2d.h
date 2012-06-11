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
// along with Hermes2D. If not, see <http://www.gnu.org/licenses/>.

// $Id$

#ifndef __HERMES_2D_H
#define __HERMES_2D_H

#include "global.h"

#include "newton_solver.h"
#include "picard_solver.h"
#include "linear_solver.h"
#include "calculation_continuity.h"

#include "boundary_conditions/essential_boundary_conditions.h"

#include "quadrature/limit_order.h"

#include "api2d.h"

#include "mesh/mesh.h"
#include "mesh/mesh_reader.h"
#include "mesh/mesh_reader_h2d.h"
#include "mesh/mesh_reader_h2d_xml.h"
#include "mesh/mesh_reader_h1d_xml.h"
#include "mesh/mesh_reader_exodusii.h"

#include "quadrature/quad.h"
#include "quadrature/quad_all.h"

#include "space/space_h1.h"
#include "space/space_hcurl.h"
#include "space/space_l2.h"
#include "space/space_hdiv.h"

#include "shapeset/shapeset_h1_all.h"
#include "shapeset/shapeset_hc_all.h"
#include "shapeset/shapeset_hd_all.h"
#include "shapeset/shapeset_l2_all.h"

#include "mesh/refmap.h"
#include "mesh/traverse.h"

#include "weakform/weakform.h"
#include "discrete_problem.h"
#include "discrete_problem_linear.h"
#include "forms.h"

#include "integrals/h1.h"
#include "integrals/hcurl.h"
#include "integrals/hdiv.h"
#include "integrals/l2.h"

#include "function/exact_solution.h"
#include "function/solution.h"
#include "function/mesh_function.h"
#include "function/filter.h"

#include "graph.h"

#include "views/view.h"
#include "views/base_view.h"
#include "views/mesh_view.h"
#include "views/order_view.h"
#include "views/scalar_view.h"
#include "views/stream_view.h"
#include "views/vector_base_view.h"
#include "views/vector_view.h"

#include "mesh/refinement_type.h"
#include "mesh/element_to_refine.h"
#include "refinement_selectors/selector.h"
#include "refinement_selectors/order_permutator.h"
#include "refinement_selectors/optimum_selector.h"
#include "refinement_selectors/proj_based_selector.h"
#include "refinement_selectors/l2_proj_based_selector.h"
#include "refinement_selectors/h1_proj_based_selector.h"
#include "refinement_selectors/hcurl_proj_based_selector.h"

#include "adapt/adapt.h"
#include "adapt/kelly_type_adapt.h"
#include "neighbor.h"
#include "projections/localprojection.h"
#include "projections/ogprojection.h"
#include "projections/ogprojection_nox.h"

#include "runge_kutta.h"
#include "spline.h"

#include "weakform_library/weakforms_elasticity.h"
#include "weakform_library/weakforms_h1.h"
#include "weakform_library/weakforms_hcurl.h"
#include "weakform_library/weakforms_maxwell.h"
#include "weakform_library/weakforms_neutronics.h"

#include "doxygen_first_page.h"
#endif
