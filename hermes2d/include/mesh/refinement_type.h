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

#ifndef __H2D_REFINEMENT_TYPE_H
#define __H2D_REFINEMENT_TYPE_H

#include "../global.h"

/// Possible refinements of an element.
enum RefinementType {
  H2D_REFINEMENT_P = -1, ///< P-refinement.
  H2D_REFINEMENT_H = 0, ///< H-refinement.
  H2D_REFINEMENT_ANISO_H = 1, ///< ANISO-refienement. The element is split along the horizontal axis. Quadrilaterals only.
  H2D_REFINEMENT_ANISO_V = 2 ///< ANISO-refienement. The element is split along the vertical axis. Quadrilaterals only.
};

/// Retuns true if a given refinement is an ANISO-refinement.
/** \param[in] refin_type A refinement type. Possible values are defined in the enum RefinementType.
 *  \return True of a given refinement is an ANISO-refinement. */
extern HERMES_API bool is_refin_aniso(const int refin_type);

/// Returns a maximum number of sons that will be generated if a given refinement is applied.
/** \param[in] refin_type A refinement type. Possible values are defined in the enum RefinementType.
 *  \return A number of possible sons. In a case of P-refinement, the function returns 1 even thought this refinement yields just a change in orders. */
extern HERMES_API int get_refin_sons(const int refin_type);

/// Returns a string representation of the refinement.
/** Used for debugging and event logging purposes.
 *  \param[in] refin_type A refinement type. Possible values are defined in the enum RefinementType.
 *  \return A string representation of a given refinement. */
extern HERMES_API const std::string get_refin_str(const int refin_type);

#endif