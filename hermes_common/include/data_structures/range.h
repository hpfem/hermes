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
/*! \file range.h
\brief File containing the class Range.
*/
#ifndef __HERMES_COMMON_RANGE_H
#define __HERMES_COMMON_RANGE_H

#include "compat.h"

namespace Hermes
{
  /// Range of values.
  class HERMES_API Range
  {
  public:
    Range();
    Range(const int& lower_bound, const int& upper_bound);
    /// True if range is empty.
    bool empty() const;
    /// Lower boundary.
    const int& lower() const;
    /// Upper boundary.
    const int& upper() const;
    bool is_in_closed(const Range& range) const;
    bool is_in_closed(const int& value) const;
    bool is_in_open(const int& value) const;
    void enlarge_to_include(const int& value);

    static Range make_envelope(const Range& a, const Range& b);
  protected:
    int lower_bound;
    int upper_bound;
    bool empty_range;
  };
}
#endif