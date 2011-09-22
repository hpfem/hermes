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

#ifndef __H2D_RANGE_H
#define __H2D_RANGE_H

#include "compat.h"

namespace Hermes
{
  /// Range of values.
  template<typename T>
  class HERMES_API Range 
  {       
    T lower_bound;    ///< Lower boundary.
    T upper_bound;    ///< Upper boundary.
    bool empty_range; ///< True if range is empty.
  public:
    Range();
    Range(const T& lower_bound, const T& upper_bound);
    bool empty() const;
    const T& lower() const;
    const T& upper() const;
    bool is_in_closed(const Range<T>& range) const;
    bool is_in_closed(const T& value) const;
    bool is_in_open(const T& value) const;
    void enlarge_to_include(const T& value);

    static Range<T> make_envelope(const Range<T>& a, const Range<T>& b);
  };
}
#endif
