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

template<typename T>
class Range {       ///< Range of values.
  T lower_bound;    ///< Lower boundary.
  T upper_bound;    ///< Upper boundary.
  bool empty_range; ///< True if range is empty.
public:
  Range() : empty_range(true) {};
  Range(const T& lower_bound, const T& upper_bound) : lower_bound(lower_bound), upper_bound(upper_bound), empty_range(lower_bound > upper_bound) {};
  bool empty() const { return empty_range; }; ///< Returns true if range is empty.
  const T& lower() const { return lower_bound; }; ///< Returns lower bound.
  const T& upper() const { return upper_bound; }; ///< Returns upper bound.
  bool is_in_closed(const Range<T>& range) const { return (range.lower_bound >= lower_bound && range.upper_bound <= upper_bound); }; ///< Returns true if value is inside the closed range.
  bool is_in_closed(const T& value) const { return (value >= lower_bound && value <= upper_bound); }; ///< Returns true if value is inside the closed range.
  bool is_in_open(const T& value) const { return (value > lower_bound && value < upper_bound); }; ///< Returns true if value is inside the open range.
  void enlarge_to_include(const T& value) { ///< Enlarge range to include a value.
    if (empty_range) {
      lower_bound = upper_bound = value;
      empty_range = false;
    }
    else {
      if (lower_bound > value)
        lower_bound = value;
      if (upper_bound < value)
        upper_bound = value;
    }
  };

  static Range<T> make_envelope(const Range<T>& a, const Range<T>& b) { ///< Create an envelope which contains both ranges.
    if (a.empty())
      return b;
    else if (b.empty())
      return a;
    else
      return Range(std::min(a.lower(), b.lower()), std::max(a.upper(), b.upper()));
  };
};

#endif
