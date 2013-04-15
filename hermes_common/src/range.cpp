// This file is part of HermesCommon
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes2D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
/*! \file range.cpp
\brief Range implementation.
*/
#include "range.h"
#include <algorithm>
namespace Hermes
{
  Range::Range() : empty_range(true) {}

  Range::Range(const int& lower_bound, const int& upper_bound) : lower_bound(lower_bound), upper_bound(upper_bound), empty_range(lower_bound > upper_bound)
  {
  }

  bool Range::empty() const
  {
    return empty_range;
  }

  const int& Range::lower() const
  {
    return lower_bound;
  }

  const int& Range::upper() const
  {
    return upper_bound;
  }

  bool Range::is_in_closed(const Range& range) const
  {
    return (range.lower_bound >= lower_bound && range.upper_bound <= upper_bound);
  }

  bool Range::is_in_closed(const int& value) const
  {
    return (value >= lower_bound && value <= upper_bound);
  }

  bool Range::is_in_open(const int& value) const
  {
    return (value > lower_bound && value < upper_bound);
  }

  void Range::enlarge_to_include(const int& value)
  {
    if(empty_range) {
      lower_bound = upper_bound = value;
      empty_range = false;
    }
    else {
      if(lower_bound > value)
        lower_bound = value;
      if(upper_bound < value)
        upper_bound = value;
    }
  }

  Range Range::make_envelope(const Range& a, const Range& b)
  {
    if(a.empty())
      return b;
    else if(b.empty())
      return a;
    else
      return Range(std::min(a.lower(), b.lower()), std::max(a.upper(), b.upper()));
  }
}