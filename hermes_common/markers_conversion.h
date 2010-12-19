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


#ifndef __MARKERS_CONVERSION_H
#define __MARKERS_CONVERSION_H

#include "error.h"
#include "compat.h"
#include <string>
#include <map>

class HERMES_API MarkersConversion
{
public:
  MarkersConversion();
  ~MarkersConversion();

  // Info about the maximum markers used so far, used in determining
  // of the internal marker for a user-supplied std::string identification for 
  // the purpose of disambiguity.
  // 
  int min_boundary_marker_unused;
  int min_element_marker_unused;
  
  // Function inserting a marker into conversion_table_for_element_markers.
  // This function controls if this user_marker x internal_marker is already 
  // present, and if not, it inserts the std::pair.
  void insert_element_marker(int internal_marker, std::string user_marker);
  // An analogy for boundary markers.
  void insert_boundary_marker(int internal_marker, std::string user_marker);

  // Lookup functions.
  // Find a user marker for this internal marker.
  std::string get_user_element_marker(int internal_marker);
  // An analogy for boundary markers.
  std::string get_user_boundary_marker(int internal_marker);

  // Find an internal marker for this user_marker.
  int get_internal_element_marker(std::string user_marker);
  // An analogy for boundary markers.
  int get_internal_boundary_marker(std::string user_marker);

  // Make sure that the internal markers do not collide.
  // This function is called whenever user supplies his own integral label.
  // This is done with respect to the task to preserve user-supplied integral markers.
  void check_boundary_marker(int marker);
  // An analogy for element markers.
  void check_element_marker(int marker);

private:
  // Conversion tables between the std::string markers the user sets and
  // the markers used internally as members of Elements, Nodes.
  std::map<int, std::string>* conversion_table_for_element_markers;
  std::map<int, std::string>* conversion_table_for_boundary_markers;

  // Inverse tables, so that it is possible to search using either
  // the internal representation, or the user std::string value.
  std::map<std::string, int>* conversion_table_for_element_markers_inverse;
  std::map<std::string, int>* conversion_table_for_boundary_markers_inverse;

};

#endif

