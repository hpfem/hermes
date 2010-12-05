// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef __BCVALUES_H
#define __BCVALUES_H

#include "common.h"
#include "tuple.h"
#include "error.h"
#include "bctypes.h"

/// Class to link boundary markers where Dirichlet BC are prescribed with 
/// the function representing them.
class HERMES_API BCValues {
public:
  BCValues() {};
  ~BCValues() {};

protected:
  /// Type of the function representing the BC.
  typedef scalar (*value_callback)(double, double);

  /// Storage of the functions. 
  std::map<int, value_callback> value_callbacks;

  /// Default (zero) function.
  static scalar default_callback(double x, double y) { return 0; };

public:
  /// Adds the function callback to represent the Dirichlet BC on the parts of
  /// the boundary marked with markers.
  void add_function(value_callback callback, Hermes::Tuple<int> markers) 
  {
    if(markers.size() == 0)
      error("BCValues::add_function() called without any boundary markers specified.");
    for(int i = 0; i < markers.size(); i++)
      /// If we find out that there is already a function present describing the Dirichlet
      /// BC on this part of the boundary, return an error, otherwise store the function.
      if(!(value_callbacks.insert(std::pair<int, value_callback>(markers[i], callback))).second)
        error("Attempt to define more than one function representing the Dirichlet BC \
               on the same part of the boundary.");
  };

  /// The same as add_function(), only supplies the default (zero) functions.
  void add_zero_function(Hermes::Tuple<int> markers) 
  {
    if(markers.size() == 0)
      error("BCValues::add_zero_function() called without any boundary markers specified.");
    for(int i = 0; i < markers.size(); i++) {
      if(value_callbacks[markers[i]] != NULL)
        error("Attempt to define more than one function representing the Dirichlet BC \
              on the same part of the boundary.");
      value_callbacks.insert(std::pair<int, value_callback>(markers[i], BCValues::default_callback));
    }
  };

  /// Checks whether all functions representing Dirichlet BCs are really set on a Dirichlet
  /// part of the boundary.
  void check_consistency(BCTypes* bc_types) 
  {
    std::map<int, value_callback>::iterator it;
    for(it = value_callbacks.begin(); it != value_callbacks.end(); it++)
      if(!bc_types->is_essential(it->first))
        error("BCTypes and BCValues incompatible.");
  };

  /// According to the bc_types variable, this function supplies the default function.
  void update(BCTypes* bc_types) 
  {
    std::vector<int>::iterator it;
    for(it = bc_types->markers_dirichlet.begin(); it != bc_types->markers_dirichlet.end(); it++)
      if(value_callbacks[*it] == NULL)
        add_zero_function(*it);
  };

  /// Main function that returns the value.
  scalar calculate(int marker, double x, double y)
  {
    if(value_callbacks[marker] == NULL)
      error("Attempt to retrieve a value of a function representing the Dirichlet BC without \
            this being set up for the current Space.");
    return value_callbacks[marker](x, y);
  }
};

#endif
