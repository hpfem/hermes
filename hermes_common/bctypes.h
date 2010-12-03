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

#ifndef __BCTYPES_H
#define __BCTYPES_H

#include "common.h"
#include "tuple.h"
#include "error.h"

// Types of boundary conditions.
enum BCType
{
  BC_ESSENTIAL, ///< Essential (Dirichlet) BC.
  BC_NATURAL,   ///< Natural (Neumann, Newton) BC.
  BC_NONE       ///< Hermes will not attempt to evaluate any boundary
                ///< integrals on this part of the boundary.
};

/// Class to link boundary markers with boundary condition types.
class HERMES_API BCTypes {
public:
 
  void add_bc_natural(Tuple<int> markers) 
  {
    int n = markers.size();
    if (n <= 0) error("BCTypes::add_bc_natural() expects at least one marker.");
    for (int i = 0; i < n; i++) this->markers_natural.push_back(markers[i]);
    return;
  };

  void add_bc_essential(Tuple<int> markers) 
  {
    int n = markers.size();
    if (n <= 0) error("BCTypes::add_bc_essential() expects at least one marker.");
    for (int i = 0; i < n; i++) this->markers_essential.push_back(markers[i]);    
    return;
  };

  void add_bc_none(Tuple<int> markers) 
  {
    int n = markers.size();
    if (n <= 0) error("BCTypes::add_bc_none() expects at least one marker.");
    for (int i = 0; i < n; i++) this->markers_none.push_back(markers[i]);    
    return;
  };

  bool is_natural(int marker) {
    int n = this->markers_natural.size();
    for (int i = 0; i < n; i++) if (this->markers_natural[i] == marker) return true;
    return false;
  }

  bool is_essential(int marker) {
    int n = this->markers_essential.size();
    for (int i = 0; i < n; i++) if (this->markers_essential[i] == marker) return true;
    return false;
  }

  bool is_none(int marker) {
    int n = this->markers_none.size();
    for (int i = 0; i < n; i++) if (this->markers_none[i] == marker) return true;
    return false;
  }

  // default (if not found) is BC_NATURAL
  BCType get_type(int marker) {
    if (this->is_essential(marker)) return BC_ESSENTIAL;
    if (this->is_none(marker)) return BC_NONE;
    return BC_NATURAL;
  }

  BCTypes() {};
  ~BCTypes() {};

  virtual BCTypes *dup() {
      BCTypes *bc = new BCTypes();
      bc->markers_natural = this->markers_natural;
      bc->markers_essential = this->markers_essential;
      bc->markers_none = this->markers_none;
      return bc;
  }

  protected:
    Tuple<int> markers_natural;
    Tuple<int> markers_essential;
    Tuple<int> markers_none;
};

static BCType default_bc_type(int marker)
{
  return BC_NATURAL;
}

class HERMES_API BCTypesCallback: public BCTypes {
public:

  BCTypesCallback():BCTypes() {
      this->bc_type_callback = NULL;
  };

  void register_callback(BCType (*bc_type_callback)(int)) {
      if (bc_type_callback == NULL) bc_type_callback = default_bc_type;
      this->bc_type_callback = bc_type_callback;
  }

  // default (if not found) is BC_NATURAL
  BCType get_type(int marker) {
    if (this->bc_type_callback == NULL) error("No callback was registered");
    return this->bc_type_callback(marker);
  }

  virtual BCTypes *dup() {
      BCTypesCallback *bc = new BCTypesCallback();
      bc->bc_type_callback = this->bc_type_callback;
      return bc;
  }

private:
    BCType (*bc_type_callback)(int);
};

#endif
