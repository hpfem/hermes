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
#include "vector.h"
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
 
  void add_bc_dirichlet(Hermes::vector<int> markers) 
  {
    unsigned int n = markers.size();
    if (n <= 0) error("BCTypes::add_bc_dirichlet() expects at least one marker.");
    for (unsigned int i = 0; i < n; i++) this->markers_dirichlet.push_back(markers[i]);    
    return;
  };

  // A wrapper utilizing the Mesh::MarkersConversion class.
  void add_bc_dirichlet(Hermes::vector<std::string> markers)
  {
    unsigned int n = markers.size();
    if (n <= 0) error("BCTypes::add_bc_dirichlet() expects at least one marker.");
    for (unsigned int i = 0; i < n; i++)
      this->markers_dirichlet_string_temp.push_back(markers[i]);
    return;
  }

  void add_bc_neumann(Hermes::vector<int> markers) 
  {
    unsigned int n = markers.size();
    if (n <= 0) error("BCTypes::add_bc_neumann() expects at least one marker.");
    for (unsigned int i = 0; i < n; i++)
        this->markers_neumann.push_back(markers[i]);
    return;
  };

  // A wrapper utilizing the Mesh::MarkersConversion class.
  void add_bc_neumann(Hermes::vector<std::string> markers) 
  {
    unsigned int n = markers.size();
    if (n <= 0) error("BCTypes::add_bc_dirichlet() expects at least one marker.");
    for (unsigned int i = 0; i < n; i++)
      this->markers_neumann_string_temp.push_back(markers[i]);
    return;
  };

  void add_bc_newton(Hermes::vector<int> markers) 
  {
    unsigned int n = markers.size();
    if (n <= 0) error("BCTypes::add_bc_newton() expects at least one marker.");
    for (unsigned int i = 0; i < n; i++) this->markers_newton.push_back(markers[i]);
    return;
  };

  // A wrapper utilizing the Mesh::MarkersConversion class.
  void add_bc_newton(Hermes::vector<std::string> markers) 
  {
    unsigned int n = markers.size();
    if (n <= 0) error("BCTypes::add_bc_dirichlet() expects at least one marker.");
    for (unsigned int i = 0; i < n; i++)
      this->markers_newton_string_temp.push_back(markers[i]);
    return;
  };

  void add_bc_none(Hermes::vector<int> markers) 
  {
    unsigned int n = markers.size();
    if (n <= 0) error("BCTypes::add_bc_none() expects at least one marker.");
    for (unsigned int i = 0; i < n; i++) this->markers_none.push_back(markers[i]);    
    return;
  };

  // A wrapper utilizing the Mesh::MarkersConversion class.
  void add_bc_none(Hermes::vector<std::string> markers) 
  {
    unsigned int n = markers.size();
    if (n <= 0) error("BCTypes::add_bc_dirichlet() expects at least one marker.");
    for (unsigned int i = 0; i < n; i++)
      this->markers_none_string_temp.push_back(markers[i]);
    return;
  };

  bool is_natural(int marker) {
    bool found = false;
    for (unsigned int i = 0; i < this->markers_neumann.size(); i++) 
      if (this->markers_neumann[i] == marker) found = true;
    for (unsigned int i = 0; i < this->markers_newton.size(); i++) 
      if (this->markers_newton[i] == marker) found = true;
    return found;
  }

  bool is_essential(int marker) {
    for (unsigned int i = 0; i < this->markers_dirichlet.size(); i++) 
      if (this->markers_dirichlet[i] == marker) return true;
    return false;
  }

  bool is_none(int marker) {
    int n = this->markers_none.size();
    for (int i = 0; i < n; i++) if (this->markers_none[i] == marker) return true;
    return false;
  }

  // Default (if not found) is BC_NATURAL.
  virtual BCType get_type(int marker) {
    if (this->is_essential(marker)) return BC_ESSENTIAL;
    if (this->is_none(marker)) return BC_NONE;
    return BC_NATURAL;
  }

  int find_index_dirichlet(int marker) {
      return this->markers_dirichlet.find_index(marker);
  }

  int find_index_neumann(int marker) {
      return this->markers_neumann.find_index(marker);
  }

  int find_index_newton(int marker) {
      return this->markers_newton.find_index(marker);
  }

  int find_index_none(int marker) {
      return this->markers_none.find_index(marker);
  }

  void check_consistency() {
      // Check whether Dirichlet boundary markers are 
      // all nonnegative and mutually distinct.
      int n_bc_dirichlet = this->markers_dirichlet.size();
      for (int i=0; i < n_bc_dirichlet; i++) {
        // Making sure that they are positive.
        if (this->markers_dirichlet[i] <= 0) error("Boundary markers need to be positive.");
        // Making sure that Dirichlet markers are mutually distinct.
        for (int j=i+1; j < n_bc_dirichlet; j++) {
          if(this->markers_dirichlet[i] == this->markers_dirichlet[j]) 
            error("Duplicated Dirichlet boundary marker %d.",
                    this->markers_dirichlet[i]);
        }
        /*
        // Debugging print statements:
        this->markers_dirichlet.print();
        this->markers_neumann.print();
        this->markers_newton.print();
        this->markers_none.print();
        */
        // Cross-checking with the array of Neumann, Newton, and None markers
        int dummy_idx = this->markers_neumann.find_index(this->markers_dirichlet[i], false);
        if (dummy_idx != -1)
            error("Mismatched boundary markers: %d.", this->markers_dirichlet[i]);
        dummy_idx = this->markers_newton.find_index(this->markers_dirichlet[i], false);
        if (dummy_idx != -1)
            error("Mismatched boundary markers: %d.", this->markers_dirichlet[i]);
        dummy_idx = this->markers_none.find_index(this->markers_dirichlet[i], false);
        if (dummy_idx != -1)
            error("Mismatched boundary markers: %d.", this->markers_dirichlet[i]);
      }
      // Check whether Neumann boundary markers are all nonnegative and mutually distinct.
      int n_bc_neumann = this->markers_neumann.size(); 
      for (int i=0; i < n_bc_neumann; i++) {
        // Making sure that they are positive.
        if(this->markers_neumann[i] <= 0) error("Boundary markers need to be positive.");
        // Making sure that Neumann markers are mutually distinct.
        for (int j=i+1; j < n_bc_neumann; j++) {
          if(this->markers_neumann[i] == this->markers_neumann[j]) 
            error("Duplicated Neumann boundary marker %d.",
                    this->markers_neumann[i]);
        }
        // Cross-checking with the array of Dirichlet, Newton, and None markers.
        int dummy_idx = this->markers_dirichlet.find_index(this->markers_neumann[i], false);
        if (dummy_idx != -1)
            error("Mismatched boundary markers: %d.", this->markers_neumann[i]);
        dummy_idx = this->markers_newton.find_index(this->markers_neumann[i], false);
        if (dummy_idx != -1)
            error("Mismatched boundary markers: %d.", this->markers_neumann[i]);
        dummy_idx = this->markers_none.find_index(this->markers_neumann[i], false);
        if (dummy_idx != -1)
            error("Mismatched boundary markers: %d.", this->markers_neumann[i]);
      }
      // Check whether Newton boundary markers are all nonnegative and mutually distinct.
      int n_bc_newton = this->markers_newton.size(); 
      for (int i=0; i < n_bc_newton; i++) {
        // Making sure that they are positive.
        if(this->markers_newton[i] <= 0) error("Boundary markers need to be positive.");
        // Making sure that Newton markers are mutually distinct.
        for (int j=i+1; j < n_bc_newton; j++) {
          if(this->markers_newton[i] == this->markers_newton[j]) 
            error("Duplicated Newton boundary marker %d.",
                    this->markers_newton[i]);
        }
        // Cross-checking with the array of Dirichlet, Neumann, and None markers.
        int dummy_idx = this->markers_dirichlet.find_index(this->markers_newton[i], false);
        if (dummy_idx != -1)
            error("Mismatched boundary markers: %d.", this->markers_newton[i]);
        dummy_idx = this->markers_neumann.find_index(this->markers_newton[i], false);
        if (dummy_idx != -1)
            error("Mismatched boundary markers: %d.", this->markers_newton[i]);
        dummy_idx = this->markers_none.find_index(this->markers_newton[i], false);
        if (dummy_idx != -1)
            error("Mismatched boundary markers: %d.", this->markers_newton[i]);
      }


  }


  BCTypes() {};
  ~BCTypes() {};

  virtual BCTypes *dup() {
      BCTypes *bc = new BCTypes();
      bc->markers_neumann = this->markers_neumann;
      bc->markers_newton = this->markers_newton;
      bc->markers_dirichlet = this->markers_dirichlet;
      bc->markers_none = this->markers_none;
      return bc;
  }

  protected:
    Hermes::vector<int> markers_neumann;
    Hermes::vector<int> markers_newton;
    Hermes::vector<int> markers_dirichlet;
    Hermes::vector<int> markers_none;

    // These members are used temporary for storing markers defined by user-supplied strings.
    Hermes::vector<std::string> markers_neumann_string_temp;
    Hermes::vector<std::string> markers_newton_string_temp;
    Hermes::vector<std::string> markers_dirichlet_string_temp;
    Hermes::vector<std::string> markers_none_string_temp;

    friend class Space;
    friend class BCValues;
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
      if (bc_type_callback == NULL)
          this->bc_type_callback = default_bc_type;
      else
          this->bc_type_callback = bc_type_callback;
  }

  // default (if not found) is BC_NATURAL
  virtual BCType get_type(int marker) {
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
