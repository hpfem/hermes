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

#ifndef __HERMES_COMMON_BCVALUES_H
#define __HERMES_COMMON_BCVALUES_H

#include "common.h"
#include "vector.h"
#include "error.h"
#include "bctypes.h"

/// Class to link boundary markers where Dirichlet BC are prescribed with 
/// the function representing them.
class HERMES_API BCValues {
public:
  BCValues(double* t) : t(t) { 
    t_set = true;
  };
  BCValues() 
  { 
    t_set = false; 
  };
  ~BCValues() {};

protected:
  /// Type of the function representing the BC. Non-time_dependent.
  typedef scalar (*value_callback)(double, double);

  /// Type of the function representing the BC. time_dependent.
  typedef scalar (*value_callback_time)(double, double, double);

  /// Storage of functions. 
  std::map<int, value_callback> value_callbacks;
  // This member is used temporary for storing markers defined by user-supplied strings.
  std::map<std::string, value_callback> value_callbacks_string_temp;

  std::map<int, value_callback_time> value_callbacks_time;
  // This member is used temporary for storing markers defined by user-supplied strings.
  std::map<std::string, value_callback_time> value_callbacks_time_string_temp;

  /// Registering which condition is time-dependent.
  std::map<int, bool> is_time_dep;

  /// Storage of 1D constants. 
  std::map<int, scalar> value_constants;
  // This member is used temporary for storing markers defined by user-supplied strings.
  std::map<std::string, scalar> value_constants_string_temp;

  /// Storage of zeroes.
  /// This is needed, so that we are able to check, if the zero has been set by user, or
  /// if it is put there by a call to map<>::operator[].
  std::map<int, bool> value_zeroes;
  // This member is used temporary for storing markers defined by user-supplied strings.
  std::map<std::string, bool> value_zeroes_string_temp;

  /// Current time. In the case of time_dependency.
  double* t;

  /// Info flag that t has been set.
  bool t_set;
  
public:
  /// This function checks that there is either a function, or a value defined on one part
  /// of the boudnary (one marker), if there are both, it gives an error.
  void check(int marker)
  {
    bool have_function = is_time_dep[marker] ? value_callbacks_time[marker] != NULL : value_callbacks[marker] != NULL;
    bool have_both_function_types = !is_time_dep[marker] ? value_callbacks_time[marker] != NULL : value_callbacks[marker] != NULL;
    if(have_function && have_both_function_types)
      error("Attempt to define more than one description of the Dirichlet BC \
               on the same part of the boundary.");

#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    if(have_function && (value_constants[marker] != 0 || value_zeroes[marker]))
#else
    if(have_function && (value_constants[marker] != std::complex<double>(0, 0) || value_zeroes[marker]))
#endif
      error("Attempt to define more than one description of the Dirichlet BC \
               on the same part of the boundary.");
    // Cleanup if everything is okay.
    if(value_callbacks[marker] == NULL)
      value_callbacks.erase(marker);
    if(value_callbacks_time[marker] == NULL)
      value_callbacks_time.erase(marker);
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    if(value_constants[marker] == 0 && !value_zeroes[marker])
#else
    if(value_constants[marker] == std::complex<double>(0, 0) && !value_zeroes[marker])
#endif
      value_constants.erase(marker);
  }

  /// Adds the function callback to represent the Dirichlet BC on the parts of
  /// the boundary marked with markers.
  void add_function(Hermes::vector<int> markers, value_callback callback) 
  {
    if(markers.size() == 0)
      error("BCValues::add_function() called without any boundary markers specified.");
    for(unsigned int i = 0; i < markers.size(); i++) {
      /// If we find out that there is already a function present describing the Dirichlet
      /// BC on this part of the boundary, return an error, otherwise store the function.
      if(!(value_callbacks.insert(std::pair<int, value_callback>(markers[i], callback))).second)
        error("Attempt to define more than one function representing the Dirichlet BC \
               on the same part of the boundary.");
      check(markers[i]);
    }
  };

  // A wrapper utilizing the Mesh::MarkersConversion class.
  void add_function(Hermes::vector<std::string> markers, value_callback callback)
  {
    if(markers.size() == 0)
      error("BCValues::add_function() called without any boundary markers specified.");
    for(unsigned int i = 0; i < markers.size(); i++)
      this->value_callbacks_string_temp.insert(std::pair<std::string, value_callback>(markers[i], callback));
  };
    
  /// Adds the function callback to represent the Dirichlet BC on the parts of
  /// the boundary marked with markers.
  void add_timedep_function(Hermes::vector<int> markers, value_callback_time callback) 
  {
    if(markers.size() == 0)
      error("BCValues::add_timedep_function() called without any boundary markers specified.");
    for(unsigned int i = 0; i < markers.size(); i++) {
      /// If we find out that there is already a function present describing the Dirichlet
      /// BC on this part of the boundary, return an error, otherwise store the function.
      if(!(value_callbacks_time.insert(std::pair<int, value_callback_time>(markers[i], callback))).second)
        error("Attempt to define more than one function representing the Dirichlet BC \
               on the same part of the boundary.");
      check(markers[i]);
      is_time_dep[markers[i]] = true;
    }
  };

  // A wrapper utilizing the Mesh::MarkersConversion class.
  void add_timedep_function(Hermes::vector<std::string> markers, value_callback_time callback)
  {
     if(markers.size() == 0)
      error("BCValues::add_timedep_function() called without any boundary markers specified.");
    for(unsigned int i = 0; i < markers.size(); i++)
      this->value_callbacks_time_string_temp.insert(std::pair<std::string, value_callback_time>(markers[i], callback));
  };

  /// The same as add_function(), only supplies a 1D constant.
  void add_const(Hermes::vector<int> markers, scalar value) 
  {
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    if(value == 0) {
#else
    if(value == std::complex<double>(0, 0)) {
#endif
      this->add_zero(markers);
      return;
    }
    if(markers.size() == 0)
      error("BCValues::add_const() called without any boundary markers specified.");
    for(unsigned int i = 0; i < markers.size(); i++) {
      /// If we find out that there is already a value present describing the Dirichlet
      /// BC on this part of the boundary, return an error, otherwise store the value.
      if(!(value_constants.insert(std::pair<int, scalar>(markers[i], value))).second)
        error("Attempt to define more than one value representing the Dirichlet BC \
               on the same part of the boundary.");
      check(markers[i]);
    }
  };

  // A wrapper utilizing the Mesh::MarkersConversion class.
  void add_const(Hermes::vector<std::string> markers, scalar value)
  {
    #if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    if(value == 0) {
    #else
        if(value == std::complex<double>(0, 0)) {
    #endif
          this->add_zero(markers);
          return;
        }
    if(markers.size() == 0)
      error("BCValues::add_const() called without any boundary markers specified.");
    for(unsigned int i = 0; i < markers.size(); i++)
      this->value_constants_string_temp.insert(std::pair<std::string, scalar>(markers[i], value));
  };

  /// The same as add_const(), only supplies a zero.
  void add_zero(Hermes::vector<int> markers) 
  {
    if(markers.size() == 0)
      error("BCValues::add_zero() called without any boundary markers specified.");
    for(unsigned int i = 0; i < markers.size(); i++) {
      /// If we find out that there is already a value present describing the Dirichlet
      /// BC on this part of the boundary, return an error, otherwise store the value.
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    if(!(value_constants.insert(std::pair<int, scalar>(markers[i], 0))).second)
#else
      if(!(value_constants.insert(std::pair<int, scalar>(markers[i], std::complex<double>(0, 0)))).second)
#endif
        error("Attempt to define more than one value representing the Dirichlet BC \
               on the same part of the boundary.");
      value_zeroes.insert(std::pair<int, bool>(markers[i], true));

      check(markers[i]);
    }
  };

  // A wrapper utilizing the Mesh::MarkersConversion class.
  void add_zero(Hermes::vector<std::string> markers)
  {
    if(markers.size() == 0)
      error("BCValues::add_zero() called without any boundary markers specified.");
    for(unsigned int i = 0; i < markers.size(); i++)
      value_zeroes_string_temp.insert(std::pair<std::string, bool>(markers[i], true));
  };

  /// Checks whether all functions representing Dirichlet BCs are really set on a Dirichlet
  /// part of the boundary.
  void check_consistency(BCTypes* bc_types) 
  {
    std::map<int, value_callback>::iterator it_func;
    for(it_func = value_callbacks.begin(); it_func != value_callbacks.end(); it_func++)
      if(!bc_types->is_essential(it_func->first))
        error("BCTypes and BCValues incompatible.");

    std::map<int, scalar>::iterator it_value;
    for(it_value = value_constants.begin(); it_value != value_constants.end(); it_value++)
      if(!bc_types->is_essential(it_value->first))
        error("BCTypes and BCValues incompatible.");
  };

  /// According to the bc_types variable, this function supplies the default function.
  void update(BCTypes* bc_types) 
  {
    std::vector<int>::iterator it;
    for(it = bc_types->markers_dirichlet.begin(); it != bc_types->markers_dirichlet.end(); it++) {
      bool have_function = is_time_dep[*it] ? value_callbacks_time[*it] != NULL : value_callbacks[*it] != NULL;
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
      if(!have_function && (value_constants[*it] == 0 && !value_zeroes[*it]))
#else
      if(!have_function && (value_constants[*it] == std::complex<double>(0, 0) && !value_zeroes[*it]))
#endif
      {
        if(is_time_dep[*it])
          value_callbacks_time.erase(*it);
        else
           value_callbacks.erase(*it);
        value_constants.erase(*it);
        value_zeroes.erase(*it);
        add_zero(*it);
      }
    }
  };

  /// This function returns a boolean value meaning that on this part of the boundary, the prescribed
  /// Dirichlet condition is constant.
  /// ! This function, as the rest of the class, does not check that this marker is really used in the
  /// problem description the user solves.
  bool is_const(int marker)
  {
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    if(value_constants[marker] == 0 && !value_zeroes[marker]) {
#else
    if(value_constants[marker] == std::complex<double>(0, 0) && !value_zeroes[marker]) {
#endif
      value_constants.erase(marker);
      return false;
    }
    return true;
  }

  /// Main function that returns the value.
  scalar calculate(int marker, double x, double y)
  {
    if(is_time_dep[marker]) {
      if(!t_set)
        error("Attempt to retrieve a value of a time-dependent function representing the Dirichlet BC without \
              time set up for the current Space.");
      if(value_callbacks_time[marker] == NULL)
        error("Attempt to retrieve a value of a function representing the Dirichlet BC without \
              this being set up for the current Space.");
      return value_callbacks_time[marker](x, y, *this->t);
    }
    else {
      if(value_callbacks[marker] == NULL)
        error("Attempt to retrieve a value of a function representing the Dirichlet BC without \
              this being set up for the current Space.");
      return value_callbacks[marker](x, y);
    }
  }

  /// An overloaded member for constant Dirichlet BCs.
  scalar calculate(int marker)
  {
#if !defined(H2D_COMPLEX) && !defined(H3D_COMPLEX)
    if(value_constants[marker] == 0 && !value_zeroes[marker])
#else
    if(value_constants[marker] == std::complex<double>(0, 0) && !value_zeroes[marker])
#endif
      error("Attempt to retrieve a value of a value representing the Dirichlet BC without \
            this being set up for the current Space.");
    return value_constants[marker]; 
  }

  /// Duplicates this instance.
  virtual BCValues *dup() {
      BCValues *bv = new BCValues();
      bv->is_time_dep = this->is_time_dep;
      bv->t = this->t;
      bv->t_set = this->t_set;
      bv->value_callbacks = this->value_callbacks;
      bv->value_callbacks_time = this->value_callbacks_time;
      bv->value_constants = this->value_constants;
      bv->value_zeroes = this->value_zeroes;
      return bv;
  }

  friend class Space;
};

#endif

