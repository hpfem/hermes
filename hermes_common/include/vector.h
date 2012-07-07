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
/*! \file vector.h
\brief Base class for representations of vectors for different solvers.
*/
#ifndef __HERMES_COMMON_TUPLE_H
#define __HERMES_COMMON_TUPLE_H

#include <vector>
#include <stdexcept>   // for exception, runtime_error, out_of_range
#include <typeinfo>
#include <sstream>
#include <stdio.h>
#include "exceptions.h"

namespace Hermes
{
  /// A vector of values.
  /** This class is used to pass a variable number of parameters in a type-safe fashion.
  *  \par Suggested Use
  *  Let us assume a function foo(Hermes::vector<Solution*>&) and instances sln1-sln3 of a class Solution. Then,
  *  - 2 up to 10 parameters: foo(Hermes::vector<Solution*>(&sln1, &sln2, &sln3));
  *  - more than 15 parameters: Fill the instance similarly to STL vector (std::vector).
  *  If needed, the one-parameter version of foo must be created separately, without using Hermes::vector.
  */
  template<typename T>
  class vector : public std::vector<T>
  {
  public:
    /// A default constructor. Creates an empty vector.
    vector() { };
    /// Default std::vector constructor.
    vector(int size) { this->reserve(size); };
    /// 1 parameters constructor.
    /// Problematic when passing as an argument, not for use.
    //vector(const T& a) { this->reserve(1); this->push_back(a);};
    /// 2 parameters constructor.
    vector(const T& a, const T& b) { this->reserve(2); this->push_back(a); this->push_back(b); };
    /// 3 parameters constructor.
    vector(const T& a, const T& b, const T& c) { this->reserve(3); this->push_back(a); this->push_back(b); this->push_back(c); };
    /// 4 parameters constructor.
    vector(const T& a, const T& b, const T& c, const T& d) { this->reserve(4); this->push_back(a); this->push_back(b); this->push_back(c); this->push_back(d); };
    /// 5 parameters constructor.
    vector(const T& a, const T& b, const T& c, const T& d, const T& e) { this->reserve(5); this->push_back(a); this->push_back(b); this->push_back(c); this->push_back(d); this->push_back(e); };
    /// 6 parameters constructor.
    vector(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f) { this->reserve(6); this->push_back(a); this->push_back(b); this->push_back(c); this->push_back(d); this->push_back(e); this->push_back(f); };
    /// 7 parameters constructor.
    vector(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g) { this->reserve(7); this->push_back(a); this->push_back(b); this->push_back(c); this->push_back(d); this->push_back(e); this->push_back(f); this->push_back(g); };
    /// 8 parameters constructor.
    vector(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g, const T& h) { this->reserve(8); this->push_back(a); this->push_back(b); this->push_back(c); this->push_back(d); this->push_back(e); this->push_back(f); this->push_back(g); this->push_back(h); };
    /// 9 parameters constructor.
    vector(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g, const T& h, const T& i) { this->reserve(9); this->push_back(a); this->push_back(b); this->push_back(c); this->push_back(d); this->push_back(e); this->push_back(f); this->push_back(g); this->push_back(h); this->push_back(i); };
    /// 10 parameters constructor.
    vector(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g, const T& h, const T& i, const T& j) { this->reserve(10); this->push_back(a); this->push_back(b); this->push_back(c); this->push_back(d); this->push_back(e); this->push_back(f); this->push_back(g); this->push_back(h); this->push_back(i); this->push_back(j); };
    /// 11 parameters constructor.
    vector(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g, const T& h, const T& i, const T& j, const T& k) { this->reserve(11); this->push_back(a); this->push_back(b); this->push_back(c); this->push_back(d); this->push_back(e); this->push_back(f); this->push_back(g); this->push_back(h); this->push_back(i); this->push_back(j); this->push_back(k);};
    /// 12 parameters constructor.
    vector(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g, const T& h, const T& i, const T& j, const T& k, const T& l) { this->reserve(12); this->push_back(a); this->push_back(b); this->push_back(c); this->push_back(d); this->push_back(e); this->push_back(f); this->push_back(g); this->push_back(h); this->push_back(i); this->push_back(j); this->push_back(k); this->push_back(l);};
    /// 13 parameters constructor.
    vector(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g, const T& h, const T& i, const T& j, const T& k, const T& l, const T& m) { this->reserve(13); this->push_back(a); this->push_back(b); this->push_back(c); this->push_back(d); this->push_back(e); this->push_back(f); this->push_back(g); this->push_back(h); this->push_back(i); this->push_back(j); this->push_back(k); this->push_back(l); this->push_back(m);};
    /// 14 parameters constructor.
    vector(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g, const T& h, const T& i, const T& j, const T& k, const T& l, const T& m, const T& n) { this->reserve(14); this->push_back(a); this->push_back(b); this->push_back(c); this->push_back(d); this->push_back(e); this->push_back(f); this->push_back(g); this->push_back(h); this->push_back(i); this->push_back(j); this->push_back(k); this->push_back(l); this->push_back(m); this->push_back(n);};
    /// 15 parameters constructor.
    vector(const T& a, const T& b, const T& c, const T& d, const T& e, const T& f, const T& g, const T& h, const T& i, const T& j, const T& k, const T& l, const T& m, const T& n, const T& o) { this->reserve(15); this->push_back(a); this->push_back(b); this->push_back(c); this->push_back(d); this->push_back(e); this->push_back(f); this->push_back(g); this->push_back(h); this->push_back(i); this->push_back(j); this->push_back(k); this->push_back(l); this->push_back(m); this->push_back(n); this->push_back(o);};

    // Look up an integer number in an array.
    int find_index_slow(const T& x) {
      for (int i=0; i < this->size(); i++) {
        if((*this)[i] == x)
          return i;
      }
      throw Hermes::Exceptions::Exception("Index not found");
    }

    // Returns maximum of the vector<T> in case of T == int.
    int max() {
      if(this->size() == 0)
        throw Hermes::Exceptions::Exception("Empty vector");
      int m;
      if(typeid((*this)[0]) != typeid(m))
        throw Hermes::Exceptions::Exception("vector<T>::max() called and T != int.");
      m = (int)(*this)[0];
      for (unsigned int i=1; i < this->size(); i++)
        if((int)(*this)[i] > m)
          m = (int)(*this)[i];
      return m;
    }

    // Returns minimum of the vector<T> in case of T == int.
    int min() {
      if(this->size() == 0)
        throw Hermes::Exceptions::Exception("Empty vector");
      int m;
      if(typeid((*this)[0]) != typeid(m))
        throw Hermes::Exceptions::Exception("vector<T>::max() called and T != int.");
      m = (int)(*this)[0];
      for (unsigned int i=1; i < this->size(); i++)
        if((int)(*this)[i] < m)
          m = (int)(*this)[i];
      return m;
    }

    // Look up an integer number in an array.
    // This prepares a permut array, so subsequent calls are very fast
    int find_index(int x, bool throw_exception=true) {
      if(this->size() == 0) {
        if(throw_exception) {
          throw Hermes::Exceptions::Exception("Empty vector");
        }
        else return -1;
      }
      int idx;
      if(typeid((*this)[0]) != typeid(idx))
        throw Hermes::Exceptions::Exception("vector<T>::find_index() called and T != int.");

      if(this->_permut.size() == 0) {
        // Initialize the permut array
        this->_min = this->min();
        this->_max = this->max();
        for (int i=0; i < (int)this->_max+1; i++) this->_permut.push_back(-1);
        for (unsigned int i=0; i < this->size(); i++) this->_permut[(int)(*this)[i]] = i;
      }
      if(((int)this->_min <= x) && (x <= (int)this->_max))
        idx = this->_permut[x];
      else
        idx = -1;
      if(idx == -1) {
        if(throw_exception)
          throw Hermes::Exceptions::Exception("Index in the vector not found");
        else
          return -1;
      }
      return idx;
    }

    void print() {
      printf("[");
      for (int i=0; i < this->size(); i++) printf("%d ", (*this)[i]);
      printf("]\n");
    }

  private:
    std::vector<int> _permut;
    int _min, _max;
  };
} // namespace Hermes

#endif