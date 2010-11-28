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

#ifndef _ADAPT_HCURL_PROJECTION_H_
#define _ADAPT_HCURL_PROJECTION_H_

#include "proj.h"

/// HCurl projection
///
/// FIXME: hex specific
///
/// @ingroup hp-adaptivity
class HERMES_API HCurlProjection : public Projection {
public:
	HCurlProjection(Solution *afn, Element *e, Shapeset *ss);

	virtual double get_error(int split, int son, const Ord3 &order);

protected:
  struct Key
    {
      int i;
      int j;
      
      Key(int i, int j)
      {
        this->i = i;
        this->j = j;
      }
    };
    
    struct Compare
    {
      bool operator()(Key a, Key b) const
      {
        if (a.i < b.i) return true;
        else if (a.i > b.i) return false;
        else
        {
          if (a.j < b.j) return true;
          else if (a.j > b.j) return false;
          else return false;
        }
      }
    };

  static std::map<Key, double*, Compare> proj_mat_entries;

	virtual void calc_projection(int split, int son, Ord3 &order);
};

#endif
