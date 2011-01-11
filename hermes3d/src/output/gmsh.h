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

#ifndef _GMSH_OUTPUT_ENGINE_H_
#define _GMSH_OUTPUT_ENGINE_H_

#include "../output.h"
#include "../../../hermes_common/matrix.h"

/// GMSH output engine.
///
///
///
/// @ingroup visualization
class HERMES_API GmshOutputEngine : public OutputEngine {
public:
	GmshOutputEngine(FILE *file);
	virtual ~GmshOutputEngine();

	/// Run the output with specified output engine
	///
	/// @return true if ok
	/// @param[in] fn A function that will be visualized
	virtual void out(MeshFunction *fn, const char *name, int item = FN_VAL);
	virtual void out(MeshFunction *fn1, MeshFunction *fn2, MeshFunction *fn3, const char *name, int item = FN_VAL_0);
	virtual void out(Mesh *mesh);
	virtual void out_bc_gmsh(Mesh *mesh, const char *name = "BCs");

	virtual void out_orders_gmsh(Space *space, const char *name = "orders");

	virtual void out(Matrix *mat);

protected:

  struct PtsKey
  {
    unsigned int * vtcs;
    unsigned int size;
    PtsKey()
    {
      vtcs = NULL;
      size = 0;
    }
    PtsKey(unsigned int vtcs_ [], unsigned int size_)
    {
      this->size = size_;
      if(size > 0)
        this->vtcs = new unsigned int [size];
      for(unsigned int i = 0; i < size; i++) {
        unsigned int temp_place = i;
        for(unsigned int j = i + 1; j < size; j++)
          if(vtcs_[j] < vtcs_[temp_place])
            temp_place = j;
        this->vtcs[i] = vtcs_[temp_place];
        vtcs_[temp_place] = vtcs_[i];
      }
    };
    ~PtsKey()
    {
      if(size > 0)
        delete [] vtcs;
    };
    PtsKey(const PtsKey &b)
    {
      size = b.size;
      if(size > 0)
        this->vtcs = new unsigned int [size];
      for(unsigned int i = 0; i < size; i++)
        vtcs[i] = b.vtcs[i];
    };
    PtsKey & operator =(const PtsKey &b)
    {
      if(size > 0)
        delete [] vtcs;
      size = b.size;
      if(size > 0)
        this->vtcs = new unsigned int [size];
      for(unsigned int i = 0; i < size; i++)
        vtcs[i] = b.vtcs[i];
      return *this;
    };
    bool operator <(const PtsKey & other) const
    {
      if(this->size < other.size)
        return true;
      else if(this->size > other.size)
        return false;
      else
        for(unsigned int i = 0; i < this->size; i++)
          if(this->vtcs[i] < other.vtcs[i])
            return true;
          else if(this->vtcs[i] > other.vtcs[i])
            return false;

      return false;
    };
    bool operator ==(const PtsKey & other) const
    {
      if(this->size < other.size)
        return false;
      else if(this->size > other.size)
        return false;
      else
        for(unsigned int i = 0; i < this->size; i++)
          if(this->vtcs[i] < other.vtcs[i])
            return false;
          else if(this->vtcs[i] > other.vtcs[i])
            return false;
      return true;
    };

    bool operator !=(const PtsKey & other) const
    {
      return (!(*this == other));
    };
  };

	/// file into which the output is done
	FILE *out_file;

	void dump_scalars(int mode, int num_pts, Point3D *pts, double *value);
	void dump_vectors(int mode, int num_pts, Point3D *pts, double *v0, double *v1, double *v2);
	void dump_mesh(Mesh *mesh);
};

#endif
