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

#ifndef __H2D_VECTORIZER_H
#define __H2D_VECTORIZER_H

#include "linearizer.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      /// \brief "Vectorizer" is a Linearizer<Scalar> for vector solutions. 
      /// The only difference is that linearized vertices are vector-valued. Also, regularization of the
      /// resulting mesh is not attempted. The class can handle different meshes in
      /// both X and Y components.
      ///
      template<typename Scalar>
      class HERMES_API Vectorizer: public Linearizer<Scalar>
      {
      public:

        Vectorizer();
        ~Vectorizer();

        void process_solution(MeshFunction<Scalar>* xsln, int xitem, MeshFunction<Scalar>* ysln, int yitem, double eps);

      public: //accessors
        double4* get_vertices() const { return verts; }

        int2* get_dashes() const { return dashes; }
        int get_num_dashes() const { return nd; }

        virtual void save_data(const char* filename);
        virtual void load_data(const char* filename);
        virtual void calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const; ///< Returns axis aligned bounding box (AABB) of vertices. Assumes lock.

        void free() {}

      protected:

        MeshFunction<Scalar>*xsln, *ysln;
        int xitem, yitem;
        int xia, xib, yia, yib;
        double4* verts;  ///< vertices: (x, y, xvalue, yvalue) quadruples
        int2* dashes;
        int nd, ed, cd;

        int get_vertex(int p1, int p2, double x, double y, double xvalue, double yvalue);
        int create_vertex(double x, double y, double xvalue, double yvalue);
        void process_dash(int iv1, int iv2);

        int add_vertex()
        {
          if (this->nv >= this->cv)
          {
            this->cv *= 2;
            verts = (double4*) realloc(verts, sizeof(double4) * this->cv);
            this->info = (int4*) realloc(this->info, sizeof(int4) * this->cv);
            verbose("Vectorizer::add_vertex(): realloc to %d", this->cv);
          }
          return this->nv++;
        }

        void add_dash(int iv1, int iv2)
        {
          if (nd >= cd) dashes = (int2*) realloc(dashes, sizeof(int2) * (cd = cd * 3 / 2));
          dashes[nd][0] = iv1;
          dashes[nd++][1] = iv2;
        }

        void push_transform(int son)
        {
          xsln->push_transform(son);
          if (ysln != xsln) ysln->push_transform(son);
        }

        void pop_transform()
        {
          xsln->pop_transform();
          if (ysln != xsln) ysln->pop_transform();
        }

        void process_triangle(int iv0, int iv1, int iv2, int level,
          Scalar* xval, Scalar* yval, double* phx, double* phy, int* indices);

        void process_quad(int iv0, int iv1, int iv2, int iv3, int level,
          Scalar* xval, Scalar* yval, double* phx, double* phy, int* indices);

        void find_min_max();

      };
    }
  }
}
#endif
