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

#include "orderizer.h"
#include "space.h"
#include "refmap.h"
#include "linear_data.cpp"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      // vertices
      static int*      ord_np[2]          = { num_vert_tri, num_vert_quad };
      static double3*  ord_tables_tri[2]  = { vert_tri0, vert_tri1 };
      static double3*  ord_tables_quad[2] = { vert_quad0, vert_quad1 };
      static double3** ord_tables[2]      = { ord_tables_tri, ord_tables_quad };

      // triangles
      static int*      num_elem[2]        = { num_elem_tri, num_elem_quad};
      static int3*     ord_elem_tri[2]    = { elem_tri0, elem_tri1 };
      static int3*     ord_elem_quad[2]   = { elem_quad0,  elem_quad1 };
      static int3**    ord_elem[2]        = { ord_elem_tri, ord_elem_quad };

      // edges
      static int*      num_edge[2]        = { num_edge_tri, num_edge_quad};
      static int3*     ord_edge_tri[2]    = { edge_tri0, edge_tri1 };
      static int3*     ord_edge_quad[2]   = { edge_quad0,  edge_quad1 };
      static int3**    ord_edge[2]        = { ord_edge_tri, ord_edge_quad };

      static class Quad2DOrd : public Quad2D
      {
      public:

        Quad2DOrd()
        {
          max_order[0] = max_order[1] = 1;
          num_tables[0] = num_tables[1] = 2;
          tables = ord_tables;
          np = ord_np;
        };

        virtual void dummy_fn() {}
      } quad_ord;

      Orderizer::Orderizer() : LinearizerBase()
      {
        verts = NULL;
        this->label_size = 0;
        ltext = NULL;
        lvert = NULL;
        lbox = NULL;
        tris_orders = NULL;

        label_count = cl1 = cl2 = cl3 = 0;

        for (int i = 0, p = 0; i <= 10; i++)
        {
          for (int j = 0; j <= 10; j++)
          {
            assert((unsigned) p < sizeof(buffer)-5);
            if(i == j)
              sprintf(buffer + p, "%d", i);
            else
              sprintf(buffer + p, "%d|%d", i, j);
            labels[i][j] = buffer + p;
            p += strlen(buffer + p) + 1;
          }
        }
      }

      int Orderizer::add_vertex()
      {
        if(this->vertex_count >= this->vertex_size)
        {
          this->vertex_size *= 2;
          verts = (double3*) realloc(verts, sizeof(double3) * vertex_size);
          this->info = (int4*) realloc(info, sizeof(int4) * vertex_size);
        }
        return this->vertex_count++;
      }

      void Orderizer::make_vert(int & index, double x, double y, double val)
      {
        index = add_vertex();
        verts[index][0] = x;
        verts[index][1] = y;
        verts[index][2] = val;
      }

      template<typename Scalar>
      void Orderizer::process_space(const Space<Scalar>* space)
      {
        // sanity check
        if(space == NULL) throw Hermes::Exceptions::Exception("Space is NULL in Orderizer:process_space().");

        if(!space->is_up_to_date())
          throw Hermes::Exceptions::Exception("The space is not up to date.");

        int type = 1;
        label_count = 0;
        vertex_count = 0;
        triangle_count = 0;
        edges_count = 0;

        // estimate the required number of vertices and triangles
        Mesh* mesh = space->get_mesh();
        if(mesh == NULL)
        {
          throw Hermes::Exceptions::Exception("Mesh is NULL in Orderizer:process_space().");
        }
        int nn = mesh->get_num_active_elements();
        this->vertex_size = std::max(this->vertex_size, 77 * nn);
        this->triangle_size = std::max(this->triangle_size, 64 * nn);
        this->edges_size = std::max(this->edges_size, 16 * nn);
        this->label_size = std::max(this->label_size, nn + 10);
        cl1 = cl2 = cl3 = label_size;

        // reuse or allocate vertex, triangle and edge arrays
        verts = (double3*) realloc(verts, sizeof(double3) * vertex_size);
				this->tris = (int3*) realloc(this->tris, sizeof(int3) * this->triangle_size);
				this->tri_markers = (int*) realloc(this->tri_markers, sizeof(int) * this->triangle_size);
				this->edges = (int2*) realloc(this->edges, sizeof(int2) * this->edges_size);
				this->edge_markers = (int*) realloc(this->edge_markers, sizeof(int) * this->edges_size);
        tris_orders = (int*) realloc(tris_orders, sizeof(int) * triangle_size);
        info = NULL;
        this->empty = false;
        lvert = (int*) realloc(lvert, sizeof(int) * label_size);
        ltext = (char**) realloc(ltext, sizeof(char*) * label_size);
        lbox = (double2*) realloc(lbox, sizeof(double2) * label_size);

        int oo, o[6];

        RefMap refmap;
        refmap.set_quad_2d(&quad_ord);

        // make a mesh illustrating the distribution of polynomial orders over the space
        Element* e;
        for_all_active_elements(e, mesh)
        {
          oo = o[4] = o[5] = space->get_element_order(e->id);
          for (unsigned int k = 0; k < e->get_nvert(); k++)
            o[k] = space->get_edge_order(e, k);

          refmap.set_active_element(e);
          double* x = refmap.get_phys_x(type);
          double* y = refmap.get_phys_y(type);

          double3* pt = quad_ord.get_points(type, e->get_mode());
          int np = quad_ord.get_num_points(type, e->get_mode());
          int id[80];
          assert(np <= 80);

          int mode = e->get_mode();
          if(e->is_quad())
          {
            o[4] = H2D_GET_H_ORDER(oo);
            o[5] = H2D_GET_V_ORDER(oo);
          }
          make_vert(lvert[label_count], x[0], y[0], o[4]);

          for (int i = 1; i < np; i++)
            make_vert(id[i-1], x[i], y[i], o[(int) pt[i][2]]);

          for (int i = 0; i < num_elem[mode][type]; i++)
            this->add_triangle(id[ord_elem[mode][type][i][0]], id[ord_elem[mode][type][i][1]], id[ord_elem[mode][type][i][2]], o[4], e->marker);

          for (int i = 0; i < num_edge[mode][type]; i++)
          {
            if(e->en[ord_edge[mode][type][i][2]]->bnd || (y[ord_edge[mode][type][i][0] + 1] < y[ord_edge[mode][type][i][1] + 1]) ||
              ((y[ord_edge[mode][type][i][0] + 1] == y[ord_edge[mode][type][i][1] + 1]) &&
              (x[ord_edge[mode][type][i][0] + 1] < x[ord_edge[mode][type][i][1] + 1])))
            {
              LinearizerBase::add_edge(id[ord_edge[mode][type][i][0]], id[ord_edge[mode][type][i][1]], e->en[ord_edge[mode][type][i][2]]->marker);
            }
          }

          double xmin = 1e100, ymin = 1e100, xmax = -1e100, ymax = -1e100;
          for (unsigned int k = 0; k < e->get_nvert(); k++)
          {
            if(e->vn[k]->x < xmin) xmin = e->vn[k]->x;
            if(e->vn[k]->x > xmax) xmax = e->vn[k]->x;
            if(e->vn[k]->y < ymin) ymin = e->vn[k]->y;
            if(e->vn[k]->y > ymax) ymax = e->vn[k]->y;
          }
          lbox[label_count][0] = xmax - xmin;
          lbox[label_count][1] = ymax - ymin;
          ltext[label_count++] = labels[o[4]][o[5]];
        }

        refmap.set_quad_2d(&g_quad_2d_std);
      }

      void Orderizer::add_triangle(int iv0, int iv1, int iv2, int order, int marker)
      {
        int index;
#pragma omp critical(realloc_triangles)
        {
          if(this->del_slot >= 0) // reuse a slot after a deleted triangle
          {
            index = this->del_slot;
            del_slot = -1;
          }
          {
            if(triangle_count >= triangle_size)
            {
							tri_markers = (int*) realloc(tri_markers, sizeof(int) * (triangle_size * 2));
              tris = (int3*) realloc(tris, sizeof(int3) * (triangle_size * 2));
              tris_orders = (int*) realloc(tris_orders, sizeof(int) * (triangle_size = triangle_size * 2));
            }
            index = triangle_count++;

            tris[index][0] = iv0;
            tris[index][1] = iv1;
            tris[index][2] = iv2;
						tris_orders[index] = order;
            tri_markers[index] = marker;
          }
        }
      }
      
      void Orderizer::free()
      {
        if(verts != NULL)
        {
          ::free(verts);
          verts = NULL;
        }
        if(lvert != NULL)
        {
          ::free(lvert);
          lvert = NULL;
        }
        if(ltext != NULL)
        {
          ::free(ltext);
          ltext = NULL;
        }
        if(lbox != NULL)
        {
          ::free(lbox);
          lbox = NULL;
        }
        if(tris_orders != NULL)
        {
          ::free(tris_orders);
          tris_orders = NULL;
        }

        LinearizerBase::free();
      }

      Orderizer::~Orderizer()
      {
        free();
      }

      template<typename Scalar>
      void Orderizer::save_orders_vtk(const Space<Scalar>* space, const char* file_name)
      {
        process_space(space);

        FILE* f = fopen(file_name, "wb");
        if(f == NULL) throw Hermes::Exceptions::Exception("Could not open %s for writing.", file_name);
        lock_data();

        // Output header for vertices.
        fprintf(f, "# vtk DataFile Version 2.0\n");
        fprintf(f, "\n");
        fprintf(f, "ASCII\n\n");
        fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

        // Output vertices.
        fprintf(f, "POINTS %d %s\n", this->vertex_count, "float");
        for (int i = 0; i < this->vertex_count; i++)
        {
          fprintf(f, "%g %g %g\n", this->verts[i][0], this->verts[i][1], 0.0);
        }

        // Output elements.
        fprintf(f, "\n");
        fprintf(f, "CELLS %d %d\n", this->triangle_count, 4 * this->triangle_count);
        for (int i = 0; i < this->triangle_count; i++)
        {
          fprintf(f, "3 %d %d %d\n", this->tris[i][0], this->tris[i][1], this->tris[i][2]);
        }

        // Output cell types.
        fprintf(f, "\n");
        fprintf(f, "CELL_TYPES %d\n", this->triangle_count);
        for (int i = 0; i < this->triangle_count; i++)
        {
          fprintf(f, "5\n");    // The "5" means triangle in VTK.
        }

        // This outputs double solution values. Look into Hermes2D/src/output/vtk.cpp
        // for how it is done for vectors.
        fprintf(f, "\n");
        fprintf(f, "CELL_DATA %d\n", this->triangle_count);
        fprintf(f, "SCALARS %s %s %d\n", "Mesh", "float", 1);
        fprintf(f, "LOOKUP_TABLE %s\n", "default");
        for (int i = 0; i < this->triangle_count; i++)
        {
          fprintf(f, "%i\n", this->tris_orders[i]);
        }

        unlock_data();
        fclose(f);
      }

      template<typename Scalar>
      void Orderizer::save_mesh_vtk(const Space<Scalar>* space, const char* file_name)
      {
        process_space(space);

        FILE* f = fopen(file_name, "wb");
        if(f == NULL) throw Hermes::Exceptions::Exception("Could not open %s for writing.", file_name);
        lock_data();

        // Output header for vertices.
        fprintf(f, "# vtk DataFile Version 2.0\n");
        fprintf(f, "\n");
        fprintf(f, "ASCII\n\n");
        fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

        // Output vertices.
        fprintf(f, "POINTS %d %s\n", this->vertex_count, "float");
        for (int i = 0; i < this->vertex_count; i++)
          fprintf(f, "%g %g %g\n", this->verts[i][0], this->verts[i][1], 0.0);

        // Output elements.
        fprintf(f, "\n");
        fprintf(f, "CELLS %d %d\n", this->edges_count, + 3 * this->edges_count);
        for (int i = 0; i < this->edges_count; i++)
          fprintf(f, "2 %d %d\n", this->edges[i][0], this->edges[i][1]);

        // Output cell types.
        fprintf(f, "\n");
        fprintf(f, "CELL_TYPES %d\n", this->edges_count);

        for (int i = 0; i < this->edges_count; i++)
          fprintf(f, "3\n");    // The "5" means triangle in VTK.

        // This outputs double solution values. Look into Hermes2D/src/output/vtk.cpp
        // for how it is done for vectors.
        fprintf(f, "\n");
        fprintf(f, "CELL_DATA %d\n", this->edges_count);
        fprintf(f, "SCALARS %s %s %d\n", "Mesh", "float", 1);
        fprintf(f, "LOOKUP_TABLE %s\n", "default");
        for (int i = 0; i < this->edges_count; i++)
          fprintf(f, "0 \n");
        unlock_data();
        fclose(f);
      }

      int Orderizer::get_labels(int*& lvert, char**& ltext, double2*& lbox) const
      {
        lvert = this->lvert;
        ltext = this->ltext;
        lbox = this->lbox;
        return label_count;
      }

      void Orderizer::calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const
      {
        if(verts == NULL)
          throw Exceptions::Exception("Cannot calculate AABB from NULL vertices");
        calc_aabb(&verts[0][0], &verts[0][1], sizeof(double3), vertex_count, min_x, max_x, min_y, max_y);
      }

      double3* Orderizer::get_vertices()
      {
        return this->verts;
      }
      int Orderizer::get_num_vertices()
      {
        return this->vertex_count;
      }

      template HERMES_API void Orderizer::save_orders_vtk<double>(const Space<double>* space, const char* file_name);
      template HERMES_API void Orderizer::save_orders_vtk<std::complex<double> >(const Space<std::complex<double> >* space, const char* file_name);
      template HERMES_API void Orderizer::save_mesh_vtk<double>(const Space<double>* space, const char* file_name);
      template HERMES_API void Orderizer::save_mesh_vtk<std::complex<double> >(const Space<std::complex<double> >* space, const char* file_name);
      template HERMES_API void Orderizer::process_space<double>(const Space<double>* space);
      template HERMES_API void Orderizer::process_space<std::complex<double> >(const Space<std::complex<double> >* space);
    }
  }
}