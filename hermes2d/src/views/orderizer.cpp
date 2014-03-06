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
#include "orderizer_quad.cpp"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      // vertices
      static int*      ord_np[2] = { num_vert_tri, num_vert_quad };
      static double3*  ord_tables_tri[2] = { vert_tri0, vert_tri1 };
      static double3*  ord_tables_quad[2] = { vert_quad0, vert_quad1 };
      static double3** ord_tables[2] = { ord_tables_tri, ord_tables_quad };

      // triangles
      static int*      num_elem[2] = { num_elem_tri, num_elem_quad };
      static int3*     ord_elem_tri[2] = { elem_tri0, elem_tri1 };
      static int3*     ord_elem_quad[2] = { elem_quad0, elem_quad1 };
      static int3**    ord_elem[2] = { ord_elem_tri, ord_elem_quad };

      // edges
      static int*      num_edge[2] = { num_edge_tri, num_edge_quad };
      static int3*     ord_edge_tri[2] = { edge_tri0, edge_tri1 };
      static int3*     ord_edge_quad[2] = { edge_quad0, edge_quad1 };
      static int3**    ord_edge[2] = { ord_edge_tri, ord_edge_quad };

      // vertices_simple
      static int*      ord_np_simple[2] = { num_vert_tri_simple, num_vert_quad_simple };
      static double3*  ord_tables_tri_simple[2] = { vert_tri_simple, vert_tri_simple };
      static double3*  ord_tables_quad_simple[2] = { vert_quad_simple, vert_quad_simple };
      static double3** ord_tables_simple[2] = { ord_tables_tri_simple, ord_tables_quad_simple };

      // triangles
      static int*      num_elem_simple[2] = { num_elem_tri_simple, num_elem_quad_simple };
      static int3*     ord_elem_tri_simple[2] = { elem_tri_simple, elem_tri_simple };
      static int3*     ord_elem_quad_simple[2] = { elem_quad_simple, elem_quad_simple };
      static int3**    ord_elem_simple[2] = { ord_elem_tri_simple, ord_elem_quad_simple };

      // edges
      static int*      num_edge_simple[2] = { num_edge_tri_simple, num_edge_quad_simple };
      static int3*     ord_edge_tri_simple[2] = { edge_tri_simple, edge_tri_simple };
      static int3*     ord_edge_quad_simple[2] = { edge_quad_simple, edge_quad_simple };
      static int3**    ord_edge_simple[2] = { ord_edge_tri_simple, ord_edge_quad_simple };

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

        virtual int get_id()
        {
          return 5;
        };
      } quad_ord;

      static class Quad2DOrdSimple : public Quad2D
      {
      public:

        Quad2DOrdSimple()
        {
          max_order[0] = max_order[1] = 1;
          num_tables[0] = num_tables[1] = 2;
          tables = ord_tables_simple;
          np = ord_np_simple;
        };

        virtual int get_id()
        {
          return 6;
        };
      } quad_ord_simple;

      Orderizer::Orderizer()
      {
        verts = nullptr;
        edges = nullptr;
        this->label_size = 0;
        ltext = nullptr;
        lvert = nullptr;
        lbox = nullptr;

        tris = nullptr;
        tri_markers = nullptr;
        edges = nullptr;
        edge_markers = nullptr;

        vertex_count = triangle_count = edges_count = this->vertex_size = this->triangle_size = this->edges_size = 0;

        label_count = 0;

        for (int i = 0, p = 0; i <= 10; i++)
        {
          for (int j = 0; j <= 10; j++)
          {
            assert((unsigned)p < sizeof(buffer)-5);
            if (i == j)
              sprintf(buffer + p, "%d", i);
            else
              sprintf(buffer + p, "%d|%d", i, j);
            labels[i][j] = buffer + p;
            p += strlen(buffer + p) + 1;
          }
        }
#ifndef NOGLUT
        pthread_mutexattr_t attr;
        pthread_mutexattr_init(&attr);
        pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
        pthread_mutex_init(&data_mutex, &attr);
        pthread_mutexattr_destroy(&attr);
#endif
      }

      int Orderizer::add_vertex()
      {
        if (this->vertex_count >= this->vertex_size)
        {
          this->vertex_size *= 2;
          this->verts = (double3*)realloc(verts, sizeof(double3)* vertex_size);
          if ((!this->verts))
          {
            free();
            throw Exceptions::Exception("Orderizer out of memory!");
          }
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

      static const int default_allocation_multiplier_vertices = 10;
      static const int default_allocation_multiplier_triangles = 15;
      static const int default_allocation_multiplier_edges = 10;

      static const int default_allocation_minsize_vertices = 10000;
      static const int default_allocation_minsize_triangles = 15000;
      static const int default_allocation_minsize_edges = 10000;

      void Orderizer::reallocate(MeshSharedPtr mesh)
      {
        int number_of_elements = mesh->get_num_elements();

        this->vertex_size = std::max(default_allocation_multiplier_vertices * number_of_elements, std::max(this->vertex_size, default_allocation_minsize_vertices));
        this->triangle_size = std::max(default_allocation_multiplier_triangles * number_of_elements, std::max(this->triangle_size, default_allocation_minsize_triangles));
        this->edges_size = std::max(default_allocation_multiplier_edges * number_of_elements, std::max(this->edges_size, default_allocation_minsize_edges));

        // Set count.
        this->vertex_count = 0;
        this->triangle_count = 0;
        this->edges_count = 0;

        this->tris = (int3*)realloc(this->tris, sizeof(int3)* this->triangle_size);
        this->tri_markers = (int*)realloc(this->tri_markers, sizeof(int)* this->triangle_size);
        this->edges = (int2*)realloc(this->edges, sizeof(int2)* this->edges_size);
        this->edge_markers = (int*)realloc(this->edge_markers, sizeof(int)* this->edges_size);

        if (!this->tris || !this->tri_markers || !this->edges || !this->edge_markers)
          throw Exceptions::Exception("Orderizer out of memory!");

        this->label_size = std::max(this->label_size, number_of_elements + 10);
        this->label_count = 0;

        this->verts = realloc_with_check<Orderizer, double3>(this->verts, this->vertex_size, this);

        this->lvert = realloc_with_check<Orderizer, int>(this->lvert, label_size, this);

        ltext = realloc_with_check<Orderizer, char *>(this->ltext, label_size, this);

        lbox = realloc_with_check<Orderizer, double2>(this->lbox, label_size, this);
      }

      Orderizer::~Orderizer()
      {
        free();
#ifndef NOGLUT
        pthread_mutex_destroy(&data_mutex);
#endif
      }

      void Orderizer::lock_data() const
      {
#ifndef NOGLUT
        pthread_mutex_lock(&data_mutex);
#endif
      }

      void Orderizer::unlock_data() const
      {
#ifndef NOGLUT
        pthread_mutex_unlock(&data_mutex);
#endif
      }
      template<typename Scalar>
      void Orderizer::process_space(SpaceSharedPtr<Scalar> space, bool show_edge_orders)
      {
        // sanity check
        if (space == nullptr)
          throw Hermes::Exceptions::Exception("Space is nullptr in Orderizer:process_space().");

        if (!space->is_up_to_date())
          throw Hermes::Exceptions::Exception("The space is not up to date.");

        MeshSharedPtr mesh = space->get_mesh();

        // Reallocate.
        this->reallocate(mesh);

        RefMap refmap;

        int oo, o[6];

        // make a mesh illustrating the distribution of polynomial orders over the space
        Element* e;
        for_all_active_elements(e, mesh)
        {
          oo = o[4] = o[5] = space->get_element_order(e->id);
          if (show_edge_orders)
          for (unsigned int k = 0; k < e->get_nvert(); k++)
            o[k] = space->get_edge_order(e, k);
          else if (e->is_curved())
          {
            if (e->is_triangle())
            for (unsigned int k = 0; k < e->get_nvert(); k++)
              o[k] = oo;
            else
            for (unsigned int k = 0; k < e->get_nvert(); k++)
              o[k] = H2D_GET_H_ORDER(oo);
          }

          double3* pt;
          int np;
          double* x;
          double* y;
          if (show_edge_orders || e->is_curved())
          {
            refmap.set_quad_2d(&quad_ord);
            refmap.set_active_element(e);
            x = refmap.get_phys_x(1);
            y = refmap.get_phys_y(1);

            pt = quad_ord.get_points(1, e->get_mode());
            np = quad_ord.get_num_points(1, e->get_mode());
          }
          else
          {
            refmap.set_quad_2d(&quad_ord_simple);
            refmap.set_active_element(e);
            x = refmap.get_phys_x(1);
            y = refmap.get_phys_y(1);

            pt = quad_ord_simple.get_points(1, e->get_mode());
            np = quad_ord_simple.get_num_points(1, e->get_mode());
          }

          int id[80];
          assert(np <= 80);

          int mode = e->get_mode();
          if (e->is_quad())
          {
            o[4] = H2D_GET_H_ORDER(oo);
            o[5] = H2D_GET_V_ORDER(oo);
          }
          if (show_edge_orders || e->is_curved())
          {
            make_vert(lvert[label_count], x[0], y[0], o[4]);

            for (int i = 1; i < np; i++)
              make_vert(id[i - 1], x[i], y[i], o[(int)pt[i][2]]);

            for (int i = 0; i < num_elem[mode][1]; i++)
              this->add_triangle(id[ord_elem[mode][1][i][0]], id[ord_elem[mode][1][i][1]], id[ord_elem[mode][1][i][2]], e->marker);

            for (int i = 0; i < num_edge[mode][1]; i++)
            {
              if (e->en[ord_edge[mode][1][i][2]]->bnd || (y[ord_edge[mode][1][i][0] + 1] < y[ord_edge[mode][1][i][1] + 1]) ||
                ((y[ord_edge[mode][1][i][0] + 1] == y[ord_edge[mode][1][i][1] + 1]) &&
                (x[ord_edge[mode][1][i][0] + 1] < x[ord_edge[mode][1][i][1] + 1])))
              {
                add_edge(id[ord_edge[mode][1][i][0]], id[ord_edge[mode][1][i][1]], e->en[ord_edge[mode][1][i][2]]->marker);
              }
            }
          }
          else
          {
            make_vert(lvert[label_count], x[0], y[0], o[4]);

            for (int i = 1; i < np; i++)
              make_vert(id[i - 1], x[i], y[i], o[(int)pt[i][2]]);

            for (int i = 0; i < num_elem_simple[mode][1]; i++)
              this->add_triangle(id[ord_elem_simple[mode][1][i][0]], id[ord_elem_simple[mode][1][i][1]], id[ord_elem_simple[mode][1][i][2]], e->marker);

            for (int i = 0; i < num_edge_simple[mode][1]; i++)
              add_edge(id[ord_edge_simple[mode][1][i][0]], id[ord_edge_simple[mode][1][i][1]], e->en[ord_edge_simple[mode][1][i][2]]->marker);
          }

          double xmin = 1e100, ymin = 1e100, xmax = -1e100, ymax = -1e100;
          for (unsigned int k = 0; k < e->get_nvert(); k++)
          {
            if (e->vn[k]->x < xmin) xmin = e->vn[k]->x;
            if (e->vn[k]->x > xmax) xmax = e->vn[k]->x;
            if (e->vn[k]->y < ymin) ymin = e->vn[k]->y;
            if (e->vn[k]->y > ymax) ymax = e->vn[k]->y;
          }
          lbox[label_count][0] = xmax - xmin;
          lbox[label_count][1] = ymax - ymin;
          ltext[label_count++] = labels[o[4]][o[5]];
        }

        refmap.set_quad_2d(&g_quad_2d_std);
      }

      void Orderizer::add_edge(int iv1, int iv2, int marker)
      {
#pragma omp critical(realloc_edges)
        {
          if (edges_count >= edges_size)
          {
            edges = (int2*)realloc(edges, sizeof(int2)* (edges_size * 1.5));
            edge_markers = (int*)realloc(edge_markers, sizeof(int)* (edges_size = edges_size * 1.5));
          }
          edges[edges_count][0] = iv1;
          edges[edges_count][1] = iv2;
          edge_markers[edges_count++] = marker;
        }
      }

      void Orderizer::add_triangle(int iv0, int iv1, int iv2, int marker)
      {
        int index;
#pragma omp critical(realloc_triangles)
        {
          if (triangle_count >= triangle_size)
          {
            this->triangle_size *= 2;
            tri_markers = (int*)realloc(tri_markers, sizeof(int)* triangle_size);
            tris = (int3*)realloc(tris, sizeof(int3)* triangle_size);
            if ((!tri_markers) || (!this->tris))
            {
              free();
              throw Exceptions::Exception("Orderizer out of memory!");
            }
          }
          index = triangle_count++;
        }
        tris[index][0] = iv0;
        tris[index][1] = iv1;
        tris[index][2] = iv2;
        tri_markers[index] = marker;
      }

      void Orderizer::free()
      {
          free_with_check(verts, true);
          free_with_check(lvert, true);
          free_with_check(ltext, true);
          free_with_check(lbox, true);
          free_with_check(tris, true);
          free_with_check(tri_markers, true);
          free_with_check(edges, true);
          free_with_check(edge_markers, true);
      }

      template<typename Scalar>
      void Orderizer::save_orders_vtk(SpaceSharedPtr<Scalar> space, const char* file_name)
      {
        process_space(space);

        FILE* f = fopen(file_name, "wb");
        if (f == nullptr) throw Hermes::Exceptions::Exception("Could not open %s for writing.", file_name);

        // Output header for vertices.
        fprintf(f, "# vtk DataFile Version 2.0\n");
        fprintf(f, "\n");
        fprintf(f, "ASCII\n\n");
        fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

        // Output vertices.
        fprintf(f, "POINTS %d %s\n", this->vertex_count, "float");
        for (int i = 0; i < this->vertex_count; i++)
        {
          fprintf(f, "%g %g %g\n", this->verts[i][0], this->verts[i][1], 0.);
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
          fprintf(f, "5\n");    // The "5" means triangle in VTK.

        // This outputs double solution values. Look into Hermes2D/src/output/vtk.cpp
        // for how it is done for vectors.
        fprintf(f, "\n");
        fprintf(f, "POINT_DATA %d\n", this->vertex_count);
        fprintf(f, "SCALARS %s %s %d\n", "Mesh", "float", 1);
        fprintf(f, "LOOKUP_TABLE %s\n", "default");
        for (int i = 0; i < this->vertex_count; i++)
          fprintf(f, "%g \n", this->verts[i][2]);
        fclose(f);
      }

      template<typename Scalar>
      void Orderizer::save_markers_vtk(SpaceSharedPtr<Scalar> space, const char* file_name)
      {
        process_space(space);

        FILE* f = fopen(file_name, "wb");
        if (f == nullptr) throw Hermes::Exceptions::Exception("Could not open %s for writing.", file_name);

        // Output header for vertices.
        fprintf(f, "# vtk DataFile Version 2.0\n");
        fprintf(f, "\n");
        fprintf(f, "ASCII\n\n");
        fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

        // Output vertices.
        fprintf(f, "POINTS %d %s\n", this->vertex_count, "float");
        for (int i = 0; i < this->vertex_count; i++)
        {
          fprintf(f, "%g %g %g\n", this->verts[i][0], this->verts[i][1], 0.);
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
          fprintf(f, "5\n");    // The "5" means triangle in VTK.

        // This outputs double solution values. Look into Hermes2D/src/output/vtk.cpp
        // for how it is done for vectors.
        fprintf(f, "\n");
        fprintf(f, "CELL_DATA %d\n", this->triangle_count);
        fprintf(f, "SCALARS %s %s %d\n", "Mesh", "float", 1);
        fprintf(f, "LOOKUP_TABLE %s\n", "default");
        for (int i = 0; i < this->triangle_count; i++)
          fprintf(f, "%d \n", this->tri_markers[i]);
        fclose(f);
      }

      template<typename Scalar>
      void Orderizer::save_mesh_vtk(SpaceSharedPtr<Scalar> space, const char* file_name)
      {
        process_space(space);

        FILE* f = fopen(file_name, "wb");
        if (f == nullptr) throw Hermes::Exceptions::Exception("Could not open %s for writing.", file_name);

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
        fprintf(f, "CELLS %d %d\n", this->edges_count, +3 * this->edges_count);
        for (int i = 0; i < this->edges_count; i++)
          fprintf(f, "2 %d %d\n", this->edges[i][0], this->edges[i][1]);

        // Output cell types.
        fprintf(f, "\n");
        fprintf(f, "CELL_TYPES %d\n", this->edges_count);

        for (int i = 0; i < this->edges_count; i++)
          fprintf(f, "3\n");    // The "3" means line in VTK.

        // This outputs double solution values. Look into Hermes2D/src/output/vtk.cpp
        // for how it is done for vectors.
        fprintf(f, "\n");
        fprintf(f, "CELL_DATA %d\n", this->edges_count);
        fprintf(f, "SCALARS %s %s %d\n", "Mesh", "float", 1);
        fprintf(f, "LOOKUP_TABLE %s\n", "default");
        for (int i = 0; i < this->edges_count; i++)
          fprintf(f, "0 \n");
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
        if (verts == nullptr)
          throw Exceptions::Exception("Cannot calculate AABB from nullptr vertices");
        calc_aabb(&verts[0][0], &verts[0][1], sizeof(double3), vertex_count, min_x, max_x, min_y, max_y);
      }

      void Orderizer::calc_aabb(double* x, double* y, int stride, int num, double* min_x, double* max_x, double* min_y, double* max_y)
      {
        *min_x = *max_x = *x;
        *min_y = *max_y = *y;

        uint8_t* ptr_x = (uint8_t*)x;
        uint8_t* ptr_y = (uint8_t*)y;
        for (int i = 0; i < num; i++, ptr_x += stride, ptr_y += stride)
        {
          *min_x = std::min(*min_x, *((double*)ptr_x));
          *min_y = std::min(*min_y, *((double*)ptr_y));
          *max_x = std::max(*max_x, *((double*)ptr_x));
          *max_y = std::max(*max_y, *((double*)ptr_y));
        }
      }

      double3* Orderizer::get_vertices()
      {
        return this->verts;
      }

      int Orderizer::get_num_vertices()
      {
        return this->vertex_count;
      }

      int3* Orderizer::get_triangles()
      {
        return this->tris;
      }

      int* Orderizer::get_triangle_markers()
      {
        return this->tri_markers;
      }

      int Orderizer::get_num_triangles()
      {
        return this->triangle_count;
      }

      int2* Orderizer::get_edges()
      {
        return this->edges;
      }

      int* Orderizer::get_edge_markers()
      {
        return this->edge_markers;
      }

      int Orderizer::get_num_edges()
      {
        return this->edges_count;
      }

      template HERMES_API void Orderizer::save_orders_vtk<double>(const SpaceSharedPtr<double> space, const char* file_name);
      template HERMES_API void Orderizer::save_orders_vtk<std::complex<double> >(const SpaceSharedPtr<std::complex<double> > space, const char* file_name);
      template HERMES_API void Orderizer::save_markers_vtk<double>(const SpaceSharedPtr<double> space, const char* file_name);
      template HERMES_API void Orderizer::save_markers_vtk<std::complex<double> >(const SpaceSharedPtr<std::complex<double> > space, const char* file_name);
      template HERMES_API void Orderizer::save_mesh_vtk<double>(const SpaceSharedPtr<double> space, const char* file_name);
      template HERMES_API void Orderizer::save_mesh_vtk<std::complex<double> >(const SpaceSharedPtr<std::complex<double> > space, const char* file_name);
      template HERMES_API void Orderizer::process_space<double>(const SpaceSharedPtr<double> space, bool);
      template HERMES_API void Orderizer::process_space<std::complex<double> >(const SpaceSharedPtr<std::complex<double> > space, bool);
    }
  }
}