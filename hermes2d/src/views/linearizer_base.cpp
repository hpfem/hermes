// This file is part of Hermes2D.
//
// Copyright 2009 Ivo Hanak <hanak@byte.cz>
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@unr.edu>
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

#include "linearizer_base.h"
namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      Quad2DLin g_quad_lin;

      Quad2DLin::Quad2DLin()
      {
        max_order[0]  = max_order[1]  = 1;
        num_tables[0] = num_tables[1] = 2;
        tables = lin_tables;
        np = lin_np;
      };

      LinearizerBase::LinearizerBase(bool auto_max) : auto_max(auto_max)
      {
        tris = NULL;
        edges = NULL;
        max = -1e100;

        vertex_count = triangle_count = edges_count = 0;

        pthread_mutexattr_t attr;
        pthread_mutexattr_init(&attr);
        pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
        pthread_mutex_init(&data_mutex, &attr);
        pthread_mutexattr_destroy(&attr);
      }
      LinearizerBase::~LinearizerBase()
      {
        if (tris != NULL)
        {
          ::free(tris);
          tris = NULL;
        }
        if (edges != NULL)
        {
          ::free(edges);
          edges = NULL;
        }
        pthread_mutex_destroy(&data_mutex);
      }

      void LinearizerBase::lock_data() const
      { 
        pthread_mutex_lock(&data_mutex);
      }

      void LinearizerBase::unlock_data() const 
      {
        pthread_mutex_unlock(&data_mutex); 
      }

      void LinearizerBase::process_edge(int iv1, int iv2, int marker)
      {
        int mid = peek_vertex(iv1, iv2);
        if (mid != -1)
        {
          process_edge(iv1, mid, marker);
          process_edge(mid, iv2, marker);
        }
        else
          add_edge(iv1, iv2, marker);
      }


      void LinearizerBase::add_edge(int iv1, int iv2, int marker)
      {
#pragma omp critical(realloc_edges)
        { 
          if (edges_count >= edges_size)
            edges = (int3*) realloc(edges, sizeof(int3) * (edges_size = edges_size * 3 / 2));
        edges[edges_count][0] = iv1;
        edges[edges_count][1] = iv2;
        edges[edges_count++][2] = marker;
        }
      }

      int LinearizerBase::peek_vertex(int p1, int p2)
      {
        // search for a vertex with parents p1, p2
        if (p1 > p2) std::swap(p1, p2);
        int index = hash(p1, p2);
        int i = hash_table[index];
        while (i >= 0)
        {
          if (info[i][0] == p1 && info[i][1] == p2) return i;
          i = info[i][2];
        }
        return -1;
      }

      void LinearizerBase::add_triangle(int iv0, int iv1, int iv2)
      {
        int index;
#pragma omp critical(realloc_triangles)
        {
          if (triangle_count >= triangle_size)
          {
            tris = (int3*) realloc(tris, sizeof(int3) * (triangle_size = triangle_size * 2));
            verbose("Linearizer::add_triangle(): realloc to %d", triangle_size);
          }
        index = triangle_count++;

        tris[index][0] = iv0;
        tris[index][1] = iv1;
        tris[index][2] = iv2;
        }
      }

      int LinearizerBase::hash(int p1, int p2)
      {
        return (984120265*p1 + 125965121*p2) & (vertex_size - 1);
      }

      void LinearizerBase::set_max_absolute_value(double max_abs)
      {
        if(max_abs < 0.0)
          warning("Setting of maximum absolute value in Linearizer with a negative value");
        else
        {
          this->auto_max = false;
          this->max = max_abs;
        }
        return;
      }

      double LinearizerBase::get_min_value() const
      {
        return min_val;
      }

      double LinearizerBase::get_max_value() const
      {
        return max_val;
      }

      void LinearizerBase::calc_aabb(double* x, double* y, int stride, int num, double* min_x, double* max_x, double* min_y, double* max_y)
      {
        *min_x = *max_x = *x;
        *min_y = *max_y = *y;

        uint8_t* ptr_x = (uint8_t*)x;
        uint8_t* ptr_y = (uint8_t*)y;
        for(int i = 0; i < num; i++, ptr_x += stride, ptr_y += stride)
        {
          *min_x = std::min(*min_x, *((double*)ptr_x));
          *min_y = std::min(*min_y, *((double*)ptr_y));
          *max_x = std::max(*max_x, *((double*)ptr_x));
          *max_y = std::max(*max_y, *((double*)ptr_y));
        }
      }

      int3* LinearizerBase::get_triangles()
      {
        return this->tris;
      }
      int LinearizerBase::get_num_triangles()
      {
        return this->triangle_count;
      }
      int3* LinearizerBase::get_edges()
      {
        return this->edges;
      }
      int LinearizerBase::get_num_edges()
      {
        return this->edges_count;
      }
    }
  }
}