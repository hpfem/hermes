// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "thread_linearizer.h"
#include "refmap.h"
#include "traverse.h"
#include "exact_solution.h"
#include "api2d.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      template<typename LinearizerDataDimensions>
      LinearizerNew<LinearizerDataDimensions>::LinearizerNew(LinearizerOutputType linearizerOutputType, bool auto_max) : states(nullptr), num_states(0), dmult(1.0), component(0), value_type(0), curvature_epsilon(1e-3), linearizerOutputType(linearizerOutputType)
      {
        xdisp = nullptr;
        user_xdisp = false;
        ydisp = nullptr;
        user_ydisp = false;

        // Threads
        // Local number of threads - to avoid calling it over and over again, and against faults caused by the
        // value being changed while assembling.
        this->threadLinearizerNew = new ThreadLinearizerNew*[this->num_threads_used];
        for (int i = 0; i < this->num_threads_used; i++)
          this->threadLinearizerNew[i] = new ThreadLinearizerNew(this);

#ifndef NOGLUT
        pthread_mutexattr_t attr;
        pthread_mutexattr_init(&attr);
        pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
        pthread_mutex_init(&data_mutex, &attr);
        pthread_mutexattr_destroy(&attr);
#endif
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::set_curvature_epsilon(double curvature_epsilon)
      {
        this->curvature_epsilon = curvature_epsilon;
      }

      template<typename LinearizerDataDimensions>
      double LinearizerNew<LinearizerDataDimensions>::get_curvature_epsilon() const
      {
        return this->curvature_epsilon;
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::set_displacement(MeshFunctionSharedPtr<double> xdisp, MeshFunctionSharedPtr<double> ydisp, double dmult)
      {
        if (xdisp)
        {
          this->xdisp = MeshFunctionSharedPtr<double>(xdisp);
          this->user_xdisp = true;
        }
        if (ydisp)
        {
          this->ydisp = MeshFunctionSharedPtr<double>(ydisp);
          this->user_ydisp = true;
        }
        if (xdisp || ydisp)
          this->dmult = dmult;
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::init(MeshFunctionSharedPtr<double> sln, int item_, double eps)
      {
        // Basic storage.
        this->item = item_;
        this->epsilon = eps;
        this->component = 0;
        this->value_type = 0;

        // Get the component and desired value from item.
        if (item >= 0x40)
        {
          component = 1;
          this->item >>= 6;
        }
        while (!(item & 1))
        {
          this->item >>= 1;
          value_type++;
        }
        //   reset the item to the value before the circus with component, value_type.
        this->item = item_;

        // Store quads & handle meshes
        meshes.clear();
        old_quad = sln->get_quad_2d();
        meshes.push_back(sln->get_mesh());
        if (this->user_xdisp)
        {
          old_quad_x = xdisp->get_quad_2d();
          meshes.push_back(xdisp->get_mesh());
        }
        if (this->user_ydisp)
        {
          old_quad_y = ydisp->get_quad_2d();
          meshes.push_back(ydisp->get_mesh());
        }
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::process_solution(MeshFunctionSharedPtr<double> sln, int item_, double eps)
      {
        // Init the caught parallel exception message.
        this->exceptionMessageCaughtInParallelBlock.clear();

        // lock data.
        lock_data();

        // Initialization of 'global' stuff.
        this->init(sln, item_, eps);

        // Parallelization.
        Traverse trav_master(ydisp == nullptr ? (xdisp == nullptr ? 1 : 2) : (xdisp == nullptr ? 2 : 3));
        states = trav_master.get_states(this->meshes, this->num_states);

#pragma omp parallel shared(trav_master) num_threads(num_threads_used)
        {
          int thread_number = omp_get_thread_num();
          int start = (this->num_states / num_threads_used) * thread_number;
          int end = (this->num_states / num_threads_used) * (thread_number + 1);
          if (thread_number == num_threads_used - 1)
            end = this->num_states;

          try
          {
            this->threadLinearizerNew[thread_number]->init_processing(sln, this);

            for (int state_i = start; state_i < end; state_i++)
            {
              // Exception already thrown -> exit the loop.
              if (!this->exceptionMessageCaughtInParallelBlock.empty())
                break;

              Traverse::State* current_state = states[state_i];

              this->threadLinearizerNew[thread_number]->process_state(current_state);
            }
            this->threadLinearizerNew[thread_number]->deinit_processing();
          }
          catch (std::exception& exception)
          {
#pragma omp critical (exceptionMessageCaughtInParallelBlock)
            this->exceptionMessageCaughtInParallelBlock = exception.what();
          }
        }

        // Free states.
        if (this->states)
        {
          for (int i = 0; i < this->num_states; i++)
            delete states[i];

          ::free(this->states);
          this->states = nullptr;
          this->num_states = 0;
        }

        // Finish.
        this->finish(sln);

        if (!this->exceptionMessageCaughtInParallelBlock.empty())
          throw Hermes::Exceptions::Exception(this->exceptionMessageCaughtInParallelBlock.c_str());
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::lock_data() const
      {
#ifndef NOGLUT
        pthread_mutex_lock(&data_mutex);
#endif
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::unlock_data() const
      {
#ifndef NOGLUT
        pthread_mutex_unlock(&data_mutex);
#endif
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::finish(MeshFunctionSharedPtr<double> sln)
      {
        // regularize the linear mesh
        if (this->exceptionMessageCaughtInParallelBlock.empty())
          find_min_max();

        // select old quadratrues
        sln->set_quad_2d(old_quad);

        // Unlock data.
        this->unlock_data();
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::find_min_max()
      {
        // find min & max vertex values
        this->min_val = 1e100;
        this->max_val = -1e100;
        for (Iterator<vertex_t> it = this->vertices_begin(); !it.end; it++)
        {
          vertex_t& vertex = it.get();

          double magnitude = 0.;
          for (int j = 0; j < LinearizerDataDimensions::dimension; j++)
            magnitude += verts[i][2 + j] * verts[i][2 + j];

          if (finite(magnitude) && magnitude < min_val)
            min_val = magnitude;
          if (finite(magnitude) && magnitude > max_val)
            max_val = magnitude;
        }
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::free()
      {
        for (int i = 0; i < this->num_threads_used; i++)
          this->threadLinearizerNew[i]->free();
      }

      template<typename LinearizerDataDimensions>
      LinearizerNew<LinearizerDataDimensions>::~LinearizerNew()
      {
#ifndef NOGLUT
        pthread_mutex_destroy(&data_mutex);
#endif
        free();
        for (int i = 0; i < this->num_threads_used; i++)
          delete this->threadLinearizerNew[i];
        delete[] this->threadLinearizerNew;
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::set_max_absolute_value(double max_abs)
      {
        if (max_abs < 0.0)
          this->warn("Setting of maximum absolute value in LinearizerNew with a negative value");
        else
        {
        }
        return;
      }

      template<typename LinearizerDataDimensions>
      double LinearizerNew<LinearizerDataDimensions>::get_min_value() const
      {
        return min_val;
      }

      template<typename LinearizerDataDimensions>
      double LinearizerNew<LinearizerDataDimensions>::get_max_value() const
      {
        return max_val;
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::save_solution_vtk(MeshFunctionSharedPtr<double> sln, const char* filename, const char *quantity_name,
        bool mode_3D, int item, double eps)
      {
        if (this->linearizerOutputType != FileExport)
          throw Exceptions::Exception("This LinearizerNew is not meant to be used for file export, create a new one with appropriate linearizerOutputType.");

        process_solution(sln, item, eps);

        FILE* f = fopen(filename, "wb");
        if (f == nullptr) throw Hermes::Exceptions::Exception("Could not open %s for writing.", filename);

        // Output header for vertices.
        fprintf(f, "# vtk DataFile Version 2.0\n");
        fprintf(f, "\n");
        fprintf(f, "ASCII\n\n");
        fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

        // Output vertices.
        fprintf(f, "POINTS %d %s\n", this->get_vertex_count(), "float");
        for (Iterator<vertex_t> it = this->vertices_begin(); !it.end; it++)
        {
          vertex_t& vertex = it.get();
          if (mode_3D == true)
            fprintf(f, "%g %g %g\n", vertex[0], vertex[1], vertex[2]);
          else
            fprintf(f, "%g %g %g\n", vertex[0], vertex[1], 0.0);
        }

        // Output elements.
        fprintf(f, "\n");
        fprintf(f, "CELLS %d %d\n", this->get_triangle_count(), 4 * this->get_triangle_count());
        for (Iterator<triangle_indices_t> it = this->triangle_indices_begin(); !it.end; it++)
        {
          triangle_indices_t& triangle_indices = it.get();
          fprintf(f, "3 %d %d %d\n", triangle_indices[0], triangle_indices[1], triangle_indices[2]);
        }

        // Output cell types.
        fprintf(f, "\n");
        fprintf(f, "CELL_TYPES %d\n", this->get_triangle_count());
        for (int i = 0; i < this->get_triangle_count(); i++)
        {
          fprintf(f, "5\n");    // The "5" means triangle in VTK.
        }

        // This outputs double solution values.
        fprintf(f, "\n");
        fprintf(f, "POINT_DATA %d\n", this->get_vertex_count());
        fprintf(f, "SCALARS %s %s %d\n", quantity_name, "float", 1);
        fprintf(f, "LOOKUP_TABLE %s\n", "default");
        for (Iterator<vertex_t> it = this->vertices_begin(); !it.end; it++)
        {
          vertex_t& vertex = it.get();
          fprintf(f, "%g\n", vertex[2]);
        }

        fclose(f);
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::save_solution_tecplot(MeshFunctionSharedPtr<double> sln, const char* filename, const char *quantity_name,
        int item, double eps)
      {
        if (this->linearizerOutputType != FileExport)
          throw Exceptions::Exception("This LinearizerNew is not meant to be used for file export, create a new one with appropriate linearizerOutputType.");

        process_solution(sln, item, eps);

        FILE* f = fopen(filename, "wb");
        if (f == nullptr) throw Hermes::Exceptions::Exception("Could not open %s for writing.", filename);

        // Output header for vertices.
        fprintf(f, "TITLE = \"%s created by Hermes.\"\n", filename);
        fprintf(f, "VARIABLES = \"X\", \"Y\", \"%s\"\n", quantity_name);
        fprintf(f, "ZONE N = %d, E = %d, DATAPACKING = POINT, ZONETYPE = FETRIANGLE\n", this->get_vertex_count(), this->get_triangle_count());

        // Output vertices.
        for (Iterator<vertex_t> it = this->vertices_begin(); !it.end; it++)
        {
          vertex_t& vertex = it.get();
          fprintf(f, "%g %g %g\n", vertex[0], vertex[1], vertex[2]);
        }

        // Output elements.
        for (Iterator<triangle_indices_t> it = this->triangle_indices_begin(); !it.end; it++)
        {
          triangle_indices_t& triangle_indices = it.get();
          fprintf(f, "%d %d %d\n", triangle_indices[0] + 1, triangle_indices[1] + 1, triangle_indices[2] + 1);
        }

        fclose(f);
      }

      template<typename LinearizerDataDimensions>
      void LinearizerNew<LinearizerDataDimensions>::calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const
      {
        *max_x = *max_y = std::numeric_limits<double>::min();
        *min_x = *min_y = std::numeric_limits<double>::max();

        for (Iterator<vertex_t> it = this->vertices_begin(); !it.end; it++)
        {
          vertex_t& vertex = it.get();
          if (vertex[0] > *max_x)
            *max_x = vertex[0];
          else if (vertex[0] < *min_x)
            *min_x = vertex[0];

          if (vertex[1] > *max_y)
            *max_y = vertex[1];
          else if (vertex[1] < *min_y)
            *min_y = vertex[1];
        }
      }

      template<typename LinearizerDataDimensions>
      template<typename T>
      LinearizerNew<LinearizerDataDimensions>::Iterator<T>::Iterator(const LinearizerNew* linearizer) : current_thread_index(0), current_thread(0), end(false), linearizer(linearizer)
      {
      }

      template<typename LinearizerDataDimensions>
      template<typename T>
      LinearizerNew<LinearizerDataDimensions>::Iterator<T>& LinearizerNew::Iterator<T>::operator++()
      {
        if (this->current_thread_index >= this->current_thread_size - 1)
        {
          if (this->current_thread == this->thread_sizes.size() - 1)
            this->end = true;
          else
          {
            this->current_thread_index = 0;
            this->current_thread_size = this->thread_sizes[++this->current_thread];
          }
        }
        else
          this->current_thread_index++;

        return *this;
      }


      template<typename LinearizerDataDimensions>
      template<>
      triangle_t& LinearizerNew<LinearizerDataDimensions>::Iterator<triangle_t>::get() const
      {
        return this->linearizer->threadLinearizerNew[this->current_thread]->triangles[this->current_thread_index];
      }

      template<typename LinearizerDataDimensions>
      template<>
      int& LinearizerNew<LinearizerDataDimensions>::Iterator<triangle_t>::get_marker() const
      {
        return this->linearizer->threadLinearizerNew[this->current_thread]->triangle_markers[this->current_thread_index];
      }

      template<typename LinearizerDataDimensions>
      template<>
      edge_t& LinearizerNew<LinearizerDataDimensions>::Iterator<edge_t>::get() const
      {
        return this->linearizer->threadLinearizerNew[this->current_thread]->edges[this->current_thread_index];
      }

      template<typename LinearizerDataDimensions>
      template<>
      int& LinearizerNew<LinearizerDataDimensions>::Iterator<edge_t>::get_marker() const
      {
        return this->linearizer->threadLinearizerNew[this->current_thread]->edge_markers[this->current_thread_index];
      }

      template<typename LinearizerDataDimensions>
      template<>
      vertex_t& LinearizerNew<LinearizerDataDimensions>::Iterator<vertex_t>::get() const
      {
        return this->linearizer->threadLinearizerNew[this->current_thread]->vertices[this->current_thread_index];
      }

      template<typename LinearizerDataDimensions>
      template<>
      triangle_indices_t& LinearizerNew<LinearizerDataDimensions>::Iterator<triangle_indices_t>::get() const
      {
        return this->linearizer->threadLinearizerNew[this->current_thread]->triangle_indices[this->current_thread_index];
      }

      template<typename LinearizerDataDimensions>
      template<>
      int& LinearizerNew<LinearizerDataDimensions>::Iterator<triangle_indices_t>::get_marker() const
      {
        return this->linearizer->threadLinearizerNew[this->current_thread]->triangle_markers[this->current_thread_index];
      }

      /// Begin - iterators.
      template<typename LinearizerDataDimensions>
      LinearizerNew<LinearizerDataDimensions>::Iterator<vertex_t> LinearizerNew<LinearizerDataDimensions>::vertices_begin() const
      {
        LinearizerNew<LinearizerDataDimensions>::Iterator<vertex_t> iterator(this);
        for (int i = 0; i < this->num_threads_used; i++)
          iterator.thread_sizes.push_back(this->threadLinearizerNew[i]->vertex_count);
        iterator.current_thread_size = iterator.thread_sizes[0];
        return iterator;
      }
      template<typename LinearizerDataDimensions>
      LinearizerNew<LinearizerDataDimensions>::Iterator<triangle_t> LinearizerNew<LinearizerDataDimensions>::triangles_begin() const
      {
        LinearizerNew<LinearizerDataDimensions>::Iterator<triangle_t> iterator(this);
        for (int i = 0; i < this->num_threads_used; i++)
          iterator.thread_sizes.push_back(this->threadLinearizerNew[i]->triangle_count);
        iterator.current_thread_size = iterator.thread_sizes[0];
        return iterator;
      }
      template<typename LinearizerDataDimensions>
      LinearizerNew<LinearizerDataDimensions>::Iterator<edge_t> LinearizerNew<LinearizerDataDimensions>::edges_begin() const
      {
        LinearizerNew<LinearizerDataDimensions>::Iterator<edge_t> iterator(this);
        for (int i = 0; i < this->num_threads_used; i++)
          iterator.thread_sizes.push_back(this->threadLinearizerNew[i]->edges_count);
        iterator.current_thread_size = iterator.thread_sizes[0];
        return iterator;
      }
      template<typename LinearizerDataDimensions>
      LinearizerNew<LinearizerDataDimensions>::Iterator<triangle_indices_t> LinearizerNew<LinearizerDataDimensions>::triangle_indices_begin() const
      {
        LinearizerNew<LinearizerDataDimensions>::Iterator<triangle_indices_t> iterator(this);
        for (int i = 0; i < this->num_threads_used; i++)
          iterator.thread_sizes.push_back(this->threadLinearizerNew[i]->triangle_count);
        iterator.current_thread_size = iterator.thread_sizes[0];
        return iterator;
      }

      template<typename LinearizerDataDimensions>
      int LinearizerNew<LinearizerDataDimensions>::get_vertex_count() const
      {
        int count = 0;
        for (int i = 0; i < this->num_threads_used; i++)
          count += this->threadLinearizerNew[i]->vertex_count;
        return count;
      }

      template<typename LinearizerDataDimensions>
      int LinearizerNew<LinearizerDataDimensions>::get_triangle_count() const
      {
        int count = 0;
        for (int i = 0; i < this->num_threads_used; i++)
          count += this->threadLinearizerNew[i]->triangle_count;
        return count;
      }

      template<typename LinearizerDataDimensions>
      int LinearizerNew<LinearizerDataDimensions>::get_edge_count() const
      {
        int count = 0;
        for (int i = 0; i < this->num_threads_used; i++)
          count += this->threadLinearizerNew[i]->edges_count;
        return count;
      }

      template<typename LinearizerDataDimensions>
      int LinearizerNew<LinearizerDataDimensions>::get_triangle_index_count() const
      {
        int count = 0;
        for (int i = 0; i < this->num_threads_used; i++)
          count += this->threadLinearizerNew[i]->triangle_count;
        return count;
      }
    }
  }
}
