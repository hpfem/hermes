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
      LinearizerCriterion::LinearizerCriterion(bool adaptive) : adaptive(adaptive)
      {
      }

      LinearizerCriterionAdaptive::LinearizerCriterionAdaptive(double error_tolerance) : LinearizerCriterion(true)
      {
        this->error_tolerance = error_tolerance;
      }

      LinearizerCriterionFixed::LinearizerCriterionFixed(int refinement_level) : LinearizerCriterion(false)
      {
        this->refinement_level = refinement_level;
      }

      template<typename LinearizerDataDimensions>
      LinearizerMultidimensional<LinearizerDataDimensions>::LinearizerMultidimensional(LinearizerOutputType linearizerOutputType, bool auto_max) :
        states(nullptr), num_states(0), dmult(1.0), curvature_epsilon(1e-5), linearizerOutputType(linearizerOutputType), criterion(LinearizerCriterionFixed(1))
      {
        xdisp = nullptr;
        user_xdisp = false;
        ydisp = nullptr;
        user_ydisp = false;

        // Threads
        // Local number of threads - to avoid calling it over and over again, and against faults caused by the
        // value being changed while assembling.
        this->threadLinearizerMultidimensional = malloc_with_check<LinearizerMultidimensional<LinearizerDataDimensions>, ThreadLinearizerMultidimensional<LinearizerDataDimensions>*>(this->num_threads_used, this);
        for (int i = 0; i < this->num_threads_used; i++)
          this->threadLinearizerMultidimensional[i] = new ThreadLinearizerMultidimensional<LinearizerDataDimensions>(this);

#ifndef NOGLUT
        pthread_mutexattr_t attr;
        pthread_mutexattr_init(&attr);
        pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
        pthread_mutex_init(&data_mutex, &attr);
        pthread_mutexattr_destroy(&attr);
#endif
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::set_criterion(LinearizerCriterion criterion)
      {
        this->criterion = criterion;
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::set_curvature_epsilon(double curvature_epsilon)
      {
        this->curvature_epsilon = curvature_epsilon;
      }

      template<typename LinearizerDataDimensions>
      double LinearizerMultidimensional<LinearizerDataDimensions>::get_curvature_epsilon() const
      {
        return this->curvature_epsilon;
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::set_displacement(MeshFunctionSharedPtr<double> xdisp, MeshFunctionSharedPtr<double> ydisp, double dmult)
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
      void LinearizerMultidimensional<LinearizerDataDimensions>::check_data(MeshFunctionSharedPtr<double>* sln)
      {
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
        if (!sln[k])
          throw Exceptions::Exception("Linearizer: too few solutions, probably wrong template argument.");
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::init(MeshFunctionSharedPtr<double>* sln, int* item_)
      {
        // Check.
        this->check_data(sln);

        // Basic storage.
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
        {
          this->item[k] = item_[k];
          this->component[k] = 0;
          this->value_type[k] = 0;

          // Get the component and desired value from item.
          if (item[k] >= 0x40)
          {
            component[k] = 1;
            this->item[k] >>= 6;
          }
          while (!(item[k] & 1))
          {
            this->item[k] >>= 1;
            value_type[k]++;
          }
          //   reset the item to the value before the circus with component, value_type.
          this->item[k] = item_[k];
        }

        // Store quads & handle meshes
        meshes.clear();
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
        {
          old_quad[k] = sln[k]->get_quad_2d();
          meshes.push_back(sln[k]->get_mesh());
        }
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
      void LinearizerMultidimensional<LinearizerDataDimensions>::process_solution(MeshFunctionSharedPtr<double>* sln, int* item_)
      {
        // Init the caught parallel exception message.
        this->exceptionMessageCaughtInParallelBlock.clear();

        // lock data.
        lock_data();

        // Initialization of 'global' stuff.
        this->init(sln, item_);

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
            double max_value_for_adaptive_refinements = 0.;

            this->threadLinearizerMultidimensional[thread_number]->init_processing(sln, this);

            for (int state_i = start; state_i < end; state_i++)
              max_value_for_adaptive_refinements = std::max(max_value_for_adaptive_refinements, this->threadLinearizerMultidimensional[thread_number]->get_max_value(states[state_i]));

            this->threadLinearizerMultidimensional[thread_number]->max_value_approx = max_value_for_adaptive_refinements;

            for (int state_i = start; state_i < end; state_i++)
            {
              // Exception already thrown -> exit the loop.
              if (!this->exceptionMessageCaughtInParallelBlock.empty())
                break;

              Traverse::State* current_state = states[state_i];

              this->threadLinearizerMultidimensional[thread_number]->process_state(current_state);
            }
            this->threadLinearizerMultidimensional[thread_number]->deinit_processing();
          }
          catch (Hermes::Exceptions::Exception& e)
          {
#pragma omp critical (exceptionMessageCaughtInParallelBlock)
            this->exceptionMessageCaughtInParallelBlock = e.info();
          }
          catch (std::exception& e)
          {
#pragma omp critical (exceptionMessageCaughtInParallelBlock)
            this->exceptionMessageCaughtInParallelBlock = e.what();
          }
        }

        // Free states.
        if (this->states)
        {
          for (int i = 0; i < this->num_states; i++)
            delete states[i];

          free_with_check(this->states);
          this->states = nullptr;
          this->num_states = 0;
        }

        // Finish.
        this->finish(sln);

        if (!this->exceptionMessageCaughtInParallelBlock.empty())
          throw Hermes::Exceptions::Exception(this->exceptionMessageCaughtInParallelBlock.c_str());
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::process_solution(MeshFunctionSharedPtr<double> sln, int item_)
      {
        MeshFunctionSharedPtr<double> slns[LinearizerDataDimensions::dimension];
        int items[LinearizerDataDimensions::dimension];
        slns[0] = sln;
        items[0] = item_;
        for (int k = 1; k < LinearizerDataDimensions::dimension; k++)
        {
          slns[k] = nullptr;
          items[k] = 0;
        }

        this->process_solution(slns, items);
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::lock_data() const
      {
#ifndef NOGLUT
        pthread_mutex_lock(&data_mutex);
#endif
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::unlock_data() const
      {
#ifndef NOGLUT
        pthread_mutex_unlock(&data_mutex);
#endif
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::finish(MeshFunctionSharedPtr<double>* sln)
      {
        // regularize the linear mesh
        if (this->exceptionMessageCaughtInParallelBlock.empty())
        {
          find_min_max();
          // Polish triangle vertex indices for FileExport case.
          if (this->linearizerOutputType == FileExport)
          {
            int running_count = 0;
            for (int i = 0; i < this->num_threads_used; i++)
            {
              if (i > 0)
              {
                for (int j = 0; j < this->threadLinearizerMultidimensional[i]->triangle_count; j++)
                {
                  for (int k = 0; k < 3; k++)
                    this->threadLinearizerMultidimensional[i]->triangle_indices[j][k] += running_count;
                }
              }
              running_count += this->threadLinearizerMultidimensional[i]->vertex_count;
            }
          }
        }

        // select old quadratrues
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
          sln[k]->set_quad_2d(old_quad[k]);

        // Unlock data.
        this->unlock_data();
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::find_min_max()
      {
        // find min & max vertex values
        this->min_val = 1e100;
        this->max_val = -1e100;
        for (Iterator<typename LinearizerDataDimensions::vertex_t> it = this->vertices_begin(); !it.end; ++it)
        {
          typename LinearizerDataDimensions::vertex_t& vertex = it.get();

          double magnitude = 0.;
          if (LinearizerDataDimensions::dimension == 1)
            magnitude = vertex[2];
          else
          {
            for (int j = 0; j < LinearizerDataDimensions::dimension; j++)
              magnitude += vertex[2 + j] * vertex[2 + j];

            magnitude = std::sqrt(magnitude);
          }

          if (finite(magnitude) && magnitude < min_val)
            min_val = magnitude;
          if (finite(magnitude) && magnitude > max_val)
            max_val = magnitude;
        }
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::free()
      {
        for (int i = 0; i < this->num_threads_used; i++)
          this->threadLinearizerMultidimensional[i]->free();
      }

      template<typename LinearizerDataDimensions>
      LinearizerMultidimensional<LinearizerDataDimensions>::~LinearizerMultidimensional()
      {
#ifndef NOGLUT
        pthread_mutex_destroy(&data_mutex);
#endif
        free();
        for (int i = 0; i < this->num_threads_used; i++)
          delete this->threadLinearizerMultidimensional[i];
        free_with_check(this->threadLinearizerMultidimensional);
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::set_max_absolute_value(double max_abs)
      {
        if (max_abs < 0.0)
          this->warn("Setting of maximum absolute value in LinearizerMultidimensional with a negative value");
        else
        {
          this->max_val = max_abs;
        }
        return;
      }

      template<typename LinearizerDataDimensions>
      double LinearizerMultidimensional<LinearizerDataDimensions>::get_min_value() const
      {
        return min_val;
      }

      template<typename LinearizerDataDimensions>
      double LinearizerMultidimensional<LinearizerDataDimensions>::get_max_value() const
      {
        return max_val;
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::save_solution_vtk(Hermes::vector<MeshFunctionSharedPtr<double> > slns, Hermes::vector<int> items, const char* filename, const char *quantity_name,
        bool mode_3D)
      {
        if (this->linearizerOutputType != FileExport)
          throw Exceptions::Exception("This LinearizerMultidimensional is not meant to be used for file export, create a new one with appropriate linearizerOutputType.");

        process_solution(&slns[0], &items[0]);

        FILE* f = fopen(filename, "wb");
        if (f == nullptr) throw Hermes::Exceptions::Exception("Could not open %s for writing.", filename);

        // Output header for vertices.
        fprintf(f, "# vtk DataFile Version 2.0\n");
        fprintf(f, "\n");
        fprintf(f, "ASCII\n\n");
        fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

        // Output vertices.
        fprintf(f, "POINTS %d %s\n", this->get_vertex_count(), "float");
        for (Iterator<typename LinearizerDataDimensions::vertex_t> it = this->vertices_begin(); !it.end; ++it)
        {
          typename LinearizerDataDimensions::vertex_t& vertex = it.get();
          if (mode_3D == true)
          {
            std::stringstream ss;
            ss << vertex[0] << " " << vertex[1];
            for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
              ss << " " << vertex[2 + k];
            fprintf(f, "%s\n", ss.str().c_str());
          }
          else
            fprintf(f, "%g %g %g\n", vertex[0], vertex[1], 0.0);
        }

        // Output elements.
        fprintf(f, "\n");
        fprintf(f, "CELLS %d %d\n", this->get_triangle_count(), 4 * this->get_triangle_count());
        for (Iterator<triangle_indices_t> it = this->triangle_indices_begin(); !it.end; ++it)
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
        fprintf(f, "SCALARS %s %s %d\n", quantity_name, "float", LinearizerDataDimensions::dimension);
        fprintf(f, "LOOKUP_TABLE %s\n", "default");
        for (Iterator<typename LinearizerDataDimensions::vertex_t> it = this->vertices_begin(); !it.end; ++it)
        {
          typename LinearizerDataDimensions::vertex_t& vertex = it.get();

          std::stringstream ssi;
          for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
          {
            ssi << vertex[2 + k];
            if (k < LinearizerDataDimensions::dimension - 1)
              ssi << " ";
          }
          fprintf(f, "%s\n", ssi.str().c_str());
        }

        fclose(f);
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::save_solution_vtk(MeshFunctionSharedPtr<double> sln, const char* filename, const char* quantity_name, bool mode_3D, int item)
      {
        Hermes::vector<MeshFunctionSharedPtr<double> > slns;
        Hermes::vector<int> items;
        slns.push_back(sln);
        items.push_back(item);
        LinearizerMultidimensional<LinearizerDataDimensions>::save_solution_vtk(slns, items, filename, quantity_name, mode_3D);
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::save_solution_tecplot(Hermes::vector<MeshFunctionSharedPtr<double> > slns, Hermes::vector<int> items, const char* filename, Hermes::vector<std::string> quantity_names)
      {
        if (this->linearizerOutputType != FileExport)
          throw Exceptions::Exception("This LinearizerMultidimensional is not meant to be used for file export, create a new one with appropriate linearizerOutputType.");

        process_solution(&slns[0], &items[0]);

        FILE* f = fopen(filename, "wb");
        if (f == nullptr) throw Hermes::Exceptions::Exception("Could not open %s for writing.", filename);

        // Output header for vertices.
        fprintf(f, "TITLE = \"%s created by Hermes.\"\n", filename);
        std::stringstream ss;
        ss << "VARIABLES = \"X\", \"Y\",";
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
          ss << " \"" << quantity_names[k] << "\"";
        fprintf(f, "%s\n", ss.str().c_str());
        fprintf(f, "ZONE N = %d, E = %d, DATAPACKING = POINT, ZONETYPE = FETRIANGLE\n", this->get_vertex_count(), this->get_triangle_count());

        // Output vertices.
        for (Iterator<typename LinearizerDataDimensions::vertex_t> it = this->vertices_begin(); !it.end; ++it)
        {
          typename LinearizerDataDimensions::vertex_t& vertex = it.get();

          std::stringstream ss;
          ss << vertex[0] << " " << vertex[1];
          for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
            ss << " " << vertex[2 + k];
          fprintf(f, "%s\n", ss.str().c_str());
        }

        // Output elements.
        for (Iterator<triangle_indices_t> it = this->triangle_indices_begin(); !it.end; ++it)
        {
          triangle_indices_t& triangle_indices = it.get();
          fprintf(f, "%d %d %d\n", triangle_indices[0] + 1, triangle_indices[1] + 1, triangle_indices[2] + 1);
        }

        fclose(f);
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::save_solution_tecplot(MeshFunctionSharedPtr<double> sln, const char* filename, const char* quantity_name, int item)
      {
        Hermes::vector<MeshFunctionSharedPtr<double> > slns;
        Hermes::vector<int> items;
        slns.push_back(sln);
        items.push_back(item);
        Hermes::vector<std::string> quantity_names;
        quantity_names.push_back(quantity_name);
        LinearizerMultidimensional<LinearizerDataDimensions>::save_solution_tecplot(slns, items, filename, quantity_names);
      }

      template<typename LinearizerDataDimensions>
      void LinearizerMultidimensional<LinearizerDataDimensions>::calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const
      {
        *max_x = *max_y = std::numeric_limits<double>::min();
        *min_x = *min_y = std::numeric_limits<double>::max();

        for (Iterator<typename LinearizerDataDimensions::vertex_t> it = this->vertices_begin(); !it.end; ++it)
        {
          typename LinearizerDataDimensions::vertex_t& vertex = it.get();
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
      LinearizerMultidimensional<LinearizerDataDimensions>::Iterator<T>::Iterator(const LinearizerMultidimensional* linearizer) : current_thread_index(0), current_thread(0), end(false), linearizer(linearizer)
      {
      }

      template<typename LinearizerDataDimensions>
      template<typename T>
      void LinearizerMultidimensional<LinearizerDataDimensions>::Iterator<T>::operator++()
      {
        if (this->current_thread_index >= this->current_thread_size - 1)
        {
          if (this->current_thread == this->thread_sizes.size() - 1)
            this->end = true;
          else
          {
            this->current_thread_index = 0;
            this->current_thread_size = this->thread_sizes[++this->current_thread];
            this->check_zero_lengths();
          }
        }
        else
          this->current_thread_index++;
      }

      template<typename LinearizerDataDimensions>
      template<typename T>
      void LinearizerMultidimensional<LinearizerDataDimensions>::Iterator<T>::check_zero_lengths()
      {
        while (this->current_thread_size == 0)
        {
          if (this->current_thread_size == this->thread_sizes[this->thread_sizes.size() - 1])
          {
            this->end = true;
            break;
          }
          this->current_thread_size = this->thread_sizes[++this->current_thread];
        }
      }

      template<>
      template<>
      typename ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t& LinearizerMultidimensional<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<typename ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t>::get() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->triangles[this->current_thread_index];
      }

      template<>
      template<>
      int& LinearizerMultidimensional<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<typename ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t>::get_marker() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->triangle_markers[this->current_thread_index];
      }

      template<>
      template<>
      typename ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t& LinearizerMultidimensional<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<typename ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t>::get() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->edges[this->current_thread_index];
      }

      template<>
      template<>
      int& LinearizerMultidimensional<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<typename ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t>::get_marker() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->edge_markers[this->current_thread_index];
      }

      template<>
      template<>
      typename ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::vertex_t& LinearizerMultidimensional<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<typename ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::vertex_t>::get() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->vertices[this->current_thread_index];
      }

      template<>
      template<>
      triangle_indices_t& LinearizerMultidimensional<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<triangle_indices_t>::get() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->triangle_indices[this->current_thread_index];
      }

      template<>
      template<>
      int& LinearizerMultidimensional<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<triangle_indices_t>::get_marker() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->triangle_markers[this->current_thread_index];
      }

      template<>
      template<>
      typename VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t& LinearizerMultidimensional<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<typename VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t>::get() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->triangles[this->current_thread_index];
      }

      template<>
      template<>
      int& LinearizerMultidimensional<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<typename VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t>::get_marker() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->triangle_markers[this->current_thread_index];
      }

      template<>
      template<>
      typename VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t& LinearizerMultidimensional<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<typename VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t>::get() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->edges[this->current_thread_index];
      }

      template<>
      template<>
      int& LinearizerMultidimensional<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<typename VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t>::get_marker() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->edge_markers[this->current_thread_index];
      }

      template<>
      template<>
      typename VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::vertex_t& LinearizerMultidimensional<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<typename VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::vertex_t>::get() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->vertices[this->current_thread_index];
      }

      template<>
      template<>
      triangle_indices_t& LinearizerMultidimensional<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<triangle_indices_t>::get() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->triangle_indices[this->current_thread_index];
      }

      template<>
      template<>
      int& LinearizerMultidimensional<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<triangle_indices_t>::get_marker() const
      {
        return this->linearizer->threadLinearizerMultidimensional[this->current_thread]->triangle_markers[this->current_thread_index];
      }


      template<typename LinearizerDataDimensions>
#ifdef _MSC_VER
      typename LinearizerMultidimensional<LinearizerDataDimensions>::Iterator<typename LinearizerDataDimensions::vertex_t> LinearizerMultidimensional<LinearizerDataDimensions>::vertices_begin() const
#else
      typename LinearizerMultidimensional<LinearizerDataDimensions>::template Iterator<typename LinearizerDataDimensions::vertex_t> LinearizerMultidimensional<LinearizerDataDimensions>::vertices_begin() const
#endif
      {
        LinearizerMultidimensional<LinearizerDataDimensions>::Iterator<typename LinearizerDataDimensions::vertex_t> iterator(this);
        for (int i = 0; i < this->num_threads_used; i++)
          iterator.thread_sizes.push_back(this->threadLinearizerMultidimensional[i]->vertex_count);
        iterator.current_thread_size = iterator.thread_sizes[0];
        iterator.check_zero_lengths();
        return iterator;
      }
      template<typename LinearizerDataDimensions>
#ifdef _MSC_VER
      typename LinearizerMultidimensional<LinearizerDataDimensions>::Iterator<typename LinearizerDataDimensions::triangle_t> LinearizerMultidimensional<LinearizerDataDimensions>::triangles_begin() const
#else
      typename LinearizerMultidimensional<LinearizerDataDimensions>::template Iterator<typename LinearizerDataDimensions::triangle_t> LinearizerMultidimensional<LinearizerDataDimensions>::triangles_begin() const
#endif
      {
        LinearizerMultidimensional<LinearizerDataDimensions>::Iterator<typename LinearizerDataDimensions::triangle_t> iterator(this);
        for (int i = 0; i < this->num_threads_used; i++)
          iterator.thread_sizes.push_back(this->threadLinearizerMultidimensional[i]->triangle_count);
        iterator.current_thread_size = iterator.thread_sizes[0];
        iterator.check_zero_lengths();
        return iterator;
      }
      template<typename LinearizerDataDimensions>
#ifdef _MSC_VER
      typename LinearizerMultidimensional<LinearizerDataDimensions>::Iterator<typename LinearizerDataDimensions::edge_t> LinearizerMultidimensional<LinearizerDataDimensions>::edges_begin() const
#else
      typename LinearizerMultidimensional<LinearizerDataDimensions>::template Iterator<typename LinearizerDataDimensions::edge_t> LinearizerMultidimensional<LinearizerDataDimensions>::edges_begin() const
#endif
      {
        LinearizerMultidimensional<LinearizerDataDimensions>::Iterator<typename LinearizerDataDimensions::edge_t> iterator(this);
        for (int i = 0; i < this->num_threads_used; i++)
          iterator.thread_sizes.push_back(this->threadLinearizerMultidimensional[i]->edges_count);
        iterator.current_thread_size = iterator.thread_sizes[0];
        iterator.check_zero_lengths();
        return iterator;
      }
      template<typename LinearizerDataDimensions>
#ifdef _MSC_VER
      typename LinearizerMultidimensional<LinearizerDataDimensions>::Iterator<triangle_indices_t> LinearizerMultidimensional<LinearizerDataDimensions>::triangle_indices_begin() const
#else
      typename LinearizerMultidimensional<LinearizerDataDimensions>::template Iterator<triangle_indices_t> LinearizerMultidimensional<LinearizerDataDimensions>::triangle_indices_begin() const
#endif
      {
        LinearizerMultidimensional<LinearizerDataDimensions>::Iterator<triangle_indices_t> iterator(this);
        for (int i = 0; i < this->num_threads_used; i++)
          iterator.thread_sizes.push_back(this->threadLinearizerMultidimensional[i]->triangle_count);
        iterator.current_thread_size = iterator.thread_sizes[0];
        iterator.check_zero_lengths();
        return iterator;
      }

      template<typename LinearizerDataDimensions>
      int LinearizerMultidimensional<LinearizerDataDimensions>::get_vertex_count() const
      {
        int count = 0;
        for (int i = 0; i < this->num_threads_used; i++)
          count += this->threadLinearizerMultidimensional[i]->vertex_count;
        return count;
      }

      template<typename LinearizerDataDimensions>
      int LinearizerMultidimensional<LinearizerDataDimensions>::get_triangle_count() const
      {
        int count = 0;
        for (int i = 0; i < this->num_threads_used; i++)
          count += this->threadLinearizerMultidimensional[i]->triangle_count;
        return count;
      }

      template<typename LinearizerDataDimensions>
      int LinearizerMultidimensional<LinearizerDataDimensions>::get_edge_count() const
      {
        int count = 0;
        for (int i = 0; i < this->num_threads_used; i++)
          count += this->threadLinearizerMultidimensional[i]->edges_count;
        return count;
      }

      template<typename LinearizerDataDimensions>
      int LinearizerMultidimensional<LinearizerDataDimensions>::get_triangle_index_count() const
      {
        int count = 0;
        for (int i = 0; i < this->num_threads_used; i++)
          count += this->threadLinearizerMultidimensional[i]->triangle_count;
        return count;
      }

      template class HERMES_API LinearizerMultidimensional<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >;
      template class HERMES_API LinearizerMultidimensional<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t>;
      template class HERMES_API LinearizerMultidimensional<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::vertex_t>;
      template class HERMES_API LinearizerMultidimensional<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t>;

      template class HERMES_API LinearizerMultidimensional<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >;
      template class HERMES_API LinearizerMultidimensional<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t>;
      template class HERMES_API LinearizerMultidimensional<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE> >::Iterator<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::vertex_t>;
      template class HERMES_API LinearizerMultidimensional<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE > >::Iterator<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t>;
    }
  }
}
