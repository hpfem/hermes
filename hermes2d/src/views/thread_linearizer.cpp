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

//#define DEBUG_LINEARIZER

static const int default_allocation_multiplier_vertices = 8;
static const int default_allocation_multiplier_triangles = 15;
static const int default_allocation_multiplier_edges = 8;

static const int default_allocation_maxsize_vertices = 500000;
static const int default_allocation_maxsize_triangles = 1000000;
static const int default_allocation_maxsize_edges = 500000;

static const int default_allocation_minsize_vertices = 10000;
static const int default_allocation_minsize_triangles = 15000;
static const int default_allocation_minsize_edges = 10000;

static const double vertex_relative_tolerance = 0.01;

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      template<typename LinearizerDataDimensions>
      ThreadLinearizerMultidimensional<LinearizerDataDimensions>::ThreadLinearizerMultidimensional(LinearizerMultidimensional<LinearizerDataDimensions>* linearizer) : criterion(linearizer->criterion)
      {
        vertex_size = 0;
        triangle_size = 0;
        edges_size = 0;

        // OpenGL part.
        triangles = nullptr;
        edges = nullptr;
        edge_markers = nullptr;

        // FileExport part.
        triangle_indices = nullptr;

        // Common part.
        vertices = nullptr;
        triangle_markers = nullptr;
        hash_table = nullptr;
        info = nullptr;
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::init_linearizer_data(LinearizerMultidimensional<LinearizerDataDimensions>* linearizer)
      {
        this->criterion = linearizer->criterion;
        this->curvature_epsilon = linearizer->curvature_epsilon;
        this->user_xdisp = linearizer->user_xdisp;
        this->user_ydisp = linearizer->user_ydisp;
        this->dmult = linearizer->dmult;
        this->linearizerOutputType = linearizer->linearizerOutputType;

        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
        {
          this->component[k] = linearizer->component[k];
          this->item[k] = linearizer->item[k];
          this->value_type[k] = linearizer->value_type[k];
        }
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::free()
      {
        // OpenGL part.
        free_with_check(this->triangles, true);
        free_with_check(this->edges, true);
        free_with_check(this->edge_markers, true);
        // FileExport part.
        free_with_check(this->triangle_indices, true);
        // Common part.
        free_with_check(this->vertices, true);
        free_with_check(this->triangle_markers, true);
        free_with_check(this->hash_table, true);
        free_with_check(this->info, true);
      }

      template<typename LinearizerDataDimensions>
      ThreadLinearizerMultidimensional<LinearizerDataDimensions>::~ThreadLinearizerMultidimensional()
      {
        free();
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::init_processing(MeshFunctionSharedPtr<double>* sln, LinearizerMultidimensional<LinearizerDataDimensions>* linearizer)
      {
        this->init_linearizer_data(linearizer);

        // Functions.
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
        {
          Solution<double>* solution = dynamic_cast<Solution<double>*>(sln[k].get());
          if (solution && solution->get_type() == HERMES_SLN)
          {
            fns[k] = new Solution<double>();
            fns[k]->copy(solution);
          }
          else
            fns[k] = sln[k]->clone();
          fns[k]->set_quad_2d(&g_quad_lin);
        }
        if (user_xdisp)
        {
          Solution<double>* xdisp_solution = dynamic_cast<Solution<double>*>(linearizer->xdisp.get());
          if (xdisp_solution && xdisp_solution->get_type() == HERMES_SLN)
          {
            fns[LinearizerDataDimensions::dimension] = new Solution<double>();
            fns[LinearizerDataDimensions::dimension]->copy(linearizer->xdisp);
          }
          else
            fns[LinearizerDataDimensions::dimension] = linearizer->xdisp->clone();

          fns[LinearizerDataDimensions::dimension]->set_quad_2d(&g_quad_lin);
        }
        if (user_ydisp)
        {
          Solution<double>* ydisp_solution = dynamic_cast<Solution<double>*>(linearizer->ydisp.get());
          if (ydisp_solution && ydisp_solution->get_type() == HERMES_SLN)
          {
            fns[LinearizerDataDimensions::dimension + (user_xdisp ? 1 : 0)] = new Solution<double>();
            fns[LinearizerDataDimensions::dimension + (user_xdisp ? 1 : 0)]->copy(linearizer->ydisp);
          }
          else
            fns[LinearizerDataDimensions::dimension + (user_xdisp ? 1 : 0)] = linearizer->ydisp->clone();

          fns[LinearizerDataDimensions::dimension + (user_xdisp ? 1 : 0)]->set_quad_2d(&g_quad_lin);
        }

        // Init storage data & counts.
        this->reallocate(sln[0]->get_mesh());
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::reallocate(MeshSharedPtr mesh)
      {
        int number_of_elements = mesh->get_num_elements();

        this->vertex_size = std::min(default_allocation_maxsize_vertices, std::max(default_allocation_multiplier_vertices * number_of_elements, std::max(this->vertex_size, default_allocation_minsize_vertices)));
        this->triangle_size = std::min(default_allocation_maxsize_triangles, std::max(default_allocation_multiplier_triangles * number_of_elements, std::max(this->triangle_size, default_allocation_minsize_triangles)));
        this->edges_size = std::min(default_allocation_maxsize_edges, std::max(default_allocation_multiplier_edges * number_of_elements, std::max(this->edges_size, default_allocation_minsize_edges)));

        // Set counts.
        this->vertex_count = 0;
        this->triangle_count = 0;
        this->edges_count = 0;

        if (this->linearizerOutputType == OpenGL)
        {
          this->triangles = realloc_with_check<ThreadLinearizerMultidimensional, typename LinearizerDataDimensions::triangle_t>(this->triangles, this->triangle_size, this);
          this->edges = realloc_with_check<ThreadLinearizerMultidimensional, typename LinearizerDataDimensions::edge_t>(this->edges, this->edges_size, this);
          this->edge_markers = realloc_with_check<ThreadLinearizerMultidimensional, int>(this->edge_markers, this->edges_size, this);
        }
        else
          this->triangle_indices = realloc_with_check<ThreadLinearizerMultidimensional, triangle_indices_t>(this->triangle_indices, this->triangle_size, this);

        this->vertices = realloc_with_check<ThreadLinearizerMultidimensional, typename LinearizerDataDimensions::vertex_t>(this->vertices, this->vertex_size, this);
        this->triangle_markers = realloc_with_check<ThreadLinearizerMultidimensional, int>(this->triangle_markers, this->triangle_size, this);

        this->hash_table = malloc_with_check<ThreadLinearizerMultidimensional<LinearizerDataDimensions>, int>(this->vertex_size, this, true);
        memset(this->hash_table, 0xff, sizeof(int)* this->vertex_size);

        this->info = malloc_with_check<ThreadLinearizerMultidimensional<LinearizerDataDimensions>, internal_vertex_info_t>(this->vertex_size, this, true);
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::deinit_processing()
      {
        for (unsigned int j = 0; j < (LinearizerDataDimensions::dimension + (this->user_xdisp ? 1 : 0) + (this->user_ydisp ? 1 : 0)); j++)
          delete fns[j];

        free_with_check(this->hash_table, true);
        free_with_check(this->info, true);
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::process_state(Traverse::State* current_state)
      {
        this->rep_element = current_state->e[0];
        this->curved = this->rep_element->is_curved();
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
        {
          fns[k]->set_active_element(current_state->e[k]);
          fns[k]->set_transform(current_state->sub_idx[k]);
          fns[k]->set_quad_order(0, this->item[k]);
          val[k] = fns[k]->get_values(component[k], value_type[k]);
        }

        const double *dx = nullptr;
        const double *dy = nullptr;

        if (user_xdisp)
        {
          fns[LinearizerDataDimensions::dimension]->set_active_element(current_state->e[LinearizerDataDimensions::dimension]);
          fns[LinearizerDataDimensions::dimension]->set_transform(current_state->sub_idx[LinearizerDataDimensions::dimension]);
          fns[LinearizerDataDimensions::dimension]->set_quad_order(0, H2D_FN_VAL);
          dx = fns[LinearizerDataDimensions::dimension]->get_fn_values();
        }

        if (user_ydisp)
        {
          fns[LinearizerDataDimensions::dimension + (user_xdisp ? 1 : 0)]->set_active_element(current_state->e[LinearizerDataDimensions::dimension + (user_xdisp ? 1 : 0)]);
          fns[LinearizerDataDimensions::dimension + (user_xdisp ? 1 : 0)]->set_transform(current_state->sub_idx[LinearizerDataDimensions::dimension + (user_xdisp ? 1 : 0)]);
          fns[LinearizerDataDimensions::dimension + (user_xdisp ? 1 : 0)]->set_quad_order(0, H2D_FN_VAL);
          dy = fns[LinearizerDataDimensions::dimension + (user_xdisp ? 1 : 0)]->get_fn_values();
        }

        int iv[H2D_MAX_NUMBER_VERTICES];
        for (unsigned int i = 0; i < this->rep_element->get_nvert(); i++)
        {
          double x_disp = fns[0]->get_refmap()->get_phys_x(0)[i];
          double y_disp = fns[0]->get_refmap()->get_phys_y(0)[i];
          if (user_xdisp)
            x_disp += dmult * dx[i];
          if (user_ydisp)
            y_disp += dmult * dy[i];

          double value[LinearizerDataDimensions::dimension];
          for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
            value[k] = val[k][i];
          iv[i] = this->get_vertex(-this->rep_element->vn[i]->id, -this->rep_element->vn[i]->id, x_disp, y_disp, value);
        }

        // recur to sub-elements
        if (current_state->e[0]->is_triangle())
          process_triangle(iv[0], iv[1], iv[2], 0);
        else
          process_quad(iv[0], iv[1], iv[2], iv[3], 0);

#ifndef DEBUG_LINEARIZER
        if (this->linearizerOutputType == OpenGL)
        {
          for (unsigned int i = 0; i < this->rep_element->get_nvert(); i++)
            process_edge(iv[i], iv[this->rep_element->next_vert(i)], this->rep_element->en[i]->marker);
        }
#endif
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::process_triangle(int iv0, int iv1, int iv2, int level)
      {
        const double* values[LinearizerDataDimensions::dimension];
        double* physical_x;
        double* physical_y;
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
        {
          fns[k]->set_quad_order(1, item[k]);
          values[k] = fns[k]->get_values(component[k], value_type[k]);
        }
        unsigned short* vertex_indices = tri_indices[0];

        if (curved)
        {
          // obtain physical element coordinates
          RefMap* refmap = fns[0]->get_refmap();
          physical_x = refmap->get_phys_x(1);
          physical_y = refmap->get_phys_y(1);

          if (user_xdisp)
          {
            fns[LinearizerDataDimensions::dimension]->set_quad_order(1, H2D_FN_VAL);
            const double* dx = fns[LinearizerDataDimensions::dimension]->get_fn_values();
            for (int i = 0; i < lin_np_tri[1]; i++)
              physical_x[i] += dmult*dx[i];
          }
          if (user_ydisp)
          {
            fns[LinearizerDataDimensions::dimension + (this->user_xdisp ? 1 : 0)]->set_quad_order(1, H2D_FN_VAL);
            const double* dy = fns[LinearizerDataDimensions::dimension + (this->user_xdisp ? 1 : 0)]->get_fn_values();
            for (int i = 0; i < lin_np_tri[1]; i++)
              physical_y[i] += dmult*dy[i];
          }
        }

        // obtain linearized values and coordinates at the midpoints
        for (int i = 0; i < 2 + LinearizerDataDimensions::dimension; i++)
        {
          midval[i][0] = (this->vertices[iv0][i] + this->vertices[iv1][i])*0.5;
          midval[i][1] = (this->vertices[iv1][i] + this->vertices[iv2][i])*0.5;
          midval[i][2] = (this->vertices[iv2][i] + this->vertices[iv0][i])*0.5;
        };

        // determine whether or not to split the element
        int split;
        if (level == MAX_LINEARIZER_DIVISION_LEVEL)
          split = 0;
        else
        {
          if (this->criterion.adaptive)
            this->split_decision(split, iv0, iv1, iv2, 0, rep_element->get_mode(), values, physical_x, physical_y, vertex_indices);
          else
            split = (level < this->criterion.refinement_level);
        }

        // split the triangle if the error is too large, otherwise produce a linear triangle
        if (split)
        {
          if (curved)
          {
            for (int i = 0; i < 3; i++)
            {
              midval[0][i] = physical_x[vertex_indices[i]];
              midval[1][i] = physical_y[vertex_indices[i]];
            }
          }

          double values_vertices[5][LinearizerDataDimensions::dimension];
          for (int v = 0; v < 3; v++)
          {
            for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
              values_vertices[v][k] = values[k][vertex_indices[v]];
          }

          // obtain mid-edge vertices
          int mid0 = get_vertex(iv0, iv1, midval[0][0], midval[1][0], values_vertices[0]);
          int mid1 = get_vertex(iv1, iv2, midval[0][1], midval[1][1], values_vertices[1]);
          int mid2 = get_vertex(iv2, iv0, midval[0][2], midval[1][2], values_vertices[2]);

          // recur to sub-elements
          this->push_transforms(0);
          process_triangle(iv0, mid0, mid2, level + 1);
          this->pop_transforms();

          this->push_transforms(1);
          process_triangle(mid0, iv1, mid1, level + 1);
          this->pop_transforms();

          this->push_transforms(2);
          process_triangle(mid2, mid1, iv2, level + 1);
          this->pop_transforms();

          this->push_transforms(3);
          process_triangle(mid1, mid2, mid0, level + 1);
          this->pop_transforms();
        }
        else
          add_triangle(iv0, iv1, iv2, this->rep_element->marker);

#ifdef DEBUG_LINEARIZER
        if (this->linearizerOutputType == OpenGL)
        {
          process_edge(iv0, iv1, this->rep_element->en[0]->marker);
          process_edge(iv1, iv2, this->rep_element->en[1]->marker);
          process_edge(iv2, iv0, this->rep_element->en[2]->marker);
        }
#endif
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::process_quad(int iv0, int iv1, int iv2, int iv3, int level)
      {
        const double* values[LinearizerDataDimensions::dimension];
        double* physical_x;
        double* physical_y;
        unsigned short* vertex_indices = quad_indices[0];
        bool flip = this->quad_flip(iv0, iv1, iv2, iv3);
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
        {
          fns[k]->set_quad_order(1, item[k]);
          values[k] = fns[k]->get_values(component[k], value_type[k]);
        }

        if (curved)
        {
          // obtain physical element coordinates
          RefMap* refmap = fns[0]->get_refmap();
          physical_x = refmap->get_phys_x(1);
          physical_y = refmap->get_phys_y(1);

          if (user_xdisp)
          {
            fns[LinearizerDataDimensions::dimension]->set_quad_order(1, H2D_FN_VAL);
            const double* dx = fns[LinearizerDataDimensions::dimension]->get_fn_values();
            for (int i = 0; i < lin_np_tri[1]; i++)
              physical_x[i] += dmult*dx[i];
          }
          if (user_ydisp)
          {
            fns[LinearizerDataDimensions::dimension + (this->user_xdisp ? 1 : 0)]->set_quad_order(1, H2D_FN_VAL);
            const double* dy = fns[LinearizerDataDimensions::dimension + (this->user_xdisp ? 1 : 0)]->get_fn_values();
            for (int i = 0; i < lin_np_tri[1]; i++)
              physical_y[i] += dmult*dy[i];
          }
        }

        // obtain linearized values and coordinates at the midpoints
        for (int i = 0; i < 2 + LinearizerDataDimensions::dimension; i++)
        {
          midval[i][0] = (this->vertices[iv0][i] + this->vertices[iv1][i]) * 0.5;
          midval[i][1] = (this->vertices[iv1][i] + this->vertices[iv2][i]) * 0.5;
          midval[i][2] = (this->vertices[iv2][i] + this->vertices[iv3][i]) * 0.5;
          midval[i][3] = (this->vertices[iv3][i] + this->vertices[iv0][i]) * 0.5;
          midval[i][4] = (midval[i][0] + midval[i][2])  * 0.5;
        };

        // the value of the middle point is not the average of the four vertex values, since quad == 2 triangles
        midval[LinearizerDataDimensions::dimension + 1][4] = flip ?
          (this->vertices[iv0][LinearizerDataDimensions::dimension + 1] + this->vertices[iv2][LinearizerDataDimensions::dimension + 1]) * 0.5
          :
          (this->vertices[iv1][LinearizerDataDimensions::dimension + 1] + this->vertices[iv3][LinearizerDataDimensions::dimension + 1]) * 0.5;

        // determine whether or not to split the element
        int split;
        if (level == MAX_LINEARIZER_DIVISION_LEVEL)
          split = 0;
        else
        {
          if (this->criterion.adaptive)
            this->split_decision(split, iv0, iv1, iv2, 0, rep_element->get_mode(), values, physical_x, physical_y, vertex_indices);
          else
            split = (level < this->criterion.refinement_level ? 3 : 0);
        }

        // split the quad if the error is too large, otherwise produce two linear triangles
        if (split)
        {
          if (curved)
          {
            for (int i = 0; i < 5; i++)
            {
              midval[0][i] = physical_x[vertex_indices[i]];
              midval[1][i] = physical_y[vertex_indices[i]];
            }
          }

          // obtain mid-edge and mid-element vertices
          int mid0, mid1, mid2, mid3, mid4;
          double values_vertices[5][LinearizerDataDimensions::dimension];
          for (int v = 0; v < 5; v++)
          {
            for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
              values_vertices[v][k] = values[k][vertex_indices[v]];
          }
          if (split != 1) mid0 = get_vertex(iv0, iv1, midval[0][0], midval[1][0], values_vertices[0]);
          if (split != 2) mid1 = get_vertex(iv1, iv2, midval[0][1], midval[1][1], values_vertices[1]);
          if (split != 1) mid2 = get_vertex(iv2, iv3, midval[0][2], midval[1][2], values_vertices[2]);
          if (split != 2) mid3 = get_vertex(iv3, iv0, midval[0][3], midval[1][3], values_vertices[3]);
          if (split == 3) mid4 = get_vertex(mid0, mid2, midval[0][4], midval[1][4], values_vertices[4]);

          // recur to sub-elements
          if (split == 3)
          {
            this->push_transforms(0);
            process_quad(iv0, mid0, mid4, mid3, level + 1);
            this->pop_transforms();

            this->push_transforms(1);
            process_quad(mid0, iv1, mid1, mid4, level + 1);
            this->pop_transforms();

            this->push_transforms(2);
            process_quad(mid4, mid1, iv2, mid2, level + 1);
            this->pop_transforms();

            this->push_transforms(3);
            process_quad(mid3, mid4, mid2, iv3, level + 1);
            this->pop_transforms();
          }
          else
            if (split == 1) // h-split
            {
            this->push_transforms(4);
            process_quad(iv0, iv1, mid1, mid3, level + 1);
            this->pop_transforms();

            this->push_transforms(5);
            process_quad(mid3, mid1, iv2, iv3, level + 1);
            this->pop_transforms();
            }
            else // v-split
            {
              this->push_transforms(6);
              process_quad(iv0, mid0, mid2, iv3, level + 1);
              this->pop_transforms();

              this->push_transforms(7);
              process_quad(mid0, iv1, iv2, mid2, level + 1);
              this->pop_transforms();
            }
        }
        else
        {
          // output two linear triangles,
          if (!flip)
          {
            add_triangle(iv3, iv0, iv1, rep_element->marker);
            add_triangle(iv1, iv2, iv3, rep_element->marker);
          }
          else
          {
            add_triangle(iv0, iv1, iv2, rep_element->marker);
            add_triangle(iv2, iv3, iv0, rep_element->marker);
          }
#ifdef DEBUG_LINEARIZER
          if (this->linearizerOutputType == OpenGL)
          {
            process_edge(iv0, iv1, this->rep_element->en[0]->marker);
            process_edge(iv1, iv2, this->rep_element->en[1]->marker);
            process_edge(iv2, iv3, this->rep_element->en[2]->marker);
            process_edge(iv3, iv0, this->rep_element->en[3]->marker);
          }
#endif
        }
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::split_decision(int& split, int iv0, int iv1, int iv2, int iv3, ElementMode2D mode, const double** values, double* physical_x, double* physical_y, unsigned short* vertex_indices) const
      {
        // Initialization.
        split = 0;
        bool done = false;

        // Core of the decision - calculate the approximate error of linearizing the normalized solution
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
        {
          // Errors (split element - non-split element values) in edge midpoints summed up.
          double error = fabs(values[k][vertex_indices[0]] - midval[2 + k][0])
            + fabs(values[k][vertex_indices[1]] - midval[2 + k][1])
            + fabs(values[k][vertex_indices[2]] - midval[2 + k][2]);

          // For the quad we have one more midpoint.
          if (mode == HERMES_MODE_QUAD)
            error += fabs(values[k][vertex_indices[3]] - midval[2 + k][3]);

          // Divide by the edge count.
          error /= (3 + mode);

          // Relative error.
          // max_value_approx here is only an approximation - only taking into account the elements being processed by this thread.
          double relative_error = error / this->max_value_approx;

          // Split ?
          // See the header of this method (split_decision) for explanation.
          // We put 3 here so that it is easier to test 'full split' both for quads && triangles.
          split = (relative_error > this->criterion.error_tolerance) ? 3 : 0;

          // Quads - division type
          if (mode == HERMES_MODE_QUAD && split)
          {
            double horizontal_error = fabs(values[k][vertex_indices[1]] - midval[2 + k][1]) + fabs(values[k][vertex_indices[3]] - midval[2 + k][3]);
            double vertical_error = fabs(values[k][vertex_indices[0]] - midval[2 + k][0]) + fabs(values[k][vertex_indices[2]] - midval[2 + k][2]);

            // Decide whether to split horizontally or vertically only
            // If one error is LINEARIZER_DIRECTIONAL_QUAD_REFINEMENT_REQUIREMENT larger than the other.
            if (horizontal_error > LINEARIZER_DIRECTIONAL_QUAD_REFINEMENT_REQUIREMENT * vertical_error)
              // h-split
              split = 1;
            else if (vertical_error > LINEARIZER_DIRECTIONAL_QUAD_REFINEMENT_REQUIREMENT * horizontal_error)
              // v-split
              split = 2;
            else
              split = 3;
          }
        }

        // If we are not splitting into four elements alreasdy and we have a curved element, check if we have to split because of the curvature.
        if (curved && split != 3)
        {
          for (int i = 0; i < 3 + mode; i++)
          {
            double error = sqr(physical_x[vertex_indices[i]] - midval[0][i])
              + sqr(physical_y[vertex_indices[i]] - midval[1][i]);

            double diameter = sqr(fns[0]->get_active_element()->diameter);

            split = (error / diameter) > this->curvature_epsilon ? 3 : split;
          }
        }
      }

      template<typename LinearizerDataDimensions>
      double ThreadLinearizerMultidimensional<LinearizerDataDimensions>::get_max_value(Traverse::State* current_state)
      {
        double local_max = 0.;
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
        {
          fns[k]->set_active_element(current_state->e[k]);
          fns[k]->set_transform(current_state->sub_idx[k]);
          fns[k]->set_quad_order(0, this->item[k]);
          const double* val = fns[k]->get_values(component[k], value_type[k]);

          for (unsigned int i = 0; i < current_state->e[k]->get_nvert(); i++)
          {
            double f = fabs(val[i]);
            if (f > local_max)
              local_max = f;
          }
        }
        return local_max;
      }

      template<typename LinearizerDataDimensions>
      bool ThreadLinearizerMultidimensional<LinearizerDataDimensions>::quad_flip(int iv0, int iv1, int iv2, int iv3) const
      {
        int a = (this->vertices[iv0][2] > this->vertices[iv1][2]) ? iv0 : iv1;
        int b = (this->vertices[iv2][2] > this->vertices[iv3][2]) ? iv2 : iv3;
        a = (this->vertices[a][2] > this->vertices[b][2]) ? a : b;
        return (a == iv1 || a == iv3);
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::push_transforms(int transform)
      {
        fns[0]->push_transform(transform);

        if (user_xdisp)
        {
          if (fns[1] != fns[0])
            fns[1]->push_transform(transform);
        }
        if (user_ydisp)
        {
          if (fns[this->user_xdisp ? 1 : 2] != fns[0])
          {
            if (user_xdisp && fns[2] == fns[1])
              return;
            fns[this->user_xdisp ? 1 : 2]->push_transform(transform);
          }
        }
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::pop_transforms()
      {
        fns[0]->pop_transform();

        if (user_xdisp)
        {
          if (fns[1] != fns[0])
            fns[1]->pop_transform();
        }
        if (user_ydisp)
        {
          if (fns[this->user_xdisp ? 1 : 2] != fns[0])
          {
            if (user_xdisp && fns[2] == fns[1])
              return;
            fns[this->user_xdisp ? 1 : 2]->pop_transform();
          }
        }
      }

      template<typename LinearizerDataDimensions>
      int ThreadLinearizerMultidimensional<LinearizerDataDimensions>::peek_vertex(int p1, int p2)
      {
        // search for a vertex with parents p1, p2
        int index = hash(p1 > p2 ? p2 : p1, p1 > p2 ? p1 : p2);
        int i = hash_table[index];
        while (i >= 0)
        {
          if (info[i][0] == p1 && info[i][1] == p2)
            return i;
          i = info[i][2];
        }
        return -1;
      }

      template<typename LinearizerDataDimensions>
      int ThreadLinearizerMultidimensional<LinearizerDataDimensions>::hash(int p1, int p2)
      {
        /// \todo Find out if this is good.
        return (984120265 * p1 + 125965121 * p2) & (vertex_size - 1);
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::process_edge(int iv1, int iv2, int marker)
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

      template<typename LinearizerDataDimensions>
      int ThreadLinearizerMultidimensional<LinearizerDataDimensions>::get_vertex(int p1, int p2, double x, double y, double* value)
      {
        // search for an existing vertex
        if (p1 > p2)
          std::swap(p1, p2);
        int index = this->hash(p1, p2);
        int i = 0;
        if (index < this->vertex_count)
        {
          i = this->hash_table[index];
          while (i >= 0 && i < this->vertex_count)
          {
            if ((this->info[i][0] == p1 && this->info[i][1] == p2)
              && (fabs(x - this->vertices[i][0]) < Hermes::HermesEpsilon)
              && (fabs(y - this->vertices[i][1]) < Hermes::HermesEpsilon))
            {
              bool check_value = true;
              for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
              {
                if (fabs(value[k] - this->vertices[i][2 + k] / value[k]) > vertex_relative_tolerance)
                  check_value = false;
              }
              if (check_value)
                return i;
            }
            // note that we won't return a vertex with a different value than the required one;
            // this takes care for discontinuities in the solution, where more vertices
            // with different values will be created
            i = info[i][2];
          }
        }

        // if not found, create a new_ one
        i = add_vertex();

        this->vertices[i][0] = x;
        this->vertices[i][1] = y;
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
          this->vertices[i][2 + k] = value[k];
        this->info[i][0] = p1;
        this->info[i][1] = p2;
        this->info[i][2] = hash_table[index];
        this->hash_table[index] = i;
        return i;
      }

      template<typename LinearizerDataDimensions>
      int ThreadLinearizerMultidimensional<LinearizerDataDimensions>::add_vertex()
      {
        if (this->vertex_count >= this->vertex_size)
        {
          int new_vertex_size = std::ceil(this->vertex_size * 1.5);

          this->vertices = realloc_with_check<ThreadLinearizerMultidimensional, typename LinearizerDataDimensions::vertex_t>(this->vertices, new_vertex_size, this);
          this->info = realloc_with_check<ThreadLinearizerMultidimensional, internal_vertex_info_t>(this->info, new_vertex_size, this);
          this->hash_table = realloc_with_check<ThreadLinearizerMultidimensional, int>(this->hash_table, new_vertex_size, this);
          memset(this->hash_table + this->vertex_size, 0xff, sizeof(int)* (new_vertex_size - this->vertex_size));

          this->vertex_size = new_vertex_size;
        }
        return this->vertex_count++;
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::add_edge(int iv1, int iv2, int marker)
      {
        if (edges_count >= edges_size)
        {
          this->edges_size = std::ceil(this->edges_size * 1.5);
          this->edges = realloc_with_check<ThreadLinearizerMultidimensional, typename LinearizerDataDimensions::edge_t>(this->edges, this->edges_size, this);
          this->edge_markers = realloc_with_check<ThreadLinearizerMultidimensional, int>(this->edge_markers, this->edges_size, this);
        }

        this->edges[edges_count][0][0] = this->vertices[iv1][0];
        this->edges[edges_count][0][1] = this->vertices[iv1][1];
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
          this->edges[edges_count][0][2 + k] = this->vertices[iv1][2 + k];
        this->edges[edges_count][1][0] = this->vertices[iv2][0];
        this->edges[edges_count][1][1] = this->vertices[iv2][1];
        for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
          this->edges[edges_count][1][2 + k] = this->vertices[iv2][2 + k];
        this->edge_markers[edges_count++] = marker;
      }

      template<typename LinearizerDataDimensions>
      void ThreadLinearizerMultidimensional<LinearizerDataDimensions>::add_triangle(int iv0, int iv1, int iv2, int marker)
      {
        if (triangle_count >= triangle_size)
        {
          this->triangle_size = std::ceil(this->triangle_size * 1.5);
          if (this->linearizerOutputType == OpenGL)
            this->triangles = realloc_with_check<ThreadLinearizerMultidimensional, typename LinearizerDataDimensions::triangle_t>(this->triangles, this->triangle_size, this);
          else
            this->triangle_indices = realloc_with_check<ThreadLinearizerMultidimensional, triangle_indices_t>(this->triangle_indices, this->triangle_size, this);
          this->triangle_markers = realloc_with_check<ThreadLinearizerMultidimensional, int>(this->triangle_markers, this->triangle_size, this);
        }

        if (this->linearizerOutputType == OpenGL)
        {
          this->triangles[triangle_count][0][0] = this->vertices[iv0][0];
          this->triangles[triangle_count][0][1] = this->vertices[iv0][1];
          for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
            this->triangles[triangle_count][0][2 + k] = this->vertices[iv0][2 + k];
          this->triangles[triangle_count][1][0] = this->vertices[iv1][0];
          this->triangles[triangle_count][1][1] = this->vertices[iv1][1];
          for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
            this->triangles[triangle_count][1][2 + k] = this->vertices[iv1][2 + k];
          this->triangles[triangle_count][2][0] = this->vertices[iv2][0];
          this->triangles[triangle_count][2][1] = this->vertices[iv2][1];
          for (int k = 0; k < LinearizerDataDimensions::dimension; k++)
            this->triangles[triangle_count][2][2 + k] = this->vertices[iv2][2 + k];
        }
        else
        {
          this->triangle_indices[triangle_count][0] = iv0;
          this->triangle_indices[triangle_count][1] = iv1;
          this->triangle_indices[triangle_count][2] = iv2;
        }

        this->triangle_markers[triangle_count++] = marker;
      }

      template class HERMES_API ThreadLinearizerMultidimensional < ScalarLinearizerDataDimensions<double> > ;
      template class HERMES_API ThreadLinearizerMultidimensional < VectorLinearizerDataDimensions<double> > ;
      template class HERMES_API ThreadLinearizerMultidimensional < ScalarLinearizerDataDimensions<float> > ;
      template class HERMES_API ThreadLinearizerMultidimensional < VectorLinearizerDataDimensions<float> > ;
    }
  }
}
