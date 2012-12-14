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

#include "vectorizer.h"
#include "refmap.h"
#include "traverse.h"
#include "mesh.h"
#include "api2d.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      extern int tri_indices[5][3];
      extern int quad_indices[9][5];

      extern int lin_np_quad[2];

      class Quad2DLin;

      Vectorizer::Vectorizer() : LinearizerBase(), component_x(0), value_type_x(0), component_y(0), value_type_y(0), curvature_epsilon(1e-3)
      {
        verts = NULL;
        dashes = NULL;
        dashes_count = dashes_size = 0;
      }

      int Vectorizer::create_vertex(double x, double y, double xvalue, double yvalue)
      {
        int i;
#pragma omp critical(realloc_vertices)
        i = add_vertex();
        verts[i][0] = x;
        verts[i][1] = y;
        verts[i][2] = xvalue;
        verts[i][3] = yvalue;
        return i;
      }

      int Vectorizer::get_vertex(int p1, int p2, double x, double y, double xvalue, double yvalue)
      {
        // search for an existing vertex
        if(p1 > p2) std::swap(p1, p2);
        int index = this->hash(p1, p2);
        int i = 0;
        if(index < this->vertex_count)
        {
          i = this->hash_table[index];
          while (i >= 0 && i < this->vertex_count)
          {
            if(this->info[i][0] == p1 && this->info[i][1] == p2)
              return i;
            // note that we won't return a vertex with a different value than the required one;
            // this takes care for discontinuities in the solution, where more vertices
            // with different values will be created
            i = info[i][2];
          }
        }

        // if not found, create a new one
#pragma omp critical(realloc_vertices)
        i = add_vertex();
        verts[i][0] = x;
        verts[i][1] = y;
        verts[i][2] = xvalue;
        verts[i][3] = yvalue;
        this->info[i][0] = p1;
        this->info[i][1] = p2;
        this->info[i][2] = this->hash_table[index];
        this->hash_table[index] = i;
        return i;
      }

      void Vectorizer::process_triangle(MeshFunction<double>** fns, int iv0, int iv1, int iv2, int level,
        double* xval, double* yval, double* phx, double* phy, int* idx, bool curved)
      {
        double midval[4][3];

        if(level < LIN_MAX_LEVEL)
        {
          int i;
          if(!(level & 1))
          {
            // obtain solution values and physical element coordinates
            fns[0]->set_quad_order(1, xitem);
            fns[1]->set_quad_order(1, yitem);
            xval = fns[0]->get_values(component_x, value_type_x);
            yval = fns[1]->get_values(component_y, value_type_y);
            for (i = 0; i < lin_np_tri[1]; i++)
            {
              double m = (sqrt(sqr(xval[i]) + sqr(yval[i])));
#pragma omp critical(max)
                if(finite(m) && fabs(m) > max)
                  max = fabs(m);
            }

            idx = tri_indices[0];

            if(curved)
            {
              RefMap* refmap = fns[0]->get_refmap();
              phx = refmap->get_phys_x(1);
              phy = refmap->get_phys_y(1);
            }
          }

          // obtain linearized values and coordinates at the midpoints
          for (i = 0; i < 4; i++)
          {
            midval[i][0] = (verts[iv0][i] + verts[iv1][i])*0.5;
            midval[i][1] = (verts[iv1][i] + verts[iv2][i])*0.5;
            midval[i][2] = (verts[iv2][i] + verts[iv0][i])*0.5;
          };

          // determine whether or not to split the element
          bool split;
          if(eps >= 1.0)
          {
            //if eps > 1, the user wants a fixed number of refinements (no adaptivity)
            split = (level < eps);
          }
          else
          {
            // calculate the approximate error of linearizing the normalized solution
            double err = fabs(sqrt(sqr(xval[idx[0]]) + sqr(yval[idx[0]])) - sqrt(sqr(midval[2][0]) + sqr(midval[3][0]))) +
              fabs(sqrt(sqr(xval[idx[1]]) + sqr(yval[idx[1]])) - sqrt(sqr(midval[2][1]) + sqr(midval[3][1]))) +
              fabs(sqrt(sqr(xval[idx[2]]) + sqr(yval[idx[2]])) - sqrt(sqr(midval[2][2]) + sqr(midval[3][2])));
            split = !finite(err) || err > max*3*eps;

            // do the same for the curvature
            if(curved)
            {
              for (i = 0; i < 4; i++)
                if(sqr(phx[idx[i]] - midval[0][i]) + sqr(phy[idx[i]] - midval[1][i]) > sqr(fns[0]->get_active_element()->get_diameter()*this->get_curvature_epsilon()))
                {
                  split = true;
                  break;
                }
            }

            // do extra tests at level 0, so as not to miss some functions with zero error at edge midpoints
            if(level == 0 && !split)
            {
              split = (fabs(sqrt(sqr(xval[8]) + sqr(yval[8])) - 0.5*(sqrt(sqr(midval[2][0]) + sqr(midval[3][0])) + sqrt(sqr(midval[2][1]) + sqr(midval[3][1])))) +
                fabs(sqrt(sqr(xval[9]) + sqr(yval[9])) - 0.5*(sqrt(sqr(midval[2][1]) + sqr(midval[3][1])) + sqrt(sqr(midval[2][2]) + sqr(midval[3][2])))) +
                fabs(sqrt(sqr(xval[4]) + sqr(yval[4])) - 0.5*(sqrt(sqr(midval[2][2]) + sqr(midval[3][2])) + sqrt(sqr(midval[2][0]) + sqr(midval[3][0]))))) > max*3*eps;
            }
          }

          // split the triangle if the error is too large, otherwise produce a linear triangle
          if(split)
          {
            if(curved)
              for (i = 0; i < 3; i++)
              {
                midval[0][i] = phx[idx[i]];
                midval[1][i] = phy[idx[i]];
              }

              // obtain mid-edge vertices
              int mid0 = get_vertex(iv0, iv1, midval[0][0], midval[1][0], xval[idx[0]], yval[idx[0]]);
              int mid1 = get_vertex(iv1, iv2, midval[0][1], midval[1][1], xval[idx[1]], yval[idx[1]]);
              int mid2 = get_vertex(iv2, iv0, midval[0][2], midval[1][2], xval[idx[2]], yval[idx[2]]);

              // recur to sub-elements
              fns[0]->push_transform(0); if(fns[1] != fns[0]) fns[1]->push_transform(0);
              process_triangle(fns, iv0, mid0, mid2,  level + 1, xval, yval, phx, phy, tri_indices[1], curved);
              fns[0]->pop_transform(); if(fns[1] != fns[0]) fns[1]->pop_transform();

              fns[0]->push_transform(1); if(fns[1] != fns[0]) fns[1]->push_transform(1);
              process_triangle(fns, mid0, iv1, mid1,  level + 1, xval, yval, phx, phy, tri_indices[2], curved);
              fns[0]->pop_transform(); if(fns[1] != fns[0]) fns[1]->pop_transform();

              fns[0]->push_transform(2); if(fns[1] != fns[0]) fns[1]->push_transform(2);
              process_triangle(fns, mid2, mid1, iv2,  level + 1, xval, yval, phx, phy, tri_indices[3], curved);
              fns[0]->pop_transform(); if(fns[1] != fns[0]) fns[1]->pop_transform();

              fns[0]->push_transform(3); if(fns[1] != fns[0]) fns[1]->push_transform(3);
              process_triangle(fns, mid1, mid2, mid0, level + 1, xval, yval, phx, phy, tri_indices[4], curved);
              fns[0]->pop_transform(); if(fns[1] != fns[0]) fns[1]->pop_transform();
              return;
          }
        }

        // no splitting: output a linear triangle
        add_triangle(iv0, iv1, iv2, fns[0]->get_active_element()->marker);
      }

      void Vectorizer::process_quad(MeshFunction<double>** fns, int iv0, int iv1, int iv2, int iv3, int level,
        double* xval, double* yval, double* phx, double* phy, int* idx, bool curved)
      {
        double midval[4][5];

        // try not to split through the vertex with the largest value
        int a = (sqr(verts[iv1][2]) + sqr(verts[iv1][3]) > sqr(verts[iv0][2]) + sqr(verts[iv0][3])) ? iv0 : iv1;
        int b = (sqr(verts[iv2][2]) + sqr(verts[iv2][3]) > sqr(verts[iv3][2]) + sqr(verts[iv3][3])) ? iv2 : iv3;
        a = (sqr(verts[a][2]) + sqr(verts[a][3]) > sqr(verts[b][2]) + sqr(verts[b][3])) ? a : b;
        int flip = (a == iv1 || a == iv3) ? 1 : 0;

        if(level < LIN_MAX_LEVEL)
        {
          int i;
          if(!(level & 1))
          {
            // obtain solution values and physical element coordinates
            fns[0]->set_quad_order(1, xitem);
            fns[1]->set_quad_order(1, yitem);
            xval = fns[0]->get_values(component_x, value_type_x);
            yval = fns[1]->get_values(component_y, value_type_y);
            for (i = 0; i < lin_np_quad[1]; i++)
            {
              double m = sqrt(sqr(xval[i]) + sqr(yval[i]));
              if(finite(m) && fabs(m) > max)
#pragma omp critical(max)
                if(finite(m) && fabs(m) > max)
                  max = fabs(m);
            }

            // This is just to make some sense.
            if(fabs(max) < 1E-10)
              max = 1E-10;

            idx = quad_indices[0];
        
            if(curved)
            {
              RefMap* refmap = fns[0]->get_refmap();
              phx = refmap->get_phys_x(1);
              phy = refmap->get_phys_y(1);
            }
          }

          // obtain linearized values and coordinates at the midpoints
          for (i = 0; i < 4; i++)
          {
            midval[i][0] = (verts[iv0][i] + verts[iv1][i]) * 0.5;
            midval[i][1] = (verts[iv1][i] + verts[iv2][i]) * 0.5;
            midval[i][2] = (verts[iv2][i] + verts[iv3][i]) * 0.5;
            midval[i][3] = (verts[iv3][i] + verts[iv0][i]) * 0.5;
            midval[i][4] = (midval[i][0]  + midval[i][2])  * 0.5;
          };

          // determine whether or not to split the element
          bool split;

          //if eps > 1, the user wants a fixed number of refinements (no adaptivity)
          if(eps >= 1.0)
            split = (level < eps);
          else
          {
            // the value of the middle point is not the average of the four vertex values, since quad == 2 triangles
            midval[2][4] = flip ? (verts[iv0][2] + verts[iv2][2]) * 0.5 : (verts[iv1][2] + verts[iv3][2]) * 0.5;
            midval[3][4] = flip ? (verts[iv0][3] + verts[iv2][3]) * 0.5 : (verts[iv1][3] + verts[iv3][3]) * 0.5;

            // calculate the approximate error of linearizing the normalized solution
            double err = fabs(sqrt(sqr(xval[idx[0]]) + sqr(yval[idx[0]])) - sqrt(sqr(midval[2][0]) + sqr(midval[3][0]))) +
              fabs(sqrt(sqr(xval[idx[1]]) + sqr(yval[idx[1]])) - sqrt(sqr(midval[2][1]) + sqr(midval[3][1]))) +
              fabs(sqrt(sqr(xval[idx[2]]) + sqr(yval[idx[2]])) - sqrt(sqr(midval[2][2]) + sqr(midval[3][2]))) +
              fabs(sqrt(sqr(xval[idx[3]]) + sqr(yval[idx[3]])) - sqrt(sqr(midval[2][3]) + sqr(midval[3][3]))) +
              fabs(sqrt(sqr(xval[idx[4]]) + sqr(yval[idx[4]])) - sqrt(sqr(midval[2][4]) + sqr(midval[3][4])));
            split = !finite(err) || err > max*40*eps;

            // do the same for the curvature
            if(curved && !split)
            {
              double cm2 = sqr(fns[0]->get_active_element()->get_diameter()*this->get_curvature_epsilon());
              if(sqr(phx[idx[1]] - midval[0][1]) + sqr(phy[idx[1]] - midval[1][1]) > cm2 ||
                sqr(phx[idx[3]] - midval[0][3]) + sqr(phy[idx[3]] - midval[1][3]) > cm2) split = true;
              
              if(sqr(phx[idx[0]] - midval[0][0]) + sqr(phy[idx[0]] - midval[1][0]) > cm2 ||
                sqr(phx[idx[2]] - midval[0][2]) + sqr(phy[idx[2]] - midval[1][2]) > cm2) split = true;
            }

            // do extra tests at level 0, so as not to miss functions with zero error at edge midpoints
            if(level == 0 && !split)
            {
              err = fabs(sqrt(sqr(xval[13]) + sqr(yval[13])) - 0.5*(sqrt(sqr(midval[2][0]) + sqr(midval[3][0])) + sqrt(sqr(midval[2][1]) + sqr(midval[3][1])))) +
                fabs(sqrt(sqr(xval[17]) + sqr(yval[17])) - 0.5*(sqrt(sqr(midval[2][1]) + sqr(midval[3][1])) + sqrt(sqr(midval[2][2]) + sqr(midval[3][2])))) +
                fabs(sqrt(sqr(xval[20]) + sqr(yval[20])) - 0.5*(sqrt(sqr(midval[2][2]) + sqr(midval[3][2])) + sqrt(sqr(midval[2][3]) + sqr(midval[3][3])))) +
                fabs(sqrt(sqr(xval[9]) + sqr(yval[9]))  - 0.5*(sqrt(sqr(midval[2][3]) + sqr(midval[3][3])) + sqrt(sqr(midval[2][0]) + sqr(midval[3][0]))));
              split = !finite(err) || (err) > max*4*eps;
            }
          }

          // split the quad if the error is too large, otherwise produce two linear triangles
          if(split)
          {
            if(curved)
            {
              for (i = 0; i < 5; i++)
              {
                midval[0][i] = phx[idx[i]];
                midval[1][i] = phy[idx[i]];
              }
            }

            // obtain mid-edge and mid-element vertices
            int mid0 = get_vertex(iv0,  iv1,  midval[0][0], midval[1][0], xval[idx[0]], yval[idx[0]]);
            int mid1 = get_vertex(iv1,  iv2,  midval[0][1], midval[1][1], xval[idx[1]], yval[idx[1]]);
            int mid2 = get_vertex(iv2,  iv3,  midval[0][2], midval[1][2], xval[idx[2]], yval[idx[2]]);
            int mid3 = get_vertex(iv3,  iv0,  midval[0][3], midval[1][3], xval[idx[3]], yval[idx[3]]);
            int mid4 = get_vertex(mid0, mid2, midval[0][4], midval[1][4], xval[idx[4]], yval[idx[4]]);

            // recur to sub-elements
            fns[0]->push_transform(0); if(fns[1] != fns[0]) fns[1]->push_transform(0); process_quad(fns, iv0, mid0, mid4, mid3, level + 1, xval, yval, phx, phy, quad_indices[1], curved);  fns[0]->pop_transform(); if(fns[1] != fns[0]) fns[1]->pop_transform();
            fns[0]->push_transform(1); if(fns[1] != fns[0]) fns[1]->push_transform(1); process_quad(fns, mid0, iv1, mid1, mid4, level + 1, xval, yval, phx, phy, quad_indices[2], curved);  fns[0]->pop_transform(); if(fns[1] != fns[0]) fns[1]->pop_transform();
            fns[0]->push_transform(2); if(fns[1] != fns[0]) fns[1]->push_transform(2); process_quad(fns, mid4, mid1, iv2, mid2, level + 1, xval, yval, phx, phy, quad_indices[3], curved);  fns[0]->pop_transform(); if(fns[1] != fns[0]) fns[1]->pop_transform();
            fns[0]->push_transform(3); if(fns[1] != fns[0]) fns[1]->push_transform(3); process_quad(fns, mid3, mid4, mid2, iv3, level + 1, xval, yval, phx, phy, quad_indices[4], curved);  fns[0]->pop_transform(); if(fns[1] != fns[0]) fns[1]->pop_transform();
            return;
          }
        }

        // output two linear triangles,
        if(!flip)
        {
          add_triangle(iv3, iv0, iv1, fns[0]->get_active_element()->marker);
          add_triangle(iv1, iv2, iv3, fns[0]->get_active_element()->marker);
        }
        else
        {
          add_triangle(iv0, iv1, iv2, fns[0]->get_active_element()->marker);
          add_triangle(iv2, iv3, iv0, fns[0]->get_active_element()->marker);
        }
      }

      void Vectorizer::process_dash(int iv1, int iv2)
      {
        int mid = this->peek_vertex(iv1, iv2);
        if(mid != -1)
        {
          process_dash(iv1, mid);
          process_dash(mid, iv2);
        }
        else
          add_dash(iv1, iv2);
      }

      void Vectorizer::find_min_max()
      {
        // find min & max vertex values
        this->min_val =  1e100;
        this->max_val = -1e100;
        for (int i = 0; i < this->vertex_count; i++)
        {
          double mag = (verts[i][2]*verts[i][2] + verts[i][3]*verts[i][3]);
          if(finite(mag) && mag < this->min_val) this->min_val = mag;
          if(finite(mag) && mag > this->max_val) this->max_val = mag;
        }
        this->max_val = sqrt(this->max_val);
        this->min_val = sqrt(this->min_val);
      }

      void Vectorizer::set_curvature_epsilon(double curvature_epsilon)
      {
        this->curvature_epsilon = curvature_epsilon;
      }

      double Vectorizer::get_curvature_epsilon()
      {
        return this->curvature_epsilon;
      }

      void Vectorizer::process_solution(MeshFunction<double>* xsln, MeshFunction<double>* ysln, int xitem_orig, int yitem_orig, double eps)
      {
        // Important, sets the current caughtException to NULL.
        this->caughtException = NULL;

        // sanity check
        if(xsln == NULL || ysln == NULL) 
          throw Hermes::Exceptions::Exception("One of the solutions is NULL in Vectorizer:process_solution().");

        lock_data();
        this->tick();

        // initialization
        this->xitem = xitem_orig;
        this->yitem = yitem_orig;
        this->eps = eps;

        // estimate the required number of vertices and triangles
        // (based on the assumption that the linear mesh will be
        // about four-times finer than the original mesh).
        int nn = xsln->get_mesh()->get_num_elements() + ysln->get_mesh()->get_num_elements();

        this->vertex_size = std::max(100 * nn, std::max(this->vertex_size, 50000));
        this->triangle_size = std::max(150 * nn, std::max(this->triangle_size, 75000));
        this->edges_size = std::max(100 * nn, std::max(this->edges_size, 50000));
        //dashes_size = edges_size;

        vertex_count = 0;
        triangle_count = 0;
        edges_count = 0;
        dashes_count = 0;

        // reuse or allocate vertex, triangle and edge arrays
        verts = (double4*) malloc(sizeof(double4) * vertex_size);
				this->tris = (int3*) realloc(this->tris, sizeof(int3) * this->triangle_size);
				this->tri_markers = (int*) realloc(this->tri_markers, sizeof(int) * this->triangle_size);
				this->edges = (int2*) realloc(this->edges, sizeof(int2) * this->edges_size);
				this->edge_markers = (int*) realloc(this->edge_markers, sizeof(int) * this->edges_size);
				dashes = (int2*) malloc(sizeof(int2) * dashes_size);
        info = (int4*) malloc(sizeof(int4) * vertex_size);
        this->empty = false;

        // initialize the hash table
        hash_table = (int*) malloc(sizeof(int) * vertex_size);
        memset(hash_table, 0xff, sizeof(int) * vertex_size);

        // select the linearization quadrature
        Quad2D *old_quad_x, *old_quad_y;
        old_quad_x = xsln->get_quad_2d();
        old_quad_y = ysln->get_quad_2d();

        Hermes::vector<const Mesh*> meshes;
        meshes.push_back(xsln->get_mesh());
        meshes.push_back(ysln->get_mesh());

        // Parallelization
        MeshFunction<double>*** fns = new MeshFunction<double>**[Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads)];
        for(unsigned int i = 0; i < Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads); i++)
        {
          fns[i] = new MeshFunction<double>*[2];
          fns[i][0] = xsln->clone();
          fns[i][0]->set_refmap(new RefMap);
          fns[i][0]->set_quad_2d(&g_quad_lin);
          fns[i][1] = ysln->clone();
          fns[i][1]->set_refmap(new RefMap);
          fns[i][1]->set_quad_2d(&g_quad_lin);
        }

        Transformable*** trfs = new Transformable**[Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads)];
        for(unsigned int i = 0; i < Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads); i++)
        {
          trfs[i] = new Transformable*[2];
          trfs[i][0] = fns[i][0];
          trfs[i][1] = fns[i][1];
        }

        // get the component and desired value from item.
        if(xitem >= 0x40)
        {
          component_x = 1;
          xitem >>= 6;
        }
        while (!(xitem & 1))
        {
          xitem >>= 1;
          value_type_x++;
        }
        // get the component and desired value from item.
        if(yitem >= 0x40)
        {
          component_y = 1;
          yitem >>= 6;
        }
        while (!(yitem & 1))
        {
          yitem >>= 1;
          value_type_y++;
        }
        xitem = xitem_orig;
        yitem = yitem_orig;

        Traverse trav_masterMax(true);
        unsigned int num_states = trav_masterMax.get_num_states(meshes);

        trav_masterMax.begin(meshes.size(), &(meshes.front()));

        Traverse* trav = new Traverse[Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads)];

        for(unsigned int i = 0; i < Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads); i++)
        {
          trav[i].begin(meshes.size(), &(meshes.front()), trfs[i]);
          trav[i].stack = trav_masterMax.stack;
        }

        int state_i;

#define CHUNKSIZE 1
int num_threads_used = Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads);
#pragma omp parallel shared(trav_masterMax) private(state_i) num_threads(num_threads_used)
        {
#pragma omp for schedule(dynamic, CHUNKSIZE)
          for(state_i = 0; state_i < num_states; state_i++)
          {
            try
            {
              Traverse::State current_state;
#pragma omp critical(get_next_state)
              current_state = trav[omp_get_thread_num()].get_next_state(&trav_masterMax.top, &trav_masterMax.id);

              fns[omp_get_thread_num()][0]->set_quad_order(0, xitem);
              fns[omp_get_thread_num()][1]->set_quad_order(0, yitem);
              double* xval = fns[omp_get_thread_num()][0]->get_values(component_x, value_type_x);
              double* yval = fns[omp_get_thread_num()][1]->get_values(component_y, value_type_y);

              for (unsigned int i = 0; i < current_state.e[0]->get_nvert(); i++)
              {
                double fx = xval[i];
                double fy = yval[i];
#pragma omp critical(max)
                  if(fabs(sqrt(fx*fx + fy*fy)) > max)
                    max = fabs(sqrt(fx*fx + fy*fy));
              }
            }
            catch(Hermes::Exceptions::Exception& e)
            {
              if(this->caughtException == NULL)
                this->caughtException = e.clone();
            }
            catch(std::exception& e)
            {
              if(this->caughtException == NULL)
                this->caughtException = new Hermes::Exceptions::Exception(e.what());
            }
          }
        }

        trav_masterMax.finish();
        for(unsigned int i = 0; i < Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads); i++)
          trav[i].finish();
        delete [] trav;

        Traverse trav_master(true);
        num_states = trav_master.get_num_states(meshes);

        trav_master.begin(meshes.size(), &(meshes.front()));

        trav = new Traverse[Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads)];

        for(unsigned int i = 0; i < Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads); i++)
        {
          trav[i].begin(meshes.size(), &(meshes.front()), trfs[i]);
          trav[i].stack = trav_master.stack;
        }

        trav_master.begin(meshes.size(), &(meshes.front()));

        for(unsigned int i = 0; i < Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads); i++)
          trav[i].begin(meshes.size(), &(meshes.front()), trfs[i]);

#pragma omp parallel shared(trav_master) private(state_i) num_threads(num_threads_used)
        {
#pragma omp for schedule(dynamic, CHUNKSIZE)
          for(state_i = 0; state_i < num_states; state_i++)
          {
            try
            {
              Traverse::State current_state;

#pragma omp critical (get_next_state)
              current_state = trav[omp_get_thread_num()].get_next_state(&trav_master.top, &trav_master.id);

              fns[omp_get_thread_num()][0]->set_quad_order(0, xitem);
              fns[omp_get_thread_num()][1]->set_quad_order(0, yitem);
              double* xval = fns[omp_get_thread_num()][0]->get_values(component_x, value_type_x);
              double* yval = fns[omp_get_thread_num()][1]->get_values(component_y, value_type_y);

              double* x = fns[omp_get_thread_num()][0]->get_refmap()->get_phys_x(0);
              double* y = fns[omp_get_thread_num()][0]->get_refmap()->get_phys_y(0);

              int iv[4];
              for (unsigned int i = 0; i < current_state.e[0]->get_nvert(); i++)
              {
                double fx = xval[i];
                double fy = yval[i];

                iv[i] = this->get_vertex(-fns[omp_get_thread_num()][0]->get_active_element()->vn[i]->id, -fns[omp_get_thread_num()][0]->get_active_element()->vn[i]->id, x[i], y[i], fx, fy);
              }

              // recur to sub-elements
              if(current_state.e[0]->is_triangle())
                process_triangle(fns[omp_get_thread_num()], iv[0], iv[1], iv[2], 0, NULL, NULL, NULL, NULL, NULL, current_state.e[0]->is_curved());
              else
                process_quad(fns[omp_get_thread_num()], iv[0], iv[1], iv[2], iv[3], 0, NULL, NULL, NULL, NULL, NULL, current_state.e[0]->is_curved());

              for (unsigned int i = 0; i < current_state.e[0]->get_nvert(); i++)
                process_edge(iv[i], iv[current_state.e[0]->next_vert(i)], current_state.e[0]->en[i]->marker);
            }
            catch(Hermes::Exceptions::Exception& e)
            {
              if(this->caughtException == NULL)
                this->caughtException = e.clone();
            }
            catch(std::exception& e)
            {
              if(this->caughtException == NULL)
                this->caughtException = new Hermes::Exceptions::Exception(e.what());
            }
          }
        }

        trav_master.finish();
        for(unsigned int i = 0; i < Hermes2DApi.get_integral_param_value(Hermes::Hermes2D::numThreads); i++)
        {
          trav[i].finish();
          for(unsigned int j = 0; j < 2; j++)
            delete fns[i][j];
          delete [] fns[i];
          delete [] trfs[i];
        }
        delete [] fns;
        delete [] trfs;
        delete [] trav;

        // regularize the linear mesh
        for (int i = 0; i < this->triangle_count; i++)
        {
          int iv0 = tris[i][0], iv1 = tris[i][1], iv2 = tris[i][2];

          int mid0 = peek_vertex(iv0, iv1);
          int mid1 = peek_vertex(iv1, iv2);
          int mid2 = peek_vertex(iv2, iv0);
          if(mid0 >= 0 || mid1 >= 0 || mid2 >= 0)
          {
            this->del_slot = i;
            regularize_triangle(iv0, iv1, iv2, mid0, mid1, mid2, tri_markers[i]);
          }
        }

        find_min_max();

        //if(verbose_mode) print_hash_stats();
        unlock_data();

        // select old quadratrues
        xsln->set_quad_2d(old_quad_x);
        ysln->set_quad_2d(old_quad_y);

        // clean up
        ::free(this->hash_table);
        ::free(this->info);
      }

      void Vectorizer::free()
      {
        if(verts != NULL)
        {
          ::free(verts);
          verts = NULL;
        }
        if(dashes != NULL)
        {
          ::free(dashes);
          dashes = NULL;
        }

        LinearizerBase::free();
      }

      Vectorizer::~Vectorizer()
      {
        free();
      }

      void Vectorizer::calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const
      {
        if(verts == NULL)
          throw Exceptions::Exception("Cannot calculate AABB from NULL vertices");

        LinearizerBase::calc_aabb(&verts[0][0], &verts[0][1], sizeof(double4), this->vertex_count, min_x, max_x, min_y, max_y);
      }

      int2* Vectorizer::get_dashes()
      {
        return this->dashes;
      }
      int Vectorizer::get_num_dashes()
      {
        return this->dashes_size;
      }

      double4* Vectorizer::get_vertices()
      {
        return this->verts;
      }
      int Vectorizer::get_num_vertices()
      {
        return this->vertex_count;
      }

      int Vectorizer::add_vertex()
      {
        if(this->vertex_count >= this->vertex_size)
        {
          this->vertex_size *= 2;
          verts = (double4*) realloc(verts, sizeof(double4) * vertex_size);
          this->info = (int4*) realloc(info, sizeof(int4) * vertex_size);
          this->hash_table = (int*) realloc(hash_table, sizeof(int) * vertex_size);
          memset(this->hash_table + this->vertex_size / 2, 0xff, sizeof(int) * this->vertex_size / 2);
        }
        return this->vertex_count++;
      }

      void Vectorizer::add_dash(int iv1, int iv2)
      {
        if(this->dashes_count >= this->dashes_size)
        {
          this->dashes_size *= 2;
          dashes = (int2*) realloc(dashes, sizeof(int2) * dashes_size);
        }
        dashes[dashes_count][0] = iv1;
        dashes[dashes_count++][1] = iv2;
      }
    }
  }
}
