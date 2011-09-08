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

#include "linearizer.h"
#include "refmap.h"
#include "traverse.h"
#include "exact_solution.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      double3 lin_pts_0_tri[] =
      {
        { -1.0, -1.0, 0.0 },
        {  1.0, -1.0, 0.0 },
        { -1.0,  1.0, 0.0 }
      };

      double3 lin_pts_0_quad[] =
      {
        { -1.0, -1.0, 0.0 },
        {  1.0, -1.0, 0.0 },
        {  1.0,  1.0, 0.0 },
        { -1.0,  1.0, 0.0 }
      };

      double3 lin_pts_1_tri[12] =
      {
        {  0.0, -1.0, 0.0 }, // 0
        {  0.0,  0.0, 0.0 }, // 1
        { -1.0,  0.0, 0.0 }, // 2
        { -0.5, -1.0, 0.0 }, // 3
        { -0.5, -0.5, 0.0 }, // 4
        { -1.0, -0.5, 0.0 }, // 5
        {  0.5, -1.0, 0.0 }, // 6
        {  0.5, -0.5, 0.0 }, // 7
        {  0.0, -0.5, 0.0 }, // 8
        { -0.5,  0.0, 0.0 }, // 9
        { -0.5,  0.5, 0.0 }, // 10
        { -1.0,  0.5, 0.0 }  // 11
      };

      double3 lin_pts_1_quad[21] =
      {
        {  0.0, -1.0, 0.0 }, // 0
        {  1.0,  0.0, 0.0 }, // 1
        {  0.0,  1.0, 0.0 }, // 2
        { -1.0,  0.0, 0.0 }, // 3
        {  0.0,  0.0, 0.0 }, // 4
        { -0.5, -1.0, 0.0 }, // 5
        {  0.0, -0.5, 0.0 }, // 6
        { -0.5,  0.0, 0.0 }, // 7
        { -1.0, -0.5, 0.0 }, // 8
        { -0.5, -0.5, 0.0 }, // 9
        {  0.5, -1.0, 0.0 }, // 10
        {  1.0, -0.5, 0.0 }, // 11
        {  0.5,  0.0, 0.0 }, // 12
        {  0.5, -0.5, 0.0 }, // 13
        {  1.0,  0.5, 0.0 }, // 14
        {  0.5,  1.0, 0.0 }, // 15
        {  0.0,  0.5, 0.0 }, // 16
        {  0.5,  0.5, 0.0 }, // 17
        { -0.5,  1.0, 0.0 }, // 18
        { -1.0,  0.5, 0.0 }, // 19
        { -0.5,  0.5, 0.0 }  // 20
      };

      int quad_indices[9][5] =
      {
        { 0, 1, 2, 3, 4 },
        { 5, 6, 7, 8, 9 }, { 10, 11, 12, 6, 13 },
        { 12, 14, 15, 16, 17 }, { 7, 16, 18, 19, 20 },
        { 0, 11, 4, 8, 6 }, { 4, 14, 2, 19, 16 },
        { 5, 4, 18, 3, 7 }, { 10, 1, 15, 4, 12 }
      };

      int tri_indices[5][3] =
      {
        { 0, 1, 2 }, { 3, 4, 5 }, { 6, 7, 8 }, { 9, 10, 11 }, { 9, 4, 8 }
      };

      int lin_np_tri[2]   = { 3, 12 };
      int lin_np_quad[2]  = { 4, 21 };
      int* lin_np[2] = { lin_np_tri, lin_np_quad };

      double3*  lin_tables_tri[2]  = { lin_pts_0_tri, lin_pts_1_tri };
      double3*  lin_tables_quad[2] = { lin_pts_0_quad, lin_pts_1_quad };
      double3** lin_tables[2]      = { lin_tables_tri, lin_tables_quad };

      Linearizer::Linearizer(MeshFunction<double>* sln) : LinearizerBase(auto_max), sln(sln), dmult(1.0), component(0), value_type(0)
      {
        verts = NULL;
        xdisp = new ZeroSolution(sln->get_mesh());
        ydisp = new ZeroSolution(sln->get_mesh());
      }

      void Linearizer::process_triangle(int iv0, int iv1, int iv2, int level,
        double* val, double* phx, double* phy, int* idx)
      {
        double midval[3][3];

        if (level < LIN_MAX_LEVEL)
        {
          int i;
          if (!(level & 1))
          {
            // obtain solution values
            sln->set_quad_order(1, item);
            val = sln->get_values(component, value_type);
            if (auto_max)
              for (i = 0; i < lin_np_tri[1]; i++) 
              {
                double v = val[i];
                if (finite(v) && fabs(v) > max) max = fabs(v);
              }

            // obtain physical element coordinates
            RefMap* refmap = sln->get_refmap();
            phx = refmap->get_phys_x(1);
            phy = refmap->get_phys_y(1);

            xdisp->set_quad_order(1, H2D_FN_VAL);
            ydisp->set_quad_order(1, H2D_FN_VAL);
            double* dx = xdisp->get_fn_values();
            double* dy = ydisp->get_fn_values();
            for (i = 0; i < lin_np_tri[1]; i++) 
            {
              phx[i] += dmult*dx[i];
              phy[i] += dmult*dy[i];
            }
            idx = tri_indices[0];
          }

          // obtain linearized values and coordinates at the midpoints
          for (i = 0; i < 3; i++)
          {
            midval[i][0] = (verts[iv0][i] + verts[iv1][i])*0.5;
            midval[i][1] = (verts[iv1][i] + verts[iv2][i])*0.5;
            midval[i][2] = (verts[iv2][i] + verts[iv0][i])*0.5;
          };

          // determine whether or not to split the element
          bool split;
          if (eps >= 1.0)
          {
            // if eps > 1, the user wants a fixed number of refinements (no adaptivity)
            split = (level < eps);
          }
          else
          {
            if (!auto_max && fabs(verts[iv0][2]) > max && fabs(verts[iv1][2]) > max && fabs(verts[iv2][2]) > max)
            {
              // do not split if the whole triangle is above the specified maximum value
              split = false;
            }
            else
            {
              // calculate the approximate error of linearizing the normalized solution
              double err = fabs(val[idx[0]] - midval[2][0]) +
                fabs(val[idx[1]] - midval[2][1]) +
                fabs(val[idx[2]] - midval[2][2]);
              split = !finite(err) || err > max*3*eps;
            }

            // do the same for the curvature
            for (i = 0; i < 3; i++)
              if (sqr(phx[idx[i]] - midval[0][i]) + sqr(phy[idx[i]] - midval[1][i]) > sqr(cmax*1.5e-3))
              { 
                split = true;
                break; 
              }

            // do extra tests at level 0, so as not to miss some functions with zero error at edge midpoints
            if (level == 0 && !split)
            {
              split = (fabs(val[8] - 0.5*(midval[2][0] + midval[2][1])) +
                fabs(val[9] - 0.5*(midval[2][1] + midval[2][2])) +
                fabs(val[4] - 0.5*(midval[2][2] + midval[2][0]))) > max*3*eps;
            }
          }

          // split the triangle if the error is too large, otherwise produce a linear triangle
          if (split)
          {
            for (i = 0; i < 3; i++) 
            {
              midval[0][i] = phx[idx[i]];
              midval[1][i] = phy[idx[i]];
            }

            // obtain mid-edge vertices
            int mid0 = get_vertex(iv0, iv1, midval[0][0], midval[1][0], val[idx[0]]);
            int mid1 = get_vertex(iv1, iv2, midval[0][1], midval[1][1], val[idx[1]]);
            int mid2 = get_vertex(iv2, iv0, midval[0][2], midval[1][2], val[idx[2]]);

            // recur to sub-elements
            sln->push_transform(0);
            process_triangle(iv0, mid0, mid2,  level+1, val, phx, phy, tri_indices[1]);
            sln->pop_transform();

            sln->push_transform(1);
            process_triangle(mid0, iv1, mid1,  level+1, val, phx, phy, tri_indices[2]);
            sln->pop_transform();

            sln->push_transform(2);
            process_triangle(mid2, mid1, iv2,  level+1, val, phx, phy, tri_indices[3]);
            sln->pop_transform();

            sln->push_transform(3);
            process_triangle(mid1, mid2, mid0, level+1, val, phx, phy, tri_indices[4]);
            sln->pop_transform();
            return;
          }
        }

        // no splitting: output a linear triangle
        add_triangle(iv0, iv1, iv2);
      }

      void Linearizer::process_quad(int iv0, int iv1, int iv2, int iv3, int level,
        double* val, double* phx, double* phy, int* idx)
      {
        double midval[3][5];

        // try not to split through the vertex with the largest value
        int a = (verts[iv0][2] > verts[iv1][2]) ? iv0 : iv1;
        int b = (verts[iv2][2] > verts[iv3][2]) ? iv2 : iv3;
        a = (verts[a][2] > verts[b][2]) ? a : b;
        int flip = (a == iv1 || a == iv3) ? 1 : 0;

        if (level < LIN_MAX_LEVEL)
        {
          int i;
          if (!(level & 1)) // this is an optimization: do the following only every other time
          {
            // obtain solution values
            sln->set_quad_order(1, item);
            val = sln->get_values(component, value_type);
            if (auto_max)
              for (i = 0; i < lin_np_quad[1]; i++) 
              {
                double v = val[i];
                if (finite(v) && fabs(v) > max) max = fabs(v);
              }

            // This is just to make some sense.
            if(fabs(max) < 1E-10)
              max = 1E-10;

            // obtain physical element coordinates
              
            RefMap* refmap = sln->get_refmap();
            phx = refmap->get_phys_x(1);
            phy = refmap->get_phys_y(1);

            xdisp->set_quad_order(1, H2D_FN_VAL);
            ydisp->set_quad_order(1, H2D_FN_VAL);
            double* dx = xdisp->get_fn_values();
            double* dy = ydisp->get_fn_values();
            for (i = 0; i < lin_np_quad[1]; i++) 
            {
              phx[i] += dmult*dx[i];
              phy[i] += dmult*dy[i];
            }
            idx = quad_indices[0];
          }

          // obtain linearized values and coordinates at the midpoints
          for (i = 0; i < 3; i++)
          {
            midval[i][0] = (verts[iv0][i] + verts[iv1][i]) * 0.5;
            midval[i][1] = (verts[iv1][i] + verts[iv2][i]) * 0.5;
            midval[i][2] = (verts[iv2][i] + verts[iv3][i]) * 0.5;
            midval[i][3] = (verts[iv3][i] + verts[iv0][i]) * 0.5;
            midval[i][4] = (midval[i][0]  + midval[i][2])  * 0.5;
          };

          // the value of the middle point is not the average of the four vertex values, since quad == 2 triangles
          midval[2][4] = flip ? (verts[iv0][2] + verts[iv2][2]) * 0.5 : (verts[iv1][2] + verts[iv3][2]) * 0.5;

          // determine whether or not to split the element
          int split;
          if (eps >= 1.0)
          {
            // if eps > 1, the user wants a fixed number of refinements (no adaptivity)
            split = (level < eps) ? 3 : 0;
          }
          else
          {
            if (!auto_max && fabs(verts[iv0][2]) > max && fabs(verts[iv1][2]) > max
              && fabs(verts[iv2][2]) > max && fabs(verts[iv3][2]) > max)
            {
              // do not split if the whole quad is above the specified maximum value
              split = 0;
            }
            else
            {
              // calculate the approximate error of linearizing the normalized solution
              double herr = fabs(val[idx[1]] - midval[2][1]) + fabs(val[idx[3]] - midval[2][3]);
              double verr = fabs(val[idx[0]] - midval[2][0]) + fabs(val[idx[2]] - midval[2][2]);
              double err  = fabs(val[idx[4]] - midval[2][4]) + herr + verr;
              split = (!finite(err) || err > max*4*eps) ? 3 : 0;

              // decide whether to split horizontally or vertically only
              if (level > 0 && split) 
              {
                if (herr > 5*verr)
                  split = 1; // h-split
                else if (verr > 5*herr)
                  split = 2; // v-split
              }
            }

            // also decide whether to split because of the curvature
            if (split != 3)
            {
              double cm2 = sqr(cmax*5e-4);
              if (sqr(phx[idx[1]] - midval[0][1]) + sqr(phy[idx[1]] - midval[1][1]) > cm2 ||
                sqr(phx[idx[3]] - midval[0][3]) + sqr(phy[idx[3]] - midval[1][3]) > cm2) split |= 1;
              if (sqr(phx[idx[0]] - midval[0][0]) + sqr(phy[idx[0]] - midval[1][0]) > cm2 ||
                sqr(phx[idx[2]] - midval[0][2]) + sqr(phy[idx[2]] - midval[1][2]) > cm2) split |= 2;

              /*for (i = 0; i < 5; i++)
              if (sqr(phx[idx[i]] - midval[0][i]) + sqr(phy[idx[i]] - midval[1][i]) > sqr(cmax*1e-3))
              { split = 1; break; }*/
            }

            // do extra tests at level 0, so as not to miss some functions with zero error at edge midpoints
            if (level == 0 && !split)
            {
              split = ((fabs(val[13] - 0.5*(midval[2][0] + midval[2][1])) +
                fabs(val[17] - 0.5*(midval[2][1] + midval[2][2])) +
                fabs(val[20] - 0.5*(midval[2][2] + midval[2][3])) +
                fabs(val[9]  - 0.5*(midval[2][3] + midval[2][0]))) > max*4*eps) ? 3 : 0;
            }
          }

          // split the quad if the error is too large, otherwise produce two linear triangles
          if (split)
          {
            for (i = 0; i < 5; i++) 
            {
              midval[0][i] = phx[idx[i]];
              midval[1][i] = phy[idx[i]];
            }

            // obtain mid-edge and mid-element vertices
            int mid0, mid1, mid2, mid3, mid4;
            if (split != 1) mid0 = get_vertex(iv0,  iv1,  midval[0][0], midval[1][0], val[idx[0]]);
            if (split != 2) mid1 = get_vertex(iv1,  iv2,  midval[0][1], midval[1][1], val[idx[1]]);
            if (split != 1) mid2 = get_vertex(iv2,  iv3,  midval[0][2], midval[1][2], val[idx[2]]);
            if (split != 2) mid3 = get_vertex(iv3,  iv0,  midval[0][3], midval[1][3], val[idx[3]]);
            if (split == 3) mid4 = get_vertex(mid0, mid2, midval[0][4], midval[1][4], val[idx[4]]);

            // recur to sub-elements
            if (split == 3)
            {
              sln->push_transform(0);
              process_quad(iv0, mid0, mid4, mid3, level+1, val, phx, phy, quad_indices[1]);
              sln->pop_transform();

              sln->push_transform(1);
              process_quad(mid0, iv1, mid1, mid4, level+1, val, phx, phy, quad_indices[2]);
              sln->pop_transform();

              sln->push_transform(2);
              process_quad(mid4, mid1, iv2, mid2, level+1, val, phx, phy, quad_indices[3]);
              sln->pop_transform();

              sln->push_transform(3);
              process_quad(mid3, mid4, mid2, iv3, level+1, val, phx, phy, quad_indices[4]);
              sln->pop_transform();
            }
            else 
              if (split == 1) // h-split
              {
                sln->push_transform(4);
                process_quad(iv0, iv1, mid1, mid3, level+1, val, phx, phy, quad_indices[5]);
                sln->pop_transform();

                sln->push_transform(5);
                process_quad(mid3, mid1, iv2, iv3, level+1, val, phx, phy, quad_indices[6]);
                sln->pop_transform();
              }
              else // v-split
              {
                sln->push_transform(6);
                process_quad(iv0, mid0, mid2, iv3, level+1, val, phx, phy, quad_indices[7]);
                sln->pop_transform();

                sln->push_transform(7);
                process_quad(mid0, iv1, iv2, mid2, level+1, val, phx, phy, quad_indices[8]);
                sln->pop_transform();
              }
            return;
          }
        }

        // output two linear triangles,
        if (!flip)
        {
          add_triangle(iv3, iv0, iv1);
          add_triangle(iv1, iv2, iv3);
        }
        else
        {
          add_triangle(iv0, iv1, iv2);
          add_triangle(iv2, iv3, iv0);
        }
      }
      
      void Linearizer::regularize_triangle(int iv0, int iv1, int iv2, int mid0, int mid1, int mid2)
      {
        // count the number of hanging mid-edge vertices
        int n = 0;
        if (mid0 >= 0) n++;
        if (mid1 >= 0) n++;
        if (mid2 >= 0) n++;
        if (n == 3)
        {
          // three hanging vertices: split into four triangles
          regularize_triangle(iv0, mid0, mid2, peek_vertex(iv0, mid0), -1, peek_vertex(mid2, iv0));
          regularize_triangle(mid0, iv1, mid1, peek_vertex(mid0, iv1), peek_vertex(iv1, mid1), -1);
          regularize_triangle(mid2, mid1, iv2, -1, peek_vertex(mid1, iv2), peek_vertex(iv2, mid2));
          regularize_triangle(mid0, mid1, mid2, -1, -1, -1);
        }
        else if (n == 2)
        {
          // two hanging vertices: split into three triangles
          if (mid0 < 0)
          {
            regularize_triangle(iv0, iv1, mid1, peek_vertex(iv0, iv1), peek_vertex(iv1, mid1), -1);
            regularize_triangle(mid2, iv0, mid1, peek_vertex(mid2, iv0), -1, -1);
            regularize_triangle(mid2, mid1, iv2, -1, peek_vertex(mid1, iv2), peek_vertex(iv2, mid2));
          }
          else if (mid1 < 0)
          {
            regularize_triangle(iv1, iv2, mid2, peek_vertex(iv1, iv2), peek_vertex(iv2, mid2), -1);
            regularize_triangle(mid0, iv1, mid2, peek_vertex(mid0, iv1), -1, -1);
            regularize_triangle(mid0, mid2, iv0, -1, peek_vertex(mid2, iv0), peek_vertex(iv0, mid0));
          }
          else
          {
            regularize_triangle(iv2, iv0, mid0, peek_vertex(iv2, iv0), peek_vertex(iv0, mid0), -1);
            regularize_triangle(mid1, iv2, mid0, peek_vertex(mid1, iv2), -1, -1);
            regularize_triangle(mid1, mid0, iv1, -1, peek_vertex(mid0, iv1), peek_vertex(iv1, mid1));
          }
        }
        else if (n == 1)
        {
          // one hanging vertex: split into two triangles
          if (mid0 >= 0)
          {
            regularize_triangle(iv0, mid0, iv2, peek_vertex(iv0, mid0), -1, peek_vertex(iv2, iv0));
            regularize_triangle(mid0, iv1, iv2, peek_vertex(mid0, iv1), peek_vertex(iv1, iv2), -1);
          }
          else if (mid1 >= 0)
          {
            regularize_triangle(iv1, mid1, iv0, peek_vertex(iv1, mid1), -1, peek_vertex(iv0, iv1));
            regularize_triangle(mid1, iv2, iv0, peek_vertex(mid1, iv2), peek_vertex(iv2, iv0), -1);
          }
          else
          {
            regularize_triangle(iv2, mid2, iv1, peek_vertex(iv2, mid2), -1, peek_vertex(iv1, iv2));
            regularize_triangle(mid2, iv0, iv1, peek_vertex(mid2, iv0), peek_vertex(iv0, iv1), -1);
          }
        }
        else
        {
          // no hanging vertices: produce a single triangle
          add_triangle(iv0, iv1, iv2);
        }
      }
      
      void Linearizer::set_displacement(MeshFunction<double>* xdisp, MeshFunction<double>* ydisp, double dmult)
      {
        if(xdisp != NULL)
          this->xdisp = xdisp;
        if(ydisp != NULL)
          this->ydisp = ydisp;
        this->dmult = dmult;
      }

      void Linearizer::process_solution(int item, double eps)
      {
        lock_data();
        Hermes::TimePeriod time_period;

        // initialization
        this->item = item;
        this->eps = eps;

        // get the component and desired value from item.
        if (item >= 0x40)
        {
          component = 1; 
          item >>= 6; 
        }
        while (!(item & 1)) 
        { 
          item >>= 1;
          value_type++; 
        }

        int nn = this->sln->get_mesh()->get_num_elements();

        vertex_size = std::max(32 * nn, 10000);
        triangle_size = std::max(64 * nn, 20000);
        edges_size = std::max(24 * nn, 7500);

        vertex_count=0;
        triangle_count=0;
        edges_count=0;
        // reuse or allocate vertex, triangle and edge arrays
        verts = (double3*) realloc(verts, sizeof(double3) * vertex_size);
        tris = (int3*) realloc(tris, sizeof(int3) * triangle_size);
        edges = (int3*) realloc(edges, sizeof(int3) * edges_size);
        info = (int4*) malloc(sizeof(int4) * vertex_size);

        // initialize the hash table
        hash_table = (int*) malloc(sizeof(int) * vertex_size);
        memset(hash_table, 0xff, sizeof(int) * vertex_size);

        // select the linearization quadratures
        Quad2D *old_quad, *old_quad_x = NULL, *old_quad_y = NULL;
        old_quad = sln->get_quad_2d();
        sln->set_quad_2d(&g_quad_lin);
        old_quad_x = xdisp->get_quad_2d();
        xdisp->set_quad_2d(&g_quad_lin);
        old_quad_y = ydisp->get_quad_2d();
        ydisp->set_quad_2d(&g_quad_lin); 

        // create all top-level vertices (corresponding to vertex nodes), with
        // all parent-son relations preserved; this is necessary for regularization to
        // work on irregular meshes
        nn = this->sln->get_mesh()->get_max_node_id();

        // obtain the solution in vertices, estimate the maximum solution value
        // meshes.
        Mesh** meshes = new Mesh*[3];

        meshes[0] = sln->get_mesh();
        meshes[1] = xdisp->get_mesh();
        meshes[2] = ydisp->get_mesh();

        // functions.
        Transformable** trfs = new Transformable*[3];
        trfs[0] = sln;
        trfs[1] = xdisp;
        trfs[2] = ydisp;
        
        // Init multi-mesh traversal.
        Traverse trav;
        trav.begin(3, meshes, trfs);

        // Loop through all elements.
        Element **e;
        // Loop through all elements.
        while ((e = trav.get_next_state(NULL, NULL)) != NULL) 
        {
          sln->set_quad_order(0, item);
          double* val = sln->get_values(component, value_type);
          if (val == NULL) 
            error("Item not defined in the solution.");

          xdisp->set_quad_order(0, H2D_FN_VAL);
          ydisp->set_quad_order(0, H2D_FN_VAL);
          
          double *dx = xdisp->get_fn_values();
          double *dy = ydisp->get_fn_values();

          int iv[4];
          for (unsigned int i = 0; i < e[0]->nvert; i++)
          {
            double f = val[i];
            if (this->auto_max && finite(f) && fabs(f) > this->max) 
              this->max = fabs(f);

            double x_disp = sln->get_refmap()->get_phys_x(0)[i];
            double y_disp = sln->get_refmap()->get_phys_y(0)[i];

            x_disp += dmult * dx[i];
            y_disp += dmult * dy[i];

            iv[i] = this->get_vertex(-rand(), -rand(), x_disp, y_disp, f);
          }

          // we won't bother calculating physical coordinates from the refmap if this is not a curved element
          this->curved = e[0]->is_curved();
          cmax = e[0]->get_diameter();

          // recur to sub-elements
          if (e[0]->is_triangle())
            process_triangle(iv[0], iv[1], iv[2], 0, NULL, NULL, NULL, NULL);
          else
            process_quad(iv[0], iv[1], iv[2], iv[3], 0, NULL, NULL, NULL, NULL);

          for (unsigned int i = 0; i < e[0]->nvert; i++)
            process_edge(iv[i], iv[e[0]->next_vert(i)], e[0]->en[i]->marker);
        }

        find_min_max();

        this->unlock_data();

        // select old quadratrues
        sln->set_quad_2d(old_quad);
        xdisp->set_quad_2d(old_quad_x);
        ydisp->set_quad_2d(old_quad_y);

        // clean up
        ::free(hash_table);
        ::free(info);
      }
      
      void Linearizer::find_min_max()
      {
        // find min & max vertex values
        this->min_val =  1e100;
        this->max_val = -1e100;
        for (int i = 0; i < this->vertex_count; i++)
        {
          if (finite(verts[i][2]) && verts[i][2] < min_val) min_val = verts[i][2];
          if (finite(verts[i][2]) && verts[i][2] > max_val) max_val = verts[i][2];
        }
      }

      int Linearizer::get_vertex(int p1, int p2, double x, double y, double value) 
      {
        // search for an existing vertex
        if (p1 > p2) std::swap(p1, p2);
        int index = this->hash(p1, p2);
        int i = this->hash_table[index];
        while (i >= 0)
        {
          if (this->info[i][0] == p1 && this->info[i][1] == p2 &&
            (value == verts[i][2] || fabs(value - verts[i][2]) < this->max*1e-4)) return i;
          // note that we won't return a vertex with a different value than the required one;
          // this takes care for discontinuities in the solution, where more vertices
          // with different values will be created
          i = info[i][2];
        }

        // if not found, create a new one
        i = add_vertex();
        verts[i][0] = x;
        verts[i][1] = y;
        verts[i][2] = value;
        this->info[i][0] = p1;
        this->info[i][1] = p2;
        this->info[i][2] = hash_table[index];
        this->hash_table[index] = i;
        return i;
      }

      int Linearizer::add_vertex()
      {
        if (this->vertex_count >= this->vertex_size)
        {
          this->vertex_size *= 2;
          verts = (double3*) realloc(verts, sizeof(double3) * vertex_size);
          this->info = (int4*) realloc(info, sizeof(int4) * vertex_size);
          this->hash_table = (int*) realloc(hash_table, sizeof(int) * vertex_size);
          memset(this->hash_table + this->vertex_size / 2, 0xff, sizeof(int) * this->vertex_size / 2);
        }
        return this->vertex_count++;
      }

      int Linearizer::get_top_vertex(int id, double value)
      {
        if (fabs(value - verts[id][2]) < max*1e-24) return id;
        return get_vertex(-rand(), -rand(), verts[id][0], verts[id][1], value);
      }

      void Linearizer::free()
      {
        if (verts != NULL) 
        {
          ::free(verts);
          verts = NULL;
        }
      }

      Linearizer::~Linearizer()
      {
        free();
      }

      void Linearizer::save_solution_vtk(const char* filename, const char *quantity_name,
        bool mode_3D, int item, double eps)
      {
        process_solution(item, eps);

        FILE* f = fopen(filename, "wb");
        if (f == NULL) error("Could not open %s for writing.", filename);
        lock_data();

        // Output header for vertices.
        fprintf(f, "# vtk DataFile Version 2.0\n");
        fprintf(f, "\n");
        fprintf(f, "ASCII\n\n");
        fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

        // Output vertices.
        fprintf(f, "POINTS %d %s\n", this->vertex_count, "float");
        for (int i=0; i < this->vertex_count; i++) 
        {
          if (mode_3D == true) fprintf(f, "%g %g %g\n", this->verts[i][0], this->verts[i][1], this->verts[i][2]);
          else fprintf(f, "%g %g %g\n", this->verts[i][0], this->verts[i][1], 0.0);
        }

        // Output elements.
        fprintf(f, "\n");
        fprintf(f, "CELLS %d %d\n", this->triangle_count, 4 * this->triangle_count);
        for (int i=0; i < this->triangle_count; i++) 
        {
          fprintf(f, "3 %d %d %d\n", this->tris[i][0], this->tris[i][1], this->tris[i][2]);
        }

        // Output cell types.
        fprintf(f, "\n");
        fprintf(f, "CELL_TYPES %d\n", this->triangle_count);
        for (int i=0; i < this->triangle_count; i++) 
        {
          fprintf(f, "5\n");    // The "5" means triangle in VTK.
        }

        // This outputs double solution values.
        fprintf(f, "\n");
        fprintf(f, "POINT_DATA %d\n", this->vertex_count);
        fprintf(f, "SCALARS %s %s %d\n", quantity_name, "float", 1);
        fprintf(f, "LOOKUP_TABLE %s\n", "default");
        for (int i=0; i < this->vertex_count; i++) 
        {
          fprintf(f, "%g\n", this->verts[i][2]);
        }

        unlock_data();
        fclose(f);
      }

      void Linearizer::calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const 
      {
        assert_msg(verts != NULL, "Cannot calculate AABB from NULL vertices");
        calc_aabb(&verts[0][0], &verts[0][1], sizeof(double3), vertex_count, min_x, max_x, min_y, max_y);
      }

      double3* Linearizer::get_vertices()
      {
        return this->verts;
      }
      int Linearizer::get_num_vertices()
      {
        return this->vertex_count;
      }
    }
  }
}
