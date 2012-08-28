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
#include "api2d.h"

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

      Linearizer::Linearizer(bool auto_max) : LinearizerBase(auto_max), dmult(1.0), component(0), value_type(0)
      {
        verts = NULL;
        xdisp = NULL;
        user_xdisp = false;
        ydisp = NULL;
        user_ydisp = false;
        trisContours = NULL;
      }

      void Linearizer::process_triangle(MeshFunction<double>** fns, int iv0, int iv1, int iv2, int level,
        double* val, double* phx, double* phy, int* idx, bool curved)
      {
        double midval[3][3];

        if(level < LIN_MAX_LEVEL)
        {
          int i;
          if(!(level & 1))
          {
            // obtain solution values
            fns[0]->set_quad_order(1, item);
            val = fns[0]->get_values(component, value_type);
            if(auto_max)
              for (i = 0; i < lin_np_tri[1]; i++)
              {
                double v = val[i];
#pragma omp critical(max)
                if(finite(v) && fabs(v) > max)
                  max = fabs(v);
              }
            idx = tri_indices[0];

            if(curved)
            {
              // obtain physical element coordinates
              RefMap* refmap = fns[0]->get_refmap();
              phx = refmap->get_phys_x(1);
              phy = refmap->get_phys_y(1);

              double* dx = NULL;
              double* dy = NULL;
              if(this->xdisp != NULL)
              {
                fns[1]->set_quad_order(1, H2D_FN_VAL);
                dx = fns[1]->get_fn_values();
              }
              if(this->ydisp != NULL)
              {
                fns[2]->set_quad_order(1, H2D_FN_VAL);
                dy = fns[2]->get_fn_values();
              }
              for (i = 0; i < lin_np_tri[1]; i++)
              {
                if(this->xdisp != NULL)
                  phx[i] += dmult*dx[i];
                if(this->ydisp != NULL)
                  phy[i] += dmult*dy[i];
              }

            }
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
          if(eps >= 1.0)
          {
            // if eps > 1, the user wants a fixed number of refinements (no adaptivity)
            split = (level + 5 < eps);
          }
          else
          {
            if(!auto_max && fabs(verts[iv0][2]) > max && fabs(verts[iv1][2]) > max && fabs(verts[iv2][2]) > max)
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
            if(curved)
            {
              for (i = 0; i < 3; i++)
                if(sqr(phx[idx[i]] - midval[0][i]) + sqr(phy[idx[i]] - midval[1][i]) > sqr(cmax*1.5e-3))
                {
                  split = true;
                  break;
                }
            }

            // do extra tests at level 0, so as not to miss some functions with zero error at edge midpoints
            if(level == 0 && !split)
            {
              split = (fabs(val[8] - 0.5*(midval[2][0] + midval[2][1])) +
                fabs(val[9] - 0.5*(midval[2][1] + midval[2][2])) +
                fabs(val[4] - 0.5*(midval[2][2] + midval[2][0]))) > max*3*eps;
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
            int mid0 = get_vertex(iv0, iv1, midval[0][0], midval[1][0], val[idx[0]]);
            int mid1 = get_vertex(iv1, iv2, midval[0][1], midval[1][1], val[idx[1]]);
            int mid2 = get_vertex(iv2, iv0, midval[0][2], midval[1][2], val[idx[2]]);

            // recur to sub-elements
            fns[0]->push_transform(0);
            if(this->xdisp != NULL)
              if(fns[1] != fns[0])
                fns[1]->push_transform(0);
            if(this->ydisp != NULL)
              if(fns[2] != fns[1])
                fns[2]->push_transform(0);
            process_triangle(fns, iv0, mid0, mid2,  level + 1, val, phx, phy, tri_indices[1], curved);
            fns[0]->pop_transform();
            if(this->xdisp != NULL)
              if(fns[1] != fns[0])
                fns[1]->pop_transform();
            if(this->ydisp != NULL)
              if(fns[2] != fns[1])
                fns[2]->pop_transform();

            fns[0]->push_transform(1);
            if(this->xdisp != NULL)
              if(fns[1] != fns[0])
                fns[1]->push_transform(1);
            if(this->ydisp != NULL)
              if(fns[2] != fns[1])
                fns[2]->push_transform(1);
            process_triangle(fns, mid0, iv1, mid1,  level + 1, val, phx, phy, tri_indices[2], curved);
            fns[0]->pop_transform();
            if(this->xdisp != NULL)
              if(fns[1] != fns[0])
                fns[1]->pop_transform();
            if(this->ydisp != NULL)
              if(fns[2] != fns[1])
                fns[2]->pop_transform();

            fns[0]->push_transform(2);
            if(this->xdisp != NULL)
              if(fns[1] != fns[0])
                fns[1]->push_transform(2);
            if(this->ydisp != NULL)
              if(fns[2] != fns[1])
                fns[2]->push_transform(2);
            process_triangle(fns, mid2, mid1, iv2,  level + 1, val, phx, phy, tri_indices[3], curved);
            fns[0]->pop_transform();
            if(this->xdisp != NULL)
              if(fns[1] != fns[0])
                fns[1]->pop_transform();
            if(this->ydisp != NULL)
              if(fns[2] != fns[1])
                fns[2]->pop_transform();

            fns[0]->push_transform(3);
            if(this->xdisp != NULL)
              if(fns[1] != fns[0])
                fns[1]->push_transform(3);
            if(this->ydisp != NULL)
              if(fns[2] != fns[1])
                fns[2]->push_transform(3);
            process_triangle(fns, mid1, mid2, mid0, level + 1, val, phx, phy, tri_indices[4], curved);
            fns[0]->pop_transform();
            if(this->xdisp != NULL)
              if(fns[1] != fns[0])
                fns[1]->pop_transform();
            if(this->ydisp != NULL)
              if(fns[2] != fns[1])
                fns[2]->pop_transform();
            return;
          }
        }

        // no splitting: output a linear triangle
        add_triangle(iv0, iv1, iv2);
      }

      void Linearizer::process_quad(MeshFunction<double>** fns, int iv0, int iv1, int iv2, int iv3, int level,
        double* val, double* phx, double* phy, int* idx, bool curved)
      {
        double midval[3][5];

        // try not to split through the vertex with the largest value
        int a = (verts[iv0][2] > verts[iv1][2]) ? iv0 : iv1;
        int b = (verts[iv2][2] > verts[iv3][2]) ? iv2 : iv3;
        a = (verts[a][2] > verts[b][2]) ? a : b;
        int flip = (a == iv1 || a == iv3) ? 1 : 0;

        if(level < LIN_MAX_LEVEL)
        {
          int i;
          if(!(level & 1)) // this is an optimization: do the following only every other time
          {
            // obtain solution values
            fns[0]->set_quad_order(1, item);
            val = fns[0]->get_values(component, value_type);
            if(auto_max)
              for (i = 0; i < lin_np_quad[1]; i++)
              {
                double v = val[i];
                if(finite(v) && fabs(v) > max)
#pragma omp critical(max)
                  if(finite(v) && fabs(v) > max)
                    max = fabs(v);
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

              double* dx = NULL;
              double* dy = NULL;

              if(this->xdisp != NULL)
                fns[1]->set_quad_order(1, H2D_FN_VAL);
              if(this->ydisp != NULL)
                fns[2]->set_quad_order(1, H2D_FN_VAL);
              if(this->xdisp != NULL)
                dx = fns[1]->get_fn_values();
              if(this->ydisp != NULL)
                dy = fns[2]->get_fn_values();
              for (i = 0; i < lin_np_quad[1]; i++)
              {
                if(this->xdisp != NULL)
                  phx[i] += dmult*dx[i];
                if(this->ydisp != NULL)
                  phy[i] += dmult*dy[i];
              }
            }
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
          if(eps >= 1.0)
          {
            // if eps > 1, the user wants a fixed number of refinements (no adaptivity)
            split = (level < eps) ? 3 : 0;
          }
          else
          {
            if(!auto_max && fabs(verts[iv0][2]) > max && fabs(verts[iv1][2]) > max
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
              if(level > 0 && split)
              {
                if(herr > 5*verr)
                  split = 1; // h-split
                else if(verr > 5*herr)
                  split = 2; // v-split
              }
            }

            // also decide whether to split because of the curvature
            if(split != 3 && curved)
            {
              double cm2 = sqr(this->cmax*1e-3);
              if(sqr(phx[idx[1]] - midval[0][1]) + sqr(phy[idx[1]] - midval[1][1]) > cm2 ||
                sqr(phx[idx[3]] - midval[0][3]) + sqr(phy[idx[3]] - midval[1][3]) > cm2) split |= 1;
              if(sqr(phx[idx[0]] - midval[0][0]) + sqr(phy[idx[0]] - midval[1][0]) > cm2 ||
                sqr(phx[idx[2]] - midval[0][2]) + sqr(phy[idx[2]] - midval[1][2]) > cm2) split |= 2;

              /*for (i = 0; i < 5; i++)
              if(sqr(phx[idx[i]] - midval[0][i]) + sqr(phy[idx[i]] - midval[1][i]) > sqr(cmax*1e-3))
              { split = 1; break; }*/
            }

            // do extra tests at level 0, so as not to miss some functions with zero error at edge midpoints
            if(level == 0 && !split)
            {
              split = ((fabs(val[13] - 0.5*(midval[2][0] + midval[2][1])) +
                fabs(val[17] - 0.5*(midval[2][1] + midval[2][2])) +
                fabs(val[20] - 0.5*(midval[2][2] + midval[2][3])) +
                fabs(val[9]  - 0.5*(midval[2][3] + midval[2][0]))) > max*4*eps) ? 3 : 0;
            }
          }

          // split the quad if the error is too large, otherwise produce two linear triangles
          if(split)
          {
            if(curved)
              for (i = 0; i < 5; i++)
              {
                midval[0][i] = phx[idx[i]];
                midval[1][i] = phy[idx[i]];
              }

            // obtain mid-edge and mid-element vertices
            int mid0, mid1, mid2, mid3, mid4;
            if(split != 1) mid0 = get_vertex(iv0,  iv1,  midval[0][0], midval[1][0], val[idx[0]]);
            if(split != 2) mid1 = get_vertex(iv1,  iv2,  midval[0][1], midval[1][1], val[idx[1]]);
            if(split != 1) mid2 = get_vertex(iv2,  iv3,  midval[0][2], midval[1][2], val[idx[2]]);
            if(split != 2) mid3 = get_vertex(iv3,  iv0,  midval[0][3], midval[1][3], val[idx[3]]);
            if(split == 3) mid4 = get_vertex(mid0, mid2, midval[0][4], midval[1][4], val[idx[4]]);

            // recur to sub-elements
            if(split == 3)
            {
              fns[0]->push_transform(0);
              if(this->xdisp != NULL)
                if(fns[1] != fns[0])
                  fns[1]->push_transform(0);
              if(this->ydisp != NULL)
                if(fns[2] != fns[1])
                  fns[2]->push_transform(0);
              process_quad(fns, iv0, mid0, mid4, mid3, level + 1, val, phx, phy, quad_indices[1], curved);
              fns[0]->pop_transform();
              if(this->xdisp != NULL)
                if(fns[1] != fns[0])
                  fns[1]->pop_transform();
              if(this->ydisp != NULL)
                if(fns[2] != fns[1])
                  fns[2]->pop_transform();

              fns[0]->push_transform(1);
              if(this->xdisp != NULL)
                if(fns[1] != fns[0])
                  fns[1]->push_transform(1);
              if(this->ydisp != NULL)
                if(fns[2] != fns[1])
                  fns[2]->push_transform(1);
              process_quad(fns, mid0, iv1, mid1, mid4, level + 1, val, phx, phy, quad_indices[2], curved);
              fns[0]->pop_transform();
              if(this->xdisp != NULL)
                if(fns[1] != fns[0])
                  fns[1]->pop_transform();
              if(this->ydisp != NULL)
                if(fns[2] != fns[1])
                  fns[2]->pop_transform();

              fns[0]->push_transform(2);
              if(this->xdisp != NULL)
                if(fns[1] != fns[0])
                  fns[1]->push_transform(2);
              if(this->ydisp != NULL)
                if(fns[2] != fns[1])
                  fns[2]->push_transform(2);
              process_quad(fns, mid4, mid1, iv2, mid2, level + 1, val, phx, phy, quad_indices[3], curved);
              fns[0]->pop_transform();
              if(this->xdisp != NULL)
                if(fns[1] != fns[0])
                  fns[1]->pop_transform();
              if(this->ydisp != NULL)
                if(fns[2] != fns[1])
                  fns[2]->pop_transform();

              fns[0]->push_transform(3);
              if(this->xdisp != NULL)
                if(fns[1] != fns[0])
                  fns[1]->push_transform(3);
              if(this->ydisp != NULL)
                if(fns[2] != fns[1])
                  fns[2]->push_transform(3);
              process_quad(fns, mid3, mid4, mid2, iv3, level + 1, val, phx, phy, quad_indices[4], curved);
              fns[0]->pop_transform();
              if(this->xdisp != NULL)
                if(fns[1] != fns[0])
                  fns[1]->pop_transform();
              if(this->ydisp != NULL)
                if(fns[2] != fns[1])
                  fns[2]->pop_transform();
            }
            else
              if(split == 1) // h-split
              {
                fns[0]->push_transform(4);
                if(this->ydisp != NULL)
                  if(fns[1] != fns[0])
                    fns[1]->push_transform(4);
                if(this->ydisp != NULL)
                  if(fns[2] != fns[1])
                    fns[2]->push_transform(4);
                process_quad(fns, iv0, iv1, mid1, mid3, level + 1, val, phx, phy, quad_indices[5], curved);
                fns[0]->pop_transform();
                if(this->xdisp != NULL)
                  if(fns[1] != fns[0])
                    fns[1]->pop_transform();
                if(this->ydisp != NULL)
                  if(fns[2] != fns[1])
                    fns[2]->pop_transform();

                fns[0]->push_transform(5);
                if(this->xdisp != NULL)
                  if(fns[1] != fns[0])
                    fns[1]->push_transform(5);
                if(this->ydisp != NULL)
                  if(fns[2] != fns[1])
                    fns[2]->push_transform(5);
                process_quad(fns, mid3, mid1, iv2, iv3, level + 1, val, phx, phy, quad_indices[6], curved);
                fns[0]->pop_transform();
                if(this->xdisp != NULL)
                  if(fns[1] != fns[0])
                    fns[1]->pop_transform();
                if(this->ydisp != NULL)
                  if(fns[2] != fns[1])
                    fns[2]->pop_transform();
              }
              else // v-split
              {
                fns[0]->push_transform(6);
                if(this->xdisp != NULL)
                  if(fns[1] != fns[0])
                    fns[1]->push_transform(6);
                if(this->ydisp != NULL)
                  if(fns[2] != fns[1])
                    fns[2]->push_transform(6);
                process_quad(fns, iv0, mid0, mid2, iv3, level + 1, val, phx, phy, quad_indices[7], curved);
                fns[0]->pop_transform();
                if(this->xdisp != NULL)
                  if(fns[1] != fns[0])
                    fns[1]->pop_transform();
                if(this->ydisp != NULL)
                  if(fns[2] != fns[1])
                    fns[2]->pop_transform();

                fns[0]->push_transform(7);
                if(this->xdisp != NULL)
                  if(fns[1] != fns[0])
                    fns[1]->push_transform(7);
                if(this->ydisp != NULL)
                  if(fns[2] != fns[1])
                    fns[2]->push_transform(7);
                process_quad(fns, mid0, iv1, iv2, mid2, level + 1, val, phx, phy, quad_indices[8], curved);
                fns[0]->pop_transform();
                if(this->xdisp != NULL)
                  if(fns[1] != fns[0])
                    fns[1]->pop_transform();
                if(this->ydisp != NULL)
                  if(fns[2] != fns[1])
                    fns[2]->pop_transform();
              }
              return;
          }
        }

        // output two linear triangles,
        if(!flip)
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

      void Linearizer::set_displacement(MeshFunction<double>* xdisp, MeshFunction<double>* ydisp, double dmult)
      {
        if(xdisp != NULL)
        {
          user_xdisp = true;
          this->xdisp = xdisp;
        }
        if(ydisp != NULL)
        {
          user_ydisp = true;
          this->ydisp = ydisp;
        }
        this->dmult = dmult;
      }

      void Linearizer::process_solution(MeshFunction<double>* sln, int item_, double eps)
      {
        // Important, sets the current caughtException to NULL.
        this->caughtException = NULL;

        lock_data();
        this->tick();

        // Initialization of 'global' stuff.
        this->sln = sln;
        this->item = item_;
        this->eps = eps;
        //   get the component and desired value from item.
        if(item >= 0x40)
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

        // Initialization of computation stuff.
        //    sizes.
        this->vertex_size = std::max(100 * this->sln->get_mesh()->get_num_elements(), 50000);
        this->triangle_size = std::max(150 * this->sln->get_mesh()->get_num_elements(), 75000);
        this->edges_size = std::max(100 * this->sln->get_mesh()->get_num_elements(), 50000);
        //    counts.
        this->vertex_count = 0;
        this->triangle_count = 0;
        this->edges_count = 0;
        //    reuse or allocate vertex, triangle and edge arrays.
        this->verts = (double3*) realloc(this->verts, sizeof(double3) * this->vertex_size);
        this->tris = (int3*) realloc(this->tris, sizeof(int3) * this->triangle_size);
        this->edges = (int3*) realloc(this->edges, sizeof(int3) * this->edges_size);
        this->info = (int4*) malloc(sizeof(int4) * this->vertex_size);
        //    initialize the hash table
        this->hash_table = (int*) malloc(sizeof(int) * this->vertex_size);
        memset(this->hash_table, 0xff, sizeof(int) * this->vertex_size);

        // select the linearization quadratures
        Quad2D *old_quad, *old_quad_x = NULL, *old_quad_y = NULL;
        old_quad = sln->get_quad_2d();
        if(xdisp != NULL)
          old_quad_x = xdisp->get_quad_2d();
        if(ydisp != NULL)
          old_quad_y = ydisp->get_quad_2d();

        // obtain the solution in vertices, estimate the maximum solution value
        // meshes.
        Hermes::vector<const Mesh*> meshes;
        meshes.push_back(sln->get_mesh());
        if(xdisp != NULL)
          meshes.push_back(xdisp->get_mesh());
        if(ydisp != NULL)
          meshes.push_back(ydisp->get_mesh());

        // Parallelization
        MeshFunction<double>*** fns = new MeshFunction<double>**[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
        for(unsigned int i = 0; i < Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads); i++)
        {
          fns[i] = new MeshFunction<double>*[3];
          fns[i][0] = sln->clone();
          fns[i][0]->set_refmap(new RefMap);
          fns[i][0]->set_quad_2d(&g_quad_lin);
          if(xdisp != NULL)
          {
            fns[i][1] = xdisp->clone();
            //fns[i][1]->set_refmap(new RefMap);
            fns[i][1]->set_quad_2d(&g_quad_lin);
          }
          if(ydisp != NULL)
          {
            fns[i][2] = ydisp->clone();
            //fns[i][2]->set_refmap(new RefMap);
            fns[i][2]->set_quad_2d(&g_quad_lin);
          }
        }

        Transformable*** trfs = new Transformable**[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];
        for(unsigned int i = 0; i < Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads); i++)
        {
          trfs[i] = new Transformable*[3];
          trfs[i][0] = fns[i][0];
          if(xdisp != NULL)
            trfs[i][1] = fns[i][1];
          if(ydisp != NULL)
            trfs[i][2] = fns[i][2];
        }

        Traverse trav_master(true);
        unsigned int num_states = trav_master.get_num_states(meshes);

        trav_master.begin(meshes.size(), &(meshes.front()));

        Traverse* trav = new Traverse[Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads)];

        for(unsigned int i = 0; i < Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads); i++)
        {
          trav[i].begin(meshes.size(), &(meshes.front()), trfs[i]);
          trav[i].stack = trav_master.stack;
        }

        int state_i;

#define CHUNKSIZE 1
        int num_threads_used = Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads);
#pragma omp parallel shared(trav_master) private(state_i) num_threads(num_threads_used)
        {
#pragma omp for schedule(dynamic, CHUNKSIZE)
          for(state_i = 0; state_i < num_states; state_i++)
          {
            try
            {
              Traverse::State current_state;
#pragma omp critical(get_next_state)
              current_state = trav[omp_get_thread_num()].get_next_state(&trav_master.top, &trav_master.id);

              fns[omp_get_thread_num()][0]->set_quad_order(0, this->item);
              double* val = fns[omp_get_thread_num()][0]->get_values(component, value_type);

              for (unsigned int i = 0; i < current_state.e[0]->get_num_surf(); i++)
              {
                double f = val[i];
#pragma omp critical (max)
                if(this->auto_max && finite(f) && fabs(f) > this->max)
                  this->max = fabs(f);
                double c = current_state.e[0]->get_diameter();
                if(c > this->cmax)
#pragma omp critical(vectorizer_get_cmax)
                  if(c > this->cmax)
                    this->cmax = c;
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
        
        trav_master.finish();
        for(unsigned int i = 0; i < Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads); i++)
          trav[i].finish();

        trav_master.begin(meshes.size(), &(meshes.front()));

        for(unsigned int i = 0; i < Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads); i++)
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

              fns[omp_get_thread_num()][0]->set_quad_order(0, this->item);
              double* val = fns[omp_get_thread_num()][0]->get_values(component, value_type);
              if(val == NULL)
              {
                delete [] trav;
                throw Hermes::Exceptions::Exception("Item not defined in the solution.");
              }

              if(xdisp != NULL)
                fns[omp_get_thread_num()][1]->set_quad_order(0, H2D_FN_VAL);
              if(ydisp != NULL)
                fns[omp_get_thread_num()][2]->set_quad_order(0, H2D_FN_VAL);

              double *dx = NULL;
              double *dy = NULL;
              if(xdisp != NULL)
                dx = fns[omp_get_thread_num()][1]->get_fn_values();
              if(ydisp != NULL)
                dy = fns[omp_get_thread_num()][2]->get_fn_values();

              int iv[4];
              for (unsigned int i = 0; i < current_state.e[0]->get_num_surf(); i++)
              {
                double f = val[i];
                double x_disp = fns[omp_get_thread_num()][0]->get_refmap()->get_phys_x(0)[i];
                double y_disp = fns[omp_get_thread_num()][0]->get_refmap()->get_phys_y(0)[i];
                if(this->xdisp != NULL)
                  x_disp += dmult * dx[i];
                if(this->ydisp != NULL)
                  y_disp += dmult * dy[i];

                iv[i] = this->get_vertex(-fns[omp_get_thread_num()][0]->get_active_element()->vn[i]->id, -fns[omp_get_thread_num()][0]->get_active_element()->vn[i]->id, x_disp, y_disp, f);
              }

              // recur to sub-elements
              if(current_state.e[0]->is_triangle())
                process_triangle(fns[omp_get_thread_num()], iv[0], iv[1], iv[2], 0, NULL, NULL, NULL, NULL, current_state.e[0]->is_curved());
              else
                process_quad(fns[omp_get_thread_num()], iv[0], iv[1], iv[2], iv[3], 0, NULL, NULL, NULL, NULL, current_state.e[0]->is_curved());

              for (unsigned int i = 0; i < current_state.e[0]->get_num_surf(); i++)
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
        for(unsigned int i = 0; i < Hermes2DApi.getParamValue(Hermes::Hermes2D::numThreads); i++)
        {
          trav[i].finish();
          for(unsigned int j = 0; j < (1 + (xdisp != NULL? 1 : 0) + (ydisp != NULL ? 1 : 0)); j++)
            delete fns[i][j];
          delete [] fns[i];
          delete [] trfs[i];
        }
        delete [] fns;
        delete [] trfs;
        delete [] trav;

        // for contours, without regularization.
        if(this->trisContours != NULL)
          delete [] this->trisContours;
        this->trisContours = new int3[this->triangle_count];
        memcpy(this->trisContours, this->tris, this->triangle_count * sizeof(int3));
        triangleContours_count = this->triangle_count;

        if(this->caughtException != NULL)
        {
          this->unlock_data();
          ::free(hash_table);
          ::free(info);
          throw *(this->caughtException);
        }

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
            regularize_triangle(iv0, iv1, iv2, mid0, mid1, mid2);
          }
        }

        find_min_max();

        this->unlock_data();

        // select old quadratrues
        sln->set_quad_2d(old_quad);
        if(user_xdisp)
          xdisp->set_quad_2d(old_quad_x);
        else
          delete xdisp;
        if(user_ydisp)
          ydisp->set_quad_2d(old_quad_y);
        else
          delete ydisp;

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
          if(finite(verts[i][2]) && verts[i][2] < min_val) min_val = verts[i][2];
          if(finite(verts[i][2]) && verts[i][2] > max_val) max_val = verts[i][2];
        }
      }

      int Linearizer::get_vertex(int p1, int p2, double x, double y, double value)
      {
        // search for an existing vertex
        if(p1 > p2) std::swap(p1, p2);
        int index = this->hash(p1, p2);
        int i = this->hash_table[index];
        while (i >= 0 && i < this->vertex_count)
        {
          if(
            this->info[i][0] == p1 && this->info[i][1] == p2 &&
            (value == verts[i][2] || fabs(value - verts[i][2]) < this->max*1e-8) &&
            (fabs(x - verts[i][0]) < 1e-8) &&
            (fabs(y - verts[i][1]) < 1e-8)
            )
            return i;
          // note that we won't return a vertex with a different value than the required one;
          // this takes care for discontinuities in the solution, where more vertices
          // with different values will be created
          i = info[i][2];
        }

        // if not found, create a new one
#pragma omp critical(realloc_vertices)
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
        if(this->vertex_count >= this->vertex_size)
        {
          this->vertex_size *= 2;
          verts = (double3*) realloc(verts, sizeof(double3) * vertex_size);
          this->info = (int4*) realloc(info, sizeof(int4) * vertex_size);
          this->hash_table = (int*) realloc(hash_table, sizeof(int) * vertex_size);
          memset(this->hash_table + this->vertex_size / 2, 0xff, sizeof(int) * this->vertex_size / 2);
        }
        return this->vertex_count++;
      }

      void Linearizer::free()
      {
        if(verts != NULL)
        {
          ::free(verts);
          verts = NULL;
        }
        if(trisContours != NULL)
          delete [] this->trisContours;
      }

      Linearizer::~Linearizer()
      {
        free();
      }

      void Linearizer::save_solution_vtk(MeshFunction<double>* sln, const char* filename, const char *quantity_name,
        bool mode_3D, int item, double eps)
      {
        this->sln = sln;
        process_solution(sln, item, eps);

        FILE* f = fopen(filename, "wb");
        if(f == NULL) throw Hermes::Exceptions::Exception("Could not open %s for writing.", filename);
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
          if(mode_3D == true) fprintf(f, "%g %g %g\n", this->verts[i][0], this->verts[i][1], this->verts[i][2]);
          else fprintf(f, "%g %g %g\n", this->verts[i][0], this->verts[i][1], 0.0);
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

        // This outputs double solution values.
        fprintf(f, "\n");
        fprintf(f, "POINT_DATA %d\n", this->vertex_count);
        fprintf(f, "SCALARS %s %s %d\n", quantity_name, "float", 1);
        fprintf(f, "LOOKUP_TABLE %s\n", "default");
        for (int i = 0; i < this->vertex_count; i++)
        {
          fprintf(f, "%g\n", this->verts[i][2]);
        }

        unlock_data();
        fclose(f);
      }

      void Linearizer::calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const
      {
        if(verts == NULL)
          throw Exceptions::Exception("Cannot calculate AABB from NULL vertices");
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

      int Linearizer::get_num_contour_triangles()
      {
        return this->triangleContours_count;
      }

      int3* Linearizer::get_contour_triangles()
      {
        return this->trisContours;
      }
    }
  }
}
