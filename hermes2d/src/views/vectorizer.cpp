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

      template<typename Scalar>
      Vectorizer<Scalar>::Vectorizer() : Linearizer<Scalar>()
      {
        verts = NULL;
        dashes = NULL;
        cd = 0;
      }

      template<typename Scalar>
      int Vectorizer<Scalar>::create_vertex(double x, double y, double xvalue, double yvalue)
      {
        int i = add_vertex();
        verts[i][0] = x;
        verts[i][1] = y;
        verts[i][2] = xvalue;
        verts[i][3] = yvalue;
        return i;
      }

      template<typename Scalar>
      int Vectorizer<Scalar>::get_vertex(int p1, int p2, double x, double y, double xvalue, double yvalue)
      {
        // search for an existing vertex
        if (p1 > p2) std::swap(p1, p2);
        int index = this->hash(p1, p2);
        int i = this->hash_table[index];
        while (i >= 0)
        {
          if (this->info[i][0] == p1 && this->info[i][1] == p2 &&
            (fabs((xvalue - verts[i][2]) / xvalue) < 1e-4 && fabs((yvalue - verts[i][3]) / yvalue) < 1e-4)) return i;
          // note that we won't return a vertex with a different value than the required one;
          // this takes care for discontinuities in the solution, where more vertices
          // with different values will be created
          i = this->info[i][2];
        }

        // if not found, create a new one
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

      template<>
      void Vectorizer<double>::process_triangle(int iv0, int iv1, int iv2, int level,
        double* xval, double* yval, double* phx, double* phy, int* idx)
      {
        if (level < LIN_MAX_LEVEL)
        {
          int i;
          if (!(level & 1))
          {
            // obtain solution values and physical element coordinates
            xsln->set_quad_order(1, xitem);
            ysln->set_quad_order(1, yitem);
            xval = xsln->get_values(xia, xib);
            yval = ysln->get_values(yia, yib);
            for (i = 0; i < lin_np_tri[1]; i++) 
            {
              double m = (sqrt(sqr(xval[i]) + sqr(yval[i])));
              if (finite(m) && fabs(m) > max) max = fabs(m);
            }

            if (curved)
            {
              RefMap* refmap = xsln->get_refmap();
              phx = refmap->get_phys_x(1);
              phy = refmap->get_phys_y(1);
            }
            idx = tri_indices[0];
          }

          // obtain linearized values and coordinates at the midpoints
          double midval[4][3];
          for (i = 0; i < 4; i++)
          {
            midval[i][0] = (verts[iv0][i] + verts[iv1][i])*0.5;
            midval[i][1] = (verts[iv1][i] + verts[iv2][i])*0.5;
            midval[i][2] = (verts[iv2][i] + verts[iv0][i])*0.5;
          };

          // determine whether or not to split the element
          bool split;
          if (eps >= 1.0)
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
            if (curved && !split)
            {
              double cerr = 0.0, cden = 0.0; // fixme
              for (i = 0; i < 3; i++)
              {
                cerr += fabs(phx[idx[i]] - midval[0][i]) + fabs(phy[idx[i]] - midval[1][i]);
                cden += fabs(phx[idx[i]]) + fabs(phy[idx[i]]);
              }
              split = (cerr > cden*2.5e-4);
            }

            // do extra tests at level 0, so as not to miss some functions with zero error at edge midpoints
            if (level == 0 && !split)
            {
              split = (fabs(sqrt(sqr(xval[8]) + sqr(yval[8])) - 0.5*(sqrt(sqr(midval[2][0]) + sqr(midval[3][0])) + sqrt(sqr(midval[2][1]) + sqr(midval[3][1])))) +
                fabs(sqrt(sqr(xval[9]) + sqr(yval[9])) - 0.5*(sqrt(sqr(midval[2][1]) + sqr(midval[3][1])) + sqrt(sqr(midval[2][2]) + sqr(midval[3][2])))) +
                fabs(sqrt(sqr(xval[4]) + sqr(yval[4])) - 0.5*(sqrt(sqr(midval[2][2]) + sqr(midval[3][2])) + sqrt(sqr(midval[2][0]) + sqr(midval[3][0]))))) > max*3*eps;
            }
          }

          // split the triangle if the error is too large, otherwise produce a linear triangle
          if (split)
          {
            if (curved)
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
              push_transform(0);  process_triangle(iv0, mid0, mid2,  level+1, xval, yval, phx, phy, tri_indices[1]);  pop_transform();
              push_transform(1);  process_triangle(mid0, iv1, mid1,  level+1, xval, yval, phx, phy, tri_indices[2]);  pop_transform();
              push_transform(2);  process_triangle(mid2, mid1, iv2,  level+1, xval, yval, phx, phy, tri_indices[3]);  pop_transform();
              push_transform(3);  process_triangle(mid1, mid2, mid0, level+1, xval, yval, phx, phy, tri_indices[4]);  pop_transform();
              return;
          }
        }

        // no splitting: output a linear triangle
        add_triangle(iv0, iv1, iv2);
      }
      template<>
      void Vectorizer<std::complex<double> >::process_triangle(int iv0, int iv1, int iv2, int level,
        std::complex<double>* xval, std::complex<double>* yval, double* phx, double* phy, int* idx)
      {
        if (level < LIN_MAX_LEVEL)
        {
          int i;
          if (!(level & 1))
          {
            // obtain solution values and physical element coordinates
            xsln->set_quad_order(1, xitem);
            ysln->set_quad_order(1, yitem);
            xval = xsln->get_values(xia, xib);
            yval = ysln->get_values(yia, yib);
            for (i = 0; i < lin_np_tri[1]; i++) 
            {
              double m = sqrt(sqr(xval[i].real()) + sqr(yval[i].real()));
              if (finite(m) && fabs(m) > max) max = fabs(m);
            }

            if (curved)
            {
              RefMap* refmap = xsln->get_refmap();
              phx = refmap->get_phys_x(1);
              phy = refmap->get_phys_y(1);
            }
            idx = tri_indices[0];
          }

          // obtain linearized values and coordinates at the midpoints
          double midval[4][3];
          for (i = 0; i < 4; i++)
          {
            midval[i][0] = (verts[iv0][i] + verts[iv1][i])*0.5;
            midval[i][1] = (verts[iv1][i] + verts[iv2][i])*0.5;
            midval[i][2] = (verts[iv2][i] + verts[iv0][i])*0.5;
          };

          // determine whether or not to split the element
          bool split;
          if (eps >= 1.0)
          {
            //if eps > 1, the user wants a fixed number of refinements (no adaptivity)
            split = (level < eps);
          }
          else
          {
            // calculate the approximate error of linearizing the normalized solution
            double err = fabs(sqrt(sqr(xval[idx[0]].real()) + sqr(yval[idx[0]].real())) - sqrt(sqr(midval[2][0]) + sqr(midval[3][0]))) +
              fabs(sqrt(sqr(xval[idx[1]].real()) + sqr(yval[idx[1]].real())) - sqrt(sqr(midval[2][1]) + sqr(midval[3][1]))) +
              fabs(sqrt(sqr(xval[idx[2]].real()) + sqr(yval[idx[2]].real())) - sqrt(sqr(midval[2][2]) + sqr(midval[3][2])));
            split = !finite(err) || err > max*3*eps;

            // do the same for the curvature
            if (curved && !split)
            {
              double cerr = 0.0, cden = 0.0; // fixme
              for (i = 0; i < 3; i++)
              {
                cerr += fabs(phx[idx[i]] - midval[0][i]) + fabs(phy[idx[i]] - midval[1][i]);
                cden += fabs(phx[idx[i]]) + fabs(phy[idx[i]]);
              }
              split = (cerr > cden*2.5e-4);
            }

            // do extra tests at level 0, so as not to miss some functions with zero error at edge midpoints
            if (level == 0 && !split)
            {
              split = (fabs(sqrt(sqr(xval[8].real()) + sqr(yval[8].real())) - 0.5*(sqrt(sqr(midval[2][0]) + sqr(midval[3][0])) + sqrt(sqr(midval[2][1]) + sqr(midval[3][1])))) +
                fabs(sqrt(sqr(xval[9].real()) + sqr(yval[9].real())) - 0.5*(sqrt(sqr(midval[2][1]) + sqr(midval[3][1])) + sqrt(sqr(midval[2][2]) + sqr(midval[3][2])))) +
                fabs(sqrt(sqr(xval[4].real()) + sqr(yval[4].real())) - 0.5*(sqrt(sqr(midval[2][2]) + sqr(midval[3][2])) + sqrt(sqr(midval[2][0]) + sqr(midval[3][0]))))) > max*3*eps;
            }
          }

          // split the triangle if the error is too large, otherwise produce a linear triangle
          if (split)
          {
            if (curved)
              for (i = 0; i < 3; i++) 
              {
                midval[0][i] = phx[idx[i]];
                midval[1][i] = phy[idx[i]];
              }

              // obtain mid-edge vertices
              int mid0 = get_vertex(iv0, iv1, midval[0][0], midval[1][0], xval[idx[0]].real(), yval[idx[0]].real());
              int mid1 = get_vertex(iv1, iv2, midval[0][1], midval[1][1], xval[idx[1]].real(), yval[idx[1]].real());
              int mid2 = get_vertex(iv2, iv0, midval[0][2], midval[1][2], xval[idx[2]].real(), yval[idx[2]].real());

              // recur to sub-elements
              push_transform(0);  process_triangle(iv0, mid0, mid2,  level+1, xval, yval, phx, phy, tri_indices[1]);  pop_transform();
              push_transform(1);  process_triangle(mid0, iv1, mid1,  level+1, xval, yval, phx, phy, tri_indices[2]);  pop_transform();
              push_transform(2);  process_triangle(mid2, mid1, iv2,  level+1, xval, yval, phx, phy, tri_indices[3]);  pop_transform();
              push_transform(3);  process_triangle(mid1, mid2, mid0, level+1, xval, yval, phx, phy, tri_indices[4]);  pop_transform();
              return;
          }
        }

        // no splitting: output a linear triangle
        add_triangle(iv0, iv1, iv2);
      }

      template<>
      void Vectorizer<double>::process_quad(int iv0, int iv1, int iv2, int iv3, int level,
        double* xval, double* yval, double* phx, double* phy, int* idx)
      {
        // try not to split through the vertex with the largest value
        int a = (sqr(verts[iv1][2]) + sqr(verts[iv1][3]) > sqr(verts[iv0][2]) + sqr(verts[iv0][3])) ? iv0 : iv1;
        int b = (sqr(verts[iv2][2]) + sqr(verts[iv2][3]) > sqr(verts[iv3][2]) + sqr(verts[iv3][3])) ? iv2 : iv3;
        a = (sqr(verts[a][2]) + sqr(verts[a][3]) > sqr(verts[b][2]) + sqr(verts[b][3])) ? a : b;

        int flip = (a == iv1 || a == iv3) ? 1 : 0;

        if (level < LIN_MAX_LEVEL)
        {
          int i;
          if (!(level & 1))
          {
            // obtain solution values and physical element coordinates
            xsln->set_quad_order(1, xitem);
            ysln->set_quad_order(1, yitem);
            xval = xsln->get_values(xia, xib);
            yval = ysln->get_values(yia, yib);
            for (i = 0; i < lin_np_quad[1]; i++) 
            {
              double m = sqrt(sqr(xval[i]) + sqr(yval[i]));
              if (finite(m) && fabs(m) > max) max = fabs(m);
            }

            if (curved)
            {
              RefMap* refmap = xsln->get_refmap();
              phx = refmap->get_phys_x(1);
              phy = refmap->get_phys_y(1);
            }
            idx = quad_indices[0];
          }

          // obtain linearized values and coordinates at the midpoints
          double midval[4][5];
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
          if (eps >= 1.0)
          {
            //if eps > 1, the user wants a fixed number of refinements (no adaptivity)
            split = (level < eps);
          }

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
            split = !finite(err) || err > max*4*eps;

            // do the same for the curvature
            if (curved && !split)
            {
              double cerr = 0.0, cden = 0.0; // fixme
              for (i = 0; i < 5; i++)
              {
                cerr += fabs(phx[idx[i]] - midval[0][i]) + fabs(phy[idx[i]] - midval[1][i]);
                cden += fabs(phx[idx[i]]) + fabs(phy[idx[i]]);
              }
              split = (cerr > cden*2.5e-4);
            }

            // do extra tests at level 0, so as not to miss functions with zero error at edge midpoints
            if (level == 0 && !split)
            {
              err = fabs(sqrt(sqr(xval[13]) + sqr(yval[13])) - 0.5*(sqrt(sqr(midval[2][0]) + sqr(midval[3][0])) + sqrt(sqr(midval[2][1]) + sqr(midval[3][1])))) +
                fabs(sqrt(sqr(xval[17]) + sqr(yval[17])) - 0.5*(sqrt(sqr(midval[2][1]) + sqr(midval[3][1])) + sqrt(sqr(midval[2][2]) + sqr(midval[3][2])))) +
                fabs(sqrt(sqr(xval[20]) + sqr(yval[20])) - 0.5*(sqrt(sqr(midval[2][2]) + sqr(midval[3][2])) + sqrt(sqr(midval[2][3]) + sqr(midval[3][3])))) +
                fabs(sqrt(sqr(xval[9]) + sqr(yval[9]))  - 0.5*(sqrt(sqr(midval[2][3]) + sqr(midval[3][3])) + sqrt(sqr(midval[2][0]) + sqr(midval[3][0]))));
              split = !finite(err) || (err) > max*2*eps; //?
            }
          }

          // split the quad if the error is too large, otherwise produce two linear triangles
          if (split)
          {
            if (curved)
              for (i = 0; i < 5; i++) 
              {
                midval[0][i] = phx[idx[i]];
                midval[1][i] = phy[idx[i]];
              }

              // obtain mid-edge and mid-element vertices
              int mid0 = get_vertex(iv0,  iv1,  midval[0][0], midval[1][0], xval[idx[0]], yval[idx[0]]);
              int mid1 = get_vertex(iv1,  iv2,  midval[0][1], midval[1][1], xval[idx[1]], yval[idx[1]]);
              int mid2 = get_vertex(iv2,  iv3,  midval[0][2], midval[1][2], xval[idx[2]], yval[idx[2]]);
              int mid3 = get_vertex(iv3,  iv0,  midval[0][3], midval[1][3], xval[idx[3]], yval[idx[3]]);
              int mid4 = get_vertex(mid0, mid2, midval[0][4], midval[1][4], xval[idx[4]], yval[idx[4]]);

              // recur to sub-elements
              push_transform(0);  process_quad(iv0, mid0, mid4, mid3, level+1, xval, yval, phx, phy, quad_indices[1]);  pop_transform();
              push_transform(1);  process_quad(mid0, iv1, mid1, mid4, level+1, xval, yval, phx, phy, quad_indices[2]);  pop_transform();
              push_transform(2);  process_quad(mid4, mid1, iv2, mid2, level+1, xval, yval, phx, phy, quad_indices[3]);  pop_transform();
              push_transform(3);  process_quad(mid3, mid4, mid2, iv3, level+1, xval, yval, phx, phy, quad_indices[4]);  pop_transform();
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
      template<>
      void Vectorizer<std::complex<double> >::process_quad(int iv0, int iv1, int iv2, int iv3, int level,
        std::complex<double>* xval, std::complex<double>* yval, double* phx, double* phy, int* idx)
      {
        // try not to split through the vertex with the largest value
        int a = (sqr(verts[iv1][2]) + sqr(verts[iv1][3]) > sqr(verts[iv0][2]) + sqr(verts[iv0][3])) ? iv0 : iv1;
        int b = (sqr(verts[iv2][2]) + sqr(verts[iv2][3]) > sqr(verts[iv3][2]) + sqr(verts[iv3][3])) ? iv2 : iv3;
        a = (sqr(verts[a][2]) + sqr(verts[a][3]) > sqr(verts[b][2]) + sqr(verts[b][3])) ? a : b;

        int flip = (a == iv1 || a == iv3) ? 1 : 0;

        if (level < LIN_MAX_LEVEL)
        {
          int i;
          if (!(level & 1))
          {
            // obtain solution values and physical element coordinates
            xsln->set_quad_order(1, xitem);
            ysln->set_quad_order(1, yitem);
            xval = xsln->get_values(xia, xib);
            yval = ysln->get_values(yia, yib);
            for (i = 0; i < lin_np_quad[1]; i++) 
            {
              double m = sqrt(sqr(xval[i].real()) + sqr(yval[i].real()));
              if (finite(m) && fabs(m) > max) max = fabs(m);
            }

            if (curved)
            {
              RefMap* refmap = xsln->get_refmap();
              phx = refmap->get_phys_x(1);
              phy = refmap->get_phys_y(1);
            }
            idx = quad_indices[0];
          }

          // obtain linearized values and coordinates at the midpoints
          double midval[4][5];
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
          if (eps >= 1.0)
          {
            //if eps > 1, the user wants a fixed number of refinements (no adaptivity)
            split = (level < eps);
          }

          else
          {
            // the value of the middle point is not the average of the four vertex values, since quad == 2 triangles
            midval[2][4] = flip ? (verts[iv0][2] + verts[iv2][2]) * 0.5 : (verts[iv1][2] + verts[iv3][2]) * 0.5;
            midval[3][4] = flip ? (verts[iv0][3] + verts[iv2][3]) * 0.5 : (verts[iv1][3] + verts[iv3][3]) * 0.5;

            // calculate the approximate error of linearizing the normalized solution
            double err = fabs(sqrt(sqr(xval[idx[0]].real()) + sqr(yval[idx[0]].real())) - sqrt(sqr(midval[2][0]) + sqr(midval[3][0]))) +
              fabs(sqrt(sqr(xval[idx[1]].real()) + sqr(yval[idx[1]].real())) - sqrt(sqr(midval[2][1]) + sqr(midval[3][1]))) +
              fabs(sqrt(sqr(xval[idx[2]].real()) + sqr(yval[idx[2]].real())) - sqrt(sqr(midval[2][2]) + sqr(midval[3][2]))) +
              fabs(sqrt(sqr(xval[idx[3]].real()) + sqr(yval[idx[3]].real())) - sqrt(sqr(midval[2][3]) + sqr(midval[3][3]))) +
              fabs(sqrt(sqr(xval[idx[4]].real()) + sqr(yval[idx[4]].real())) - sqrt(sqr(midval[2][4]) + sqr(midval[3][4])));
            split = !finite(err) || err > max*4*eps;

            // do the same for the curvature
            if (curved && !split)
            {
              double cerr = 0.0, cden = 0.0; // fixme
              for (i = 0; i < 5; i++)
              {
                cerr += fabs(phx[idx[i]] - midval[0][i]) + fabs(phy[idx[i]] - midval[1][i]);
                cden += fabs(phx[idx[i]]) + fabs(phy[idx[i]]);
              }
              split = (cerr > cden*2.5e-4);
            }

            // do extra tests at level 0, so as not to miss functions with zero error at edge midpoints
            if (level == 0 && !split)
            {
              err = fabs(sqrt(sqr(xval[13]) + sqr(yval[13])) - 0.5*(sqrt(sqr(midval[2][0]) + sqr(midval[3][0])) + sqrt(sqr(midval[2][1]) + sqr(midval[3][1])))) +
                fabs(sqrt(sqr(xval[17]) + sqr(yval[17])) - 0.5*(sqrt(sqr(midval[2][1]) + sqr(midval[3][1])) + sqrt(sqr(midval[2][2]) + sqr(midval[3][2])))) +
                fabs(sqrt(sqr(xval[20]) + sqr(yval[20])) - 0.5*(sqrt(sqr(midval[2][2]) + sqr(midval[3][2])) + sqrt(sqr(midval[2][3]) + sqr(midval[3][3])))) +
                fabs(sqrt(sqr(xval[9]) + sqr(yval[9]))  - 0.5*(sqrt(sqr(midval[2][3]) + sqr(midval[3][3])) + sqrt(sqr(midval[2][0]) + sqr(midval[3][0]))));
              split = !finite(err) || (err) > max*2*eps; //?
            }
          }

          // split the quad if the error is too large, otherwise produce two linear triangles
          if (split)
          {
            if (curved)
              for (i = 0; i < 5; i++) 
              {
                midval[0][i] = phx[idx[i]];
                midval[1][i] = phy[idx[i]];
              }

              // obtain mid-edge and mid-element vertices
              int mid0 = get_vertex(iv0,  iv1,  midval[0][0], midval[1][0], xval[idx[0]].real(), yval[idx[0]].real());
              int mid1 = get_vertex(iv1,  iv2,  midval[0][1], midval[1][1], xval[idx[1]].real(), yval[idx[1]].real());
              int mid2 = get_vertex(iv2,  iv3,  midval[0][2], midval[1][2], xval[idx[2]].real(), yval[idx[2]].real());
              int mid3 = get_vertex(iv3,  iv0,  midval[0][3], midval[1][3], xval[idx[3]].real(), yval[idx[3]].real());
              int mid4 = get_vertex(mid0, mid2, midval[0][4], midval[1][4], xval[idx[4]].real(), yval[idx[4]].real());

              // recur to sub-elements
              push_transform(0);  process_quad(iv0, mid0, mid4, mid3, level+1, xval, yval, phx, phy, quad_indices[1]);  pop_transform();
              push_transform(1);  process_quad(mid0, iv1, mid1, mid4, level+1, xval, yval, phx, phy, quad_indices[2]);  pop_transform();
              push_transform(2);  process_quad(mid4, mid1, iv2, mid2, level+1, xval, yval, phx, phy, quad_indices[3]);  pop_transform();
              push_transform(3);  process_quad(mid3, mid4, mid2, iv3, level+1, xval, yval, phx, phy, quad_indices[4]);  pop_transform();
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

      template<typename Scalar>
      void Vectorizer<Scalar>::process_dash(int iv1, int iv2)
      {
        int mid = this->peek_vertex(iv1, iv2);
        if (mid != -1)
        {
          process_dash(iv1, mid);
          process_dash(mid, iv2);
        }
        else
          add_dash(iv1, iv2);
      }

      template<typename Scalar>
      void Vectorizer<Scalar>::find_min_max()
      {
        // find min & max vertex values
        this->min_val =  1e100;
        this->max_val = -1e100;
        for (int i = 0; i < this->nv; i++)
        {
          double mag = (verts[i][2]*verts[i][2] + verts[i][3]*verts[i][3]);
          if (finite(mag) && mag < this->min_val) this->min_val = mag;
          if (finite(mag) && mag > this->max_val) this->max_val = mag;
        }
        this->max_val = sqrt(this->max_val);
        this->min_val = sqrt(this->min_val);

      }

      template<>
      void Vectorizer<double>::process_solution(MeshFunction<double>* xsln, int xitem, MeshFunction<double>* ysln, int yitem, double eps)
      {
        // sanity check
        if (xsln == NULL || ysln == NULL) error("One of the solutions is NULL in Vectorizer:process_solution().");


        lock_data();
        Hermes::TimePeriod cpu_time;

        // initialization
        this->xsln = xsln;
        this->ysln = ysln;
        this->xitem = xitem;
        this->yitem = yitem;
        this->eps = eps;
        this->nv = nt = ne = nd = 0;
        del_slot = -1;

        Mesh* meshes[2] = { xsln->get_mesh(), ysln->get_mesh() };
        if (meshes[0] == NULL || meshes[1] == NULL) 
        {
          error("One of the meshes is NULL in Vectorizer:process_solution().");
        }

        Transformable* fns[2] = { xsln, ysln };
        Traverse trav;

        // estimate the required number of vertices and triangles
        // (based on the assumption that the linear mesh will be
        // about four-times finer than the original mesh).
        int nn = meshes[0]->get_num_elements() + meshes[1]->get_num_elements();
        int ev = std::max(32 * nn, 10000);
        int et = std::max(64 * nn, 20000);
        int ee = std::max(24 * nn, 7500);
        int ed = ee;

        lin_init_array(verts, double4, cv, ev);
        lin_init_array(tris, int3, ct, et);
        lin_init_array(edges, int3, ce, ee);
        lin_init_array(dashes, int2, cd, ed);

        this->info = (int4*) malloc(sizeof(int4) * cv);

        // initialize the hash table
        int size = 0x1000;
        while (size*2 < cv) size *= 2;
        this->hash_table = (int*) malloc(sizeof(int) * size);
        memset(this->hash_table, 0xff, sizeof(int) * size);
        mask = size-1;


        // select the linearization quadrature
        Quad2D *old_quad_x, *old_quad_y;
        old_quad_x = xsln->get_quad_2d();
        old_quad_y = ysln->get_quad_2d();

        xsln->set_quad_2d((Quad2D*) &quad_lin);
        ysln->set_quad_2d((Quad2D*) &quad_lin);

        if (!xitem) error("Parameter 'xitem' cannot be zero.");
        if (!yitem) error("Parameter 'yitem' cannot be zero.");
        get_gv_a_b(xitem, xia, xib);
        get_gv_a_b(yitem, yia, yib);
        if (xib >= 6) error("Ithis->nvalid value of paremeter 'xitem'.");
        if (yib >= 6) error("Ithis->nvalid value of paremeter 'yitem'.");

        max = 1e-10;
        trav.begin(2, meshes, fns);
        Element** e;
        while ((e = trav.get_next_state(NULL, NULL)) != NULL)
        {
          xsln->set_quad_order(0, xitem);
          ysln->set_quad_order(0, yitem);
          double* xval = xsln->get_values(xia, xib);
          double* yval = ysln->get_values(yia, yib);

          for (unsigned int i = 0; i < e[0]->nvert; i++)
          {
            double fx = xval[i];
            double fy = yval[i];
            if (fabs(sqrt(fx*fx + fy*fy)) > max) max = fabs(sqrt(fx*fx + fy*fy));
          }
        }
        trav.finish();

        trav.begin(2, meshes, fns);
        // process all elements of the mesh
        while ((e = trav.get_next_state(NULL, NULL)) != NULL)
        {
          xsln->set_quad_order(0, xitem);
          ysln->set_quad_order(0, yitem);
          double* xval = xsln->get_values(xia, xib);
          double* yval = ysln->get_values(yia, yib);

          double* x = xsln->get_refmap()->get_phys_x(0);
          double* y = ysln->get_refmap()->get_phys_y(0);

          int iv[4];
          for (unsigned int i = 0; i < e[0]->nvert; i++)
          {
            double fx = xval[i];
            double fy = yval[i];
            iv[i] = create_vertex(x[i], y[i], fx, fy);
          }

          // we won't bother calculating physical coordinates from the refmap if this is not a curved element
          curved = (e[0]->cm != NULL);

          // recur to sub-elements
          if (e[0]->is_triangle())
            process_triangle(iv[0], iv[1], iv[2], 0, NULL, NULL, NULL, NULL, NULL);
          else
            process_quad(iv[0], iv[1], iv[2], iv[3], 0, NULL, NULL, NULL, NULL, NULL);

          // process edges and dashes (bold line for edge in both meshes, dashed line for edge in one of the meshes)
          Trf* xctm = xsln->get_ctm();
          Trf* yctm = ysln->get_ctm();
          double r[4] = { -1.0, 1.0, 1.0, -1.0 };
          double ref[4][2] = { {-1.0,-1.0}, {1.0,-1.0}, {1.0,1.0}, {-1.0,1.0} };
          for (unsigned int i = 0; i < e[0]->nvert; i++)
          {
            bool bold = false;
            double px = ref[i][0];
            double py = ref[i][1];
            // for odd edges (1, 3) we check x coordinate after ctm transformation, if it's the same (1 or -1) in both meshes => bold
            if (i & 1) 
            {
              if ((xctm->m[0]*px + xctm->t[0] == r[i]) && (yctm->m[0]*px + yctm->t[0] == r[i]))
                bold = true;
            }
            // for even edges (0, 4) we check y coordinate after ctm transformation, if it's the same (-1 or 1) in both meshes => bold
            else 
            {
              if ((xctm->m[1]*py + xctm->t[1] == r[i]) && (yctm->m[1]*py + yctm->t[1] == r[i]))
                bold = true;
            }
            int j = e[0]->next_vert(i);
            // we draw a line only if both edges lies on the boundary or if the line is from left top to right bottom
            if (((e[0]->en[i]->bnd) && (e[1]->en[i]->bnd)) ||
              (verts[iv[i]][1] < verts[iv[j]][1]) ||
              (verts[iv[i]][1] == verts[iv[j]][1] && verts[iv[i]][0] < verts[iv[j]][0]))
            {
              if (bold)
                process_edge(iv[i], iv[j], e[0]->en[i]->marker);
              else
                process_dash(iv[i], iv[j]);
            }
          }
        }
        trav.finish();

        find_min_max();

        verbose("Vectorizer<Scalar>created %d verts and %d tris in %0.3g s", this->nv, nt, cpu_time.tick().last());
        //if (verbose_mode) print_hash_stats();
        unlock_data();

        // select old quadratrues
        xsln->set_quad_2d(old_quad_x);
        ysln->set_quad_2d(old_quad_y);

        // clean up
        ::free(this->hash_table);
        ::free(this->info);

      }
      template<>
      void Vectorizer<std::complex<double> >::process_solution(MeshFunction<std::complex<double> >* xsln, int xitem, MeshFunction<std::complex<double> >* ysln, int yitem, double eps)
      {
        // sanity check
        if (xsln == NULL || ysln == NULL) error("One of the solutions is NULL in Vectorizer:process_solution().");


        lock_data();
        Hermes::TimePeriod cpu_time;

        // initialization
        this->xsln = xsln;
        this->ysln = ysln;
        this->xitem = xitem;
        this->yitem = yitem;
        this->eps = eps;
        this->nv = nt = ne = nd = 0;
        del_slot = -1;

        Mesh* meshes[2] = { xsln->get_mesh(), ysln->get_mesh() };
        if (meshes[0] == NULL || meshes[1] == NULL) 
        {
          error("One of the meshes is NULL in Vectorizer:process_solution().");
        }

        Transformable* fns[2] = { xsln, ysln };
        Traverse trav;

        // estimate the required number of vertices and triangles
        // (based on the assumption that the linear mesh will be
        // about four-times finer than the original mesh).
        int nn = meshes[0]->get_num_elements() + meshes[1]->get_num_elements();
        int ev = std::max(32 * nn, 10000);
        int et = std::max(64 * nn, 20000);
        int ee = std::max(24 * nn, 7500);
        int ed = ee;

        lin_init_array(verts, double4, cv, ev);
        lin_init_array(tris, int3, ct, et);
        lin_init_array(edges, int3, ce, ee);
        lin_init_array(dashes, int2, cd, ed);

        this->info = (int4*) malloc(sizeof(int4) * cv);

        // initialize the hash table
        int size = 0x1000;
        while (size*2 < cv) size *= 2;
        this->hash_table = (int*) malloc(sizeof(int) * size);
        memset(this->hash_table, 0xff, sizeof(int) * size);
        mask = size-1;


        // select the linearization quadrature
        Quad2D *old_quad_x, *old_quad_y;
        old_quad_x = xsln->get_quad_2d();
        old_quad_y = ysln->get_quad_2d();

        xsln->set_quad_2d((Quad2D*) &quad_lin);
        ysln->set_quad_2d((Quad2D*) &quad_lin);

        if (!xitem) error("Parameter 'xitem' cannot be zero.");
        if (!yitem) error("Parameter 'yitem' cannot be zero.");
        get_gv_a_b(xitem, xia, xib);
        get_gv_a_b(yitem, yia, yib);
        if (xib >= 6) error("Ithis->nvalid value of paremeter 'xitem'.");
        if (yib >= 6) error("Ithis->nvalid value of paremeter 'yitem'.");

        max = 1e-10;
        trav.begin(2, meshes, fns);
        Element** e;
        while ((e = trav.get_next_state(NULL, NULL)) != NULL)
        {
          xsln->set_quad_order(0, xitem);
          ysln->set_quad_order(0, yitem);
          std::complex<double>* xval = xsln->get_values(xia, xib);
          std::complex<double>* yval = ysln->get_values(yia, yib);

          for (unsigned int i = 0; i < e[0]->nvert; i++)
          {
            double fx = xval[i].real();
            double fy = yval[i].real();
            if (fabs(sqrt(fx*fx + fy*fy)) > max) max = fabs(sqrt(fx*fx + fy*fy));
          }
        }
        trav.finish();

        trav.begin(2, meshes, fns);
        // process all elements of the mesh
        while ((e = trav.get_next_state(NULL, NULL)) != NULL)
        {
          xsln->set_quad_order(0, xitem);
          ysln->set_quad_order(0, yitem);
          std::complex<double>* xval = xsln->get_values(xia, xib);
          std::complex<double>* yval = ysln->get_values(yia, yib);

          double* x = xsln->get_refmap()->get_phys_x(0);
          double* y = ysln->get_refmap()->get_phys_y(0);

          int iv[4];
          for (unsigned int i = 0; i < e[0]->nvert; i++)
          {
            double fx = xval[i].real();
            double fy = yval[i].real();
            iv[i] = create_vertex(x[i], y[i], fx, fy);
          }

          // we won't bother calculating physical coordinates from the refmap if this is not a curved element
          curved = (e[0]->cm != NULL);

          // recur to sub-elements
          if (e[0]->is_triangle())
            process_triangle(iv[0], iv[1], iv[2], 0, NULL, NULL, NULL, NULL, NULL);
          else
            process_quad(iv[0], iv[1], iv[2], iv[3], 0, NULL, NULL, NULL, NULL, NULL);

          // process edges and dashes (bold line for edge in both meshes, dashed line for edge in one of the meshes)
          Trf* xctm = xsln->get_ctm();
          Trf* yctm = ysln->get_ctm();
          double r[4] = { -1.0, 1.0, 1.0, -1.0 };
          double ref[4][2] = { {-1.0,-1.0}, {1.0,-1.0}, {1.0,1.0}, {-1.0,1.0} };
          for (unsigned int i = 0; i < e[0]->nvert; i++)
          {
            bool bold = false;
            double px = ref[i][0];
            double py = ref[i][1];
            // for odd edges (1, 3) we check x coordinate after ctm transformation, if it's the same (1 or -1) in both meshes => bold
            if (i & 1) 
            {
              if ((xctm->m[0]*px + xctm->t[0] == r[i]) && (yctm->m[0]*px + yctm->t[0] == r[i]))
                bold = true;
            }
            // for even edges (0, 4) we check y coordinate after ctm transformation, if it's the same (-1 or 1) in both meshes => bold
            else 
            {
              if ((xctm->m[1]*py + xctm->t[1] == r[i]) && (yctm->m[1]*py + yctm->t[1] == r[i]))
                bold = true;
            }
            int j = e[0]->next_vert(i);
            // we draw a line only if both edges lies on the boundary or if the line is from left top to right bottom
            if (((e[0]->en[i]->bnd) && (e[1]->en[i]->bnd)) ||
              (verts[iv[i]][1] < verts[iv[j]][1]) ||
              (verts[iv[i]][1] == verts[iv[j]][1] && verts[iv[i]][0] < verts[iv[j]][0]))
            {
              if (bold)
                process_edge(iv[i], iv[j], e[0]->en[i]->marker);
              else
                process_dash(iv[i], iv[j]);
            }
          }
        }
        trav.finish();

        find_min_max();

        verbose("Vectorizer<Scalar>created %d verts and %d tris in %0.3g s", this->nv, nt, cpu_time.tick().last());
        //if (verbose_mode) print_hash_stats();
        unlock_data();

        // select old quadratrues
        xsln->set_quad_2d(old_quad_x);
        ysln->set_quad_2d(old_quad_y);

        // clean up
        ::free(this->hash_table);
        ::free(this->info);

      }

      template<typename Scalar>
      void Vectorizer<Scalar>::save_data(const char* filename)
      {
        FILE* f = fopen(filename, "wb");
        if (f == NULL) error("Could not open %s for writing.", filename);
        this->lock_data();

        if (fwrite("H2DV\001\000\000\000", 1, 8, f) != 8 ||
          fwrite(&this->nv, sizeof(int), 1, f) != 1 ||
          fwrite(this->verts, sizeof(double4), this->nv, f) != (unsigned) this->nv ||
          fwrite(&this->nt, sizeof(int), 1, f) != 1 ||
          fwrite(this->tris, sizeof(int3), this->nt, f) != (unsigned) this->nt ||
          fwrite(&this->ne, sizeof(int), 1, f) != 1 ||
          fwrite(this->edges, sizeof(int3), this->ne, f) != (unsigned) this->ne ||
          fwrite(&this->nd, sizeof(int), 1, f) != 1 ||
          fwrite(this->dashes, sizeof(int2), this->nd, f) != (unsigned) this->nd)
        {
          error("Error writing data to %s", filename);
        }

        this->unlock_data();
        fclose(f);
      }

      template<typename Scalar>
      void Vectorizer<Scalar>::load_data(const char* filename)
      {
        FILE* f = fopen(filename, "rb");
        if (f == NULL) error("Could not open %s for reading.", filename);
        this->lock_data();

        struct { char magic[4]; int ver; } hdr;
        if (fread(&hdr, sizeof(hdr), 1, f) != 1)
          error("Error reading %s", filename);

        if (hdr.magic[0] != 'H' || hdr.magic[1] != '2' || hdr.magic[2] != 'D' || hdr.magic[3] != 'V')
          error("File %s is not a Hermes2D Vectorizer<Scalar>file.", filename);
        if (hdr.ver > 1)
          error("File %s -- unsupported file version.", filename);

#define read_array(array, type, n, c, what) \
  if (fread(&n, sizeof(int), 1, f) != 1) \
  error("Error reading the number of " what " from %s", filename); \
  lin_init_array(array, type, c, n); \
  if (fread(array, sizeof(type), n, f) != (unsigned) n) \
  error("Error reading " what " from %s", filename);

        read_array(this->verts, double4, this->nv, this->cv, "vertices");
        read_array(this->tris,  int3,    this->nt, this->ct, "triangles");
        read_array(this->edges, int3,    this->ne, this->ce, "edges");
        read_array(this->dashes, int2,   this->nd, this->cd, "dashes");

        find_min_max();
        this->unlock_data();
        fclose(f);
      }

      template<typename Scalar>
      Vectorizer<Scalar>::~Vectorizer()
      {
        lin_free_array(verts, this->nv, this->cv);
        // fixme: free dashes?
      }

      template<typename Scalar>
      void Vectorizer<Scalar>::calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const 
      {
        assert_msg(verts != NULL, "Cannot calculate AABB from NULL vertices");
        calc_aabb(&verts[0][0], &verts[0][1], sizeof(double4), this->nv, min_x, max_x, min_y, max_y);
      }
      template class HERMES_API Vectorizer<double>;
      template class HERMES_API Vectorizer<std::complex<double> >;
    }
  }
}
