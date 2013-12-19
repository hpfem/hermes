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

#ifndef NOGLUT

#include <GL/freeglut.h>
#include "global.h"
#include "vector_view.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      VectorView::VectorView(const char* title, WinGeom* wg)
        : View(title, wg)
      {
        vec = new Vectorizer(OpenGL);
        gx = gy = 0.0;
        gs = 20.0;
        hexa = true;
        mode = 0;
        lines = false;
        pmode = false;
        length_coef = 1.0;
      }

      VectorView::VectorView(char* title, WinGeom* wg)
        : View(title, wg)
      {
        vec = new Vectorizer(OpenGL);
        gx = gy = 0.0;
        gs = 20.0;
        hexa = true;
        mode = 0;
        lines = false;
        pmode = false;
        length_coef = 1.0;
      }

      VectorView::~VectorView()
      {
        delete vec;
      }

      void VectorView::show(MeshFunctionSharedPtr<double> vsln)
      {
        if (vsln->get_num_components() < 2)
          throw Hermes::Exceptions::Exception("The single-argument version of show() is only for vector-valued solutions.");
        show(vsln, vsln, H2D_FN_VAL_0, H2D_FN_VAL_1);
      }

      void VectorView::show(MeshFunctionSharedPtr<double> xsln, MeshFunctionSharedPtr<double> ysln, int xitem, int yitem, MeshFunctionSharedPtr<double> xdisp, MeshFunctionSharedPtr<double> ydisp, double dmult)
      {
        if (xsln != nullptr && ysln != nullptr && xsln == ysln)
          this->warn("Identical solutions passed to the two-argument version of show(). Most likely this is a mistake.");

        vec->set_displacement(xdisp, ydisp, dmult);

        vec->lock_data();
        MeshFunctionSharedPtr<double> slns[2] = { xsln, ysln };
        int items[2] = { xitem, yitem };
        vec->process_solution(slns, items);

        if (range_auto)
        {
          range_min = vec->get_min_value();
          range_max = vec->get_max_value();
        }
        vec->calc_vertices_aabb(&vertices_min_x, &vertices_max_x, &vertices_min_y, &vertices_max_y);
        vec->unlock_data();

        create();
        update_layout();
        reset_view(false);

        refresh();

        wait_for_draw();
      }

      static int n_vert(int i) { return (i + 1) % 3; }
      static int p_vert(int i) { return (i + 2) % 3; }

      void VectorView::set_mode(int mode)
      {
        this->mode = mode % 3;
        refresh();
      }

      Vectorizer* VectorView::get_vectorizer()
      {
        return this->vec;
      }

      void VectorView::plot_arrow(double x, double y, double xval, double yval, double max, double min, double gs)
      {
        if (mode == 1)
          glColor3f(0.0f, 0.0f, 0.0f);
        else
          glColor3f(0.5f, 0.5f, 0.5f);

        // magnitude
        double Real_mag = sqrt(sqr(xval) + sqr(yval));
        double mag = Real_mag;
        if (Real_mag > max) mag = max;
        double length = mag / max * gs * length_coef;
        double width = 0.1 * gs;
        if (mode == 1) width *= 1.2;
        double xnew = x + gs * xval * mag / (max * Real_mag) * length_coef;
        double ynew = y - gs * yval * mag / (max * Real_mag) * length_coef;

        if ((mag) / (max - min) < 1e-5)
        {
          glTranslated(x, y, 0.0);

          glBegin(GL_QUADS);
          glVertex2d(width, width);
          glVertex2d(width, -width);
          glVertex2d(-width, -width);
          glVertex2d(-width, width);
          glEnd();
        }
        else
        {
          glBegin(GL_LINES);
          glVertex2d(x, y);
          glVertex2d(xnew, ynew);
          glEnd();

          glTranslated(x, y, 0.0);
          glRotated(atan2(-yval, xval) * 180.0 / M_PI, 0.0, 0.0, 1.0);

          glBegin(GL_TRIANGLES);
          glVertex2d(length + 3 * width, 0.0);
          glVertex2d(length - 2 * width, width);
          glVertex2d(length - 2 * width, -width);
          glEnd();
        }

        glLoadIdentity();

        if (mode == 1)
        {
          float color[3];
          get_palette_color((mag - min) / (max - min), color); //  0.0 -- 1.0
          glColor3f(color[0], color[1], color[2]);

          if (mag / (max - min) < 1e-5)
          {
            glBegin(GL_QUADS);
            glVertex2d(width, width);
            glVertex2d(width, -width);
            glVertex2d(-width, -width);
            glVertex2d(-width, width);
            glEnd();
          }
          else
          {
            glBegin(GL_LINES);
            glVertex2d(x, y);
            glVertex2d(xnew, ynew);
            glEnd();

            glTranslated(x - 1, y, 0.0);
            glRotated(atan2(-yval, xval) * 180.0 / M_PI, 0.0, 0.0, 1.0);

            glBegin(GL_TRIANGLES);
            glVertex2d(length + 3 * width, 0.0);
            glVertex2d(length - 2 * width, width);
            glVertex2d(length - 2 * width, -width);
            glEnd();

            glLoadIdentity();
          }
        }
      }

      void VectorView::draw_edges_2d()
      {
        glColor3f(0.5, 0.5, 0.5);
        glBegin(GL_LINES);
        for (Vectorizer::Iterator<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t> it = vec->edges_begin(); !it.end; ++it)
        {
          VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t& edge = it.get();
          int& edge_marker = it.get_marker();

          if (lines || edge_marker)
          {
            glVertex2d(transform_x(edge[0][0]), transform_y(edge[0][1]));
            glVertex2d(transform_x(edge[1][0]), transform_y(edge[1][1]));
          }
        }
        glEnd();
      }

      void VectorView::on_display()
      {
        set_ortho_projection();
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_TEXTURE_1D);
        glPolygonMode(GL_FRONT_AND_BACK, pmode ? GL_LINE : GL_FILL);

        // initial grid point and grid step
        double gt = gs;
        if (hexa) gt *= sqrt(3.0) / 2.0;

        double max_length = 0.0;

        // transform all vertices
        vec->lock_data();

        // value range
        double min = range_min, max = range_max;
        if (range_auto) { min = vec->get_min_value(); max = vec->get_max_value(); }
        double irange = 1.0 / (max - min);
        // special case: constant solution
        if (fabs(min - max) < Hermes::HermesEpsilon) { irange = 1.0; min -= 0.5; }

        // draw all triangles
        if (mode != 1)
          glEnable(GL_TEXTURE_1D);
        glBindTexture(GL_TEXTURE_1D, gl_pallete_tex_id);
        glBegin(GL_TRIANGLES);
        glColor3f(0.95f, 0.95f, 0.95f);

        for (Vectorizer::Iterator<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t> it = vec->triangles_begin(); !it.end; ++it)
        {
          VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t& triangle = it.get();

          double mag = sqrt(sqr(triangle[0][2]) + sqr(triangle[0][3]));
          glTexCoord2d((mag - min) * irange * tex_scale + tex_shift, 0.0);
          glVertex2d(transform_x(triangle[0][0]), transform_y(triangle[0][1]));

          double length = sqr(triangle[0][2]) + sqr(triangle[0][3]);
          if (length > max_length)
            max_length = length;

          mag = sqrt(sqr(triangle[1][2]) + sqr(triangle[1][3]));
          glTexCoord2d((mag - min) * irange * tex_scale + tex_shift, 0.0);
          glVertex2d(transform_x(triangle[1][0]), transform_y(triangle[1][1]));

          length = sqr(triangle[1][2]) + sqr(triangle[1][3]);
          if (length > max_length)
            max_length = length;

          mag = sqrt(sqr(triangle[2][2]) + sqr(triangle[2][3]));
          glTexCoord2d((mag - min) * irange * tex_scale + tex_shift, 0.0);
          glVertex2d(transform_x(triangle[2][0]), transform_y(triangle[2][1]));

          length = sqr(triangle[2][2]) + sqr(triangle[2][3]);
          if (length > max_length)
            max_length = length;
        }

        max_length = sqrt(max_length);
        glEnd();
        glDisable(GL_TEXTURE_1D);

        draw_edges_2d();

        // draw arrows
        if (mode != 2)
        {
          for (Vectorizer::Iterator<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t> it = vec->triangles_begin(); !it.end; ++it)
          {
            VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t& triangle = it.get();
            {
              double miny = 1e100;
              int idx = -1, k, l1, l2, r2, r1, s;
              double lry, x;
              double mr, ml, lx, rx, xval, yval;

              double wh = output_height + gt, ww = output_width + gs;
              if ((triangle[0][0] < -gs) && (triangle[1][0] < -gs) && (triangle[2][0] < -gs)) continue;
              if ((triangle[0][0] >  ww) && (triangle[1][0] >  ww) && (triangle[2][0] >  ww)) continue;
              if ((triangle[0][1] < -gt) && (triangle[1][1] < -gt) && (triangle[2][1] < -gt)) continue;
              if ((triangle[0][1] >  wh) && (triangle[1][1] >  wh) && (triangle[2][1] >  wh)) continue;

              // find vertex with min y-coordinate
              for (k = 0; k < 3; k++)
              if (triangle[k][1] < miny)
              {
                idx = k;
                miny = triangle[k][1];
              }
              l1 = r1 = idx;
              l2 = n_vert(idx);
              r2 = p_vert(idx);

              double tx_l1_0 = transform_x(triangle[l1][0]);
              double tx_l1_1 = transform_x(triangle[l1][1]);

              double tx_l2_0 = transform_x(triangle[l2][0]);
              double tx_l2_1 = transform_x(triangle[l2][1]);

              double tx_r1_0 = transform_x(triangle[r1][0]);
              double tx_r1_1 = transform_x(triangle[r1][1]);

              double tx_r2_0 = transform_x(triangle[r2][0]);
              double tx_r2_1 = transform_x(triangle[r2][1]);

              double ty_l1_0 = transform_y(triangle[l1][0]);
              double ty_l1_1 = transform_y(triangle[l1][1]);

              double ty_l2_0 = transform_y(triangle[l2][0]);
              double ty_l2_1 = transform_y(triangle[l2][1]);

              double ty_r1_0 = transform_y(triangle[r1][0]);
              double ty_r1_1 = transform_y(triangle[r1][1]);

              double ty_r2_0 = transform_y(triangle[r2][0]);
              double ty_r2_1 = transform_y(triangle[r2][1]);

              // plane of x and y values on triangle
              double a[2], b[2], c[2], d[2];
              for (int n = 0; n < 2; n++)
              {
                a[n] = (ty_l1_1-ty_l2_1) * (triangle[r1][2 + n] - triangle[r2][2 + n])
                  - (triangle[l1][2 + n] - triangle[l2][2 + n]) * (ty_r1_1-ty_r2_1);
                b[n] = (triangle[l1][2 + n] - triangle[l2][2 + n]) * (tx_r1_0-tx_r2_0)
                  - (tx_l1_0-tx_l2_0) * (triangle[r1][2 + n] - triangle[r2][2 + n]);
                c[n] = (tx_l1_0-tx_l2_0) * (ty_r1_1-ty_r2_1) - (ty_l1_1-ty_l2_1) * (tx_r1_0-tx_r2_0);
                d[n] = -a[n] * tx_l1_0 - b[n] * ty_l1_1 - c[n] * triangle[l1][2 + n];
                a[n] /= c[n]; b[n] /= c[n]; d[n] /= c[n];
              }

              s = (int)ceil((ty_l1_1-gy) / gt);  // first step
              lry = gy + s*gt;
              bool shift = hexa && (s & 1);

              // if there are two points with min y-coordinate, switch to the next segment
              if (ty_l1_1 == ty_l2_1 || ty_r1_1 == ty_r2_1)
              {
                if (ty_l1_1 == ty_l2_1)
                {
                  l1 = l2; l2 = r2;
                }
                else if (ty_r1_1 == ty_r2_1)
                {
                  r1 = r2; r2 = l2;
                }
              }

              // slope of the left and right segment
              ml = (tx_l1_0-tx_l2_0) / (ty_l1_1-ty_l2_1);
              mr = (tx_r1_0-tx_r2_0) / (ty_r1_1-ty_r2_1);
              // x-coordinates of the endpoints of the first line
              lx = tx_l1_0 + ml * (lry - (ty_l1_1));
              rx = tx_r1_0 + mr * (lry - (ty_r1_1));

              if (lry < -gt)
              {
                k = (int)floor(-lry / gt);
                lry += gt * k;
                lx += k * ml * gt;
                rx += k * mr * gt;
              }

              // while we are in triangle
              while (((lry < ty_l2_1) || (lry < ty_r2_1)) && (lry < wh))
              {
                // while we are in the segment
                while (((lry <= ty_l2_1) && (lry <= ty_r2_1)) && (lry < wh))
                {
                  double gz = gx;
                  if (shift) gz -= 0.5*gs;
                  s = (int)ceil((lx - gz) / gs);
                  x = gz + s*gs;
                  if (hexa) shift = !shift;

                  if (x < -gs)
                  {
                    k = (int)floor(-x / gs);
                    x += gs * k;
                  }
                  // go along the line
                  while ((x < rx) && (x < ww))
                  {
                    // plot the arrow
                    xval = -a[0] * x - b[0] * lry - d[0];
                    yval = -a[1] * x - b[1] * lry - d[1];
                    plot_arrow(x, lry, xval, yval, max, min, gs);
                    x += gs;
                  }
                  // move to the next line
                  lx += ml*gt;
                  rx += mr*gt;
                  lry += gt;
                }
                // change segment
                if (lry >= ty_l2_1)
                {
                  l1 = l2; l2 = r2;
                  ml = (tx_l1_0-tx_l2_0) / (ty_l1_1-ty_l2_1);
                  lx = tx_l1_0 + ml * (lry - (ty_l1_1));
                }
                else
                {
                  r1 = r2; r2 = l2;
                  mr = (tx_r1_0-tx_r2_0) / (ty_r1_1-ty_r2_1);
                  rx = tx_r1_0 + mr * (lry - (ty_r1_1));
                }
              }
            }
          }

          vec->unlock_data();
        }
      }

      void VectorView::on_mouse_move(int x, int y)
      {
        if (dragging)
        {
          gx += (x - mouse_x);
          gy += (y - mouse_y);
        }
        View::on_mouse_move(x, y);
      }

      void VectorView::on_key_down(unsigned char key, int x, int y)
      {
        switch (key)
        {
        case 'm':
          lines = !lines;
          refresh();
          break;

        case 'l':
          pmode = !pmode;
          refresh();
          break;

        case 'c':
          reset_view(true);
          refresh();
          break;

        case 'f':
          set_palette_filter(pal_filter != GL_LINEAR);
          break;

        case 'x':
          set_grid_type(!hexa);
          break;

        case 'b':
          mode++;
          if (mode > 2) mode = 0;
          refresh();
          break;

        case '*':
        case '/':
          if (key == '*') length_coef *= 1.1; else length_coef /= 1.1;
          refresh();
          break;

        default:
          View::on_key_down(key, x, y);
          break;
        }
      }

      const char* VectorView::get_help_text() const
      {
        return
          "VectorView\n\n"
          "Controls:\n"
          "  Left mouse - pan\n"
          "  Right mouse - zoom\n"
          "  B - toggle view mode (type of arrows x no arrows)\n"
          "  * - extend arrows\n"
          "  / - shorten arrows\n"
          "  C - center image\n"
          "  F - toggle smooth palette\n"
          "  X - toggle hexagonal grid\n"
          "  H - render high-quality frame\n"
          "  M - toggle mesh\n"
          "  P - cycle palettes\n"
          "  S - save screenshot\n"
          "  F1 - this help\n"
          "  Esc, Q - quit";
      }
    }
  }
}
#endif