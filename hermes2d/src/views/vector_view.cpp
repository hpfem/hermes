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
        : View(title, wg), vec(NULL)
      {
        gx = gy = 0.0;
        gs = 20.0;
        hexa = true;
        mode = 0;
        lines = false;
        pmode = false;
        length_coef = 1.0;
      }

      VectorView::VectorView(char* title, WinGeom* wg)
        : View(title, wg), vec(NULL)
      {
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

      void VectorView::show(MeshFunction<double>* vsln, double eps)
      {
        if(vec == NULL)
          vec = new Vectorizer;
        if(vsln->get_num_components() < 2)
          throw Hermes::Exceptions::Exception("The single-argument version of show() is only for vector-valued solutions.");
        show(vsln, vsln, eps, H2D_FN_VAL_0, H2D_FN_VAL_1);
      }

      void VectorView::show(MeshFunction<double>* xsln, MeshFunction<double>* ysln, double eps)
      {
        if(vec == NULL)
          vec = new Vectorizer;
        if(xsln == ysln)
          this->warn("Identical solutions passed to the two-argument version of show(). Most likely this is a mistake.");
        show(xsln, ysln, eps, H2D_FN_VAL_0, H2D_FN_VAL_0);
      }

      void VectorView::show(MeshFunction<double>* xsln, MeshFunction<double>* ysln, double eps, int xitem, int yitem)
      {
        if(vec == NULL)
          vec = new Vectorizer;

        vec->lock_data();
        vec->process_solution(xsln, ysln, xitem, yitem, eps);
        if(range_auto) 
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
        if(mode == 1)
          glColor3f(0.0f, 0.0f, 0.0f);
        else
          glColor3f(0.5f, 0.5f, 0.5f);

        // magnitude
        double Real_mag = sqrt(sqr(xval) + sqr(yval));
        double mag = Real_mag;
        if(Real_mag > max) mag = max;
        double length = mag/max * gs * length_coef;
        double width = 0.1 * gs;
        if(mode == 1) width *= 1.2;
        double xnew = x + gs * xval * mag / (max * Real_mag) * length_coef;
        double ynew = y - gs * yval * mag / (max * Real_mag) * length_coef;

        if((mag)/(max - min) < 1e-5)
        {
          glTranslated(x, y, 0.0);

          glBegin(GL_QUADS);
          glVertex2d( width,  width);
          glVertex2d( width, -width);
          glVertex2d(-width, -width);
          glVertex2d(-width,  width);
          glEnd();
        }
        else
        {
          glBegin(GL_LINES);
          glVertex2d(x, y);
          glVertex2d(xnew, ynew);
          glEnd();

          glTranslated(x, y, 0.0);
          glRotated(atan2(-yval, xval) * 180.0/M_PI, 0.0, 0.0, 1.0);

          glBegin(GL_TRIANGLES);
          glVertex2d(length + 3 * width,  0.0);
          glVertex2d(length - 2 * width,  width);
          glVertex2d(length - 2 * width, -width);
          glEnd();
        }

        glLoadIdentity();

        if(mode == 1)
        {
          float color[3];
          get_palette_color((mag - min)/(max - min), color); //  0.0 -- 1.0
          glColor3f(color[0], color[1], color[2]);

          if(mag/(max - min) < 1e-5)
          {
            glBegin(GL_QUADS);
            glVertex2d( width,  width);
            glVertex2d( width, -width);
            glVertex2d(-width, -width);
            glVertex2d(-width,  width);
            glEnd();
          }
          else
          {
            glBegin(GL_LINES);
            glVertex2d(x, y);
            glVertex2d(xnew, ynew);
            glEnd();

            glTranslated(x - 1, y, 0.0);
            glRotated(atan2(-yval, xval) * 180.0/M_PI, 0.0, 0.0, 1.0);

            glBegin(GL_TRIANGLES);
            glVertex2d(length + 3 * width,  0.0);
            glVertex2d(length - 2 * width,  width);
            glVertex2d(length - 2 * width, -width);
            glEnd();

            glLoadIdentity();
          }
        }
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
        if(hexa) gt *= sqrt(3.0)/2.0;

        double max_length = 0.0;

        // transform all vertices
        vec->lock_data();
        int i;
        int nv = vec->get_num_vertices();
        double4* vert = vec->get_vertices();
        double2* tvert = new double2[nv];

        for (i = 0; i < nv; i++)
        {
          tvert[i][0] = transform_x(vert[i][0]);
          tvert[i][1] = transform_y(vert[i][1]);

          // find max length of vectors
          double length = sqr(vert[i][2]) + sqr(vert[i][3]);
          if(length > max_length) max_length = length;
        }
        max_length = sqrt(max_length);

        // value range
        double min = range_min, max = range_max;
        if(range_auto) { min = vec->get_min_value(); max = vec->get_max_value(); }
        double irange = 1.0 / (max - min);
        // special case: constant solution
        if(fabs(min - max) < 1e-8) { irange = 1.0; min -= 0.5; }

        // draw all triangles
        int3* xtris = vec->get_triangles();

        if(mode != 1) glEnable(GL_TEXTURE_1D);
        glBindTexture(GL_TEXTURE_1D, gl_pallete_tex_id);
        glBegin(GL_TRIANGLES);
        glColor3f(0.95f, 0.95f, 0.95f);
        for (i = 0; i < vec->get_num_triangles(); i++)
        {
          double mag = sqrt(sqr(vert[xtris[i][0]][2]) + sqr(vert[xtris[i][0]][3]));
          glTexCoord2d((mag -min) * irange * tex_scale + tex_shift, 0.0);
          glVertex2d(tvert[xtris[i][0]][0], tvert[xtris[i][0]][1]);

          mag = sqrt(sqr(vert[xtris[i][1]][2]) + sqr(vert[xtris[i][1]][3]));
          glTexCoord2d((mag -min) * irange * tex_scale + tex_shift, 0.0);
          glVertex2d(tvert[xtris[i][1]][0], tvert[xtris[i][1]][1]);

          mag = sqrt(sqr(vert[xtris[i][2]][2]) + sqr(vert[xtris[i][2]][3]));
          glTexCoord2d((mag -min) * irange * tex_scale + tex_shift, 0.0);
          glVertex2d(tvert[xtris[i][2]][0], tvert[xtris[i][2]][1]);
        }
        glEnd();
        glDisable(GL_TEXTURE_1D);

        // draw all edges
        /*if(mode == 0) glColor3f(0.3, 0.3, 0.3);
        else*/ glColor3f(0.5, 0.5, 0.5);
        glBegin(GL_LINES);
        int2* edges = vec->get_edges();
        for (i = 0; i < vec->get_num_edges(); i++)
        {
          if(lines || edges[i][2] != 0)
          {
            glVertex2d(tvert[edges[i][0]][0], tvert[edges[i][0]][1]);
            glVertex2d(tvert[edges[i][1]][0], tvert[edges[i][1]][1]);
          }
        }
        glEnd();

        // draw dashed edges
        if(lines)
        {
          glEnable(GL_LINE_STIPPLE);
          glLineStipple(1, 0xCCCC);
          glBegin(GL_LINES);
          int2* dashes = vec->get_dashes();
          for (i = 0; i < vec->get_num_dashes(); i++)
          {
            glVertex2d(tvert[dashes[i][0]][0], tvert[dashes[i][0]][1]);
            glVertex2d(tvert[dashes[i][1]][0], tvert[dashes[i][1]][1]);
          }
          glEnd();
          glDisable(GL_LINE_STIPPLE);
        }

        // draw arrows
        if(mode != 2)
        {
          for (i = 0; i < vec->get_num_triangles(); i++)
          {
            double miny = 1e100;
            int idx = -1, k, l1, l2, r2, r1, s;
            double lry, x;
            double mr, ml, lx, rx, xval, yval;

            double wh = output_height + gt, ww = output_width + gs;
            if((tvert[xtris[i][0]][0] < -gs) && (tvert[xtris[i][1]][0] < -gs) && (tvert[xtris[i][2]][0] < -gs)) continue;
            if((tvert[xtris[i][0]][0] >  ww) && (tvert[xtris[i][1]][0] >  ww) && (tvert[xtris[i][2]][0] >  ww)) continue;
            if((tvert[xtris[i][0]][1] < -gt) && (tvert[xtris[i][1]][1] < -gt) && (tvert[xtris[i][2]][1] < -gt)) continue;
            if((tvert[xtris[i][0]][1] >  wh) && (tvert[xtris[i][1]][1] >  wh) && (tvert[xtris[i][2]][1] >  wh)) continue;

            // find vertex with min y-coordinate
            for (k = 0; k < 3; k++)
              if(tvert[xtris[i][k]][1] < miny)
                miny = tvert[xtris[i][idx = k]][1];
            l1 = r1 = xtris[i][idx];
            l2 = xtris[i][n_vert(idx)];
            r2 = xtris[i][p_vert(idx)];

            // plane of x and y values on triangle
            double a[2], b[2], c[2], d[2];
            for (int n = 0; n < 2; n++)
            {
              a[n] = (tvert[l1][1] - tvert[l2][1])*(vert[r1][2 +n] - vert[r2][2 + n]) - (vert[l1][2 + n] - vert[l2][2 + n])*(tvert[r1][1] - tvert[r2][1]);
              b[n] = (vert[l1][2 + n] - vert[l2][2 + n])*(tvert[r1][0] - tvert[r2][0]) - (tvert[l1][0] - tvert[l2][0])*(vert[r1][2 + n] - vert[r2][2 + n]);
              c[n] = (tvert[l1][0] - tvert[l2][0])*(tvert[r1][1] - tvert[r2][1]) - (tvert[l1][1] - tvert[l2][1])*(tvert[r1][0] - tvert[r2][0]);
              d[n] = -a[n] * tvert[l1][0] - b[n] * tvert[l1][1] - c[n] * vert[l1][2 + n];
              a[n] /= c[n]; b[n] /= c[n]; d[n] /= c[n];
            }

            s = (int) ceil((tvert[l1][1] - gy)/gt);  // first step
            lry = gy + s*gt;
            bool shift = hexa && (s & 1);

            // if there are two points with min y-coordinate, switch to the next segment
            if((tvert[l1][1] == tvert[l2][1]) || (tvert[r1][1] == tvert[r2][1]))
            {
              if(tvert[l1][1] == tvert[l2][1])
              {l1 = l2; l2 = r2;}
              else if(tvert[r1][1] == tvert[r2][1])
              {r1 = r2; r2 = l2;}
            }

            // slope of the left and right segment
            ml = (tvert[l1][0] - tvert[l2][0])/(tvert[l1][1] - tvert[l2][1]);
            mr = (tvert[r1][0] - tvert[r2][0])/(tvert[r1][1] - tvert[r2][1]);
            // x-coordinates of the endpoints of the first line
            lx = tvert[l1][0] + ml * (lry - (tvert[l1][1]));
            rx = tvert[r1][0] + mr * (lry - (tvert[r1][1]));

            if(lry < -gt)
            {
              k = (int) floor(-lry/gt);
              lry += gt * k;
              lx += k * ml * gt;
              rx += k * mr * gt;
            }

            // while we are in triangle
            while (((lry < tvert[l2][1]) || (lry < tvert[r2][1])) && (lry < wh))
            {
              // while we are in the segment
              while (((lry <= tvert[l2][1]) && (lry <= tvert[r2][1])) && (lry < wh))
              {
                double gz = gx;
                if(shift) gz -= 0.5*gs;
                s = (int) ceil((lx - gz)/gs);
                x = gz + s*gs;
                if(hexa) shift = !shift;

                if(x < -gs)
                {
                  k = (int) floor(-x/gs);
                  x += gs * k;
                }
                // go along the line
                while ((x < rx) && (x < ww))
                {
                  // plot the arrow
                  xval = -a[0]*x - b[0]*lry - d[0];
                  yval = -a[1]*x - b[1]*lry - d[1];
                  plot_arrow(x, lry, xval, yval, max, min, gs);
                  x += gs;
                }
                // move to the next line
                lx += ml*gt;
                rx += mr*gt;
                lry += gt;
              }
              // change segment
              if(lry >= tvert[l2][1])
              {
                l1 = l2; l2 = r2;
                ml = (tvert[l1][0] - tvert[l2][0])/(tvert[l1][1] - tvert[l2][1]);
                lx = tvert[l1][0] + ml * (lry - (tvert[l1][1]));
              }
              else
              {
                r1 = r2; r2 = l2;
                mr = (tvert[r1][0] - tvert[r2][0])/(tvert[r1][1] - tvert[r2][1]);
                rx = tvert[r1][0] + mr * (lry - (tvert[r1][1]));
              }
            }
          }
        }

        delete [] tvert;
        vec->unlock_data();
      }

      void VectorView::on_mouse_move(int x, int y)
      {
        if(dragging)
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
          if(mode > 2) mode = 0;
          refresh();
          break;

        case '*':
        case '/':
          if(key == '*') length_coef *= 1.1; else length_coef /= 1.1;
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