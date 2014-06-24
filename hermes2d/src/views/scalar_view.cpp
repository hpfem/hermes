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

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <algorithm>
#include <cmath>
#include <list>
#include "global.h"
#include "scalar_view.h"

#define GL_BUFFER_OFFSET(i) ((char *)nullptr + (i))

/* constants */
#define MIN_CONT_STEP 1.0E-2 ///< A minimal step of a contour.
#define CONT_CHANGE 1.0E-2 ///< A change of a contour in GUI
#define D3DV_SCALE_STEP_COEF 1.1 ///< A scale coefficient for changing contours and scaling along the Y-axis using keyboard.

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      const int ScalarView::fovy = 50;
      const double ScalarView::znear = 0.05;
      const double ScalarView::zfar = 10;

      void ScalarView::init()
      {
        lin = new Linearizer(OpenGL);
        pmode = mode3d = false;
        panning = false;

        contours = false;
        cont_orig = 0.0;
        cont_step = 0.2;
        cont_color[0] = 0.0f; cont_color[1] = 0.0f; cont_color[2] = 0.0f;

        show_edges = true;
        show_aabb = false;
        edges_color[0] = 0.5f; edges_color[1] = 0.4f; edges_color[2] = 0.4f;

        show_values = true;
        lin_updated = false;

        do_zoom_to_fit = true;
        is_constant = false;
      }

#ifndef _MSC_VER

      ScalarView::ScalarView(const char* title, WinGeom* wg) :
        View(title, wg),
        element_id_widget(0),
        show_element_info(false)
      {
        init();
      }
#else
      ScalarView::ScalarView(WinGeom* wg) :
        View("ScalarView", wg),
        element_id_widget(0),
        show_element_info(false)
      {
        init();
      }
#endif

      ScalarView::ScalarView(char* title, WinGeom* wg) :
        View(title, wg),
        element_id_widget(0),
        show_element_info(false)
      {
        init();
      }

      ScalarView::~ScalarView()
      {
        delete lin;
      }

      void ScalarView::on_close()
      {
        if (element_id_widget != 0)
        {
          glDeleteLists(element_id_widget, 1);
          element_id_widget = 0;
        }

        //call of parent implementation
        View::on_close();
      }

      void ScalarView::show(MeshFunctionSharedPtr<double> sln, int item, MeshFunctionSharedPtr<double> xdisp, MeshFunctionSharedPtr<double> ydisp, double dmult)
      {
        // For preservation of the sln's active element. Will be set back after the visualization.
        Element* active_element = sln->get_active_element();

        if (!range_auto)
          lin->set_max_absolute_value(std::max(fabs(range_min), fabs(range_max)));

        lin->set_displacement(xdisp, ydisp, dmult);
        lin->lock_data();

        lin->process_solution(sln, item);

        update_mesh_info();

        // Initialize element info.
        init_element_info(sln->get_mesh());

        lin->unlock_data();

        create();
        update_layout();
        wait_for_draw();

        // FIXME: find out why this has to be called after wait_for_draw in order for the view to be reset initially.
        reset_view(false); // setting true here makes the view always reset after calling 'show'; particularly in the adaptivity process,
        // it would disallow the observation of the process from a manually set viewpoint.
        refresh();

        // Now we reset the active element if it was set before the MeshFunction sln entered this method.
        // Only for Solutions. This method may fail for filters, as they may not have RefMaps correctly set.
        if (dynamic_cast<Solution<double>*>(sln.get()) != nullptr)
        if (active_element != nullptr)
          // Also when there has not been a call to set_active_element since assignment to this MeshFunction,
          // there is nothing to restore to.
        if (active_element->active)
          sln->set_active_element(active_element);
      }

      void ScalarView::show_linearizer_data(double eps, int item)
      {
        update_mesh_info();

        create();
        update_layout();
        wait_for_draw();
        // FIXME: find out why this has to be called after wait_for_draw in order for the view to be reset initially.
        reset_view(false); // setting true here makes the view always reset after calling 'show'; particularly in the adaptivity process,
        // it would disallow the observation of the process from a manually set viewpoint.
        refresh();
      }

      void ScalarView::update_mesh_info()
      {
        // Get a range of vertex values (or use the range set by the user).
        double vert_min = lin->get_min_value();
        double vert_max = lin->get_max_value();
        // Special case: constant function; offset the lower limit of range so that the domain is drawn under the
        // function and also the scale is drawn correctly.
        if ((vert_max - vert_min) < Hermes::HermesEpsilon)
        {
          is_constant = true;
          vert_min -= 0.5;
        }

        if (range_auto)
        {
          range_min = vert_min;
          range_max = vert_max;
        }

        if (fabs(range_max - range_min) < Hermes::HermesEpsilon)
          value_irange = 1.0;
        else
          value_irange = 1.0 / (range_max - range_min);

        // Calculate the axes-aligned bounding box in the xy-plane.
        lin->calc_vertices_aabb(&vertices_min_x, &vertices_max_x, &vertices_min_y, &vertices_max_y);

        // Calculate average value.
        value_range_avg = 0.0;

        for (Linearizer::Iterator<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::vertex_t> it = lin->vertices_begin(); !it.end; ++it)
        {
          ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::vertex_t& vertex = it.get();
          if (vertex[2] > range_max)
            value_range_avg += range_max;
          else if (vertex[2] < range_min)
            value_range_avg += range_min;
          else
            value_range_avg += vertex[2];
        }

        value_range_avg /= lin->get_vertex_count();
        lin_updated = true;
      }

      void ScalarView::init_element_info(MeshSharedPtr mesh)
      {
        //cleanup
        element_infos.clear();

        //check how many active elements are neccessary and prepare space for them
        element_infos.reserve(mesh->get_num_active_elements());

        //build element info
        Element *element = nullptr;
        for_all_active_elements(element, mesh)
        {
          double sum_x = 0.0, sum_y = 0.0;
          double max_x, max_y, min_x, min_y;
          max_x = min_x = element->vn[0]->x;
          max_y = min_y = element->vn[0]->y;
          for (unsigned int i = 0; i < element->get_nvert(); i++)
          {
            sum_x += element->vn[i]->x;
            sum_y += element->vn[i]->y;

            if (element->vn[i]->x > max_x)
              max_x = element->vn[i]->x;
            if (element->vn[i]->x < min_x)
              min_x = element->vn[i]->x;
            if (element->vn[i]->y > max_y)
              max_y = element->vn[i]->y;
            if (element->vn[i]->y < min_y)
              min_y = element->vn[i]->y;
          }
          element_infos.push_back(ElementInfo(element->id,
            (float)(sum_x / element->get_nvert()), (float)(sum_y / element->get_nvert()),
            (float)(max_x - min_x), (float)(max_y - min_y)));
        }
      }

      Linearizer* ScalarView::get_linearizer()
      {
        return this->lin;
      }

      void ScalarView::draw_element_infos_2d()
      {
        //create widgets
        create_element_info_widgets();

        //prepare environment
        glMatrixMode(GL_MODELVIEW);
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);
        glDisable(GL_BLEND);

        //draw element IDs
        std::vector<ElementInfo>::const_iterator iter = element_infos.begin();
        while (iter != element_infos.end())
        {
          //check element dimension in pixels
          float width = (float)(iter->width * scale);
          float height = (float)(iter->height * scale);

          //draw if AABB of element is large enough
          if (width > 40.f && height > 20.f)
          {
            //prepare environment
            glPushMatrix();
            glTranslatef(iter->x, iter->y, 0.0f);
            glScalef(1 / (float)scale, 1 / (float)scale, 1.0f);

            //prepare number
            void* font = GLUT_BITMAP_HELVETICA_10;
            unsigned char buffer[128];
            sprintf((char*)buffer, "%d", iter->id);
            int width_px = glutBitmapLength(font, buffer);
            int height_px = glutBitmapHeight(font);

            //draw background
            glCallList(element_id_widget);

            //draw number
            glTranslatef(-width_px / 2.0f, -height_px / 3.0f, 0.0f);
            glColor3f(1.0f, 1.0f, 1.0f);
            glRasterPos2f(0.0f, 0.0f);
            glutBitmapString(font, buffer);

            //clear environment
            glPopMatrix();
          }

          //advance to next
          ++iter;
        }
      }

      void ScalarView::create_element_info_widgets()
      {
        if (element_id_widget == 0)
        {
          element_id_widget = glGenLists(1);
          glNewList(element_id_widget, GL_COMPILE);
          {
            glBegin(GL_QUADS);

            //background
            float radius = 20.f;
            glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
            glVertex2f(-radius*1.1, radius*0.5);
            glVertex2f(radius*1.1, radius*0.5);
            glVertex2f(radius*1.1, -radius*0.5);
            glVertex2f(-radius*1.1, -radius*0.5);

            //foreground
            radius = 18.f;
            glColor4f(0.2f, 0.2f, 0.4f, 1.0f);
            glVertex2f(-radius*1.1, radius*0.5);
            glVertex2f(radius*1.1, radius*0.5);
            glVertex2f(radius*1.1, -radius*0.5);
            glVertex2f(-radius*1.1, -radius*0.5);

            glEnd();
          }
          glEndList();
        }
      }

      void ScalarView::show_contours(double step, double orig)
      {
        if (step == 0.0)
          throw Exceptions::ValueException("step", step, 0.0);
        contours = true;
        cont_orig = orig;
        cont_step = step;
        set_palette_filter(true);
        refresh();
      }

      void ScalarView::draw_tri_contours(ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t& triangle)
      {
        // sort the vertices by their value, keep track of the permutation sign.
        int i, idx[3] = { 0, 1, 2 }, perm = 0;
        for (i = 0; i < 2; i++)
        {
          if (triangle[idx[0]][2] > triangle[idx[1]][2])
          {
            std::swap(idx[0], idx[1]);
            perm++;
          }
          if (triangle[idx[1]][2] > triangle[idx[2]][2])
          {
            std::swap(idx[1], idx[2]);
            perm++;
          }
        }
        if (fabs(triangle[idx[0]][2] - triangle[idx[2]][2]) < 1e-3 * fabs(cont_step))
          return;

        // get the first (lowest) contour value
        double val = triangle[idx[0]][2];
        val = std::ceil((val - cont_orig) / cont_step) * cont_step + cont_orig;

        int l1 = 0, l2 = 1;
        int r1 = 0, r2 = 2;

        while (val < triangle[idx[r2]][2])
        {
          double ld = triangle[idx[l2]][2] - triangle[idx[l1]][2];
          double rd = triangle[idx[r2]][2] - triangle[idx[r1]][2];

          // draw a slice of the triangle
          while (val < triangle[idx[l2]][2])
          {
            double lt = (val - triangle[idx[l1]][2]) / ld;
            double rt = (val - triangle[idx[r1]][2]) / rd;

            double x1 = (1.0 - lt) * triangle[idx[l1]][0] + lt * triangle[idx[l2]][0];
            double y1 = (1.0 - lt) * triangle[idx[l1]][1] + lt * triangle[idx[l2]][1];
            double x2 = (1.0 - rt) * triangle[idx[r1]][0] + rt * triangle[idx[r2]][0];
            double y2 = (1.0 - rt) * triangle[idx[r1]][1] + rt * triangle[idx[r2]][1];

            if (perm & 1)
            {
              glVertex2d(x1, y1);
              glVertex2d(x2, y2);
            }
            else
            {
              glVertex2d(x2, y2);
              glVertex2d(x1, y1);
            }

            val += cont_step;
          }
          l1 = 1;
          l2 = 2;
        }
      }

      void ScalarView::prepare_gl_geometry()
      {
        if (lin_updated)
          lin_updated = false;
      }

      void ScalarView::draw_values_2d()
      {
        //set texture for coloring
        glEnable(GL_TEXTURE_1D);
        glBindTexture(GL_TEXTURE_1D, gl_pallete_tex_id);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
        glColor3f(0, 1, 0);

        //set texture transformation matrix
        glMatrixMode(GL_TEXTURE);
        glLoadIdentity();
        glTranslated(tex_shift, 0.0, 0.0);
        glScaled(tex_scale, 0.0, 0.0);

        //render triangles
        glBegin(GL_TRIANGLES);

        for (Linearizer::Iterator<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t> it = lin->triangles_begin(); !it.end; ++it)
        {
          ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t& triangle = it.get();
          glTexCoord1d((triangle[0][2] - range_min) * value_irange);
          glVertex2d(triangle[0][0], triangle[0][1]);
          glTexCoord1d((triangle[1][2] - range_min) * value_irange);
          glVertex2d(triangle[1][0], triangle[1][1]);
          glTexCoord1d((triangle[2][2] - range_min) * value_irange);
          glVertex2d(triangle[2][0], triangle[2][1]);
        }
        glEnd();

        //GL clenaup
        glMatrixMode(GL_TEXTURE); //switch-off texture transform
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
      }

      void ScalarView::draw_edges_2d()
      {
        glColor3fv(edges_color);
        glBegin(GL_LINES);
        for (Linearizer::Iterator<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t> it = lin->edges_begin(); !it.end; ++it)
        {
          ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t& edge = it.get();
          int& edge_marker = it.get_marker();

          if (show_edges || edge_marker)
          {
            glVertex2d(edge[0][0], edge[0][1]);
            glVertex2d(edge[1][0], edge[1][1]);
          }
        }
        glEnd();
      }

#define V0    vertices_min_x - xctr, range_min - yctr, -(vertices_min_y - zctr)
#define V1    vertices_max_x - xctr, range_min - yctr, -(vertices_min_y - zctr)
#define V2    vertices_max_x - xctr, range_min - yctr, -(vertices_max_y - zctr)
#define V3    vertices_min_x - xctr, range_min - yctr, -(vertices_max_y - zctr)
#define V4    vertices_min_x - xctr, range_max - yctr, -(vertices_min_y - zctr)
#define V5    vertices_max_x - xctr, range_max - yctr, -(vertices_min_y - zctr)
#define V6    vertices_max_x - xctr, range_max - yctr, -(vertices_max_y - zctr)
#define V7    vertices_min_x - xctr, range_max - yctr, -(vertices_max_y - zctr)

      void ScalarView::draw_aabb()
      {
        // Axis-aligned bounding box of the model.
        GLdouble aabb[] =
        {
          V0, V1, V2, V3,    // bottom
          V0, V1, V5, V4,    // front
          V0, V3, V7, V4,    // left
          V1, V2, V6, V5,    // right
          V2, V3, V7, V6,    // back
          V4, V5, V6, V7     // top
        };

        // Set the edge color.
        glColor3fv(edges_color);

        // Make the cube wire frame.
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        // Scale the box appropriately.
        glPushMatrix();
        glScaled(xzscale, yscale, xzscale);

        // Activate and specify pointer to vertex array.
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_DOUBLE, 0, aabb);

        // Draw the box (glDrawElements would be better, but do not work for me).
        glDrawArrays(GL_QUADS, 0, 24);

        // Deactivate vertex arrays after drawing.
        glDisableClientState(GL_VERTEX_ARRAY);

        // Recover the original model/view matrix
        glPopMatrix();
      }

      void ScalarView::on_display_2d()
      {
        set_ortho_projection();
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_TEXTURE_2D);
        glDisable(GL_BLEND);

        // setup transformation (follows View::transform_x and View::transform_y)
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        glTranslated(center_x, center_y, 0.0);
        glScaled(1.0, -1.0, 1.0);
        glTranslated(trans_x, trans_y, 0.0);
        glScaled(scale, scale, 1.0);

        // draw all triangles
        if (show_values)
          draw_values_2d();

        // draw contours
        glDisable(GL_TEXTURE_1D);
        if (contours)
        {
          glColor3fv(cont_color);
          //glLineWidth(2.0f);
          glBegin(GL_LINES);
          for (Linearizer::Iterator<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t> it = this->lin->triangles_begin(); !it.end; ++it)
          {
            ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t& triangle = it.get();
            draw_tri_contours(triangle);
          }
          glEnd();
        }

        // draw edges and boundary of mesh
        draw_edges_2d();

        //draw element IDS
        if (show_element_info)
          draw_element_infos_2d();

        //cleanup
        glPopMatrix();
      }

      void ScalarView::on_display_3d()
      {
        set_3d_projection(fovy, znear, zfar);

        glClear(GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);

        // Set camera transforamtion.
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Initialize light and material.
        init_lighting();

        if (do_zoom_to_fit)
        {
          // If the user presses 'c' to center the view automatically, calculate how far to move
          // the visualized model away from the viewer so that it is completely visible in current window.
          ztrans = calculate_ztrans_to_fit_view();
          do_zoom_to_fit = false;
        }

        // Set model transformation.
        glTranslated(xtrans, ytrans, ztrans);
        glRotated(xrot, 1, 0, 0);
        glRotated(yrot, 0, 1, 0);

        // Draw the surface.
        glEnable(GL_LIGHTING);
        glEnable(GL_TEXTURE_1D);
        glBindTexture(GL_TEXTURE_1D, gl_pallete_tex_id);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        glEnable(GL_NORMALIZE);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0, 1.0);
        glBegin(GL_TRIANGLES);
        double normal_xzscale = 1.0 / xzscale, normal_yscale = 1.0 / yscale;

        for (Linearizer::Iterator<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t> it = this->lin->triangles_begin(); !it.end; ++it)
        {
          ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t& triangle = it.get();

          for (int j = 0; j < 3; j++)
          {
            glTexCoord2d((triangle[j][2] - range_min) * value_irange * tex_scale + tex_shift, 0.0);

            glVertex3d((triangle[j][0] - xctr) * xzscale, (triangle[j][2] - yctr) * yscale, -(triangle[j][1] - zctr) * xzscale);
          }
        }

        glEnd();
        glDisable(GL_POLYGON_OFFSET_FILL);

        // Draw edges.
        glDisable(GL_LIGHTING);
        glDisable(GL_TEXTURE_1D);
        if (show_edges)
        {
          glColor3fv(edges_color);
          glBegin(GL_LINES);
          for (Linearizer::Iterator<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t> it = lin->edges_begin(); !it.end; ++it)
          {
            ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t& edge = it.get();

            for (int j = 0; j < 2; j++)
              glVertex3d((edge[j][0] - xctr) * xzscale, (edge[j][2] - yctr) * yscale, -(edge[j][1] - zctr) * xzscale);
          }
          glEnd();
        }

        // Draw the whole bounding box or only the boundary edges.
        if (show_aabb)
          draw_aabb();
        else
        {
          glColor3fv(edges_color);
          glBegin(GL_LINES);
          for (Linearizer::Iterator<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t> it = lin->edges_begin(); !it.end; ++it)
          {
            // Outline of the domain boundary at the bottom of the plot or at the bottom user-defined limit
            double y_coord = (range_min - yctr) * yscale;
            ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::edge_t& edge = it.get();
            int& edge_marker = it.get_marker();

            if (edge_marker)
            {
              glVertex3d((edge[0][0] - xctr) * xzscale, y_coord, -(edge[0][1] - zctr) * xzscale);
              glVertex3d((edge[1][0] - xctr) * xzscale, y_coord, -(edge[1][1] - zctr) * xzscale);
            }
          }
          glEnd();
        }
      }

      void ScalarView::on_display()
      {
        // lock and get data
        lin->lock_data();

        glPolygonMode(GL_FRONT_AND_BACK, pmode ? GL_LINE : GL_FILL);

        //prepare vertices
        prepare_gl_geometry();

        if (mode3d)
          this->on_display_3d();
        else
          this->on_display_2d();

        lin->unlock_data();
      }

      static inline void normalize(double& x, double& y, double& z)
      {
        double l = 1.0 / sqrt(sqr(x) + sqr(y) + sqr(z));
        x *= l; y *= l; z *= l;
      }

      void ScalarView::update_layout()
      {
        View::update_layout();
        // (x, y, -z) coordinates (in the eye coordinate system) of the point that lies at the center of solution domain
        // and vertically at the position of the average value of all vertices. This point will be the origin for all
        // drawing and also the initial position of the camera.
        xctr = (vertices_max_x + vertices_min_x) / 2.0;
        yctr = value_range_avg;
        zctr = (vertices_max_y + vertices_min_y) / 2.0;
      }

      void ScalarView::reset_view(bool force_reset)
      {
        if (force_reset || view_not_reset) { // Reset 3d view.
          xrot = 40.0; yrot = 0.0;
          xtrans = ytrans = ztrans = 0.0;
          // Ensure that the model (before applying any transformations) may be translated along the view axis completely
          // into the viewing volume, and that it is not too flat (so that considerable amount of detail is visible).
          // Later on, we will determine the translation distance so that the model occupies as much of the view space as possible.

          // Set scaling into the (-1, 1) range on the x- and z- axes
          double max_radial = std::max(vertices_max_x - vertices_min_x, vertices_max_y - vertices_min_y);
          xzscale = 2.0 / max_radial;

          // Determine the largest y-axis distance of the function surface from the line of sight
          double max_axial = std::max(range_max - yctr, fabs(range_min - yctr));

          // We use the same scaling for the y-axis as for the x- and z- axes if after this scaling,
          //  1. the model could be translated so that it lies completely below the top clipping plane of the viewing frustum and
          //     its farthest (away from the camera) corner is at least 1 unit before the far clipping plane (to allow some zooming out),
          //  2. the model's bounding box's part above (or below, whichever is bigger) the average function value (yctr) has dimensions
          //     (-1, 1)x(0, height)x(-1, 1), where height > 0.1.
          // If this is not true, the model is either too tall or too flat and will be subject to different scaling
          // along the y-axis, such that it would fit to the viewing frustum at one third of its depth (i.e. at one third of
          // the total available zooming range).
          /// \todo allow the user to decide whether he always wants equal axes scaling or prefers actually seeing something sensible...
          double tan_fovy_half = tan((double)fovy / 2.0 / 180.0 * M_PI);
          double max_allowed_height = (zfar - 3)*tan_fovy_half;
          if (max_axial * xzscale > max_allowed_height || (max_axial * xzscale < 0.1 && !is_constant))
            yscale = (znear + zfar) / 3.0 * tan_fovy_half * value_irange;
          else
            yscale = xzscale;

          do_zoom_to_fit = true;
        }
        View::reset_view(force_reset); // Reset 2d view.
      }

      double ScalarView::calculate_ztrans_to_fit_view()
      {
        // Axis-aligned bounding box of the model (homogeneous coordinates in the model space), divided into the bottom and top base.
        GLdouble aabb[2][16] =
        {
          {
            V0, 1,
            V1, 1,
            V2, 1,
            V3, 1
          },
          {
            V4, 1,
            V5, 1,
            V6, 1,
            V7, 1
          }
        };
        GLdouble aabb_base[16];

        // Tangents of the half of the vertical and horizontal field of view angles.
        double tan_fovy_half = tan((double)fovy / 2.0 / 180.0 * M_PI);
        double tan_fovx_half = tan_fovy_half * (double)output_width / output_height;

        // The viewpoint z-coordinate we are looking for.
        double optimal_viewpoint_pos = 0.0;

        // We will change the model transformation matrix, so save the current one.
        glPushMatrix();

        // Create the model transformation matrix with the default z-translation.
        glLoadIdentity();
        glTranslated(xtrans, ytrans, ztrans);
        glRotated(xrot, 1, 0, 0);
        glRotated(yrot, 0, 1, 0);
        glScaled(xzscale, yscale, xzscale);

        // As far as I know, OpenGL can only perform 4x4 matrix multiplication, so we find the 8 transformed bounding box corners by
        // first multiplying the transformation matrix with a matrix having as its 4 columns the coordinates of the bottom base corners
        // and then with a matrix created from the top base corners.
        for (int i = 0; i < 2; i++)
        {
          glPushMatrix();
          glMultMatrixd(aabb[i]);
          glGetDoublev(GL_MODELVIEW_MATRIX, aabb_base);
          glPopMatrix();

          // Go through the transformed corners and find the biggest distance of its perspective projection center to the origin.
          GLdouble *coord_ptr = &aabb_base[0];
          for (int j = 0; j < 4; j++)
          {
            double perspective_center_to_origin_dist = fabs(coord_ptr[0]) / tan_fovx_half + coord_ptr[2];
            if (perspective_center_to_origin_dist > optimal_viewpoint_pos)
              optimal_viewpoint_pos = perspective_center_to_origin_dist;

            perspective_center_to_origin_dist = fabs(coord_ptr[1]) / tan_fovy_half + coord_ptr[2];
            if (perspective_center_to_origin_dist > optimal_viewpoint_pos)
              optimal_viewpoint_pos = perspective_center_to_origin_dist;

            coord_ptr += 4;
          }
        }

        // Restore the original model transformation matrix.
        glPopMatrix();

        return -optimal_viewpoint_pos;
      }

      void ScalarView::set_vertical_scaling(double sc)
      {
        if (mode3d)
          yscale *= sc;
        else if (contours)
          cont_step *= sc;
        refresh();
      }

      void ScalarView::set_min_max_range(double min, double max)
      {
        /// \todo allow settin min = max, in which case draw the corresponding contour.
        if (fabs(max - min) < Hermes::HermesEpsilon)
        {
          this->warn("Range (%f, %f) is too narrow: adjusted to (%f, %f)", min, max, min - 0.5, max);
          min -= 0.5;
        }
        View::set_min_max_range(min, max);
      }

      void ScalarView::init_lighting()
      {
        float light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
        float light_ambient[] = { 0.1f, 0.1f, 0.1f, 1.0f };
        float light_diffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };

        glEnable(GL_LIGHT0);
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
        glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
        //  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

        float material_ambient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
        float material_diffuse[] = { 0.8f, 0.8f, 0.8f, 1.0f };
        float material_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };

        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128);
        glDisable(GL_COLOR_MATERIAL);

        glShadeModel(GL_SMOOTH);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
        if (GLEW_EXT_separate_specular_color)
        {
          glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL_EXT, GL_SEPARATE_SPECULAR_COLOR_EXT);
        }
      }

      void ScalarView::set_3d_mode(bool enable)
      {
        mode3d = enable; dragging = scaling = false;
        if (mode3d)
        {
          lin->lock_data();
          lin->unlock_data();
        }
        refresh();
      }

      void ScalarView::on_key_down(unsigned char key, int x, int y)
      {
        switch (key)
        {
        case 'm':
          show_edges = !show_edges;
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

        case 'k':
          contours = !contours;
          if (contours)
            this->cont_step = (this->range_max - this->range_min) / 50.;
          refresh();
          break;

        case 'i':
          show_element_info = !show_element_info;
          if (!mode3d)
            refresh();
          break;

        case 'b':
          show_aabb = !show_aabb;
          if (mode3d)
            refresh();
          break;

        case '3':
          mode3d = !mode3d;
          dragging = scaling = false;
          refresh();
          break;

        case '*':
          lin->lock_data();
          if (mode3d)
            yscale *= D3DV_SCALE_STEP_COEF;
          else if (contours)
            cont_step *= D3DV_SCALE_STEP_COEF;
          lin->unlock_data();
          refresh();
          break;

        case '/':
          lin->lock_data();
          if (mode3d)
            yscale /= D3DV_SCALE_STEP_COEF;
          else if (contours)
            cont_step /= D3DV_SCALE_STEP_COEF;
          lin->unlock_data();
          refresh();
          break;

        default:
          View::on_key_down(key, x, y);
          break;
        }
      }


      void ScalarView::on_mouse_move(int x, int y)
      {
        if (mode3d && (dragging || scaling || panning))
        {
          if (dragging)
          {
            yrot += 0.4 * (x - mouse_x);
            xrot += 0.4 * (y - mouse_y);

            if (xrot < -90) xrot = -90;
            else if (xrot > 90) xrot = 90;
          }
          else if (scaling)
          {
            ztrans += 0.01 * (mouse_y - y);
            if (ztrans > -0.25) ztrans = -0.25;
            else if (ztrans < -7) ztrans = -7;
          }
          else
          {
            xtrans += 0.002 * (x - mouse_x);
            ytrans += 0.002 * (mouse_y - y);
          }

          refresh();
          mouse_x = x;
          mouse_y = y;
          return;
        }
        else
        {
          if (!mode3d && show_edges && !dragging && !scaling && !panning)
            refresh();
          else
            View::on_mouse_move(x, y);
        }
      }

      void ScalarView::on_right_mouse_down(int x, int y)
      {
        View::on_right_mouse_down(x, y);
      }

      void ScalarView::on_middle_mouse_down(int x, int y)
      {
        if (!mode3d) return;
        dragging = scaling = false;
        panning = true;
      }

      void ScalarView::on_middle_mouse_up(int x, int y)
      {
        panning = false;
      }


      const char* ScalarView::get_help_text() const
      {
        return
          "ScalarView\n\n"
          "Controls:\n"
          "  Left mouse - pan\n"
          "  Right mouse - zoom\n"
          "  3 - toggle 3D mode\n"
          "  C - center image\n"
          "  F - toggle smooth palette\n"
          "  H - render high-quality frame\n"
          "  K - toggle contours\n"
          "  M - toggle mesh\n"
          "  B - toggle bounding box\n"
          "  I - toggle element IDs (2d only)\n"
          "  P - cycle palettes\n"
          "  S - save screenshot\n"
          "  * - increase contour density\n"
          "  / - decrease contour density\n"
          "  F1 - this help\n"
          "  Esc, Q - quit\n\n"
          "3D mode:\n"
          "  Left mouse - rotate\n"
          "  Right mouse - zoom\n"
          "  Middle mouse - pan\n"
          "  * - increase Z scale\n"
          "  / - decrease Z scale";
      }
    }
  }
}
#endif
