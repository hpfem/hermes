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

#define GL_BUFFER_OFFSET(i) ((char *)NULL + (i))

#ifdef ENABLE_VIEWER_GUI
# include <AntTweakBar.h>
#endif

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
        pmode = mode3d = false;
        normals = NULL;
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
        gl_coord_buffer = 0; gl_index_buffer = 0; gl_edge_inx_buffer = 0;

        do_zoom_to_fit = true;
        is_constant = false;
      }

#ifndef _MSC_VER

      ScalarView::ScalarView(const char* title, WinGeom* wg) :
      View(title, wg), lin(NULL),
        vertex_nodes(0),
        pointed_vertex_node(NULL),
        allow_node_selection(false),
        pointed_node_widget(0),
        selected_node_widget(0),
        node_pixel_radius(10),
        node_widget_vert_cnt(32),
        element_id_widget(0),
        show_element_info(false)
#ifdef ENABLE_VIEWER_GUI
        , tw_wnd_id(TW_WND_ID_NONE), tw_setup_bar(NULL)
#endif
      {
        init();
      }
#endif

      ScalarView::ScalarView(char* title, WinGeom* wg) :
      View(title, wg), lin(NULL),
        vertex_nodes(0),
        pointed_vertex_node(NULL),
        allow_node_selection(false),
        pointed_node_widget(0),
        selected_node_widget(0),
        node_pixel_radius(10),
        node_widget_vert_cnt(32),
        element_id_widget(0),
        show_element_info(false)
#ifdef ENABLE_VIEWER_GUI
        , tw_wnd_id(TW_WND_ID_NONE), tw_setup_bar(NULL)
#endif
      {
        init();
      }

      ScalarView::~ScalarView()
      {
        delete [] normals;
        vertex_nodes.clear();
        if(lin != NULL)
          delete lin;

# ifdef ENABLE_VIEWER_GUI
        if(tw_wnd_id != TW_WND_ID_NONE)
        {
          TwSetCurrentWndID(tw_wnd_id);
          TwTerminate();
          tw_wnd_id = TW_WND_ID_NONE;
        }
# endif
      }

      void ScalarView::on_create(int output_id)
      {
        View::on_create(output_id);

# ifdef ENABLE_VIEWER_GUI
        if(tw_wnd_id == TW_WND_ID_NONE) //GUI has to be initialized before any displaying occurs (on_create calls post_redisplay by default)
        {
          TwSetCurrentWndID(TW_WND_ID_NONE);
          if(TwInit(TW_OPENGL, NULL))
          {
            tw_wnd_id = TwGetCurrentWndID();
            create_setup_bar();
          }
          else
          {
            error("TW init failed: %s", TwGetLastError());
          }
        }
# endif
      }

      void ScalarView::on_close()
      {
        //GUI cleanup
# ifdef ENABLE_VIEWER_GUI
        if(tw_wnd_id != TW_WND_ID_NONE)
        {
          TwSetCurrentWndID(tw_wnd_id);
          TwTerminate();
          tw_wnd_id = TW_WND_ID_NONE;
        }
# endif

        //OpenGL clenaup
        if(pointed_node_widget != 0)
        {
          glDeleteLists(pointed_node_widget, 1);
          pointed_node_widget = 0;
        }
        if(selected_node_widget != 0)
        {
          glDeleteLists(selected_node_widget, 1);
          selected_node_widget = 0;
        }
        if(element_id_widget != 0)
        {
          glDeleteLists(element_id_widget, 1);
          element_id_widget = 0;
        }
        if(gl_coord_buffer != 0)
        {
          glDeleteBuffersARB(1, &gl_coord_buffer);
          gl_coord_buffer = 0;
        }
        if(gl_index_buffer != 0)
        {
          glDeleteBuffersARB(1, &gl_index_buffer);
          gl_index_buffer = 0;
        }

        //call of parent implementation
        View::on_close();
      }

      void ScalarView::create_setup_bar()
      {
#ifdef ENABLE_VIEWER_GUI
        char buffer[1024];

        TwBar* tw_bar = TwNewBar("View setup");

        // Contours.
        TwAddVarRW(tw_bar, "contours", TW_TYPE_BOOLCPP, &contours, " group = Contour2D label = 'Show contours'");
        sprintf(buffer, " group = Contour2D label = 'Begin' step = %g", CONT_CHANGE);
        TwAddVarRW(tw_bar, "cont_orig", TW_TYPE_DOUBLE, &cont_orig, buffer);
        sprintf(buffer, " group = Contour2D label = 'Step' min = %g step = %g", MIN_CONT_STEP, CONT_CHANGE);
        TwAddVarRW(tw_bar, "cont_step", TW_TYPE_DOUBLE, &cont_step, buffer);
        TwAddVarRW(tw_bar, "cont_color", TW_TYPE_COLOR3F, &cont_color, " group = Contour2D label = 'Color'");

        // Mesh.
        TwAddVarRW(tw_bar, "show_values", TW_TYPE_BOOLCPP, &show_values, " group = Elements2D label = 'Show Value'");
        TwAddVarRW(tw_bar, "show_edges", TW_TYPE_BOOLCPP, &show_edges, " group = Elements2D label = 'Show Edges'");
        TwAddVarRW(tw_bar, "show_aabb", TW_TYPE_BOOLCPP, &show_aabb, " group = Elements2D label = 'Show Bounding box'");
        TwAddVarRW(tw_bar, "show_element_info", TW_TYPE_BOOLCPP, &show_element_info, " group = Elements2D label = 'Show ID'");
        TwAddVarRW(tw_bar, "allow_node_selection", TW_TYPE_BOOLCPP, &allow_node_selection, " group = Elements2D label = 'Allow Node Sel.'");
        TwAddVarRW(tw_bar, "edges_color", TW_TYPE_COLOR3F, &edges_color, " group = Elements2D label = 'Edge color'");

        // Help.
        const char* help_text = get_help_text();
        TwSetParam(tw_bar, NULL, "help", TW_PARAM_CSTRING, 1, help_text);

        tw_setup_bar = tw_bar;

#endif
      }

      void ScalarView::show(MeshFunction<double>* sln, double eps, int item,
        MeshFunction<double>* xdisp, MeshFunction<double>* ydisp, double dmult)
      {
        // For preservation of the sln's active element. Will be set back after the visualization.
        Element* active_element = sln->get_active_element();

        if(lin == NULL)
          lin = new Linearizer;

        if(!range_auto)
          lin->set_max_absolute_value(std::max(fabs(range_min), fabs(range_max)));

        lin->set_displacement(xdisp, ydisp, dmult);
        lin->lock_data();

        lin->process_solution(sln, item, eps);

        update_mesh_info();

        // Initialize mesh nodes for displaying and selection.
        init_vertex_nodes(sln->get_mesh());

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
        if(dynamic_cast<Solution<double>*>(sln) != NULL)
          if(active_element != NULL)
            // Also when there has not been a call to set_active_element since assignment to this MeshFunction,
            // there is nothing to restore to.
            if(active_element->active)
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
        // Calculate normals if necessary.
        if(mode3d)
          calculate_normals(lin->get_vertices(), lin->get_num_vertices(), lin->get_triangles(), lin->get_num_triangles());
        else
        {
          delete [] normals;
          normals = NULL;
        }

        // Get a range of vertex values (or use the range set by the user).
        double vert_min = lin->get_min_value();
        double vert_max = lin->get_max_value();
        // Special case: constant function; offset the lower limit of range so that the domain is drawn under the
        // function and also the scale is drawn correctly.
        if((vert_max - vert_min) < 1e-8)
        {
          is_constant = true;
          vert_min -= 0.5;
        }

        if(range_auto)
        {
          range_min = vert_min;
          range_max = vert_max;
        }

        if(fabs(range_max - range_min) < 1e-8)
          value_irange = 1.0;
        else
          value_irange = 1.0 / (range_max - range_min);

        // Calculate the axes-aligned bounding box in the xy-plane.
        lin->calc_vertices_aabb(&vertices_min_x, &vertices_max_x, &vertices_min_y, &vertices_max_y);

        // Calculate average value.
        value_range_avg = 0.0;
        double3* verts = lin->get_vertices();
        const int num_verts = lin->get_num_vertices();
        for(int i = 0; i < num_verts; i++)
          if(verts[i][2] > range_max)
            value_range_avg += range_max;
          else if(verts[i][2] < range_min)
            value_range_avg += range_min;
          else
            value_range_avg += verts[i][2];
          value_range_avg /= num_verts;

          lin_updated = true;
      }

      bool ScalarView::compare_vertex_nodes_x(const VertexNodeInfo& a, const VertexNodeInfo& b)
      {
        return a.x < b.x;
      }

      void ScalarView::init_vertex_nodes(const Mesh* mesh)
      {
        //clear all selections
        pointed_vertex_node = NULL;
#ifdef ENABLE_VIEWER_GUI
        if(tw_wnd_id != TW_WND_ID_NONE)
        {
          TwSetCurrentWndID(tw_wnd_id);
          vector<VertexNodeInfo>::iterator iter = vertex_nodes.begin();
          while (iter != vertex_nodes.end())
          {
            if(iter->tw_bar != NULL)
              TwDeleteBar((TwBar*)iter->tw_bar);
            ++iter;
          }
        }
#endif
        vertex_nodes.clear();

        //count a number of active nodes
        int active_nodes = 0;
        int max_node_id = mesh->get_max_node_id();
        for(int i = 0; i < max_node_id; i++)
        {
          Node* node = mesh->get_node(i);
          if(node->type == 0 && node->used) //inspect only vertex nodes
            active_nodes++;
        }

        //allocate
        vertex_nodes.clear();
        vertex_nodes.reserve(active_nodes);

        //copy
        for(int i = 0; i < max_node_id; i++)
        {
          Node* node = mesh->get_node(i);
          if(node->type == 0 && node->used)
            vertex_nodes.push_back(VertexNodeInfo(node->id, (float)node->x, (float)node->y));
        }

        //sort
        std::sort(vertex_nodes.begin(), vertex_nodes.end(), compare_vertex_nodes_x);
      }

      ScalarView::VertexNodeInfo* ScalarView::find_nearest_node_in_range(float x, float y, float radius)
      {
        VertexNodeInfo node_info(-1, x - radius, y); //right side of the widget
        Hermes::vector<VertexNodeInfo>::iterator found_iter = std::lower_bound(vertex_nodes.begin(), vertex_nodes.end(), node_info, compare_vertex_nodes_x);
        Hermes::vector<VertexNodeInfo>::iterator found_nearest = vertex_nodes.end();
        float found_nearest_dist = -1;
        while (found_iter != vertex_nodes.end() && abs(found_iter->x - x) <= radius)
        {
          if(abs(found_iter->y - y) <= radius)
          {
            float dist = std::min(found_iter->x - x, found_iter->y - y);
            if(found_nearest_dist < 0 || dist < found_nearest_dist)
            {
              found_nearest_dist = dist;
              found_nearest = found_iter;
            }
          }
          ++found_iter;
        }

        if(found_nearest != vertex_nodes.end())
          return found_nearest.operator->();
        else
          return NULL;
      }

      void ScalarView::draw_single_vertex_node(const VertexNodeInfo& node)
      {
        //prepare environment
        glPushMatrix();
        glTranslatef(node.x, node.y, 0.0f);
        glScalef(1/(float)scale, 1/(float)scale, 1.0f);

        //prepare number
        void* font = GLUT_BITMAP_HELVETICA_10;
        unsigned char buffer[128];
        sprintf((char*)buffer, "%d", node.id);
        int width_px = glutBitmapLength(font, buffer);
        int height_px = glutBitmapHeight(font);

        //draw target
        if(width_px > (2*node_pixel_radius))
        {
          glPushMatrix();
          float coef = width_px / (2.0f * node_pixel_radius);
          glScalef(coef, coef, 1.0f);
          glCallList(selected_node_widget);
          glPopMatrix();
        }
        else
          glCallList(selected_node_widget);

        //draw number
        glTranslatef(-width_px/2.0f, -height_px/3.0f, 0.0f);
        glColor3f(1.0f, 1.0f, 1.0f);
        glRasterPos2f(0.0f, 0.0f);
        glutBitmapString(font, buffer);

        //clear environment
        glPopMatrix();
      }

      void ScalarView::draw_vertex_nodes()
      {
        //create widgets
        create_nodes_widgets();

        //prepare environment
        glMatrixMode(GL_MODELVIEW);
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);

        //draw selected nodes
        glDisable(GL_BLEND);
        Hermes::vector<VertexNodeInfo>::const_iterator iter = vertex_nodes.begin();
        while (iter != vertex_nodes.end())
        {
          if(iter->selected)
            draw_single_vertex_node(*iter);
          ++iter;
        }

        //draw a node under cursor
        if(pointed_vertex_node != NULL)
        {
          glEnable(GL_BLEND);
          glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
          glPushMatrix();
          glTranslatef(pointed_vertex_node->x, pointed_vertex_node->y, 0.0f);
          glScalef(1/(float)scale, 1/(float)scale, 1.0f);
          glCallList(pointed_node_widget);
          glPopMatrix();
          glDisable(GL_BLEND);
        }
      }

      void ScalarView::create_nodes_widgets()
      {
        //pointed node
        if(pointed_node_widget == 0)
        {
          pointed_node_widget  = glGenLists(1);
          glNewList(pointed_node_widget, GL_COMPILE);
          {
            //background
            glColor4f(1.0f, 1.0f, 1.0f, 0.5f);
            glBegin(GL_TRIANGLE_FAN);
            glVertex2f(0.0f, 0.0f);
            float radius = 1.3f * node_pixel_radius;
            for(int i = 0; i < node_widget_vert_cnt; i++)
            {
              float angle = (float)((i / (float)(node_widget_vert_cnt-1)) * (2.0f * M_PI));
              glVertex2f(radius * cos(angle), radius * sin(angle));
            }
            glEnd();

            //foreground
            glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
            glBegin(GL_LINE_STRIP);
            radius = (float)node_pixel_radius;
            for(int i = 0; i < node_widget_vert_cnt; i++)
            {
              float angle = (float)((i / (float)(node_widget_vert_cnt-1)) * (2.0f * M_PI));
              glVertex2f(radius * cos(angle), radius * sin(angle));
            }
            glEnd();
          }
          glEndList();
        }

        //selected mesh node
        if(selected_node_widget == 0)
        {
          selected_node_widget = glGenLists(1);
          glNewList(selected_node_widget, GL_COMPILE);
          {
            //background
            glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
            glBegin(GL_TRIANGLE_FAN);
            glVertex2f(0.0f, 0.0f);
            float radius = 1.3f * node_pixel_radius;
            for(int i = 0; i < node_widget_vert_cnt; i++)
            {
              float angle = (float)((i / (float)(node_widget_vert_cnt-1)) * (2.0f * M_PI));
              glVertex2f(radius * cos(angle), radius * sin(angle));
            }
            glEnd();

            //foreground
            glColor4f(0.4f, 0.4f, 0.2f, 1.0f);
            glBegin(GL_TRIANGLE_FAN);
            glVertex2f(0.0f, 0.0f);
            radius = (float)node_pixel_radius;
            for(int i = 0; i < node_widget_vert_cnt; i++)
            {
              float angle = (float)((i / (float)(node_widget_vert_cnt-1)) * (2.0f * M_PI));
              glVertex2f(radius * cos(angle), radius * sin(angle));
            }
            glEnd();
          }
          glEndList();
        }
      }

      void ScalarView::init_element_info(const Mesh* mesh)
      {
        //cleanup
        element_infos.clear();

        //check how many active elements are neccessary and prepare space for them
        element_infos.reserve(mesh->get_num_active_elements());

        //build element info
        Element *element = NULL;
        for_all_active_elements(element, mesh)
        {
          double sum_x = 0.0, sum_y = 0.0;
          double max_x, max_y, min_x, min_y;
          max_x = min_x = element->vn[0]->x;
          max_y = min_y = element->vn[0]->y;
          for(unsigned int i = 0; i < element->get_num_surf(); i++)
          {
            sum_x += element->vn[i]->x;
            sum_y += element->vn[i]->y;

            if(element->vn[i]->x > max_x)
              max_x = element->vn[i]->x;
            if(element->vn[i]->x < min_x)
              min_x = element->vn[i]->x;
            if(element->vn[i]->y > max_y)
              max_y = element->vn[i]->y;
            if(element->vn[i]->y < min_y)
              min_y = element->vn[i]->y;
          }
          element_infos.push_back(ElementInfo(element->id,
            (float)(sum_x / element->get_num_surf()), (float)(sum_y / element->get_num_surf()),
            (float)(max_x - min_x), (float)(max_y - min_y)));
        }
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
        Hermes::vector<ElementInfo>::const_iterator iter = element_infos.begin();
        while (iter != element_infos.end())
        {
          //check element dimension in pixels
          float width = (float)(iter->width * scale);
          float height = (float)(iter->height * scale);

          //draw if AABB of element is large enough
          if(width > 3*node_pixel_radius && height > 3*node_pixel_radius)
          {
            //prepare environment
            glPushMatrix();
            glTranslatef(iter->x, iter->y, 0.0f);
            glScalef(1/(float)scale, 1/(float)scale, 1.0f);

            //prepare number
            void* font = GLUT_BITMAP_HELVETICA_10;
            unsigned char buffer[128];
            sprintf((char*)buffer, "%d", iter->id);
            int width_px = glutBitmapLength(font, buffer);
            int height_px = glutBitmapHeight(font);

            //draw background
            glCallList(element_id_widget);

            //draw number
            glTranslatef(-width_px/2.0f, -height_px/3.0f, 0.0f);
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
        if(element_id_widget == 0)
        {
          element_id_widget = glGenLists(1);
          glNewList(element_id_widget, GL_COMPILE);
          {
            glBegin(GL_QUADS);

            //background
            float radius = 1.3f * node_pixel_radius;
            glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
            glVertex2f(-radius, radius);
            glVertex2f( radius, radius);
            glVertex2f( radius, -radius);
            glVertex2f(-radius, -radius);

            //foreground
            radius = (float)node_pixel_radius;
            glColor4f(0.2f, 0.2f, 0.4f, 1.0f);
            glVertex2f(-radius, radius);
            glVertex2f( radius, radius);
            glVertex2f( radius, -radius);
            glVertex2f(-radius, -radius);

            glEnd();
          }
          glEndList();
        }
      }

      void ScalarView::show_contours(double step, double orig)
      {
        if(step == 0.0)
          throw Exceptions::ValueException("step", step, 0.0);
        contours = true;
        cont_orig = orig;
        cont_step = step;
        set_palette_filter(true);
        refresh();
      }

      static double my_ceil(double x)
      {
        double y = ceil(x);
        if(y > x) return y;
        return y + 1.0;
      }

      void ScalarView::draw_tri_contours(double3* vert, int3* tri)
      {
        // sort the vertices by their value, keep track of the permutation sign
        int i, idx[3], perm = 0;
        memcpy(idx, tri, sizeof(idx));
        for (i = 0; i < 2; i++)
        {
          if(vert[idx[0]][2] > vert[idx[1]][2]) { std::swap(idx[0], idx[1]); perm++; }
          if(vert[idx[1]][2] > vert[idx[2]][2]) { std::swap(idx[1], idx[2]); perm++; }
        }
        if(fabs(vert[idx[0]][2] - vert[idx[2]][2]) < 1e-3 * fabs(cont_step)) return;

        // get the first (lowest) contour value
        double val = vert[idx[0]][2];
        val = my_ceil((val - cont_orig) / cont_step) * cont_step + cont_orig;

        int l1 = 0, l2 = 1;
        int r1 = 0, r2 = 2;

        // Adjustment of the contour step.
        /*
        while(std::abs(val - vert[idx[r2]][2]) > 50 * cont_step)
          cont_step *= 2;
        while(std::abs(val - vert[idx[r2]][2]) < 2E-2 * cont_step)
          cont_step /= 2;
        */

        while (val < vert[idx[r2]][2])
        {
          double ld = vert[idx[l2]][2] - vert[idx[l1]][2];
          double rd = vert[idx[r2]][2] - vert[idx[r1]][2];

          // draw a slice of the triangle
          while (val < vert[idx[l2]][2])
          {
            double lt = (val - vert[idx[l1]][2]) / ld;
            double rt = (val - vert[idx[r1]][2]) / rd;

            double x1 = (1.0 - lt) * vert[idx[l1]][0] + lt * vert[idx[l2]][0];
            double y1 = (1.0 - lt) * vert[idx[l1]][1] + lt * vert[idx[l2]][1];
            double x2 = (1.0 - rt) * vert[idx[r1]][0] + rt * vert[idx[r2]][0];
            double y2 = (1.0 - rt) * vert[idx[r1]][1] + rt * vert[idx[r2]][1];

            if(perm & 1) { glVertex2d(x1, y1); glVertex2d(x2, y2); }
            else { glVertex2d(x2, y2); glVertex2d(x1, y1); }

            val += cont_step;
          }
          l1 = 1;
          l2 = 2;
        }
      }

      void ScalarView::prepare_gl_geometry()
      {
        if(lin_updated)
        {
          lin_updated = false;

          try
          {
            //get input data
            int vert_cnt = lin->get_num_vertices();
            double3* verts = lin->get_vertices();
            int tri_cnt = lin->get_num_triangles();
            int3* tris = lin->get_triangles();

            //check if extension is supported
            if(!GLEW_ARB_vertex_buffer_object)
              throw std::runtime_error("ARB_vertex_buffer_object not supported");

            //reallocate indices
            if(gl_index_buffer == 0 || tri_cnt > max_gl_tris)
            {
              if(gl_index_buffer == 0)
                glGenBuffersARB(1, &gl_index_buffer);
              glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, gl_index_buffer);
              glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, sizeof(GLint) * tri_cnt * 3, NULL, GL_DYNAMIC_DRAW_ARB);
              GLenum err = glGetError();
              if(err != GL_NO_ERROR)
                throw std::runtime_error("unable to allocate vertex buffer: " + err);
              max_gl_tris = tri_cnt;
            }
            else
            {
              glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, gl_index_buffer);
            }

            //fill indices
            GLuint* gl_triangle = (GLuint*)glMapBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, GL_WRITE_ONLY_ARB);
            if(gl_triangle == NULL)
              throw std::runtime_error("unable to map index buffer: " + glGetError());
            gl_tri_cnt = 0;
            for(int i = 0; i < tri_cnt; i++)
            {
              const int3& triangle = tris[i];
              const double3& vert_a = verts[triangle[0]];
              const double3& vert_b = verts[triangle[1]];
              const double3& vert_c = verts[triangle[2]];
              if(finite(vert_a[2]) && finite(vert_b[2]) && finite(vert_c[2]))
              {
                gl_triangle[0] = (GLint)triangle[0];
                gl_triangle[1] = (GLint)triangle[1];
                gl_triangle[2] = (GLint)triangle[2];
                gl_tri_cnt++;
                gl_triangle += 3; //three indices per triangle
              }
            }
            glUnmapBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB);

            //reallocate vertices
            if(gl_coord_buffer == 0 || vert_cnt > max_gl_verts)
            {
              if(gl_coord_buffer == 0)
                glGenBuffersARB(1, &gl_coord_buffer);
              glBindBufferARB(GL_ARRAY_BUFFER_ARB, gl_coord_buffer);
              glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLVertex2) * vert_cnt, NULL, GL_DYNAMIC_DRAW_ARB);
              GLenum err = glGetError();
              if(err != GL_NO_ERROR)
                throw std::runtime_error("unable to allocate coord buffer: " + err);
              max_gl_verts = vert_cnt;
            }
            else
            {
              glBindBufferARB(GL_ARRAY_BUFFER_ARB, gl_coord_buffer);
            }

            //fill vertices
            GLVertex2* gl_verts = (GLVertex2*)glMapBufferARB(GL_ARRAY_BUFFER_ARB, GL_WRITE_ONLY_ARB);
            if(gl_verts == NULL)
              throw std::runtime_error("unable to map coord buffer: " + glGetError());
            for(int i = 0; i < vert_cnt; i++)
              gl_verts[i] = GLVertex2((float)verts[i][0], (float)verts[i][1], (float)((verts[i][2] - range_min) * value_irange));
            glUnmapBufferARB(GL_ARRAY_BUFFER_ARB);

            //allocate edge indices
            if(gl_edge_inx_buffer == 0)
            {
              glGenBuffersARB(1, &gl_edge_inx_buffer);
              glBindBufferARB(GL_ARRAY_BUFFER_ARB, gl_edge_inx_buffer);
              glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLuint) * H2DV_GL_MAX_EDGE_BUFFER * 2, NULL, GL_DYNAMIC_DRAW_ARB);
              GLenum err = glGetError();
              if(err != GL_NO_ERROR) { //if it fails, no problem
                glDeleteBuffersARB(1, &gl_edge_inx_buffer);
                gl_edge_inx_buffer = 0;
              }
            }
            glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
          }
          catch(std::exception &e)
          { //out-of-memory or any other failure
            if(gl_coord_buffer) { glDeleteBuffersARB(1, &gl_coord_buffer); gl_coord_buffer = 0; }
            if(gl_index_buffer) { glDeleteBuffersARB(1, &gl_index_buffer); gl_index_buffer = 0; }
            if(gl_edge_inx_buffer) { glDeleteBuffersARB(1, &gl_edge_inx_buffer); gl_edge_inx_buffer = 0; }
          }
        }
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
        if(gl_coord_buffer == 0 || gl_index_buffer == 0)
        { //render using the safe but slow methoold and slow method
          //obtain data
          int tri_cnt = lin->get_num_triangles();
          const double3* vert = lin->get_vertices();
          const int3* tris = lin->get_triangles();

          //render
          glBegin(GL_TRIANGLES);
          for (int i = 0; i < tri_cnt; i++)
          {
            const int3& triangle = tris[i];
            const double3& vert_a = vert[triangle[0]];
            const double3& vert_b = vert[triangle[1]];
            const double3& vert_c = vert[triangle[2]];

            if(finite(vert_a[2]) && finite(vert_b[2]) && finite(vert_c[2]))
            {
              glTexCoord1d((vert_a[2] - range_min) * value_irange);
              glVertex2d(vert_a[0], vert_a[1]);
              glTexCoord1d((vert_b[2] - range_min) * value_irange);
              glVertex2d(vert_b[0], vert_b[1]);
              glTexCoord1d((vert_c[2] - range_min) * value_irange);
              glVertex2d(vert_c[0], vert_c[1]);
            }
          }
          glEnd();
        }
        else { //render using vertex buffer object
          if(gl_tri_cnt > 0)
          {
            //bind vertices
            glBindBufferARB(GL_ARRAY_BUFFER_ARB, gl_coord_buffer);
            glVertexPointer(2, GL_FLOAT, sizeof(GLVertex2), GL_BUFFER_OFFSET(0));
            glTexCoordPointer(1, GL_FLOAT, sizeof(GLVertex2), GL_BUFFER_OFFSET(GLVertex2::H2D_OFFSETOF_COORD));
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_TEXTURE_COORD_ARRAY);

            //bind indices
            glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, gl_index_buffer);

            //render
            glDrawElements(GL_TRIANGLES, 3*gl_tri_cnt, GL_UNSIGNED_INT, GL_BUFFER_OFFSET(0));

            //GL cleanup
            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_TEXTURE_COORD_ARRAY);
          }
        }

        //GL clenaup
        glMatrixMode(GL_TEXTURE); //switch-off texture transform
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);
      }

      void ScalarView::draw_edges_2d()
      {
        glColor3fv(edges_color);
        bool displayed = false;
        if(gl_edge_inx_buffer != 0) {//VBO
          //prepage GL
          glEnableClientState(GL_VERTEX_ARRAY);

          //map edge buffer
          glBindBufferARB(GL_ARRAY_BUFFER_ARB, gl_edge_inx_buffer);
          GLuint *gl_inx_buffer = (GLuint*)glMapBufferARB(GL_ARRAY_BUFFER_ARB, GL_WRITE_ONLY);
          GLenum err = glGetError();
          if(err == GL_NO_ERROR) {//render edges
            unsigned int buffer_inx = 0;
            int3* edges = lin->get_edges();
            int edge_cnt = lin->get_num_edges();
            for (int i = 0; i < edge_cnt; i++) //TODO: we could draw only left-right, top-bottom ones
            {
              const int3 &edge = edges[i];
              if(show_edges || edge[2] != 0) { //select internal edges to draw
                gl_inx_buffer[buffer_inx] = (GLuint)edge[0];
                gl_inx_buffer[buffer_inx + 1] = (GLuint)edge[1];
                buffer_inx += 2;
              }

              //render buffer if it is full or if this is the last edge processed
              if(buffer_inx == (2*H2DV_GL_MAX_EDGE_BUFFER) || (buffer_inx > 0 && i == (edge_cnt-1)))
              {
                glUnmapBufferARB(GL_ARRAY_BUFFER_ARB);

                //bind vertices and buffers
                glBindBufferARB(GL_ARRAY_BUFFER_ARB, gl_coord_buffer);
                glVertexPointer(2, GL_FLOAT, sizeof(GLVertex2), GL_BUFFER_OFFSET(0));
                glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, gl_edge_inx_buffer);

                //render
                glDrawElements(GL_LINES, buffer_inx, GL_UNSIGNED_INT, GL_BUFFER_OFFSET(0));

                //map edge buffer
                glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
                glBindBufferARB(GL_ARRAY_BUFFER_ARB, gl_edge_inx_buffer);
                gl_inx_buffer = (GLuint*)glMapBufferARB(GL_ARRAY_BUFFER_ARB, GL_WRITE_ONLY);
                buffer_inx = 0;
              }
            }
            glUnmapBufferARB(GL_ARRAY_BUFFER_ARB);

            //GL cleanup
            glDisableClientState(GL_VERTEX_ARRAY);

            displayed = true;
          }
        }

        //safe fallback
        if(!displayed) { //VBO not suppored
          glBegin(GL_LINES);
          draw_edges(&draw_gl_edge, NULL, !show_edges);
          glEnd();
        }
      }

      void ScalarView::draw_normals_3d()
      {
        double normal_xzscale = 1.0 / xzscale, normal_yscale = 1.0 / yscale;

        glPushAttrib(GL_ENABLE_BIT);
        glPushMatrix();

        const int num_vert = lin->get_num_vertices();
        double3* vert = lin->get_vertices();

        glDisable(GL_LIGHTING);
        glDisable(GL_TEXTURE_1D);
        glScaled(1, 1, -1);
        glColor3f(0.8f, 0.5f, 0.5f);
        glBegin(GL_LINES);
        for(int i = 0; i < num_vert; i++)
        {
          double x = (vert[i][0] - xctr) * xzscale;
          double y = (vert[i][2] - yctr) * yscale;
          double z = (vert[i][1] - zctr) * xzscale;
          glVertex3d(x, y, z);
          glVertex3d(x + 0.1*normals[i][0]*normal_xzscale, y + 0.1*normals[i][2]*normal_yscale, z + 0.1*normals[i][1]*normal_xzscale);
        }
        glEnd();

        glPopMatrix();
        glPopAttrib();
      }

      void ScalarView::draw_gl_edge(int inx_vert_a, int inx_vert_b, ScalarView* viewer, void* param)
      {
        double3* verts = viewer->lin->get_vertices();
        glVertex2d(verts[inx_vert_a][0], verts[inx_vert_a][1]);
        glVertex2d(verts[inx_vert_b][0], verts[inx_vert_b][1]);
      }

      void ScalarView::draw_edges(DrawSingleEdgeCallback draw_single_edge, void* param, bool boundary_only)
      {
        int3* edges = lin->get_edges();
        int edge_cnt = lin->get_num_edges();
        for (int i = 0; i < edge_cnt; i++) //TODO: we could draw only left-right, top-bottom ones
        {
          const int3 &edge = edges[i];
          if(!boundary_only || edge[2] != 0) //select internal edges to draw
            draw_single_edge(edge[0], edge[1], this, param);
        }
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

      void ScalarView::on_display()
      {
        int i, j;

        // lock and get data
        lin->lock_data();
        double3* vert = lin->get_vertices();
        int3* tris = lin->get_triangles();
        int3* edges = lin->get_edges();

        glPolygonMode(GL_FRONT_AND_BACK, pmode ? GL_LINE : GL_FILL);

        //prepare vertices
        prepare_gl_geometry();

        if(!mode3d)
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
          if(show_values)
            draw_values_2d();

          // draw contours
          glDisable(GL_TEXTURE_1D);
          if(contours)
          {
            glColor3fv(cont_color);
            glBegin(GL_LINES);
            for (i = 0; i < lin->get_num_triangles(); i++)
            {
              if(finite(vert[tris[i][0]][2]) && finite(vert[tris[i][1]][2]) && finite(vert[tris[i][2]][2]))
              {
                draw_tri_contours(vert, &tris[i]);
              }
            }
            glEnd();
          }

          // draw edges and boundary of mesh
          draw_edges_2d();

          //draw element IDS
          if(show_element_info)
            draw_element_infos_2d();

          //draw nodes and node info
          if(show_edges && allow_node_selection)
            draw_vertex_nodes();

          //cleanup
          glPopMatrix();
        }
        else
        {
          set_3d_projection(fovy, znear, zfar);

          glClear(GL_DEPTH_BUFFER_BIT);
          glEnable(GL_DEPTH_TEST);

          // Set camera transforamtion.
          glMatrixMode(GL_MODELVIEW);
          glLoadIdentity();

          // Initialize light and material.
          init_lighting();

          if(do_zoom_to_fit)
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
          for (i = 0; i < lin->get_num_triangles(); i++)
          {
            for (j = 0; j < 3; j++)
            {
              glNormal3d(normals[tris[i][j]][0] * normal_xzscale, normals[tris[i][j]][2] * normal_yscale, -normals[tris[i][j]][1] * normal_xzscale);
              glTexCoord2d((vert[tris[i][j]][2] - range_min) * value_irange * tex_scale + tex_shift, 0.0);
              glVertex3d((vert[tris[i][j]][0] - xctr) * xzscale,
                (vert[tris[i][j]][2] - yctr) * yscale,
                -(vert[tris[i][j]][1] - zctr) * xzscale);
            }
          }
          glEnd();
          glDisable(GL_POLYGON_OFFSET_FILL);

          // Draw edges.
          glDisable(GL_LIGHTING);
          glDisable(GL_TEXTURE_1D);
          if(show_edges)
          {
            glColor3fv(edges_color);
            glBegin(GL_LINES);
            for (i = 0; i < lin->get_num_edges(); i++)
            {
              glVertex3d((vert[edges[i][0]][0] - xctr) * xzscale,
                (vert[edges[i][0]][2] - yctr) * yscale,
                -(vert[edges[i][0]][1] - zctr) * xzscale);
              glVertex3d((vert[edges[i][1]][0] - xctr) * xzscale,
                (vert[edges[i][1]][2] - yctr) * yscale,
                -(vert[edges[i][1]][1] - zctr) * xzscale);
            }
            glEnd();
          }

          // Draw the whole bounding box or only the boundary edges.
          if(show_aabb)
            draw_aabb();
          else
          {
            glColor3fv(edges_color);
            glBegin(GL_LINES);
            for (i = 0; i < lin->get_num_edges(); i++)
            {
              // Outline of the domain boundary at the bottom of the plot or at the bottom user-defined limit
              double y_coord = (range_min - yctr) * yscale;
              int3& edge = edges[i];
              if(edge[2])
              {
                glVertex3d((vert[edge[0]][0] - xctr) * xzscale, y_coord,
                  -(vert[edge[0]][1] - zctr) * xzscale);
                glVertex3d((vert[edge[1]][0] - xctr) * xzscale, y_coord,
                  -(vert[edge[1]][1] - zctr) * xzscale);
              }
            }
            glEnd();
          }
        }

        lin->unlock_data();

        // Draw TW.
#ifdef ENABLE_VIEWER_GUI
        if(tw_wnd_id != TW_WND_ID_NONE)
        {
          TwSetCurrentWndID(tw_wnd_id);
          TwDraw();
        }
#endif
      }

      static inline void normalize(double& x, double& y, double& z)
      {
        double l = 1.0 / sqrt(sqr(x) + sqr(y) + sqr(z));
        x *= l; y *= l; z *= l;
      }

      void ScalarView::calculate_normals(double3* vert, int num_verts, int3* tris, int num_tris)
      {
        if(normals != NULL)
          delete [] normals;
        normals = new double3[num_verts];
        memset(normals, 0, sizeof(double3) * num_verts);
        for (int i = 0; i < num_tris; i++)
        {
          int3 &tri = tris[i];
          double ax = (vert[tri[1]][0] - vert[tri[0]][0]);
          double ay = (vert[tri[1]][1] - vert[tri[0]][1]);
          double az = (vert[tri[1]][2] - vert[tri[0]][2]);

          double bx = (vert[tri[2]][0] - vert[tri[0]][0]);
          double by = (vert[tri[2]][1] - vert[tri[0]][1]);
          double bz = (vert[tri[2]][2] - vert[tri[0]][2]);

          double nx = ay * bz - az * by;
          double ny = az * bx - ax * bz;
          double nz = ax * by - ay * bx;
          normalize(nx, ny, nz);

          for (int j = 0; j < 3; j++)
          {
            normals[tri[j]][0] += nx;
            normals[tri[j]][1] += ny;
            normals[tri[j]][2] += nz;
          }
        }

        for (int i = 0; i < num_verts; i++)
          normalize(normals[i][0], normals[i][1], normals[i][2]);
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
        if(force_reset || view_not_reset) { // Reset 3d view.
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
          double tan_fovy_half = tan((double) fovy / 2.0 / 180.0 * M_PI);
          double max_allowed_height = (zfar-3)*tan_fovy_half;
          if(max_axial * xzscale > max_allowed_height || (max_axial * xzscale < 0.1 && !is_constant) )
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
        double tan_fovy_half = tan((double) fovy / 2.0 / 180.0 * M_PI);
        double tan_fovx_half = tan_fovy_half * (double) output_width / output_height;

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
          glGetDoublev( GL_MODELVIEW_MATRIX, aabb_base );
          glPopMatrix();

          // Go through the transformed corners and find the biggest distance of its perspective projection center to the origin.
          GLdouble *coord_ptr = &aabb_base[0];
          for (int j = 0; j < 4; j++)
          {
            double perspective_center_to_origin_dist = fabs(coord_ptr[0]) / tan_fovx_half + coord_ptr[2];
            if(perspective_center_to_origin_dist > optimal_viewpoint_pos)
              optimal_viewpoint_pos = perspective_center_to_origin_dist;

            perspective_center_to_origin_dist = fabs(coord_ptr[1]) / tan_fovy_half + coord_ptr[2];
            if(perspective_center_to_origin_dist > optimal_viewpoint_pos)
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
        if(mode3d)
          yscale *= sc;
        else if(contours)
          cont_step *= sc;
        refresh();
      }

      void ScalarView::set_min_max_range(double min, double max)
      {
        /// \todo allow settin min = max, in which case draw the corresponding contour.
        if(fabs(max-min) < 1e-8)
        {
          this->warn("Range (%f, %f) is too narrow: adjusted to (%f, %f)", min, max, min-0.5, max);
          min -= 0.5;
        }
        View::set_min_max_range(min, max);
      }

      void ScalarView::init_lighting()
      {
        float light_specular[] = {  1.0f, 1.0f, 1.0f, 1.0f };
        float light_ambient[]  = {  0.1f, 0.1f, 0.1f, 1.0f };
        float light_diffuse[]  = {  1.0f, 1.0f, 1.0f, 1.0f };

        glEnable(GL_LIGHT0);
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
        glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
        glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
        //  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

        float material_ambient[]  = { 0.5f, 0.5f, 0.5f, 1.0f };
        float material_diffuse[]  = { 0.8f, 0.8f, 0.8f, 1.0f };
        float material_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };

        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128);
        glDisable(GL_COLOR_MATERIAL);

        glShadeModel(GL_SMOOTH);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
        if(GLEW_EXT_separate_specular_color)
        {
          glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL_EXT, GL_SEPARATE_SPECULAR_COLOR_EXT);
        }
      }

      void ScalarView::on_key_down(unsigned char key, int x, int y)
      {
        VIEWER_GUI(TwSetCurrentWndID(tw_wnd_id));
        VIEWER_GUI_CALLBACK(TwEventKeyboardGLUT(key, x, y))
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
            {
              reset_view(true);
              refresh();
              break;
            }

          case 'f':
            set_palette_filter(pal_filter != GL_LINEAR);
            break;

          case 'k':
            contours = !contours;
            refresh();
            break;

          case 'i':
            show_element_info = !show_element_info;
            if(!mode3d)
              refresh();
            break;

          case 'b':
            show_aabb = !show_aabb;
            if(mode3d)
              refresh();
            break;

          case '3':
            {
              mode3d = !mode3d;
              dragging = scaling = false;
              if(mode3d)
              {
                lin->lock_data();
                if(normals == NULL)
                  calculate_normals(lin->get_vertices(), lin->get_num_vertices(), lin->get_triangles(), lin->get_num_triangles());
                lin->unlock_data();
              }
              refresh();
              break;
            }

          case '*':
            lin->lock_data();
            if(mode3d)
              yscale *= D3DV_SCALE_STEP_COEF;
            else if(contours)
              cont_step *= D3DV_SCALE_STEP_COEF;
            lin->unlock_data();
            refresh();
            break;

          case '/':
            lin->lock_data();
            if(mode3d)
              yscale /= D3DV_SCALE_STEP_COEF;
            else if(contours)
              cont_step /= D3DV_SCALE_STEP_COEF;
            lin->unlock_data();
            refresh();
            break;

          default:
            View::on_key_down(key, x, y);
            break;
          }
        }
      }

      void ScalarView::on_special_key(int key, int x, int y)
      {
        VIEWER_GUI(TwSetCurrentWndID(tw_wnd_id));
        VIEWER_GUI_CALLBACK(TwEventSpecialGLUT(key, x, y))
        {
          View::on_special_key(key, x, y);
        }
      }

      void ScalarView::on_mouse_move(int x, int y)
      {
        VIEWER_GUI(TwSetCurrentWndID(tw_wnd_id));
        if(mode3d && (dragging || scaling || panning))
        {
          if(dragging)
          {
            yrot += 0.4 * (x - mouse_x);
            xrot += 0.4 * (y - mouse_y);

            if(xrot < -90) xrot = -90;
            else if(xrot > 90) xrot = 90;
          }
          else if(scaling)
          {
            ztrans += 0.01 * (mouse_y - y);
            if(ztrans > -0.25) ztrans = -0.25;
            else if(ztrans < -7) ztrans = -7;
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
          VIEWER_GUI_CALLBACK(TwEventMouseMotionGLUT(x, y))
          {
            if(!mode3d && show_edges && !dragging && !scaling && !panning)
            {
              if(allow_node_selection)
              {
                VertexNodeInfo* found_node = find_nearest_node_in_range((float)untransform_x(x), (float)untransform_y(y), (float)(node_pixel_radius / scale));
                if(found_node != pointed_vertex_node)
                  pointed_vertex_node = found_node;
              }
              else
                pointed_vertex_node = NULL;
              refresh();
            }
            else
            {
              View::on_mouse_move(x, y);
            }
          }
        }
      }

      void ScalarView::on_left_mouse_down(int x, int y)
      {
        VIEWER_GUI(TwSetCurrentWndID(tw_wnd_id));
        VIEWER_GUI_CALLBACK(TwEventMouseButtonGLUT(GLUT_LEFT_BUTTON, GLUT_DOWN, x, y))
        {
          View::on_left_mouse_down(x, y);
        }
      }

      void ScalarView::on_left_mouse_up(int x, int y)
      {
        VIEWER_GUI(TwSetCurrentWndID(tw_wnd_id));
        VIEWER_GUI_CALLBACK(TwEventMouseButtonGLUT(GLUT_LEFT_BUTTON, GLUT_UP, x, y))
        {
          View::on_left_mouse_up(x, y);
        }
      }

      void ScalarView::on_right_mouse_down(int x, int y)
      {
        VIEWER_GUI(TwSetCurrentWndID(tw_wnd_id));
        VIEWER_GUI_CALLBACK(TwEventMouseButtonGLUT(GLUT_RIGHT_BUTTON, GLUT_DOWN, x, y))
        {
          //handle node selection
          if(allow_node_selection && pointed_vertex_node != NULL)
          {
            if(pointed_vertex_node->selected) //deselct node
            {
# ifdef ENABLE_VIEWER_GUI
              if(pointed_vertex_node->tw_bar != NULL)
              {
                TwDeleteBar((TwBar*)pointed_vertex_node->tw_bar);
                pointed_vertex_node->tw_bar = NULL;
              }
# endif
              pointed_vertex_node->selected = false;
            }
            else { //select node
# ifdef ENABLE_VIEWER_GUI
              stringstream str;
              str << "Vertex node: " << pointed_vertex_node->id;
              TwBar* bar = TwNewBar(str.str().c_str());
              TwAddVarRO(bar, "X", TW_TYPE_FLOAT, &pointed_vertex_node->x, NULL);
              TwAddVarRO(bar, "Y", TW_TYPE_FLOAT, &pointed_vertex_node->y, NULL);
              int intTmp = 0; TwSetParam(bar, NULL, "iconifiable", TW_PARAM_INT32, 1, &intTmp);
              pointed_vertex_node->tw_bar = bar;
# endif
              pointed_vertex_node->selected = true;
            }
            refresh();
          }
          else { //avoid mixing of functionality
            View::on_right_mouse_down(x, y);
          }
        }
      }

      void ScalarView::on_right_mouse_up(int x, int y)
      {
        VIEWER_GUI(TwSetCurrentWndID(tw_wnd_id));
        VIEWER_GUI_CALLBACK(TwEventMouseButtonGLUT(GLUT_RIGHT_BUTTON, GLUT_UP, x, y))
        {
          View::on_right_mouse_up(x, y);
        }
      }

      void ScalarView::on_middle_mouse_down(int x, int y)
      {
        VIEWER_GUI(TwSetCurrentWndID(tw_wnd_id));
        VIEWER_GUI_CALLBACK(TwEventMouseButtonGLUT(GLUT_MIDDLE_BUTTON, GLUT_DOWN, x, y))
        {
          if(!mode3d) return;
          dragging = scaling = false;
          panning = true;
        }
      }

      void ScalarView::on_middle_mouse_up(int x, int y)
      {
        VIEWER_GUI(TwSetCurrentWndID(tw_wnd_id));
        VIEWER_GUI_CALLBACK(TwEventMouseButtonGLUT(GLUT_MIDDLE_BUTTON, GLUT_UP, x, y))
        {
          panning = false;
        }
      }

      void ScalarView::on_reshape(int width, int height)
      {
        VIEWER_GUI(TwSetCurrentWndID(tw_wnd_id));
        VIEWER_GUI(TwWindowSize(width, height));

        View::on_reshape(width, height);
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