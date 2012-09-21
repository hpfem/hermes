// This file is part of Hermes2D.
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@unr.edu>
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

// $Id: view4.cpp 1086 2008-10-21 09:05:44Z jakub $

#ifndef NOGLUT

#include <GL/freeglut.h>
#include "global.h"
#include "mesh_view.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      MeshView::MeshView(const char* title, WinGeom* wg)
        : View(title, wg), lin(NULL)
      {
        nodes = elems = NULL;
        b_scale = false;
        b_ids = false;
        b_markers = true;
        b_elem_mrk = false;
      }

      MeshView::MeshView(char* title, WinGeom* wg)
        : View(title, wg), lin(NULL)
      {
        nodes = elems = NULL;
        b_scale = false;
        b_ids = false;
        b_markers = true;
        b_elem_mrk = false;
      }

      MeshView::~MeshView()
      {
        if(nodes != NULL) delete [] nodes;
        if(elems != NULL) delete [] elems;
        if(lin != NULL)
          delete this->lin;
      }

      void MeshView::show(Mesh* mesh)
      {
        ZeroSolution<double> sln(mesh);
        if(mesh == NULL) throw Hermes::Exceptions::Exception("mesh == NULL in MeshView::show().");
        if(mesh->get_max_element_id() == 0) throw Hermes::Exceptions::Exception("Attempt to visualize empty mesh in MeshView::show().");

        this->mesh = mesh;

        if(lin == NULL)
          lin = new Linearizer();

        lin->process_solution(&sln);
        lin->lock_data();
        lin->calc_vertices_aabb(&vertices_min_x, &vertices_max_x, &vertices_min_y, &vertices_max_y);
        lin->unlock_data();

        int i;

        if(elems != NULL) delete [] elems;
        ne = mesh->get_max_element_id() + 1;
        elems = new ObjInfo[ne];
        for (i = 0; i < ne; i++)
          elems[i].id = -1;

        int active_element_cnt = 0;
        float min_error = -1, max_error = -1;
        Element* e;
        for_all_active_elements(e, mesh)
        {
          ObjInfo* oi = elems + e->id;
          oi->id = e->id;
          oi->type = e->marker;
          oi->x = oi->y = 0.0;
          for (unsigned i = 0; i < e->get_nvert(); i++)
          {
            oi->x += e->vn[i]->x;
            oi->y += e->vn[i]->y;
          }
          oi->x /= e->get_nvert();
          oi->y /= e->get_nvert();
        }

        create();
        update_layout();
        refresh();
        reset_view(false);
        wait_for_draw();
      }

      void MeshView::set_b_elem_mrk(bool set)
      {
        if(b_ids) 
          b_ids = false;
        b_elem_mrk = !b_elem_mrk;
        refresh();
      } 

      void MeshView::on_display()
      {
        set_ortho_projection();
        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);

        // transform all vertices
        lin->lock_data();
        int i, nv = lin->get_num_vertices();
        double3* vert = lin->get_vertices();
        double2* tvert = new double2[nv];
        for (i = 0; i < nv; i++)
        {
          tvert[i][0] = transform_x(vert[i][0]);
          tvert[i][1] = transform_y(vert[i][1]);
        }

        // draw all triangles
        int3* tris = lin->get_triangles();
        glColor3f(0.9f, 0.9f, 0.9f);
        glBegin(GL_TRIANGLES);
        for (i = 0; i < lin->get_num_triangles(); i++)
        {
          glVertex2d(tvert[tris[i][0]][0], tvert[tris[i][0]][1]);
          glVertex2d(tvert[tris[i][1]][0], tvert[tris[i][1]][1]);
          glVertex2d(tvert[tris[i][2]][0], tvert[tris[i][2]][1]);
        }
        glEnd();

        // draw all edges
        glLineStipple(5, 0x5555);
        int3* edges = lin->get_edges();
        for (i = 0; i < lin->get_num_edges(); i++)
        {
          int mrk = b_markers ? edges[i][2] : 0;

          if(!edges[i][2] &&
            ((tvert[edges[i][0]][1] == tvert[edges[i][1]][1] &&
            tvert[edges[i][0]][0] < tvert[edges[i][1]][0]) ||
            tvert[edges[i][0]][1] < tvert[edges[i][1]][1])) continue;

          float* color = get_marker_color(mrk);
          glColor3f(color[0], color[1], color[2]);
          glLineWidth(mrk ? 1.5f : 1.0f);
          glBegin(GL_LINES);
          glVertex2d(tvert[edges[i][0]][0], tvert[edges[i][0]][1]);
          glVertex2d(tvert[edges[i][1]][0], tvert[edges[i][1]][1]);
          glEnd();

          if(mrk)
          {
            glEnable(GL_LINE_STIPPLE);
            glColor3f(0.4f, 0.4f, 0.4f);
            glBegin(GL_LINES);
            glVertex2d(tvert[edges[i][0]][0], tvert[edges[i][0]][1]);
            glVertex2d(tvert[edges[i][1]][0], tvert[edges[i][1]][1]);
            glEnd();
            glDisable(GL_LINE_STIPPLE);
          }
        }
        glLineWidth(1.0);

        if(b_ids)  // draw element ids
        {
          glColor3f(0, 0, 0);
          for (i = 0; i < ne; i++)
          {
            if(elems[i].id < 0) continue;
            char text[20];
            sprintf(text, "#%d", elems[i].id);
            draw_text(transform_x(elems[i].x), transform_y(elems[i].y), text, 0);
          }
        }
        else if(b_elem_mrk)  // draw element markers
        {
          glColor3f(0, 0, 0);
          for (i = 0; i < ne; i++)
          {
            if(elems[i].id < 0) continue;
            char text[2000];
            sprintf(text, "%s", mesh->get_element_markers_conversion().get_user_marker(elems[i].type).marker.c_str());
            draw_text(transform_x(elems[i].x), transform_y(elems[i].y), text, 0);
          }
        }

        delete [] tvert;
        lin->unlock_data();
      }

      void MeshView::on_key_down(unsigned char key, int x, int y)
      {
        switch (key)
        {
        case 'c':
          reset_view(true);
          refresh();
          break;

        case 'b':
          b_markers = !b_markers;
          refresh();
          break;

        case 'i':
          if(b_elem_mrk) b_elem_mrk = false;
          b_ids = !b_ids;
          refresh();
          break;

        case 'm':
          if(b_ids) b_ids = false;
          b_elem_mrk = !b_elem_mrk;
          refresh();
          break;

        default:
          View::on_key_down(key, x, y);
          break;
        }
      }

      float* MeshView::get_marker_color(int marker)
      {
        static float edgecol[3] = { 0.3f, 0.3f, 0.3f };
        static float randcol[3];
        static float mc[8][3] =
        {
          { 1.0f, 0.3f, 0.3f },
          { 0.0f, 0.9f, 0.0f },
          { 0.0f, 0.0f, 0.7f },
          { 1.0f, 1.0f, 0.2f },
          { 0.7f, 0.0f, 0.0f },
          { 0.0f, 0.5f, 0.0f },
          { 0.3f, 0.5f, 1.0f },
          { 0.8f, 0.8f, 0.0f },
        };

        if(marker == 0)
          return edgecol;
        else if(marker > 0 && marker < 8)
          return mc[marker];
        else
        {
          srand(marker + 2);
          randcol[0] = (float) rand() / RAND_MAX;
          randcol[1] = (float) rand() / RAND_MAX;
          randcol[2] = (float) rand() / RAND_MAX;
          return randcol;
        }
      }

      const char* MeshView::get_help_text() const
      {
        return
          "MeshView\n"
          "Controls:\n"
          "  Left mouse - pan\n"
          "  Right mouse - zoom\n"
          "  C - center image\n"
          "  H - render high-quality frame\n"
          "  B - toggle boundary markers\n"
          "  I - toggle element IDs\n"
          "  M - toggle element markers\n"
          "  S - save screenshot\n"
          "  F1 - this help\n"
          "  Esc, Q - quit";
      }
    }
  }
}
#endif // NOGLUT