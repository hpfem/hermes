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
#include "stream_view.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      StreamView::StreamView(const char* title, WinGeom* wg)
        : View(title, wg), vec(NULL)
      {
        lines = false;
        pmode = false;
        num_stream = 0;
        root_x_min = 1e100;
        root_y_min = 1e100;
        root_x_max = -1e100;
        root_y_max = -1e100;
        root = NULL;
      }

      StreamView::StreamView(char* title, WinGeom* wg)
        : View(title, wg), vec(NULL)
      {
        lines = false;
        pmode = false;
        num_stream = 0;
        root_x_min = 1e100;
        root_y_min = 1e100;
        root_x_max = -1e100;
        root_y_max = -1e100;
        root = NULL;
      }

      void StreamView::show(MeshFunctionSharedPtr<double> xsln, MeshFunctionSharedPtr<double> ysln, int marker, double step, double eps)
      {
        if(this->vec == NULL)
          this->vec = new Vectorizer;
        if(xsln == ysln)
          throw Hermes::Exceptions::Exception("Identical solutions passed to the two-argument version of show(). This is most likely a mistake.");
        show(xsln, ysln, marker, step, eps, H2D_FN_VAL_0, H2D_FN_VAL_0);
      }

      bool StreamView::is_in_triangle(int idx, double x, double y, double3& bar)
      {
        double4* vert = vec->get_vertices();
        int3* xtris = vec->get_triangles();
        int3& tri = xtris[idx];
        double x1 = vert[tri[0]][0], x2 = vert[tri[1]][0], x3 = vert[tri[2]][0];
        double y1 = vert[tri[0]][1], y2 = vert[tri[1]][1], y3 = vert[tri[2]][1];
        double jac = ((x1 - x3)*(y2 - y3) - (x2 - x3)*(y1 - y3));
        double eps = jac * Hermes::Epsilon;
        double a = ((y2 - y3) * (x - x3) - (x2 - x3) * (y - y3));
        double b = ((x1 - x3) * (y - y3) - (y1 - y3) * (x - x3));
        double c = jac - a - b;
        bar[0] = a / jac; bar[1] = b / jac; bar[2] = c / jac;
        if((a >= - eps && a <= jac + eps) && (b >= - eps && b <= jac + eps) && (c >= - eps && c <= jac + eps))
          return true;
        else
          return false;
      }

      void StreamView::add_element_to_tree(Node* father, int e_idx, double x_min, double x_max, double y_min, double y_max)
      {
        double4* vert = vec->get_vertices();
        int3* xtris = vec->get_triangles();
        if(father->leaf == true)
        {
          father->elements[father->num_elem++] = e_idx;
          if(father->num_elem >= 100) // too many elements
          {
            father->leaf = false;
            for (int k = 0; k < 2; k++)
            {
              father->sons[k] = new Node;
              father->sons[k]->level = father->level + 1;
              father->sons[k]->leaf = true;
              father->sons[k]->num_elem = 0;
            }
            for (int i = 0; i < father->num_elem; i++)
              add_element_to_tree(father, father->elements[i], x_min, x_max, y_min, y_max);
          }
        }
        else
        {
          double x_mid = (x_min + x_max)/2;
          double y_mid = (y_min + y_max)/2;
          int3& tri = xtris[e_idx];
          if(father->level % 2) // level = 1, 3, 5, ...
          {
            if(vert[tri[0]][1] <= y_mid || vert[tri[1]][1] <= y_mid || vert[tri[2]][1] <= y_mid)
              add_element_to_tree(father->sons[0], e_idx, x_min, x_max, y_min, y_mid);
            if(vert[tri[0]][1] >= y_mid || vert[tri[1]][1] >= y_mid || vert[tri[2]][1] >= y_mid)
              add_element_to_tree(father->sons[1], e_idx, x_min, x_max, y_mid, y_max);
          }
          else // level 0, 2, 4, ...
          {
            if(vert[tri[0]][0] <= x_mid || vert[tri[1]][0] <= x_mid || vert[tri[2]][0] <= x_mid)
              add_element_to_tree(father->sons[0], e_idx, x_min, x_mid, y_min, y_max);
            if(vert[tri[0]][0] >= x_mid || vert[tri[1]][0] >= x_mid || vert[tri[2]][0] >= x_mid)
              add_element_to_tree(father->sons[1], e_idx, x_mid, x_max, y_min, y_max);
          }
        }
      }

      void StreamView::build_tree()
      {
        root->leaf = true;
        root->level = 0;
        root->num_elem = 0;
        for (int i = 0; i < vec->get_num_triangles(); i++)
        {
          add_element_to_tree(root, i, root_x_min, root_x_max, root_y_min, root_y_max);
        }
      }

      int StreamView::find_triangle_in_tree(double x, double y, Node* father, double x_min, double x_max, double y_min, double y_max, double3& bar)
      {
        if(father->leaf == true)
        {
          double4* vert = vec->get_vertices();
          int3* xtris = vec->get_triangles();
          for (int idx = 0; idx < father->num_elem; idx++)
          {
            int i = father->elements[idx];
            if(is_in_triangle(i, x, y, bar)) return i;
          }
          return -1;
        }
        else
        {
          double x_mid = (x_max + x_min)/2;
          double y_mid = (y_max + y_min)/2;
          if(father->level % 2) // level = 1, 3, 5, ...
            if(y <= y_mid)
              return find_triangle_in_tree(x, y, father->sons[0], x_min, x_max, y_min, y_mid, bar);
            else
              return find_triangle_in_tree(x, y, father->sons[1], x_min, x_max, y_mid, y_max, bar);
          else // level 0, 2, 4, ...
            if(x <= x_mid)
              return find_triangle_in_tree(x, y, father->sons[0], x_min, x_mid, y_min, y_max, bar);
            else
              return find_triangle_in_tree(x, y, father->sons[1], x_mid, x_max, y_min, y_max, bar);
        }
      }

      bool StreamView::get_solution_values(double x, double y, double& xval, double& yval)
      {
        double4* vert = vec->get_vertices();
        int3* xtris = vec->get_triangles();
        double3 bar;
        int e_idx;
        if((e_idx = find_triangle_in_tree(x, y, root, root_x_min, root_x_max, root_y_min, root_y_max, bar)) == -1) return false;
        int3& tri = xtris[e_idx];
        xval = bar[0] * vert[tri[0]][2] + bar[1] * vert[tri[1]][2] + bar[2] * vert[tri[2]][2];
        yval = bar[0] * vert[tri[0]][3] + bar[1] * vert[tri[1]][3] + bar[2] * vert[tri[2]][3];
        return true;
      }

      void StreamView::delete_tree(Node* father)
      {
        if(father->leaf == false)
        {
          delete_tree(father->sons[0]);
          delete_tree(father->sons[1]);
        }
        delete [] father;
      }

      int StreamView::create_streamline(double x_start, double y_start, int idx)
      {
        double ODE_EPS = 1e-5;
        double tau = initial_tau;
        double x = x_start;
        double y = y_start;
        int k = 0;
        const int buffer_length = 5000;
        double2* buffer = new double2[buffer_length];
        bool tau_ok = false;
        bool end = false;
        bool almost_end_of_domain = false;

        double k1, k2, k3, k4, k5;
        double l1, l2, l3, l4, l5;
        double x1, x2, x3, x4, x5;
        double y1, y2, y3, y4, y5;

        while(1)
        {
          if(get_solution_values(x, y, k1, l1) == false) // point (x, y) does not lie in the domain
          {
            tau = tau/2;  // draw streamline to the end of the domain
            if(tau < min_tau) break;
            continue;
          }
          if(fabs(k1) / max_mag  < 1e-5 && fabs(l1) / max_mag < 1e-5) break;  // stop streamline when zero solution

          // add new point to steamline
          buffer[k][0] = x;
          buffer[k][1] = y;
          k++;
          if(k >= buffer_length) break;

          do
          {
            if(almost_end_of_domain)  // draw streamline to the end of the domain
            {
              almost_end_of_domain = false;
              tau = tau/2;
              if(tau < min_tau) { end = true; break; }
            }

            // Merson's adaptive Runge-Kutta method
            x1 = x + 1.0/3.0 * tau * k1; y1 = y + 1.0/3.0 * tau * l1;
            if(get_solution_values(x1, y1, k2, l2) == false) {  almost_end_of_domain = true;  continue;  }
            x2 = x + 1.0/6.0 * tau * (k1 + k2); y2 = y + 1.0/6.0 * tau * (l1 + l2);
            if(get_solution_values(x2, y2, k3, l3) == false) {  almost_end_of_domain = true;  continue;  }
            x3 = x + tau * (0.125 * k1 + 0.375 * k3); y3 = y + tau * (0.125 * l1 + 0.375 * l3);
            if(get_solution_values(x3, y3, k4, l4) == false) {  almost_end_of_domain = true;  continue;  }
            x4 = x + tau * (0.5 * k1 - 1.5 * k3 + 2.0 * k4); y4 = y + tau * (0.5 * l1 - 1.5 * l3 + 2.0 * l4);
            if(get_solution_values(x4, y4, k5, l5) == false) {  almost_end_of_domain = true;  continue;  }
            x5 = x + tau * 1.0/6.0 * (k1 + 4.0 * k4 + k5); y5 = y + tau * 1.0/6.0 * (l1 + 4.0 * l4 + l5);

            // error according to Merson
            double x_err = 1.0/5.0 * (x4 - x5) / (root_x_max - root_x_min);
            double y_err = 1.0/5.0 * (y4 - y5) / (root_y_max - root_y_min);
            double err = std::max(fabs(x_err), fabs(y_err));
            if(err < ODE_EPS)
            {
              tau_ok = true;  x = x5;  y = y5;
            }
            else
            {
              tau_ok = false; tau = tau/2;
              if(tau < min_tau)  tau = min_tau;
            }

            if(err < ODE_EPS/32)
            {
              // new tau according to Merson
              tau = 0.8 * tau * pow(ODE_EPS/err, 0.2);
              if(tau > max_tau)  tau = max_tau;
            }
          }
          while (!tau_ok);
          if(end) break; // get out from both while cycles
        }

        streamlines[idx] = new double2[k];
        memcpy(streamlines[idx], buffer, k*sizeof(double2));
        delete [] buffer;

        return k;
      }

      static double4* comp_vert;
      static int compare(const void* p1, const void* p2)
      {
        const int3* e1 = ((const int3*) p1);
        const int3* e2 = ((const int3*) p2);
        double x1 = comp_vert[(*e1)[0]][0];
        double y1 = comp_vert[(*e1)[0]][1];
        double x2 = comp_vert[(*e2)[0]][0];
        double y2 = comp_vert[(*e2)[0]][1];
        if(x1 < x2) return -1;
        if(x1 == x2 && y1 < y2) return -1;
        if(x1 > x2) return 1;
        if(x1 == x2 && y1 > y2) return 1;
        if(x1 == x2 && y1 == y2) return 0;
        throw Hermes::Exceptions::Exception("internal error: reached end of non-void function");
        return 0;
      }

      int StreamView::find_initial_edge(int num_edges, int3* edges)
      {
        int i, j;
        double4* vert = vec->get_vertices();
        for (i = 0; i < num_edges; i++)
        {
          if(edges[i][2] == 0) // not visited yet
          {
            for (j = 0; j < num_edges; j++)
              if(vert[edges[j][1]][0] == vert[edges[i][0]][0] && vert[edges[j][1]][1] == vert[edges[i][0]][1])
                break;
            if(j == num_edges) return i;
          }
        }
        return -1;
      }

      static int find_next_edge(int num_edges, int3* edges, int b_idx)
      {
        int3 key;
        key[0] = b_idx;
        int3* edge = (int3*) bsearch(&key, edges, num_edges, sizeof(int3), compare);
        if(edge == NULL)
          return -1; // not found
        else
          return edge - edges;
      }

      void StreamView::find_initial_points(int marker, double step, double2*& initial_points)
      {
        int k = 0;
        int ne = vec->get_num_edges();
				int2* edges = vec->get_edges();
        int* edge_markers = vec->get_edge_markers();
        double4* vert = vec->get_vertices();
        int3* bnd_edges = new int3[ne];
        for (int i = 0; i < ne; i++)
        {
          if(edge_markers[i] == marker)
          {
            bnd_edges[k][0] = edges[i][0];
            bnd_edges[k][1] = edges[i][1];
            bnd_edges[k++][2] = 0; // not visited
          }
        }
        int num_edges = k;

        // sort edges by their first vertex coordinates
        // sort first by x, then by y
        comp_vert = vert;
        qsort(bnd_edges, num_edges, sizeof(int3), compare);

        // create list of initial points
        int buffer_length = 1000;
        initial_points = new double2[buffer_length];
        k = 0;
        int idx;
        while ((idx = find_initial_edge(num_edges, bnd_edges)) != -1)
        {
          double tmp_step = step;
          do
          {
            bnd_edges[idx][2] = 1; // visited
            double ax = vert[bnd_edges[idx][0]][0]; double bx = vert[bnd_edges[idx][1]][0];
            double ay = vert[bnd_edges[idx][0]][1]; double by = vert[bnd_edges[idx][1]][1];
            double len = sqrt(sqr(bx - ax) + sqr(by - ay));
            double remaining_len = len; double init_x = ax; double init_y = ay;
            while (tmp_step < remaining_len)
            {
              initial_points[k][0] = init_x + tmp_step * ((bx - ax) / len);
              initial_points[k][1] = init_y + tmp_step * ((by - ay) / len);
              remaining_len = remaining_len - tmp_step;
              tmp_step = step;
              init_x = initial_points[k][0];
              init_y = initial_points[k][1];
              k++;
            }
            tmp_step = tmp_step - remaining_len;
          }
          while ((idx = find_next_edge(num_edges, bnd_edges, bnd_edges[idx][1])) != -1);
        }
        num_stream = k;

        delete [] bnd_edges;
      }

      void StreamView::show(MeshFunctionSharedPtr<double> xsln, MeshFunctionSharedPtr<double> ysln, int marker, double step, double eps, int xitem, int yitem)
      {
        if(vec == NULL)
          vec = new Vectorizer;
        vec->process_solution(xsln, ysln, xitem, yitem, eps);

        vec->lock_data();
        if(range_auto)
        {
          range_min = vec->get_min_value();
          range_max = vec->get_max_value();
        }
        vec->calc_vertices_aabb(&vertices_min_x, &vertices_max_x, &vertices_min_y, &vertices_max_y);

        // create streamlines
        double4* vert = vec->get_vertices();
        for (int i = 0; i < vec->get_num_vertices(); i++)
        {
          if(vert[i][0] < root_x_min) root_x_min = vert[i][0];
          if(vert[i][0] > root_x_max) root_x_max = vert[i][0];
          if(vert[i][1] < root_y_min) root_y_min = vert[i][1];
          if(vert[i][1] > root_y_max) root_y_max = vert[i][1];
        }

        initial_tau = std::max(root_x_max - root_x_min, root_y_max - root_y_min) / 100;
        max_tau = initial_tau * 10;
        min_tau = initial_tau / 50;
        max_mag = vec->get_max_value();

        this->tick();
        root = new Node;
        build_tree();

        double2* initial_points;
        find_initial_points(marker, step, initial_points);

        streamlines = (double2**) malloc(sizeof(double2*) * (num_stream));
        streamlength = (int*) malloc(sizeof(int) * (num_stream));
        for (int i = 0; i < num_stream; i++)
          streamlength[i] =  create_streamline(initial_points[i][0], initial_points[i][1], i);

        delete [] initial_points;

        vec->unlock_data();

        create();
        update_layout();
        reset_view(false);
        refresh();
        wait_for_draw();
      }

      void StreamView::add_streamline(double x, double y)
      {
        if(root == NULL)
          throw Hermes::Exceptions::Exception("Function add_streamline must be called after StreamView::show().");
        this->tick();
        streamlines = (double2**) realloc(streamlines, sizeof(double2*) * (num_stream + 1));
        streamlength = (int*) realloc(streamlength, sizeof(int) * (num_stream + 1));
        streamlength[num_stream] = create_streamline(x, y, num_stream);
        num_stream++;
        refresh();
      }

      static int n_vert(int i) { return (i + 1) % 3; }
      static int p_vert(int i) { return (i + 2) % 3; }

      void StreamView::on_display()
      {
        set_ortho_projection();
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_TEXTURE_1D);
        glPolygonMode(GL_FRONT_AND_BACK, pmode ? GL_LINE : GL_FILL);

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
        }

        // value range
        double min = range_min, max = range_max;
        if(range_auto) { min = vec->get_min_value(); max = vec->get_max_value(); }
        double irange = 1.0 / (max - min);
        // special case: constant solution
        if(fabs(min - max) < Hermes::Epsilon) { irange = 1.0; min -= 0.5; }

        // draw all triangles
        int3* xtris = vec->get_triangles();

        glEnable(GL_TEXTURE_1D);
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
        glColor3f(0.5, 0.5, 0.5);
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

        // draw streamlines
        glColor3f(0.0, 0.0, 0.0);
        for (i = 0; i < num_stream; i++)
        {
          glBegin(GL_LINE_STRIP);
          int k = 0;
          for (int j = 0; j < streamlength[i] - 1; j++)
          {
            glVertex2d(transform_x(streamlines[i][j][0]), transform_y(streamlines[i][j][1]));
            glVertex2d(transform_x(streamlines[i][j + 1][0]), transform_y(streamlines[i][j + 1][1]));
          }
          glEnd();
        }

        delete [] tvert;
        vec->unlock_data();
      }

      void StreamView::on_mouse_move(int x, int y)
      {
        View::on_mouse_move(x, y);
      }

      void StreamView::on_key_down(unsigned char key, int x, int y)
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

          // delete last streamline
        case 26: // ctrl z
          if(num_stream > 0)
          {
            num_stream--;
            delete [] streamlines[num_stream];
            refresh();
          }
          break;

        default:
          View::on_key_down(key, x, y);
          break;
        }
      }

      void StreamView::on_left_mouse_down(int x, int y)
      {
        View::on_left_mouse_down(x, y);

        // adding streamline (initial point set at (x, y))
        if(!scale_focused && glutGetModifiers() == GLUT_ACTIVE_CTRL)
        {
          this->tick();
          streamlines = (double2**) realloc(streamlines, sizeof(double2*) * (num_stream + 1));
          streamlength = (int*) realloc(streamlength, sizeof(int) * (num_stream + 1));
          streamlength[num_stream] = create_streamline(untransform_x(x), untransform_y(y), num_stream);
          num_stream++;
          refresh();
        }
      }

      const char* StreamView::get_help_text() const
      {
        return
          "StreamView\n\n"
          "Controls:\n"
          "  Left mouse - pan\n"
          "  Right mouse - zoom\n"
          "  CTRL + Left mouse click - add steamline\n"
          "  CTRL z - delete last streamline\n"
          "  C - center image\n"
          "  F - toggle smooth palette\n"
          "  H - render high-quality frame\n"
          "  M - toggle mesh\n"
          "  P - cycle palettes\n"
          "  S - save screenshot\n"
          "  F1 - this help\n"
          "  Esc, Q - quit";
      }

      StreamView::~StreamView()
      {
        delete_tree(root);
        for (int i = 0; i < num_stream; i++)
          delete [] streamlines[i];
        delete [] streamlines;
        delete [] streamlength;
        delete vec;
      }
    }
  }
}

#endif