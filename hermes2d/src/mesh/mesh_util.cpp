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

#include "mesh.h"
#include "refmap.h"

namespace Hermes
{
  namespace Hermes2D
  {
    void MeshUtil::assign_curve(Node* en, Curve* curve, int p1, int p2)
    {
      // assign the arc to the elements sharing the edge node
      for (unsigned int node_i = 0; node_i < 2; node_i++)
      {
        Element* e = en->elem[node_i];
        if (e == nullptr)
          continue;

        if (e->cm == nullptr)
        {
          e->cm = new CurvMap;
          e->cm->toplevel = true;
          e->cm->order = 4;
        }

        int idx = -1;
        for (unsigned j = 0; j < e->get_nvert(); j++)
        {
          if (e->en[j] == en)
          {
            idx = j;
            break;
          }
        }
        assert(idx >= 0);

        if (e->vn[idx]->id == p1)
          e->cm->curves[idx] = curve;
        else
        {
          Curve* curve_rev = MeshUtil::reverse_curve(curve);
          e->cm->curves[idx] = curve_rev;
        }
      }
    }

    Curve* MeshUtil::reverse_curve(Curve* curve)
    {
      if (curve->type == ArcType)
      {
        Arc* toReturn = new Arc((Arc*)curve);
        ((Arc*)toReturn)->angle = -((Arc*)toReturn)->angle;
        for (int i = 0; i < toReturn->np; i++)
        {
          toReturn->pt[((Arc*)curve)->np - 1 - i][0] = ((Arc*)curve)->pt[i][0];
          toReturn->pt[((Arc*)curve)->np - 1 - i][1] = ((Arc*)curve)->pt[i][1];
          toReturn->pt[((Arc*)curve)->np - 1 - i][2] = ((Arc*)curve)->pt[i][2];
        }

        return toReturn;
      }
      else
      {
        Nurbs* toReturn = new Nurbs((Nurbs*)curve);

        for (int i = 0; i < ((Nurbs*)curve)->np; i++)
        {
          toReturn->pt[((Nurbs*)curve)->np - 1 - i][0] = ((Nurbs*)curve)->pt[i][0];
          toReturn->pt[((Nurbs*)curve)->np - 1 - i][1] = ((Nurbs*)curve)->pt[i][1];
          toReturn->pt[((Nurbs*)curve)->np - 1 - i][2] = ((Nurbs*)curve)->pt[i][2];
        }

        for (int i = 0; i < ((Nurbs*)curve)->nk; i++)
          toReturn->kv[i] = ((Nurbs*)curve)->kv[i];
        for (int i = ((Nurbs*)curve)->degree + 1; i < ((Nurbs*)curve)->nk - ((Nurbs*)curve)->degree - 1; i++)
          toReturn->kv[((Nurbs*)curve)->nk - 1 - i] = 1.0 - ((Nurbs*)curve)->kv[i];

          return toReturn;
      }
    }

    Node* MeshUtil::get_base_edge_node(Element* base, int edge)
    {
      while (!base->active) // we need to go down to an active element
      {
        int son1, son2;
        get_edge_sons(base, edge, son1, son2);
        base = base->sons[son1];
      }
      return base->en[edge];
    }

    int MeshUtil::get_edge_sons(Element* e, int edge, int& son1, int& son2)
    {
      assert(!e->active);

      if (!e->is_triangle())
      {
        if (e->sons[2] == nullptr) // horz quad
        {
          if (edge == 0 || edge == 2) { son1 = edge >> 1;   return 1; }
          else if (edge == 1) { son1 = 0; son2 = 1; return 2; }
          else { son1 = 1; son2 = 0; return 2; }
        }
        else if (e->sons[0] == nullptr) // vert quad
        {
          if (edge == 1 || edge == 3) { son1 = (edge == 1) ? 3 : 2; return 1; }
          else if (edge == 0) { son1 = 2; son2 = 3; return 2; }
          else { son1 = 3; son2 = 2; return 2; }
        }
      }

      // triangle or 4-son quad
      son1 = edge;
      son2 = e->next_vert(edge);
      return 2;
    }

    Arc* MeshUtil::load_arc(MeshSharedPtr mesh, int id, Node** en, int p1, int p2, double angle, bool skip_check)
    {
      Arc* curve = new Arc(angle);

      *en = mesh->peek_edge_node(p1, p2);

      if(*en == nullptr)
      {
        if(!skip_check)
          throw Hermes::Exceptions::MeshLoadFailureException("Curve #%d: edge %d-%d does not exist.", id, p1, p2);
        else
          return nullptr;
      }
    
      // edge endpoints control points.
      curve->pt[0][0] = mesh->nodes[p1].x;
      curve->pt[0][1] = mesh->nodes[p1].y;
      curve->pt[0][2] = 1.0;
      curve->pt[2][0] = mesh->nodes[p2].x;
      curve->pt[2][1] = mesh->nodes[p2].y;
      curve->pt[2][2] = 1.0;

      // read the arc angle
      double a = (180.0 - angle) / 180.0 * M_PI;

      // generate one inner control point
      double x = 1.0 / std::tan(a * 0.5);
      curve->pt[1][0] = 0.5*((curve->pt[2][0] + curve->pt[0][0]) + (curve->pt[2][1] - curve->pt[0][1]) * x);
      curve->pt[1][1] = 0.5*((curve->pt[2][1] + curve->pt[0][1]) - (curve->pt[2][0] - curve->pt[0][0]) * x);
      curve->pt[1][2] = Hermes::cos((M_PI - a) * 0.5);

      return curve;
    }

    MeshHashGrid::MeshHashGrid(Mesh* mesh) : mesh_seq(mesh->get_seq())
    {
      mesh->calc_bounding_box();

      // create grid
      double interval_len_x = (mesh->top_right_x - mesh->bottom_left_x) / GRID_SIZE;
      double interval_len_y = (mesh->top_right_y - mesh->bottom_left_y) / GRID_SIZE;

      intervals_x[0] = mesh->bottom_left_x;
      intervals_y[0] = mesh->bottom_left_y;

      for (int i = 1; i < GRID_SIZE; i++)
      {
        intervals_x[i] = intervals_x[i - 1] + interval_len_x;
        intervals_y[i] = intervals_y[i - 1] + interval_len_y;
      }

      intervals_x[GRID_SIZE] = mesh->top_right_x;
      intervals_y[GRID_SIZE] = mesh->top_right_y;

      for (int i = 0; i < GRID_SIZE; i++)
      {
        for (int j = 0; j < GRID_SIZE; j++)
        {
          m_grid[i][j] = new MeshHashGridElement(intervals_x[i], intervals_y[j], intervals_x[i + 1], intervals_y[j + 1]);
        }
      }

      // assign elements
      Element *element;
      double2 p1, p2;
      int x_min, x_max, y_min, y_max;
      for_all_active_elements(element, mesh)
      {
        elementBoundingBox(element, p1, p2);

        x_min = 0;
        while (intervals_x[x_min + 1] < p1[0])
          x_min++;

        x_max = GRID_SIZE - 1;
        while (intervals_x[x_max] > p2[0])
          x_max--;

        y_min = 0;
        while (intervals_y[y_min + 1] < p1[1])
          y_min++;

        y_max = GRID_SIZE - 1;
        while (intervals_y[y_max] > p2[1])
          y_max--;

        assert((x_max >= x_min) && (y_max >= y_min));

        for (int i = x_min; i <= x_max; i++)
        {
          for (int j = y_min; j <= y_max; j++)
          {
            assert(m_grid[i][j]->belongs(element));
            m_grid[i][j]->insert(element);
          }
        }
      }
    }

    MeshHashGrid::~MeshHashGrid()
    {
      for (int i = 0; i < GRID_SIZE; i++)
      {
        for (int j = 0; j < GRID_SIZE; j++)
        {
          delete m_grid[i][j];
        }
      }
    }

    MeshHashGridElement::MeshHashGridElement(double lower_left_x, double lower_left_y, double upper_right_x, double upper_right_y, int depth) : lower_left_x(lower_left_x), lower_left_y(lower_left_y), upper_right_x(upper_right_x), upper_right_y(upper_right_y), m_depth(depth), m_active(true), element_count(0)
    {
      this->elements = new Element*[MAX_ELEMENTS];
      for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
        m_sons[i][j] = nullptr;
    }

    MeshHashGridElement::~MeshHashGridElement()
    {
      if (elements)
        delete[] elements;
      for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
      if (m_sons[i][j])
        delete m_sons[i][j];
    }

    bool MeshHashGridElement::belongs(Element *element)
    {
      double2 p1, p2;
      MeshHashGrid::elementBoundingBox(element, p1, p2);
      return ((p1[0] <= upper_right_x) && (p2[0] >= lower_left_x) && (p1[1] <= upper_right_y) && (p2[1] >= lower_left_y));
    }

    void MeshHashGridElement::insert(Element *element)
    {
      if (m_active)
      {
        elements[element_count++] = element;

        if ((element_count >= MAX_ELEMENTS) && (m_depth < MAX_DEPTH))
        {
          m_active = false;
          double xx[3] = { lower_left_x, (lower_left_x + upper_right_x) / 2., upper_right_x };
          double yy[3] = { lower_left_y, (lower_left_y + upper_right_y) / 2., upper_right_y };
          for (int i = 0; i < 2; i++)
          {
            double x0 = xx[i];
            double x1 = xx[i + 1];
            for (int j = 0; j < 2; j++)
            {
              double y0 = yy[j];
              double y1 = yy[j + 1];

              assert(m_sons[i][j] == nullptr);

              m_sons[i][j] = new MeshHashGridElement(x0, y0, x1, y1, m_depth + 1);

              for (int elem_i = 0; elem_i < this->element_count; elem_i++)
              {
                if (m_sons[i][j]->belongs(elements[elem_i]))
                  m_sons[i][j]->insert(elements[elem_i]);
              }
            }
          }

          delete[] elements;
          elements = nullptr;
        }
      }
      else
      {
        for (int i = 0; i < 2; i++)
        {
          for (int j = 0; j < 2; j++)
          {
            assert(m_sons[i][j]);
            if (m_sons[i][j]->belongs(element))
              m_sons[i][j]->insert(element);
          }
        }
      }
    }

    bool MeshHashGridElement::belongs(double x, double y)
    {
      return (x >= lower_left_x) && (x <= upper_right_x) && (y <= lower_left_y) && (y >= upper_right_y);
    }

    Element* MeshHashGridElement::getElement(double x, double y)
    {
      if (m_active)
      {
        for (int elem_i = 0; elem_i < this->element_count; elem_i++)
        if (RefMap::is_element_on_physical_coordinates(elements[elem_i], x, y))
          return elements[elem_i];

        return nullptr;
      }
      else
      {
        Element* element;
        for (int i = 0; i < 2; i++)
        {
          for (int j = 0; j < 2; j++)
          {
            element = m_sons[i][j]->getElement(x, y);
            if (element)
              return element;
          }
        }
        return nullptr;
      }
    }

    void MeshHashGrid::elementBoundingBox(Element *element, double2 &p1, double2 &p2)
    {
      p1[0] = p2[0] = element->vn[0]->x;
      p1[1] = p2[1] = element->vn[0]->y;

      for (int i = 1; i < element->get_nvert(); i++)
      {
        double xx = element->vn[i]->x;
        double yy = element->vn[i]->y;
        if (xx > p2[0])
          p2[0] = xx;
        if (xx < p1[0])
          p1[0] = xx;
        if (yy > p2[1])
          p2[1] = yy;
        if (yy < p1[1])
          p1[1] = yy;
      }

      if (element->is_curved())
      {
        // todo: should be improved
        double diameter_x = p2[0] - p1[0];
        p2[0] += diameter_x;
        p1[0] -= diameter_x;

        double diameter_y = p2[1] - p1[1];
        p2[1] += diameter_y;
        p1[1] -= diameter_y;
      }
    }

    Element* MeshHashGrid::getElement(double x, double y)
    {
      int i = 0;
      while ((i < GRID_SIZE) && (intervals_x[i + 1] < x))
        i++;

      int j = 0;
      while ((j < GRID_SIZE) && (intervals_y[j + 1] < y))
        j++;

      // this means that x or y is outside mesh, but it can hapen
      if ((i >= GRID_SIZE) || (j >= GRID_SIZE))
        return nullptr;
      else
        return m_grid[i][j]->getElement(x, y);
    }

    int MeshHashGrid::get_mesh_seq() const
    {
      return this->mesh_seq;
    }

    MarkerArea::MarkerArea(Mesh *mesh, int marker) : mesh_seq(mesh->get_seq())
    {
      area = 0;
      Element* elem;
      for_all_active_elements(elem, mesh)
      {
        if (elem->marker == marker)
        {
          elem->calc_area(true);
          area += elem->area;
        }
      }
    }

    int MarkerArea::get_mesh_seq() const
    {
      return this->mesh_seq;
    }

    double MarkerArea::get_area() const
    {
      return this->area;
    }
  }
}
