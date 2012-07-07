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

// $Id: view.h 1086 2008-10-21 09:05:44Z jakub $

#ifndef __H2D_MESH_VIEW_H
#define __H2D_MESH_VIEW_H
#include "view.h"
#include "../function/exact_solution.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      // you can define NOGLUT to turn off all OpenGL stuff in Hermes2D
#ifndef NOGLUT

      /// \brief Displays a mesh.
      ///
      /// MeshView is a debugging tool for displaying meshes together with its element
      /// id numbers, boundary markers etc.
      ///
      class HERMES_API MeshView : public View
      {
      public:

        MeshView(const char* title = "MeshView", WinGeom* wg = NULL);
        MeshView(char* title, WinGeom* wg = NULL);
        virtual ~MeshView();

        void show(Mesh* mesh);

      protected:

        Linearizer* lin;

        bool b_ids, b_markers, b_elem_mrk;

        struct ObjInfo
        {
          double x, y;
          int id, type;
        };

        ObjInfo* nodes;
        ObjInfo* elems;
        int nn, ne;

        float* get_marker_color(int marker);

        virtual void on_display();
        virtual void on_key_down(unsigned char key, int x, int y);
        virtual void scale_dispatch() {}
        virtual const char* get_help_text() const;

        Mesh* mesh;
      };
#else
      class HERMES_API MeshView : public View
      {
      public:

        MeshView(const char* title = "MeshView", WinGeom* wg = NULL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        MeshView(char* title, WinGeom* wg = NULL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void show(Mesh* mesh) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
      };
#endif
    }
  }
}
#endif