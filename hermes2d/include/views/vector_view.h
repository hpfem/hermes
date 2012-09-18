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

#ifndef __H2D_VECTOR_VIEW_H
#define __H2D_VECTOR_VIEW_H
#include "view.h"
namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      // you can define NOGLUT to turn off all OpenGL stuff in Hermes2D
#ifndef NOGLUT

      /// \brief Visualizes a vector PDE solution.
      ///
      /// VectorView is a visualization window for all vector-valued PDE solutions.
      ///
      class HERMES_API VectorView : public View
      {
      public:

        VectorView(const char* title = "VectorView", WinGeom* wg = NULL);
        VectorView(char* title, WinGeom* wg = NULL);
        ~VectorView();

        void show(MeshFunction<double>* vsln, double eps = HERMES_EPS_NORMAL);
        void show(MeshFunction<double>* xsln, MeshFunction<double>* ysln, double eps = HERMES_EPS_NORMAL);
        void show(MeshFunction<double>* xsln, MeshFunction<double>* ysln, double eps, int xitem, int yitem);

        inline void set_grid_type(bool hexa) { this->hexa = hexa; refresh(); };
        void set_mode(int mode);
      protected:
        Vectorizer* vec;

        double gx, gy, gs;
        bool hexa; ///< false - quad grid, true - hexa grid
        int mode;  ///< 0 - magnitude is on the background, 1 - arrows are colored, 2 - no arrows, just magnitude on the background
        bool lines, pmode;
        double length_coef; ///< for extending or shortening arrows

        void plot_arrow(double x, double y, double xval, double yval, double max, double min, double gs);

        virtual void on_display();
        virtual void on_mouse_move(int x, int y);
        virtual void on_key_down(unsigned char key, int x, int y);
        virtual const char* get_help_text() const;
      };
#else
      class HERMES_API VectorView : public View
      {
      public:
        VectorView(const char* title = "VectorView", WinGeom* wg = NULL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        VectorView(char* title, WinGeom* wg = NULL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        void show(MeshFunction<double>* vsln, double eps = HERMES_EPS_NORMAL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void show(MeshFunction<double>* xsln, MeshFunction<double>* ysln, double eps = HERMES_EPS_NORMAL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void show(MeshFunction<double>* xsln, MeshFunction<double>* ysln, double eps, int xitem, int yitem) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        inline void set_grid_type(bool hexa) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
      };
#endif
    }
  }
}
#endif