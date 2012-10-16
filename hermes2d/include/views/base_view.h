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

#ifndef __H2D_BASE_VIEW_H
#define __H2D_BASE_VIEW_H

#include "../global.h"
#include "scalar_view.h"
#include "../space/space.h"
#include "../shapeset/precalc.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      // you can define NOGLUT to turn off all OpenGL stuff in Hermes2D
#ifndef NOGLUT

      /// \brief Visualizes the basis functions of a space.
      ///
      /// BaseView is a debugging tool for the visualization of the basis functions
      /// of a given space.
      ///
      template<typename Scalar>
      class HERMES_API BaseView : public ScalarView
      {
      public:

        BaseView(const char* title = "BaseView", WinGeom* wg = NULL);
        BaseView(char* title, WinGeom* wg = NULL);

        void show(const Space<Scalar>* space, double eps = HERMES_EPS_LOW, int item = H2D_FN_VAL_0);

        ~BaseView() { free(); }

      protected:

        const Space<Scalar>* space;
        PrecalcShapeset* pss;
        Solution<Scalar>* sln;

        double eps;
        int ndof, item;
        int base_index;

        std::string basic_title;

        void free();
        void update_solution();
        void update_title();

        virtual void on_special_key(int key, int x, int y);
        virtual const char* get_help_text() const;
      };
#else
      template<typename Scalar>
      class HERMES_API BaseView : public ScalarView
      {
      public:
        BaseView(const char* title = "BaseView", WinGeom* wg = NULL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        BaseView(char* title, WinGeom* wg = NULL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        void show(const Space<Scalar>* space, double eps = HERMES_EPS_LOW, int item = H2D_FN_VAL_0) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
      };

#endif
    }
  }
}
#endif