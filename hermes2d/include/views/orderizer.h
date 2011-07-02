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

#ifndef __H2D_ORDERIZER_H
#define __H2D_ORDERIZER_H

#include "linearizer.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      /// Like the Linearizer, but generates a triangular mesh showing polynomial
      /// orders in a space, hence the funky name.
      ///
      class HERMES_API Orderizer : public Linearizer<double>
      {
      public:

        Orderizer();
        ~Orderizer();

        void process_space(Space<double>* space);
        void process_space(Space<std::complex<double> >* space);

        int get_labels(int*& lvert, char**& ltext, double2*& lbox) const
        { lvert = this->lvert; ltext = this->ltext; lbox = this->lbox; return nl; };

        virtual void save_data(const char* filename);
        virtual void load_data(const char* filename);
        /// Saves a MeshFunction (Solution, Filter) in VTK format.
        virtual void save_orders_vtk(Space<double>* space, const char* file_name);
        /// This function is used by save_solution_vtk().
        virtual void save_data_vtk(const char* file_name);

      protected:

        char  buffer[1000];
        char* labels[11][11];

        int  nl, cl1, cl2, cl3;
        int* lvert;
        char** ltext;
        double2* lbox;

      };
    }
  }
}
#endif
