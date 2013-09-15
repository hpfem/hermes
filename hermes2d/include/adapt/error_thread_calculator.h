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

#ifndef __H2D_ERROR_THREAD_CALCULATOR_H
#define __H2D_ERROR_THREAD_CALCULATOR_H

#include "error_calculator.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar> class ErrorCalculator;
    template<typename Scalar> class NeighborSearch;

    template<typename Scalar>
    class ErrorThreadCalculator
    {
    public:
      ErrorThreadCalculator(ErrorCalculator<Scalar>* errorCalculator);
      ~ErrorThreadCalculator();
      void evaluate_one_state(Traverse::State* current_state);

      class DGErrorCalculator
      {
      public:
        DGErrorCalculator(ErrorThreadCalculator* errorThreadCalculator);

        /// Assemble DG forms.
        void assemble_one_edge();

        /// Initialize neighbors.
        bool init_neighbors();
        void deinit_neighbors();

        void assemble_one_neighbor(unsigned int neighbor_i);
      private:
        void initialize_error_and_norm_functions(NormFormDG<Scalar>* mfs, DiscontinuousFunc<Scalar>* error_func[2], DiscontinuousFunc<Scalar>* norm_func[2]);

        ErrorThreadCalculator* errorThreadCalculator;
        Traverse::State* current_state;

        NeighborSearch<Scalar>** neighbor_searches;
        int num_neighbors;
      };
    private:
      void evaluate_volumetric_forms(Traverse::State* current_state, int order);
      void evaluate_surface_forms_one_edge(Traverse::State* current_state, int order);
      
      void evaluate_volumetric_form(NormFormVol<Scalar>* form, Func<Scalar>* difference_func_i, Func<Scalar>* difference_func_j, Func<Scalar>* rsln_i, Func<Scalar>* rsln_j, double* error, double* norm);
      
      void evaluate_surface_form(NormFormSurf<Scalar>* form, Func<Scalar>* difference_func_i, Func<Scalar>* difference_func_j, Func<Scalar>* rsln_i, Func<Scalar>* rsln_j, double* error, double* norm);
      void initialize_error_and_norm_functions_surf(NormFormSurf<Scalar>* mfs, Func<Scalar>* error_func[2], Func<Scalar>* norm_func[2]);
      
      template<typename NormFormType>
      void initialize_error_and_norm_functions(NormFormType* mf, Func<Scalar>* error_func[2], Func<Scalar>* norm_func[2], int order);
      template<typename NormFormType, typename FuncType>
      void deinitialize_error_and_norm_functions(NormFormType* mf, FuncType* error_func[2], FuncType* norm_func[2]);
      
      void evaluate_DG_form(NormFormDG<Scalar>* form, DiscontinuousFunc<Scalar>* difference_func_i, DiscontinuousFunc<Scalar>* difference_func_j, DiscontinuousFunc<Scalar>* rsln_i, DiscontinuousFunc<Scalar>* rsln_j, double* error, double* norm);

      int n_quadrature_points;
      Geom<double>* geometry;
      double* jacobian_x_weights;
      Solution<Scalar>** slns;
      Solution<Scalar>** rslns;

      Traverse::State* current_state;

      ErrorCalculator<Scalar>* errorCalculator;
    };
  }
}

#endif