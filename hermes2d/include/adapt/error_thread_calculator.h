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
    template<typename Scalar>
    class ErrorThreadCalculator
    {
    public:
      ErrorThreadCalculator(Hermes::vector<MeshFunctionSharedPtr<Scalar> >& coarse_solutions, Hermes::vector<MeshFunctionSharedPtr<Scalar> >& fine_solutions, ErrorCalculator<Scalar>* errorCalculator);
      ~ErrorThreadCalculator();
      void evaluate_one_state(Traverse::State* current_state);

      void evaluate_volumetric_forms(Traverse::State* current_state, int order);
      void evaluate_surface_forms_one_edge(Traverse::State* current_state, int order);
      void evaluate_DG_forms_one_edge(Traverse::State* current_state, int order);
      void evaluate_form(MatrixForm<Scalar>* form, Func<Scalar>* difference_func_i, Func<Scalar>* difference_func_j, Func<Scalar>* rsln_i, Func<Scalar>* rsln_j, double* error, double* norm);

      int n_quadrature_points;
      Geom<double>* geometry;
      double* jacobian_x_weights;
      Solution<Scalar>** slns;
      Solution<Scalar>** rslns;

      ErrorCalculator<Scalar>* errorCalculator;

      class DGErrorCalculator : public DiscreteProblemDGAssembler<Scalar>
      {
      public:
        DGErrorCalculator(ErrorThreadCalculator* errorThreadCalculator);
        void assemble_one_neighbor(bool edge_processed, unsigned int neighbor_i, NeighborSearch<Scalar>** neighbor_searches);
      private:
        ErrorThreadCalculator* errorThreadCalculator;
      };
    };
  }
}

#endif