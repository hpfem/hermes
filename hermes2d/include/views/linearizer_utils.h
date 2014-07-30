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
/*! \file linearizer_utils.h
\brief File containing utilities for class.
*/
#ifndef __H2D_LINEARIZER_UTILS_H
#define __H2D_LINEARIZER_UTILS_H

#include "hermes_common.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      /// Standard "quality" defining constants.
      const double HERMES_EPS_VERYLOW = 0.25;
      /// Standard "quality" defining constants.
      const double HERMES_EPS_LOW = 0.05;
      /// Standard "quality" defining constants.
      const double HERMES_EPS_NORMAL = 0.01;
      /// Standard "quality" defining constants.
      const double HERMES_EPS_HIGH = 0.005;
      /// Standard "quality" defining constants.
      const double HERMES_EPS_VERYHIGH = 0.001;

#ifndef LINEARIZER_DATA_TYPE
#define LINEARIZER_DATA_TYPE double
#endif

      /// We refine a quad directionally (horizontally, vertically) only if the error in one direction is this much larger than in the other.
#ifndef LINEARIZER_DIRECTIONAL_QUAD_REFINEMENT_REQUIREMENT
#define LINEARIZER_DIRECTIONAL_QUAD_REFINEMENT_REQUIREMENT 5.0
#endif

      /// Very important constant putting an upper bound on the maximum number of successive element division (when dealing with a higher-order FEM solution).
#define MAX_LINEARIZER_DIVISION_LEVEL 6

      /// Typedefs used throughout the Linearizer functionality.
      template<typename Scalar>
      struct ScalarLinearizerDataDimensions
      {
      };

      /// Typedefs used throughout the Linearizer functionality.
      template<>
      struct ScalarLinearizerDataDimensions < float >
      {
        static const int dimension = 1;

        typedef float3x3 triangle_t;
        typedef float2x3 edge_t;
        typedef float3 vertex_t;
      };

      /// Typedefs used throughout the Linearizer functionality.
      template<>
      struct ScalarLinearizerDataDimensions < double >
      {
        static const int dimension = 1;

        typedef double3x3 triangle_t;
        typedef double2x3 edge_t;
        typedef double3 vertex_t;
      };

      /// Typedefs used throughout the Linearizer functionality.
      template<typename Scalar>
      struct VectorLinearizerDataDimensions
      {
      };

      /// Typedefs used throughout the Linearizer functionality.
      template<>
      struct VectorLinearizerDataDimensions < float >
      {
        static const int dimension = 2;

        typedef float3x4 triangle_t;
        typedef float2x4 edge_t;
        typedef float4 vertex_t;
      };

      /// Typedefs used throughout the Linearizer functionality.
      template<>
      struct VectorLinearizerDataDimensions < double >
      {
        static const int dimension = 2;

        typedef double3x4 triangle_t;
        typedef double2x4 edge_t;
        typedef double4 vertex_t;
      };

      /// Typedefs used throughout the Linearizer functionality.
      typedef int3 internal_vertex_info_t;
      /// Typedefs used throughout the Linearizer functionality.
      typedef int3 triangle_indices_t;

      template<typename LinearizerDataDimensions>
      class HERMES_API ThreadLinearizerMultidimensional;

      /// \brief Abstract class for criterion according to which the linearizer stops dividing elements at some point
      /// Class is not abstract per say, but works as a base class for the following classes.
      class HERMES_API LinearizerCriterion
      {
      public:
        LinearizerCriterion(bool adaptive);
        double error_tolerance;
        int refinement_level;
        bool adaptive;
      };

      /// \brief Adaptive Linearizer criterion - error tolerance (see further) where the element division stops
      /// Error tolerance here is the relative improvement of quality that the currently proposed element division would bring.
      /// If this quantity is below the specified tolerance, the currently proposed division is not made and the division algorithm for the current element stops.
      class HERMES_API LinearizerCriterionAdaptive : public LinearizerCriterion
      {
      public:
        LinearizerCriterionAdaptive(double error_tolerance);
      };

      /// \brief Simple Linearizer criterion - every element is refined exactly the same number of times.
      /// This number is specified in the constructor.
      class HERMES_API LinearizerCriterionFixed : public LinearizerCriterion
      {
      public:
        LinearizerCriterionFixed(int refinement_level);
      };
    }
  }
}
#endif
