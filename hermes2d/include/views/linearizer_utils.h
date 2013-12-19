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

#ifndef __H2D_LINEARIZER_UTILS_H
#define __H2D_LINEARIZER_UTILS_H

#include "common.h"
#include "compat.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      /// Standard "quality" defining constants.
      const double HERMES_EPS_VERYLOW = 0.25;
      const double HERMES_EPS_LOW = 0.05;
      const double HERMES_EPS_NORMAL = 0.01;
      const double HERMES_EPS_HIGH = 0.005;
      const double HERMES_EPS_VERYHIGH = 0.001;

#ifndef LINEARIZER_DATA_TYPE
#define LINEARIZER_DATA_TYPE double
#endif

      /// Typedefs used throughout the Linearizer functionality.
      template<typename Scalar>
      struct ScalarLinearizerDataDimensions
      {
      };

      template<>
      struct ScalarLinearizerDataDimensions<float>
      {
        static const int dimension = 1;

        typedef float3x3 triangle_t;
        typedef float2x3 edge_t;
        typedef float3 vertex_t;
      };

      template<>
      struct ScalarLinearizerDataDimensions<double>
      {
        static const int dimension = 1;

        typedef double3x3 triangle_t;
        typedef double2x3 edge_t;
        typedef double3 vertex_t;
      };

      template<typename Scalar>
      struct VectorLinearizerDataDimensions
      {
      };

      template<>
      struct VectorLinearizerDataDimensions<float>
      {
        static const int dimension = 2;

        typedef float3x4 triangle_t;
        typedef float2x4 edge_t;
        typedef float4 vertex_t;
      };

      template<>
      struct VectorLinearizerDataDimensions<double>
      {
        static const int dimension = 2;

        typedef double3x4 triangle_t;
        typedef double2x4 edge_t;
        typedef double4 vertex_t;
      };

      typedef int3 internal_vertex_info_t;
      typedef int3 triangle_indices_t;

      template<typename LinearizerDataDimensions>
      class HERMES_API ThreadLinearizerMultidimensional;

#define MAX_LINEARIZER_DIVISION_LEVEL 6

      class HERMES_API LinearizerCriterion
      {
      public:
        LinearizerCriterion(bool adaptive);
        double error_tolerance;
        int refinement_level;
        bool adaptive;
      };

      class HERMES_API LinearizerCriterionAdaptive : public LinearizerCriterion
      {
      public:
        LinearizerCriterionAdaptive(double error_tolerance);
      };

      class HERMES_API LinearizerCriterionFixed : public LinearizerCriterion
      {
      public:
        LinearizerCriterionFixed(int refinement_level);
      };
    }
  }
}
#endif
