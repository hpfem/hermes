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

#ifndef __H2D_POSTPROCESSING_H
#define __H2D_POSTPROCESSING_H

#include "../function/mesh_function.h"
#include "../space/space.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// \brief Solution (mesh function) Postprocessing capabilities.
    namespace PostProcessing
    {
      template<typename Scalar>
      class HERMES_API Limiter
        : public Hermes::Mixins::TimeMeasurable,
        public Hermes::Mixins::Loggable,
        public Hermes::Hermes2D::Mixins::Parallel,
        public Hermes::Hermes2D::Mixins::StateQueryable
      {
      public:
        Limiter(SpaceSharedPtr<Scalar> space, Scalar* solution_vector);
        Limiter(Hermes::vector<SpaceSharedPtr<Scalar> > spaces, Scalar* solution_vector);
        virtual ~Limiter();

        /// Get the zero-th solution.
        MeshFunctionSharedPtr<Scalar> get_solution();
        /// Get all solutions.
        void get_solutions(Hermes::vector<MeshFunctionSharedPtr<Scalar> > solutions);
        /// Get changed element ids.
        Hermes::vector<int> get_changed_element_ids() const;

        /// Helpers for state querying.
        virtual bool isOkay() const;
        inline std::string getClassName() const { return "Limiter"; }

        /// Set the solution vector.
        void set_solution_vector(Scalar* sln);

        /// Get the solution vector.
        Scalar* get_solution_vector();

      protected:
        int component_count;
        Hermes::vector<SpaceSharedPtr<Scalar> > spaces;
        Scalar* solution_vector;
        Hermes::vector<MeshFunctionSharedPtr<Scalar> > limited_solutions;
        Hermes::vector<int> changed_element_ids;

        virtual void process() = 0;

      private:
        void init(Scalar* solution_vector_);
      };

      class HERMES_API VertexBasedLimiter
        : public Limiter<double>
      {
      public:
        VertexBasedLimiter(SpaceSharedPtr<double> space, double* solution_vector, int maximum_polynomial_order);
        VertexBasedLimiter(Hermes::vector<SpaceSharedPtr<double> > spaces, double* solution_vector, int maximum_polynomial_order);
        virtual ~VertexBasedLimiter();
        Hermes::vector<std::pair<int, double> > get_correction_factors() const;
        void print_detailed_info(bool print_details = true);
        int maximum_polynomial_order;
        void set_p_coarsening_only();
        static bool wider_bounds_on_boundary;

      private:
        bool p_coarsening_only;

        void init(int maximum_polynomial_order);

        void process();

        void prepare_min_max_vertex_values(bool quadratic);

        /// Get mean value of the mixed derivative (mixed_derivative_index) on element e, of the "component" - component
        /// of the solution.
        double get_centroid_value_multiplied(Element* e, int component, int mixed_derivative_index);
        
        double get_edge_midpoint_value_multiplied(Element* e, int component, int mixed_derivative_index, int edge);

        void impose_linear_correction_factor(Element* e, int component);

        /// Return if there was a need to limit the second derivatives.
        bool impose_quadratic_correction_factor(Element* e, int component);

        double*** vertex_min_values;
        double*** vertex_max_values;
        void allocate_vertex_values();
        void deallocate_vertex_values();

        int mixed_derivatives_count;
        Hermes::vector<std::pair<int, double> > correction_factors;
        bool print_details;
      };
    }
  }
}

#endif
