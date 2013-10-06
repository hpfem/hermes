#include "postprocessing.h"
#include "../space/space.h"
#include "../function/solution.h"
#include "../views/scalar_view.h"
#include "refmap.h"
#include "forms.h"
#include "limit_order.h"
#include "discrete_problem/discrete_problem_helpers.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace PostProcessing
    {
      bool VertexBasedLimiter::wider_bounds_on_boundary = false;

      template<typename Scalar>
      Limiter<Scalar>::Limiter(SpaceSharedPtr<Scalar> space, Scalar* solution_vector) : component_count(1)
      {
        spaces.push_back(space);
        this->init(solution_vector);
      }

      template<typename Scalar>
      Limiter<Scalar>::Limiter(Hermes::vector<SpaceSharedPtr<Scalar> > spaces, Scalar* solution_vector) : spaces(spaces), component_count(spaces.size())
      {
        this->init(solution_vector);
      }

      template<typename Scalar>
      void Limiter<Scalar>::init(Scalar* solution_vector_)
      {
        if (solution_vector_)
        {
          try
          {
            int ndof = Space<Scalar>::get_num_dofs(this->spaces);
            Scalar value = solution_vector_[ndof - 1];

            this->solution_vector = new Scalar[Space<Scalar>::get_num_dofs(this->spaces)];
            memcpy(this->solution_vector, solution_vector_, sizeof(Scalar)* Space<Scalar>::get_num_dofs(this->spaces));
          }
          catch (...)
          {
            throw Exceptions::Exception("Wrong combination of space(s) and solution_vector passed to Limiter().");
          }
        }
      }

      template<typename Scalar>
      Limiter<Scalar>::~Limiter()
      {
        if (this->solution_vector)
          delete[] this->solution_vector;
      }

      template<typename Scalar>
      MeshFunctionSharedPtr<Scalar> Limiter<Scalar>::get_solution()
      {
        // A check.
        warn_if(this->component_count > 1, "One solution asked from a Limiter, but multiple solutions exist for limiting.");

        MeshFunctionSharedPtr<Scalar> solution(new Solution<Scalar>());
        Hermes::vector<MeshFunctionSharedPtr<Scalar> > solutions;
        solutions.push_back(solution);
        this->get_solutions(solutions);
        return solutions.back();
      }

      template<typename Scalar>
      void Limiter<Scalar>::get_solutions(Hermes::vector<MeshFunctionSharedPtr<Scalar> > solutions)
      {
        if (solutions.size() != this->component_count)
          throw Exceptions::Exception("Limiter does not have correct number of spaces, solutions.");

        this->limited_solutions = solutions;

        // Processing.
        this->process();

        if (this->limited_solutions.empty() || this->limited_solutions.size() != this->component_count)
          throw Exceptions::Exception("Limiter failed for unknown reason.");
      }

      template<typename Scalar>
      Hermes::vector<int> Limiter<Scalar>::get_changed_element_ids() const
      {
        return this->changed_element_ids;
      }

      template<typename Scalar>
      bool Limiter<Scalar>::isOkay() const
      {
        bool okay = true;

        return true;
      }

      template<typename Scalar>
      Scalar* Limiter<Scalar>::get_solution_vector()
      {
        return this->solution_vector;
      }

      template<typename Scalar>
      void Limiter<Scalar>::set_solution_vector(Scalar* solution_vector_)
      {
        if (solution_vector_)
        {
          try
          {
            int ndof = Space<Scalar>::get_num_dofs(this->spaces);
            Scalar value = solution_vector_[ndof - 1];

            this->solution_vector = new Scalar[Space<Scalar>::get_num_dofs(this->spaces)];
            memcpy(this->solution_vector, solution_vector_, sizeof(Scalar)* Space<Scalar>::get_num_dofs(this->spaces));
          }
          catch (...)
          {
            throw Exceptions::Exception("Wrong combination of space(s) and solution_vector passed to Limiter().");
          }
        }
      }

      VertexBasedLimiter::VertexBasedLimiter(SpaceSharedPtr<double> space, double* solution_vector, int maximum_polynomial_order)
        : Limiter<double>(space, solution_vector)
      {
        this->init(maximum_polynomial_order);
      }

      VertexBasedLimiter::VertexBasedLimiter(Hermes::vector<SpaceSharedPtr<double> > spaces, double* solution_vector, int maximum_polynomial_order)
        : Limiter<double>(spaces, solution_vector)
      {
        this->init(maximum_polynomial_order);
      }

      void VertexBasedLimiter::set_p_coarsening_only()
      {
        this->p_coarsening_only = true;
      }

      void VertexBasedLimiter::init(int maximum_polynomial_order)
      {
        this->maximum_polynomial_order = maximum_polynomial_order;
        this->p_coarsening_only = false;

        // Checking that this is the Taylor shapeset.
        for (int i = 0; i < this->component_count; i++)
        {
          if (this->spaces[i]->get_shapeset()->get_id() != 31)
            throw Exceptions::Exception("VertexBasedLimiter designed for L2ShapesetTaylor. Ignore this exception for unforeseen problems.");
        }

        vertex_min_values = nullptr;
        vertex_max_values = nullptr;

        // This is what is the key aspect of the necessity to use L2ShapesetTaylor (or any other one that uses P_{} also for quads).
        this->mixed_derivatives_count = (maximum_polynomial_order)*(maximum_polynomial_order + 1) / 2;

        this->print_details = false;
      }

      VertexBasedLimiter::~VertexBasedLimiter()
      {
        deallocate_vertex_values();
      }

      void VertexBasedLimiter::print_detailed_info(bool print_details_)
      {
        this->print_details = print_details_;
      }

      Hermes::vector<std::pair<int, double> > VertexBasedLimiter::get_correction_factors() const
      {
        return this->correction_factors;
      }

      void VertexBasedLimiter::process()
      {
        // 0. Preparation.
        // Start by creating temporary solutions and states for paralelism.
        Solution<double>::vector_to_solutions(this->solution_vector, this->spaces, this->limited_solutions);

        // 1. Quadratic
        // Prepare the vertex values for the quadratic part.
        prepare_min_max_vertex_values(true);

        // Use those to incorporate the correction factor.
        Element* e;

        // Vector to remember if there was limiting of the second derivatives.
        Hermes::vector<bool> quadratic_correction_done;

        if (this->get_verbose_output())
          std::cout << "Quadratic correction" << std::endl;

        for (int component = 0; component < this->component_count; component++)
        {
          MeshSharedPtr mesh = this->spaces[component]->get_mesh();
          if (this->get_verbose_output() && this->component_count > 1)
            std::cout << "Component: " << component << std::endl;
          for_all_active_elements(e, mesh)
          {
            bool second_order = H2D_GET_H_ORDER(this->spaces[component]->get_element_order(e->id)) >= 2 || H2D_GET_V_ORDER(this->spaces[component]->get_element_order(e->id)) >= 2;
            if (!second_order)
            {
              quadratic_correction_done.push_back(false);
              continue;
            }

            if (this->get_verbose_output())
              std::cout << "Element: " << e->id << std::endl;

            quadratic_correction_done.push_back(this->impose_quadratic_correction_factor(e, component));

            if (this->get_verbose_output())
              std::cout << std::endl;
          }
        }

        // Adjust the solutions according to the quadratic terms handling.
        Solution<double>::vector_to_solutions(this->solution_vector, this->spaces, this->limited_solutions);


        // 2. Linear
        // Prepare the vertex values for the linear part.
        prepare_min_max_vertex_values(false);

        if (this->get_verbose_output())
          std::cout << "Linear correction" << std::endl;

        int running_i = 0;
        for (int component = 0; component < this->component_count; component++)
        {
          MeshSharedPtr mesh = this->spaces[component]->get_mesh();
          if (this->get_verbose_output() && this->component_count > 1)
            std::cout << "Component: " << component << std::endl;
          for_all_active_elements(e, mesh)
          {
            if (this->get_verbose_output())
              std::cout << "Element: " << e->id << std::endl;

            bool second_order = H2D_GET_H_ORDER(this->spaces[component]->get_element_order(e->id)) >= 2 || H2D_GET_V_ORDER(this->spaces[component]->get_element_order(e->id)) >= 2;
            if (quadratic_correction_done[running_i++] || !second_order)
              this->impose_linear_correction_factor(e, component);

            if (this->get_verbose_output())
              std::cout << std::endl;
          }
        }

        // Create the final solutions.
        Solution<double>::vector_to_solutions(this->solution_vector, this->spaces, this->limited_solutions);
      }

      void VertexBasedLimiter::impose_linear_correction_factor(Element* e, int component)
      {
        double correction_factor = std::numeric_limits<double>::infinity();

        double centroid_value_multiplied = this->get_centroid_value_multiplied(e, component, 0);
        if (this->get_verbose_output())
          std::cout << std::endl << "center-value: " << centroid_value_multiplied << " (" << 0 << ") ";

        Solution<double>* sln = dynamic_cast<Solution<double>*>(this->limited_solutions[component].get());

        for (int i_vertex = 0; i_vertex < e->get_nvert(); i_vertex++)
        {
          if (this->get_verbose_output())
            std::cout << std::endl << "vertex: " << i_vertex;

          Node* vertex = e->vn[i_vertex];
          double x, y;
          RefMap::untransform(e, vertex->x, vertex->y, x, y);

          double vertex_value = sln->get_ref_value_transformed(e, x, y, 0, 0);

          if (this->get_verbose_output())
            std::cout << "\tvalue: " << vertex_value;

          double fraction;
          if (std::abs(vertex_value - centroid_value_multiplied) < Hermes::epsilon)
          {
            fraction = 1.;
            if (this->get_verbose_output())
              std::cout << "\tcenter_value";

          }
          else
            if (vertex_value > centroid_value_multiplied)
            {
              fraction = std::min(1., (this->vertex_max_values[component][vertex->id][0] - centroid_value_multiplied) / (vertex_value - centroid_value_multiplied));
              if (this->get_verbose_output())
                std::cout << "\tmax_value: " << this->vertex_max_values[component][vertex->id][0];
            }
            else
            {
              fraction = std::min(1., (this->vertex_min_values[component][vertex->id][0] - centroid_value_multiplied) / (vertex_value - centroid_value_multiplied));
              if (this->get_verbose_output())
                std::cout << "\tmin_value: " << this->vertex_min_values[component][vertex->id][0];
            }

            correction_factor = std::min(correction_factor, fraction);
        }
        if (this->get_verbose_output())
          std::cout << std::endl << "correction_factor " << correction_factor << std::endl;

        if (correction_factor < (1 - 1e-3))
        {
          AsmList<double> al;
          this->spaces[component]->get_element_assembly_list(e, &al);

          for (int i_basis_fn = 0; i_basis_fn < al.cnt; i_basis_fn++)
          {
            int order = this->spaces[component]->get_shapeset()->get_order(al.idx[i_basis_fn], e->get_mode());
            if (H2D_GET_H_ORDER(order) == 1 || H2D_GET_V_ORDER(order) == 1)
            {
              if (this->p_coarsening_only)
                this->solution_vector[al.dof[i_basis_fn]] = 0.;
              else
                this->solution_vector[al.dof[i_basis_fn]] *= correction_factor;
            }
          }

          this->changed_element_ids.push_back(e->id);
          this->correction_factors.push_back(std::pair<int, double>(1, correction_factor));
        }
      }

      bool VertexBasedLimiter::impose_quadratic_correction_factor(Element* e, int component)
      {
        if (this->get_verbose_output())
          std::cout << "quadratic: ";

        double correction_factor = std::numeric_limits<double>::infinity();

        Solution<double>* sln = dynamic_cast<Solution<double>*>(this->limited_solutions[component].get());

        for (int i_derivative = 1; i_derivative <= 2; i_derivative++)
        {
          double centroid_value_multiplied = this->get_centroid_value_multiplied(e, component, i_derivative);
          if (this->get_verbose_output())
            std::cout << std::endl << "center-value: " << centroid_value_multiplied << " (" << i_derivative << ") ";

          for (int i_vertex = 0; i_vertex < e->get_nvert(); i_vertex++)
          {
            if (this->get_verbose_output())
              std::cout << std::endl << "vertex: " << i_vertex;

            Node* vertex = e->vn[i_vertex];
            double x, y;
            RefMap::untransform(e, vertex->x, vertex->y, x, y);

            double vertex_value = sln->get_ref_value_transformed(e, x, y, 0, i_derivative);

            if (this->get_verbose_output())
              std::cout << "\tvalue: " << vertex_value;

            double fraction;
            if (std::abs(vertex_value - centroid_value_multiplied) < Hermes::epsilon)
            {
              if (this->get_verbose_output())
                std::cout << "\tcenter_value";

              fraction = 1.;
            }
            else
              if (vertex_value > centroid_value_multiplied)
              {
                fraction = std::min(1., (this->vertex_max_values[component][vertex->id][i_derivative] - centroid_value_multiplied) / (vertex_value - centroid_value_multiplied));
                if (this->get_verbose_output())
                  std::cout << "\tmax_value: " << this->vertex_max_values[component][vertex->id][i_derivative];

              }
              else
              {
                fraction = std::min(1., (this->vertex_min_values[component][vertex->id][i_derivative] - centroid_value_multiplied) / (vertex_value - centroid_value_multiplied));
                if (this->get_verbose_output())
                  std::cout << "\tmin_value: " << this->vertex_min_values[component][vertex->id][i_derivative];
              }

              correction_factor = std::min(correction_factor, fraction);
          }
        }

        if (this->get_verbose_output())
          std::cout << std::endl << "correction_factor " << correction_factor << std::endl;

        if (correction_factor < (1 - 1e-3))
        {
          AsmList<double> al;
          this->spaces[component]->get_element_assembly_list(e, &al);

          for (int i_basis_fn = 0; i_basis_fn < al.cnt; i_basis_fn++)
          {
            int order = this->spaces[component]->get_shapeset()->get_order(al.idx[i_basis_fn], e->get_mode());
            if (H2D_GET_H_ORDER(order) == 2 || H2D_GET_V_ORDER(order) == 2)
            {
              if (this->p_coarsening_only)
                this->solution_vector[al.dof[i_basis_fn]] = 0.;
              else
                this->solution_vector[al.dof[i_basis_fn]] *= correction_factor;
            }
          }

          this->changed_element_ids.push_back(e->id);
          this->correction_factors.push_back(std::pair<int, double>(2, correction_factor));
          return true;
        }
        else
          return false;
      }

      void VertexBasedLimiter::prepare_min_max_vertex_values(bool quadratic)
      {
        // Reallocate if calculating quadratic part (first of the two).
        if (quadratic)
        {
          deallocate_vertex_values();
          allocate_vertex_values();
        }

        // Calculate min/max vertex values.
        Element* e;
        for (int component = 0; component < this->component_count; component++)
        {
          MeshSharedPtr mesh = this->spaces[component]->get_mesh();

          for_all_active_elements(e, mesh)
          {
            for (int i_vertex = 0; i_vertex < e->get_nvert(); i_vertex++)
            {
              Node* vertex = e->vn[i_vertex];
              for (int i_derivative = (quadratic ? 1 : 0); i_derivative < (quadratic ? this->mixed_derivatives_count : 1); i_derivative++)
              {
                double element_centroid_value_multiplied = this->get_centroid_value_multiplied(e, component, i_derivative);
                this->vertex_min_values[component][vertex->id][i_derivative] = std::min(this->vertex_min_values[component][vertex->id][i_derivative], element_centroid_value_multiplied);
                this->vertex_max_values[component][vertex->id][i_derivative] = std::max(this->vertex_max_values[component][vertex->id][i_derivative], element_centroid_value_multiplied);
                if (e->en[i_vertex]->bnd && this->wider_bounds_on_boundary)
                {
                  double element_mid_edge_value_multiplied = this->get_edge_midpoint_value_multiplied(e, component, i_derivative, i_vertex);

                  this->vertex_min_values[component][vertex->id][i_derivative] = std::min(this->vertex_min_values[component][vertex->id][i_derivative], element_mid_edge_value_multiplied);
                  this->vertex_max_values[component][vertex->id][i_derivative] = std::max(this->vertex_max_values[component][vertex->id][i_derivative], element_mid_edge_value_multiplied);

                  Node* next_vertex = e->vn[(i_vertex + 1) % e->get_nvert()];
                  this->vertex_min_values[component][next_vertex->id][i_derivative] = std::min(this->vertex_min_values[component][next_vertex->id][i_derivative], element_mid_edge_value_multiplied);
                  this->vertex_max_values[component][next_vertex->id][i_derivative] = std::max(this->vertex_max_values[component][next_vertex->id][i_derivative], element_mid_edge_value_multiplied);
                }
              }
            }
          }
        }
      }

      double VertexBasedLimiter::get_centroid_value_multiplied(Element* e, int component, int mixed_derivative_index)
      {
        if (mixed_derivative_index > 5)
        {
          throw Exceptions::MethodNotImplementedException("VertexBasedLimiter::get_centroid_value_multiplied only works for first and second derivatives.");
          return 0.;
        }

        Solution<double>* sln = dynamic_cast<Solution<double>*>(this->limited_solutions[component].get());
        double result;
        if (e->get_mode() == HERMES_MODE_TRIANGLE)
          result = sln->get_ref_value_transformed(e, CENTROID_TRI_X, CENTROID_TRI_Y, 0, mixed_derivative_index);
        else
          result = sln->get_ref_value_transformed(e, CENTROID_QUAD_X, CENTROID_QUAD_Y, 0, mixed_derivative_index);

        return result;
      }

      double VertexBasedLimiter::get_edge_midpoint_value_multiplied(Element* e, int component, int mixed_derivative_index, int edge)
      {
        if (mixed_derivative_index > 5)
        {
          throw Exceptions::MethodNotImplementedException("VertexBasedLimiter::get_centroid_value_multiplied only works for first and second derivatives.");
          return 0.;
        }

        Solution<double>* sln = dynamic_cast<Solution<double>*>(this->limited_solutions[component].get());
        double result;

        double x;
        double y;

        if (e->get_mode() == HERMES_MODE_TRIANGLE)
        {
          if (edge == 0)
          {
            x = 0.;
            y = -1;
          }
          else if (edge == 1)
          {
            x = 0;
            y = 0;
          }
          else if (edge == 2)
          {
            x = -1.;
            y = 0;
          }
          result = sln->get_ref_value_transformed(e, x, y, 0, mixed_derivative_index);
        }
        else
        {
          if (edge == 0)
          {
            x = 0.;
            y = -1;
          }
          else if (edge == 1)
          {
            x = 1;
            y = 0;
          }
          else if (edge == 2)
          {
            x = 0.;
            y = 1.;
          }
          else if (edge == 3)
          {
            x = -1.;
            y = 0.;
          }
          result = sln->get_ref_value_transformed(e, x, y, 0, mixed_derivative_index);
        }

        return result;
      }

      void VertexBasedLimiter::allocate_vertex_values()
      {
        this->vertex_min_values = new double**[this->component_count];
        for (int i = 0; i < this->component_count; i++)
        {
          this->vertex_min_values[i] = new double*[this->spaces[i]->get_mesh()->get_max_node_id()];

          for (int j = 0; j < this->spaces[i]->get_mesh()->get_max_node_id(); j++)
          {
            this->vertex_min_values[i][j] = new double[this->mixed_derivatives_count];
            for (int k = 0; k < this->mixed_derivatives_count; k++)
              this->vertex_min_values[i][j][k] = std::numeric_limits<double>::infinity();
          }
        }

        this->vertex_max_values = new double**[this->component_count];
        for (int i = 0; i < this->component_count; i++)
        {
          this->vertex_max_values[i] = new double*[this->spaces[i]->get_mesh()->get_max_node_id()];

          for (int j = 0; j < this->spaces[i]->get_mesh()->get_max_node_id(); j++)
          {
            this->vertex_max_values[i][j] = new double[this->mixed_derivatives_count];
            for (int k = 0; k < this->mixed_derivatives_count; k++)
              this->vertex_max_values[i][j][k] = -std::numeric_limits<double>::infinity();
          }
        }
      }

      void VertexBasedLimiter::deallocate_vertex_values()
      {
        if (this->vertex_min_values)
        {
          for (int i = 0; i < this->component_count; i++)
          {
            for (int j = 0; j < this->spaces[i]->get_mesh()->get_max_node_id(); j++)
            {
              delete[] this->vertex_min_values[i][j];
              delete[] this->vertex_max_values[i][j];
            }

            delete[] this->vertex_min_values[i];
            delete[] this->vertex_max_values[i];
          }

          delete[] this->vertex_min_values;
          delete[] this->vertex_max_values;
        }
      }

      template<typename Scalar>
      IntegralCalculator<Scalar>::IntegralCalculator(MeshFunctionSharedPtr<Scalar> source_function, int number_of_integrals) : number_of_integrals(number_of_integrals)
      {
        source_functions.push_back(source_function);
      }

      template<typename Scalar>
      IntegralCalculator<Scalar>::IntegralCalculator(Hermes::vector<MeshFunctionSharedPtr<Scalar> > source_functions, int number_of_integrals) : source_functions(source_functions), number_of_integrals(number_of_integrals)
      {
      }

      template<typename Scalar>
      Scalar* IntegralCalculator<Scalar>::calculate(std::string marker)
      {
        Hermes::vector<std::string> markers;
        markers.push_back(marker);
        return this->calculate(markers);
      }

      template<>
      void IntegralCalculator<double>::add_results(double* results_local, double* results)
      {
        for (int i = 0; i < this->number_of_integrals; i++)
        {
#pragma openmp atomic
          results[i] += results_local[i];
        }
      }

      template<>
      void IntegralCalculator<std::complex<double> >::add_results(std::complex<double> * results_local, std::complex<double> * results)
      {
        for (int i = 0; i < this->number_of_integrals; i++)
        {
#pragma openmp critical integralAddition
          results[i] += results_local[i];
        }
      }

      template<typename Scalar>
      VolumetricIntegralCalculator<Scalar>::VolumetricIntegralCalculator(MeshFunctionSharedPtr<Scalar> source_function, int number_of_integrals) : IntegralCalculator<Scalar>(source_function, number_of_integrals)
      {
      }

      template<typename Scalar>
      VolumetricIntegralCalculator<Scalar>::VolumetricIntegralCalculator(Hermes::vector<MeshFunctionSharedPtr<Scalar> > source_functions, int number_of_integrals) : IntegralCalculator<Scalar>(source_functions, number_of_integrals)
      {
      }

      template<typename Scalar>
      Scalar* VolumetricIntegralCalculator<Scalar>::calculate(Hermes::vector<std::string> markers)
      {
        Hermes::vector<int> internal_markers;
        for (int i = 0; i < markers.size(); i++)
        {
          Hermes::Hermes2D::Mesh::MarkersConversion::IntValid internalMarker = this->source_functions[0]->get_mesh()->get_element_markers_conversion().get_internal_marker(markers[i]);
          if (internalMarker.valid)
            internal_markers.push_back(internalMarker.marker);
        }

        int source_functions_size = this->source_functions.size();
        Traverse trav(source_functions_size);
        int num_states;
        Traverse::State** states = trav.get_states(this->source_functions, num_states);

        for (int i = 0; i < source_functions_size; i++)
          this->source_functions[i]->set_quad_2d(&g_quad_2d_std);

        Scalar* result = (Scalar*)calloc(this->number_of_integrals, sizeof(Scalar));

#pragma omp parallel num_threads(this->num_threads_used)
        {
          RefMap* refmap = new RefMap;
          refmap->set_quad_2d(&g_quad_2d_std);

          int thread_number = omp_get_thread_num();
          int start = (num_states / this->num_threads_used) * thread_number;
          int end = (num_states / this->num_threads_used) * (thread_number + 1);
          if (thread_number == this->num_threads_used - 1)
            end = num_states;

          Hermes::Ord* orders = new Hermes::Ord[this->number_of_integrals];
          Func<Hermes::Ord>** func_ord = (Func<Hermes::Ord>**)malloc(this->number_of_integrals * sizeof(Func<Hermes::Ord>*));

          MeshFunction<Scalar>** source_fuctions_cloned = new MeshFunction<Scalar>*[source_functions_size];
          for (int i = 0; i < source_functions_size; i++)
            source_fuctions_cloned[i] = this->source_functions[i]->clone();
          Func<Scalar>** func = (Func<Scalar>**)malloc(this->number_of_integrals * sizeof(Func<Scalar>*));

          Scalar* result_thread_local = (Scalar*)calloc(this->number_of_integrals, sizeof(Scalar));
          Scalar* result_local = (Scalar*)malloc(this->number_of_integrals * sizeof(Scalar));
          double* jacobian_x_weights;

          for (int state_i = start; state_i < end; state_i++)
          {
            Traverse::State* current_state = states[state_i];
            bool target_marker = false;
            for (int i = 0; i < internal_markers.size(); i++)
              if (current_state->e[0]->marker == internal_markers[i])
              {
                target_marker = true;
                break;
              }
              if (target_marker)
                continue;

              memset(result_local, 0, sizeof(Scalar)* this->number_of_integrals);

              // Set active element.
              for (int i = 0; i < source_functions_size; i++)
              {
                source_fuctions_cloned[i]->set_active_element(current_state->e[i]);
                source_fuctions_cloned[i]->set_transform(current_state->sub_idx[i]);
                refmap->set_active_element(current_state->e[i]);
              }

              // Integration order.
              Hermes::Ord order = Hermes::Ord(refmap->get_inv_ref_order());
              memset(orders, 0, sizeof(nullptr) * this->number_of_integrals);
              for (int i = 0; i < source_functions_size; i++)
                func_ord[i] = init_fn_ord(this->source_functions[i]->get_fn_order());

              this->order(func_ord, orders);

              for (int i = 0; i < source_functions_size; i++)
                order += orders[i];

              int order_int = order.get_order();
              limit_order(order_int, refmap->get_active_element()->get_mode());

              for (int i = 0; i < source_functions_size; i++)
                func[i] = init_fn(source_fuctions_cloned[i], order_int);

              Geom<double>* geometry = init_geom_vol(refmap, order_int);
              int n = init_geometry_points(&refmap, 1, order_int, geometry, jacobian_x_weights);

              this->integral(n, jacobian_x_weights, func, geometry, result_local);

              geometry->free();
              delete geometry;
              delete[] jacobian_x_weights;

              for (int i = 0; i < source_functions_size; i++)
              {
                func_ord[i]->free_ord();
                delete func_ord[i];
                func[i]->free_fn();
                delete func[i];
              }

              for (int i = 0; i < this->number_of_integrals; i++)
                result_thread_local[i] += result_local[i];
          }

          this->add_results(result_thread_local, result);
          ::free(result_thread_local);

          for (int i = 0; i < source_functions_size; i++)
            delete source_fuctions_cloned[i];
          delete[] source_fuctions_cloned;
          delete[] orders;
          ::free(func_ord);
          ::free(func);
          ::free(result_local);
          delete refmap;
        }

        for (int i = 0; i < num_states; i++)
          delete states[i];
        free(states);

        return result;
      }

      template<typename Scalar>
      SurfaceIntegralCalculator<Scalar>::SurfaceIntegralCalculator(MeshFunctionSharedPtr<Scalar> source_function, int number_of_integrals) : IntegralCalculator<Scalar>(source_function, number_of_integrals)
      {
      }

      template<typename Scalar>
      SurfaceIntegralCalculator<Scalar>::SurfaceIntegralCalculator(Hermes::vector<MeshFunctionSharedPtr<Scalar> > source_functions, int number_of_integrals) : IntegralCalculator<Scalar>(source_functions, number_of_integrals)
      {
      }

      template<typename Scalar>
      Scalar* SurfaceIntegralCalculator<Scalar>::calculate(Hermes::vector<std::string> markers)
      {
        Hermes::vector<int> internal_markers;
        for (int i = 0; i < markers.size(); i++)
        {
          Hermes::Hermes2D::Mesh::MarkersConversion::IntValid internalMarker = this->source_functions[0]->get_mesh()->get_element_markers_conversion().get_internal_marker(markers[i]);
          if (internalMarker.valid)
            internal_markers.push_back(internalMarker.marker);
        }

        int source_functions_size = this->source_functions.size();
        Traverse trav(source_functions_size);
        int num_states;
        Traverse::State** states = trav.get_states(this->source_functions, num_states);

        for (int i = 0; i < source_functions_size; i++)
          this->source_functions[i]->set_quad_2d(&g_quad_2d_std);

        Scalar* result = (Scalar*)calloc(this->number_of_integrals, sizeof(Scalar));

#pragma omp parallel num_threads(this->num_threads_used)
        {
          RefMap* refmap = new RefMap;
          refmap->set_quad_2d(&g_quad_2d_std);

          int thread_number = omp_get_thread_num();
          int start = (num_states / this->num_threads_used) * thread_number;
          int end = (num_states / this->num_threads_used) * (thread_number + 1);
          if (thread_number == this->num_threads_used - 1)
            end = num_states;

          Hermes::Ord* orders = new Hermes::Ord[this->number_of_integrals];
          Func<Hermes::Ord>** func_ord = (Func<Hermes::Ord>**)malloc(this->number_of_integrals * sizeof(Func<Hermes::Ord>*));

          MeshFunction<Scalar>** source_fuctions_cloned = new MeshFunction<Scalar>*[source_functions_size];
          for (int i = 0; i < source_functions_size; i++)
            source_fuctions_cloned[i] = this->source_functions[i]->clone();
          Func<Scalar>** func = (Func<Scalar>**)malloc(this->number_of_integrals * sizeof(Func<Scalar>*));

          Scalar* result_thread_local = (Scalar*)calloc(this->number_of_integrals, sizeof(Scalar));
          Scalar* result_local = (Scalar*)malloc(this->number_of_integrals * sizeof(Scalar));
          double* jacobian_x_weights;

          for (int state_i = start; state_i < end; state_i++)
          {
            Traverse::State* current_state = states[state_i];
            if (!current_state->isBnd)
              continue;

            // Set active element.
            for (int i = 0; i < source_functions_size; i++)
            {
              source_fuctions_cloned[i]->set_active_element(current_state->e[i]);
              source_fuctions_cloned[i]->set_transform(current_state->sub_idx[i]);
              refmap->set_active_element(current_state->e[i]);
            }

            Hermes::Ord order = Hermes::Ord(refmap->get_inv_ref_order());

            for (int edge = 0; edge < current_state->e[0]->nvert; edge++)
            {
              if (!current_state->bnd[edge])
                continue;

              bool target_marker = false;
              for (int i = 0; i < internal_markers.size(); i++)
                if (current_state->e[0]->en[edge]->marker == internal_markers[i])
                {
                  target_marker = true;
                  break;
                }
                if (target_marker)
                  continue;

                memset(result_local, 0, sizeof(Scalar)* this->number_of_integrals);

                // Integration order.
                memset(orders, 0, sizeof(nullptr) * this->number_of_integrals);
                for (int i = 0; i < source_functions_size; i++)
                  func_ord[i] = init_fn_ord(this->source_functions[i]->get_fn_order());

                this->order(func_ord, orders);

                for (int i = 0; i < source_functions_size; i++)
                  order += orders[i];

                int order_int = order.get_order();
                limit_order(order_int, refmap->get_active_element()->get_mode());

                double3* tan;
                Geom<double>* geometry;
                int n = init_surface_geometry_points(&refmap, 1, order_int, edge, current_state->e[0]->en[edge]->marker, geometry, jacobian_x_weights);

                for (int i = 0; i < source_functions_size; i++)
                  func[i] = init_fn(source_fuctions_cloned[i], order_int);

                this->integral(n, jacobian_x_weights, func, geometry, result_local);

                geometry->free();
                delete geometry;
                delete[] jacobian_x_weights;

                for (int i = 0; i < source_functions_size; i++)
                {
                  func_ord[i]->free_ord();
                  delete func_ord[i];
                  func[i]->free_fn();
                  delete func[i];
                }

                for (int i = 0; i < this->number_of_integrals; i++)
                  result_thread_local[i] += .5 * result_local[i];
            }
          }

          this->add_results(result_thread_local, result);
          ::free(result_thread_local);

          for (int i = 0; i < source_functions_size; i++)
            delete source_fuctions_cloned[i];
          delete[] source_fuctions_cloned;
          delete[] orders;
          ::free(func_ord);
          ::free(func);
          ::free(result_local);
          delete refmap;
        }

        for (int i = 0; i < num_states; i++)
          delete states[i];
        free(states);

        return result;
      }

      template class HERMES_API Limiter<double>;
      template class HERMES_API VolumetricIntegralCalculator<double>;
      template class HERMES_API SurfaceIntegralCalculator<double>;
      template class HERMES_API Limiter<std::complex<double> >;
      template class HERMES_API VolumetricIntegralCalculator<std::complex<double> >;
      template class HERMES_API SurfaceIntegralCalculator<std::complex<double> >;
    }
  }
}
