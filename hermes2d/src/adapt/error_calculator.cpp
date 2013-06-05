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

#include "adapt/error_calculator.h"
#include "function/exact_solution.h"
#include "adapt/error_thread_calculator.h"
#include "norm_form.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    ErrorCalculator<Scalar>::ErrorCalculator(CalculatedErrorType errorType) :
      errorType(errorType),
      elements_stored(false),
      element_references(NULL),
      errors_squared_sum(0.0),
      norms_squared_sum(0.0)
    {
      memset(errors, 0, sizeof(double*) * H2D_MAX_COMPONENTS);
      memset(norms, 0, sizeof(double*) * H2D_MAX_COMPONENTS);
    }

    template<typename Scalar>
    ErrorCalculator<Scalar>::~ErrorCalculator()
    {
      for(int i = 0; i < this->component_count; i++)
      {
        if(errors[i])
          ::free(errors[i]);
        if(norms[i])
          ::free(norms[i]);
      }

      if(this->element_references)
        ::free(this->element_references);
    }

    template<typename Scalar>
    MeshFunctionSharedPtr<double> ErrorCalculator<Scalar>::get_errorMeshFunction(int component)
    {
      if(component >= this->component_count)
        throw Exceptions::ValueException("component", component, this->component_count);

      // The value is ready to be returned if it has been initialized and no other error calculation has been
      // performed since.
      if(this->errorMeshFunction[component])
        return this->errorMeshFunction[component];
      else
      {
        this->errorMeshFunction[component].reset(new ExactSolutionConstantArray<double, double>(coarse_solutions[component]->get_mesh(), this->errors[component]));
        return this->errorMeshFunction[component];
      }
    }

    template<typename Scalar>
    void ErrorCalculator<Scalar>::init_data_storage()
    {
      this->num_act_elems = 0;
      errors_squared_sum = 0.;
      norms_squared_sum = 0.;

      // For all components, get the number of active elements,
      // make room for elements, optionally for norms, and incerementaly
      // calculate the total number of elements.
      for(int i = 0; i < this->component_count; i++)
      {
        int num_elements_i = this->coarse_solutions[i]->get_mesh()->get_max_element_id();

        if(errors[i] == NULL)
          errors[i] = (double*)calloc(num_elements_i, sizeof(double));
        else
        {
          errors[i] = (double*)realloc(errors[i], num_elements_i * sizeof(double));
          memset(errors[i], 0, sizeof(double) * num_elements_i);
        }
         
        if(norms[i] == NULL)
          norms[i] = (double*)calloc(num_elements_i, sizeof(double));
        else
        {
          norms[i] = (double*)realloc(norms[i], num_elements_i * sizeof(double));
          memset(norms[i], 0, sizeof(double) * num_elements_i);
        }

        component_errors[i] = 0.;
        component_norms[i] = 0.;
        int num_active_elements_i = this->coarse_solutions[i]->get_mesh()->get_num_active_elements();
        element_count[i] = num_active_elements_i;
        this->num_act_elems += num_active_elements_i;
      }

      // Create the array for references and initialize it.
      if(this->element_references)
        ::free(this->element_references);

      this->element_references = (ElementReference*)malloc(this->num_act_elems * sizeof(ElementReference));

      int running_count_total = 0;
      for(int i = 0; i < this->component_count; i++)
      {
        Element* e;
        for_all_active_elements(e, this->coarse_solutions[i]->get_mesh())
        {
          this->element_references[running_count_total++] = ErrorCalculator<Scalar>::ElementReference(i, e->id, &this->errors[i][e->id], &this->norms[i][e->id]);
        }
      }

      // Also handle the errorMeshFunction.
      for(int i = 0; i < this->component_count; i++)
        if(this->errorMeshFunction[i])
          this->errorMeshFunction[i].reset();
    }

    template<typename Scalar>
    void ErrorCalculator<Scalar>::calculate_errors(MeshFunctionSharedPtr<Scalar>& coarse_solution, MeshFunctionSharedPtr<Scalar>& fine_solution, bool sort_and_store)
    {
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > coarse_solutions;
      coarse_solutions.push_back(coarse_solution);
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > fine_solutions;
      fine_solutions.push_back(fine_solution);
      this->calculate_errors(coarse_solutions, fine_solutions, sort_and_store);
    }

    template<typename Scalar>
    bool ErrorCalculator<Scalar>::isOkay() const
    {
      bool okay = true;

      if(fine_solutions.size() != this->component_count)
      {
        okay = false;
        throw Exceptions::LengthException(0, fine_solutions.size(), this->component_count);
      }

      if(this->mfvol.empty() && this->mfsurf.empty() && this->mfDG.empty())
      {
        okay = false;
        throw Exceptions::Exception("No error forms - instances of NormForm<Scalar> passed to ErrorCalculator via add_error_form().");
      }

      return okay;
    }
    
    template<typename Scalar>
    void ErrorCalculator<Scalar>::calculate_errors(Hermes::vector<MeshFunctionSharedPtr<Scalar> >& coarse_solutions_, Hermes::vector<MeshFunctionSharedPtr<Scalar> >& fine_solutions_, bool sort_and_store)
    {
      this->coarse_solutions = coarse_solutions_;
      this->fine_solutions = fine_solutions_;
      this->component_count = this->coarse_solutions.size();
      
      this->check();

      this->init_data_storage();

      // Prepare multi-mesh traversal and error arrays.
      Hermes::vector<MeshSharedPtr > meshes;

      for (int i = 0; i < this->component_count; i++)
        meshes.push_back(coarse_solutions[i]->get_mesh());
      for (int i = 0; i < this->component_count; i++)
        meshes.push_back(fine_solutions[i]->get_mesh());

      int num_states;
      Traverse trav(this->component_count);
      Traverse::State** states = trav.get_states(meshes, num_states);

#pragma omp parallel num_threads(this->num_threads_used)
      {
        int thread_number = omp_get_thread_num();
        int start = (num_states / this->num_threads_used) * thread_number;
        int end = (num_states / this->num_threads_used) * (thread_number + 1);
        if(thread_number == this->num_threads_used - 1)
          end = num_states;

        // Create a calculator for this thread.
        ErrorThreadCalculator<Scalar> errorThreadCalculator(this);

        // Do the work.
        for(int state_i = start; state_i < end; state_i++)
          errorThreadCalculator.evaluate_one_state(states[state_i]);
      }

      for(int i = 0; i < num_states; i++)
        delete states[i];
      free(states);

      // Clean after ourselves.
      for(int i = 0; i < this->component_count; i++)
      {
        Element* e;
        for_all_active_elements(e, coarse_solutions[i]->get_mesh())
          e->visited = false;
        for_all_active_elements(e, fine_solutions[i]->get_mesh())
          e->visited = false;
      }

      // Sums calculation & error postprocessing.
      this->postprocess_error();

      if(sort_and_store)
      {
        std::qsort(this->element_references, this->num_act_elems, sizeof(ElementReference), &this->compareElementReference);
        elements_stored = true;
      }
      else
        elements_stored = false;
    }

    template<typename Scalar>
    void ErrorCalculator<Scalar>::postprocess_error()
    {
      // Indexer through active elements on all meshes.
      int running_indexer = 0;
      for(int i = 0; i < this->component_count; i++)
      {
        for(int j = 0; j < this->element_count[i]; j++)
        {
          component_errors[i] += *(this->element_references[running_indexer + j].error);
          component_norms[i] += *(this->element_references[running_indexer + j].norm);
          
          if(this->errorType == RelativeErrorToElementNorm)
					   *(this->element_references[running_indexer + j].error) /= *(this->element_references[running_indexer + j].norm);
        }

        if(this->errorType == RelativeErrorToGlobalNorm)
          for(int j = 0; j < this->element_count[i]; j++)
					   *(this->element_references[running_indexer + j].error) /= component_norms[i];

        norms_squared_sum += component_norms[i];
        errors_squared_sum += component_errors[i];
        
        if(this->errorType == RelativeErrorToGlobalNorm || this->errorType == RelativeErrorToElementNorm)
          component_errors[i] /= component_norms[i];

        running_indexer += this->element_count[i];
      }

      if(this->errorType == RelativeErrorToGlobalNorm || this->errorType == RelativeErrorToElementNorm)
        errors_squared_sum /= norms_squared_sum;
    }

    template<typename Scalar>
    void ErrorCalculator<Scalar>::add_error_form(NormFormVol<Scalar>* form)
    {
      this->mfvol.push_back(form);
    }

    template<typename Scalar>
    void ErrorCalculator<Scalar>::add_error_form(NormFormSurf<Scalar>* form)
    {
      this->mfsurf.push_back(form);
    }

    template<typename Scalar>
    void ErrorCalculator<Scalar>::add_error_form(NormFormDG<Scalar>* form)
    {
      this->mfDG.push_back(form);
    }
    
    template<typename Scalar>
    bool ErrorCalculator<Scalar>::data_prepared_for_querying() const
    {
      if(!this->element_references)
      {
        throw Hermes::Exceptions::Exception("Elements / norms have not been calculated so far.");
        return false;
      }

      return true;
    }

    template<typename Scalar>
    const typename ErrorCalculator<Scalar>::ElementReference& ErrorCalculator<Scalar>::get_element_reference(unsigned int id) const
    {
      return this->element_references[id];
    }

    template<typename Scalar>
    double ErrorCalculator<Scalar>::get_element_error_squared(int component, int element_id) const
    {
      if(!this->data_prepared_for_querying())
        return 0.0;
      
      if(component >= this->component_count)
        throw Hermes::Exceptions::ValueException("component", component, this->component_count);
      if(element_id >= element_count[component])
        throw Hermes::Exceptions::ValueException("element_id", element_id, this->element_count[component]);

      return this->errors[component][element_id];
    }

    template<typename Scalar>
    double ErrorCalculator<Scalar>::get_element_norm_squared(int component, int element_id) const
    {
      if(!this->data_prepared_for_querying())
        return 0.0;
      if(component >= this->component_count)
        throw Hermes::Exceptions::ValueException("component", component, this->component_count);
      if(element_id >= element_count[component])
        throw Hermes::Exceptions::ValueException("element_id", element_id, this->element_count[component]);

      return this->norms[component][element_id];
    }

    template<typename Scalar>
    double ErrorCalculator<Scalar>::get_error_squared(int component) const
    {
      if(!this->data_prepared_for_querying())
        return 0.0;
      if(component >= this->component_count)
        throw Hermes::Exceptions::ValueException("component", component, this->component_count);

      return this->component_errors[component];
    }

    template<typename Scalar>
    double ErrorCalculator<Scalar>::get_norm_squared(int component) const
    {
      if(!this->data_prepared_for_querying())
        return 0.0;
      if(component >= this->component_count)
        throw Hermes::Exceptions::ValueException("component", component, this->component_count);

      return this->component_norms[component];
    }

    template<typename Scalar>
    double ErrorCalculator<Scalar>::get_total_error_squared() const
    {
      if(!this->data_prepared_for_querying())
        return 0.0;
      return this->errors_squared_sum;
    }

    template<typename Scalar>
    double ErrorCalculator<Scalar>::get_total_norm_squared() const
    {
      if(!this->data_prepared_for_querying())
        return 0.0;
      return this->norms_squared_sum;
    }

    template<typename Scalar, NormType normType>
    DefaultErrorCalculator<Scalar, normType>::DefaultErrorCalculator(CalculatedErrorType errorType, int component_count) : ErrorCalculator<Scalar>(errorType)
    {
      for(int i = 0; i < component_count; i++)
      {
        this->add_error_form(new DefaultNormFormVol<Scalar>(i, i, normType));
      }
    }

    template<typename Scalar, NormType normType>
    DefaultNormCalculator<Scalar, normType>::DefaultNormCalculator(int component_count) : ErrorCalculator<Scalar>(AbsoluteError)
    {
      for(int i = 0; i < component_count; i++)
      {
        this->add_error_form(new DefaultNormFormVol<Scalar>(i, i, normType));
      }
    }
    
    template<typename Scalar, NormType normType>
    double DefaultNormCalculator<Scalar, normType>::calculate_norms(Hermes::vector<MeshFunctionSharedPtr<Scalar> >& solutions)
    {
      Hermes::vector<MeshFunctionSharedPtr<Scalar> > zero_fine_solutions;
      for(int i = 0; i < solutions.size(); i++)
        zero_fine_solutions.push_back(MeshFunctionSharedPtr<Scalar>(new ZeroSolution<Scalar>(solutions[i]->get_mesh())));
      this->calculate_errors(solutions, zero_fine_solutions, false);

      return this->get_total_error_squared();
    }

    template<typename Scalar, NormType normType>
    double DefaultNormCalculator<Scalar, normType>::calculate_norm(MeshFunctionSharedPtr<Scalar> & solution)
    {
      MeshFunctionSharedPtr<Scalar> zero_fine_solution(new ZeroSolution<Scalar>(solution->get_mesh()));
      this->calculate_errors(solution, zero_fine_solution, false);

      return this->get_total_error_squared();
    }

    template<typename Scalar, NormType normType>
    DefaultErrorCalculator<Scalar, normType>::~DefaultErrorCalculator()
    {
      for(int i = 0; i < this->mfvol.size(); i++)
        delete this->mfvol[i];
    }

    template<typename Scalar, NormType normType>
    DefaultNormCalculator<Scalar, normType>::~DefaultNormCalculator()
    {
      for(int i = 0; i < this->mfvol.size(); i++)
        delete this->mfvol[i];
    }

    template HERMES_API class DefaultErrorCalculator<double, HERMES_H1_NORM>;
    template HERMES_API class DefaultErrorCalculator<std::complex<double>, HERMES_H1_NORM>;
    template HERMES_API class DefaultErrorCalculator<double, HERMES_L2_NORM>;
    template HERMES_API class DefaultErrorCalculator<std::complex<double>, HERMES_L2_NORM>;
    template HERMES_API class DefaultErrorCalculator<double, HERMES_HCURL_NORM>;
    template HERMES_API class DefaultErrorCalculator<std::complex<double>, HERMES_HCURL_NORM>;
    template HERMES_API class DefaultErrorCalculator<double, HERMES_HDIV_NORM>;
    template HERMES_API class DefaultErrorCalculator<std::complex<double>, HERMES_HDIV_NORM>;

    template HERMES_API class DefaultNormCalculator<double, HERMES_H1_NORM>;
    template HERMES_API class DefaultNormCalculator<std::complex<double>, HERMES_H1_NORM>;
    template HERMES_API class DefaultNormCalculator<double, HERMES_L2_NORM>;
    template HERMES_API class DefaultNormCalculator<std::complex<double>, HERMES_L2_NORM>;
    template HERMES_API class DefaultNormCalculator<double, HERMES_HCURL_NORM>;
    template HERMES_API class DefaultNormCalculator<std::complex<double>, HERMES_HCURL_NORM>;
    template HERMES_API class DefaultNormCalculator<double, HERMES_HDIV_NORM>;
    template HERMES_API class DefaultNormCalculator<std::complex<double>, HERMES_HDIV_NORM>;
    
    template HERMES_API class ErrorCalculator<double>;
    template HERMES_API class ErrorCalculator<std::complex<double> >;
  }
}