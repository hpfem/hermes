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

#ifndef __H2D_ADAPT_H
#define __H2D_ADAPT_H

#include "../space/space.h"
#include "error_calculator.h"
#include "../refinement_selectors/selector.h"
#include "exceptions.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar> class ErrorCalculator;

    /** \defgroup g_adapt Adaptivity
    *  \brief Adaptivity provides framework for modifying elements in order to decrease errors of the solution.
    *
    *  Adaptivity classes calculate error of every element.
    *  An error of an element is calculated either by comparing a
    *  coarse solution with a reference solution or by evaluating
    *  suitable error estimate. Errors of
    *  elements define the order in which elements are examined.
    *  During examining an element, a refinement is proposed
    *  and the element is refined if applicable. The refinement
    *  is proposed through refinement selectors, see \ref g_selectors.
    *
    */

    /// Enum for selecting the stopping criterion of the loop over all elements.
    /// \ingroup g_adapt
    /// Serves for only inspect such a portion of all elements of all meshes in the system
    /// for potential refinement.
    enum AdaptivityStoppingCriterion
    {
      /// Refine elements until prescribed portion of error is processed.
      AdaptStoppingCriterionCumulative,
      AdaptStoppingCriterionSingleElement,
      AdaptStoppingCriterionLevels
    };

    /// Evaluation of an error between a (coarse) solution and a reference solution and adaptivity.
    /// \ingroup g_adapt
    /** The class provides basic functionality necessary to adaptively refine elements.
    *  Given a reference solution and a coarse solution, it calculates error estimates
    *  and it acts as a container for the calculated errors.
    *  For default values of attributes, see the method set_defaults().
    */
    template<typename Scalar>
    class HERMES_API Adapt :
      public Hermes::Mixins::TimeMeasurable,
      public Hermes::Mixins::Loggable,
      public Hermes::Hermes2D::Mixins::Parallel
    {
    public:
      /// Constructor. Suitable for problems where various solution components belong to different spaces (L2, H1, Hcurl,
      /// Hdiv). If proj_norms are not specified, they are defined according to the spaces.
      Adapt(Hermes::vector<SpaceSharedPtr<Scalar> > spaces, ErrorCalculator<Scalar>* error_calculator);
			Adapt(SpaceSharedPtr<Scalar> space, ErrorCalculator<Scalar>* error_calculator);
      virtual ~Adapt();  ///< Destructor. Deallocates allocated private data.

      /// Refines elements based on results from the ErrorCalculator class.
      /**
      *  \param[in] refinement_selectors Vector of selectors.
      *  \return True if no element was refined. In usual case, this indicates that adaptivity is not able to refine anything and the adaptivity loop should end. */
      bool adapt(Hermes::vector<RefinementSelectors::Selector<Scalar>*>& refinement_selectors);

      /// Refines elements based on results from the ErrorCalculator class.
      /**
      *  \param[in] refinement_selector A pointer to a selector which will select a refinement.
      *  \return True if no element was refined. In usual case, this indicates that adaptivity is not able to refine anything and the adaptivity loop should end. */
      bool adapt(RefinementSelectors::Selector<Scalar>* refinement_selector);

      /// Set the current strategy.
      /// \param[in] strategy The strategy, see the info for AdaptivityStoppingCriterion enum.
      /// \param[in] threshold The number representing a threshold in a meaning specific to the strategy.
      /// Default strategy : AdaptStoppingCriterionCumulative.
      /// Default threshold: 0.3.
      void set_strategy(AdaptivityStoppingCriterion strategy, double threshold = 0.3);

      /// (Experimental)
      /// Iterative improvement - at the end of the adaptation algorithm, the resulting space is used for error calculation between the
      /// Reference solution and its coarse component and if the previous error divided by error after refinements are applied 
      /// is not below iterative_improvement_factor, another adaptivity step is performed.
      /// \param[in] iterative_improvement_factor The factor.
      void set_iterative_improvement(double iterative_improvement_factor = 1e-1);

      /// Set the regularization level.
      /// See attribute regularization.
      void set_regularization_level(int regularization);

    protected:
      /// Set default values.
      void set_defaults();

      /// (Experimental)
      /// Iterative improvement - at the end of the adaptation algorithm, the resulting space is used for error calculation between the
      /// Reference solution and its coarse component and if the previous error divided by error after refinements are applied 
      /// is not below iterative_improvement_factor, another adaptivity step is performed.
      bool iterative_improvement;
      /// Prescribed factor for the iterative improvement.
      double iterative_improvement_factor;
      /// Internal.
      int iterative_improvement_iteration;

      /// Initialization.
      void init_adapt(Hermes::vector<RefinementSelectors::Selector<Scalar>*>& refinement_selectors, int** element_refinement_location, MeshSharedPtr* meshes);
      /// Return the number of element where a refinement will be sought.
      int calculate_attempted_element_refinements_count();
      /// Handle meshes and spaces at the end of the routine.
      void adapt_postprocess(MeshSharedPtr* meshes, int element_refinements_count);
      /// Deinitialization.
      void deinit_adapt(int** element_refinement_location);

      /// Common code for the constructors.
      void init();

      /// Current strategy.
      AdaptivityStoppingCriterion strategy;

      /// Current threshold.
      double threshold;

      /// Decide if the refinement at hand will be carried out.
      bool add_refinement(double processed_error_squared, double max_error_squared, int element_inspected_i);

      /// Apply a single refinement.
      /** \param[in] A refinement to apply. */
      void apply_refinement(const ElementToRefine& elem_ref);

      /// Apply a vector of refinements.
      /** \param[in] A vector of refinements to apply. */
      virtual void apply_refinements(ElementToRefine* elems_to_refine, int num_elem_to_process);

      /// Fixes refinements of a mesh which is shared among multiple components of a multimesh.
      /** If a mesh is shared among components, it has to be refined similarly in order to avoid inconsistency.
      *  \param[in] meshes An array of meshes of components.
      *  \param[in] elems_to_refine An array of refinements.
      *  \param[in] num_elem_to_process Length of the array elems_to_refine.
      *  \param[in] idx A 2D array that translates a pair (a component index, an element id) to an index of a refinement in the vector of refinements. If the index is below zero, a given element was not refined.
      *  \param[in] refinement_selectors Selectors used by the adaptivity. The selector is used to correct orders of modified refinements using RefinementSelectors::Selector::update_shared_mesh_orders(). */
      void fix_shared_mesh_refinements(MeshSharedPtr * meshes, ElementToRefine*& elems_to_refine, int num_elem_to_process, int** idx, RefinementSelectors::Selector<Scalar>** refinement_selectors);

      /// Enforces the same order to an element of a mesh which is shared among multiple components.
      /** \param[in] meshes An array of meshes of components. */
      void homogenize_shared_mesh_orders(MeshSharedPtr * meshes);

      /// Number of solution components (as in wf->neq).
      int num;

      /// Total error from the error calculator.
      /// Stored because of iterative_improvement.
      double sum_error_squared;

      /// Regularization (max. level of hanging nodes) level.
      int regularization;

      /// Spaces.
      Hermes::vector<SpaceSharedPtr<Scalar> > spaces;

      /// Error calculator.
      ErrorCalculator<Scalar>* errorCalculator;

      /// Internal.
      std::exception* caughtException;
    };
  }
}
#endif