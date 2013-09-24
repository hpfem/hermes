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

    /// Class for selecting the stopping criterion of the loop over all elements in an descending order wrt. their error.
    /// \ingroup g_adapt
    /// Serves for only inspect such a portion of all elements of all meshes in the system
    /// for potential refinement.
    template<typename Scalar>
    class AdaptivityStoppingCriterion
    {
    public:
      /// Main AdaptivityStoppingCriterion method. It is evaluated in the loop over all elements in an descending order wrt. their error.
      /// The loop ends with the first negative result of a call to this method.
      /// Decide if the refinement at hand will be carried out.
      virtual bool add_refinement(ErrorCalculator<Scalar>* error_calculator, double processed_error_squared, double max_error_squared, int element_inspected_i) = 0;
    };

    /// Stopping criterion based on cumulative processed error.
    /// The method add_refinement will return false as soon as the already processed refinements counted for AdaptStoppingCriterionCumulative::threshold
    /// of the total error.
    template<typename Scalar>
    class HERMES_API AdaptStoppingCriterionCumulative : public AdaptivityStoppingCriterion<Scalar>
    {
    public:
      /// Constructor specifying the threshold (see description of threshold).
      AdaptStoppingCriterionCumulative(double threshold);

      /// Decide if the refinement at hand will be carried out.
      /// Will return false as soon as the already processed refinements counted for AdaptStoppingCriterionCumulative::threshold
      /// of the total error.
      bool add_refinement(ErrorCalculator<Scalar>* error_calculator, double processed_error_squared, double max_error_squared, int element_inspected_i);
    private:
      /// The quantity representing the portion (fraction) of total error that is processed.
      /// See comments above.
      double threshold;
    };
    
    /// Stopping criterion based on maximum element error.
    /// The method add_refinement will return false as soon as the particular element carries lower error than AdaptStoppingCriterionSingleElement::threshold
    /// times the maximum element error.
    template<typename Scalar>
    class HERMES_API AdaptStoppingCriterionSingleElement : public AdaptivityStoppingCriterion<Scalar>
    {
    public:
      /// Constructor specifying the threshold (see description of threshold).
      AdaptStoppingCriterionSingleElement(double threshold);

      /// Decide if the refinement at hand will be carried out.
      bool add_refinement(ErrorCalculator<Scalar>* error_calculator, double processed_error_squared, double max_error_squared, int element_inspected_i);
    private:
      /// The quantity representing the portion (fraction) of maximum error that the processed elements contain.
      /// See comments above.
      double threshold;
    };
    
    /// Stopping criterion based on refining elements with similar errors.
    /// The method add_refinement will return false as soon as the particular element carries significantly less error than the previous one in the descending sequence.
    /// Useful e.g. when we are more interested in overall solution quality than resolution of a steep singularity etc.
    template<typename Scalar>
    class HERMES_API AdaptStoppingCriterionLevels : public AdaptivityStoppingCriterion<Scalar>
    {
    public:
      /// Constructor specifying the threshold (see description of threshold).
      AdaptStoppingCriterionLevels(double threshold);

      /// Decide if the refinement at hand will be carried out.
      bool add_refinement(ErrorCalculator<Scalar>* error_calculator, double processed_error_squared, double max_error_squared, int element_inspected_i);
    private:
      /// The quantity representing the portion (fraction) of the current element error to the previous one for those that will be refined..
      /// See comments above.
      double threshold;
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
      public Hermes::Hermes2D::Mixins::StateQueryable, 
      public Hermes::Hermes2D::Mixins::Parallel
    {
    public:
      /// Constructor. Suitable for problems where various solution components belong to different spaces (L2, H1, Hcurl,
      /// Hdiv). If proj_norms are not specified, they are defined according to the spaces.
      Adapt(Hermes::vector<SpaceSharedPtr<Scalar> > spaces, ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* strategy = NULL);
			Adapt(SpaceSharedPtr<Scalar> space, ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* strategy = NULL);
			Adapt(ErrorCalculator<Scalar>* error_calculator, AdaptivityStoppingCriterion<Scalar>* strategy = NULL);
      virtual ~Adapt();  ///< Destructor. Deallocates allocated private data.

      /// Refines elements based on results from the ErrorCalculator class.
      /**
      *  \param[in] refinement_selectors Vector of selectors.
      *  \return True if no element was refined. In usual case, this indicates that adaptivity is not able to refine anything and the adaptivity loop should end. */
      bool adapt(Hermes::vector<RefinementSelectors::Selector<Scalar>*> refinement_selectors);

      /// Refines elements based on results from the ErrorCalculator class.
      /**
      *  \param[in] refinement_selector A pointer to a selector which will select a refinement.
      *  \return True if no element was refined. In usual case, this indicates that adaptivity is not able to refine anything and the adaptivity loop should end. */
      bool adapt(RefinementSelectors::Selector<Scalar>* refinement_selector);

      /// Set the current strategy.
      /// \param[in] strategy The strategy, see the info for AdaptivityStoppingCriterion.
      void set_strategy(AdaptivityStoppingCriterion<Scalar>* strategy);

      /// Set the regularization level.
      /// See attribute regularization.
      void set_regularization_level(int regularization);

      /// Internal checking.
      virtual bool isOkay() const;
      inline std::string getClassName() const { return "Adapt"; }

      /// Set spaces.
      void set_spaces(Hermes::vector<SpaceSharedPtr<Scalar> > spaces);
      void set_space(SpaceSharedPtr<Scalar> space);

      /// Return the error mesh function - for postprocessing the information about which elements have been refined.
      /// \param component The component.
      MeshFunctionSharedPtr<double> get_refinementInfoMeshFunction(int component = 0);
    protected:
      /// Set default values.
      void set_defaults();

      /// Initialization.
      void init_adapt(Hermes::vector<RefinementSelectors::Selector<Scalar>*>& refinement_selectors, ElementToRefine*** element_refinement_location, MeshSharedPtr* meshes);
      /// Return the number of element where a refinement will be sought.
      int calculate_attempted_element_refinements_count();
      /// Handle meshes and spaces at the end of the routine.
      void adapt_postprocess(MeshSharedPtr* meshes, int element_refinements_count);
      /// Deinitialization.
      void deinit_adapt(ElementToRefine*** element_refinement_location);

      /// Common code for the constructors.
      void init();

      /// Current strategy.
      AdaptivityStoppingCriterion<Scalar>* strategy;

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
      void fix_shared_mesh_refinements(MeshSharedPtr * meshes, ElementToRefine*& elems_to_refine, int& num_elem_to_process, ElementToRefine*** refinement_location, RefinementSelectors::Selector<Scalar>** refinement_selectors);

      /// Enforces the same order to an element of a mesh which is shared among multiple components.
      /** \param[in] meshes An array of meshes of components. */
      void homogenize_shared_mesh_orders(MeshSharedPtr * meshes);

      /// Number of solution components (as in wf->neq).
      int num;

      /// Regularization (max. level of hanging nodes) level.
      int regularization;

      /// Spaces.
      Hermes::vector<SpaceSharedPtr<Scalar> > spaces;

      /// Error calculator.
      ErrorCalculator<Scalar>* errorCalculator;

      /// Information about performed refinements.
      ElementToRefine* elements_to_refine;
      int elements_to_refine_count;
      
      /// Mesh function for postprocessing the information about which elements have been refined.
      MeshFunctionSharedPtr<double> refinementInfoMeshFunction[H2D_MAX_COMPONENTS];
    };
  }
}
#endif