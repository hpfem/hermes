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

#include "../forms.h"
#include "../space/space.h"
#include "../weakform/weakform.h"
#include "../integrals/h1.h"
#include "../integrals/hcurl.h"
#include "../integrals/hdiv.h"
#include "../integrals/l2.h"
#include "../mesh/element_to_refine.h"
#include "../refinement_selectors/selector.h"
#include "exceptions.h"
#include "../global.h"

namespace Hermes
{
  namespace Hermes2D
  {
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

    template<typename Scalar> class Global;

    /// Evaluation of an error between a (coarse) solution and a reference solution and adaptivity. \ingroup g_adapt
    /** The class provides basic functionality necessary to adaptively refine elements.
    *  Given a reference solution and a coarse solution, it calculates error estimates
    *  and it acts as a container for the calculated errors.
    */
    template<typename Scalar>
    class HERMES_API Adapt : public Hermes::Mixins::TimeMeasurable, public Hermes::Mixins::Loggable
    {
    public:
      /// Constructor. Suitable for problems where various solution components belong to different spaces (L2, H1, Hcurl,
      /// Hdiv). If proj_norms are not specified, they are defined according to the spaces.
      Adapt(const Hermes::vector<Space<Scalar>*>& spaces, Hermes::vector<ProjNormType> proj_norms = Hermes::vector<ProjNormType>());
      Adapt(Space<Scalar>* space, ProjNormType proj_norm = HERMES_UNSET_NORM);
      virtual ~Adapt();  ///< Destructor. Deallocates allocated private data.

      /// Matrix forms for error calculation.
      class HERMES_API MatrixFormVolError : public MatrixFormVol<Scalar>
      {
      public:
        /// Empty constructor.
        MatrixFormVolError(int i, int j);
        /// Constructor that takes the norm identification.
        MatrixFormVolError(int i, int j, ProjNormType type);

        /// Error bilinear form callback function.
        virtual Scalar value(int n, double *wt, Func<Scalar> *u_ext[],
          Func<Scalar> *u, Func<Scalar> *v, Geom<double> *e,
          Func<Scalar> **ext) const;

        /// Error bilinear form to estimate order of a function.
        virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[],
          Func<Hermes::Ord> *u, Func<Hermes::Ord> *v, Geom<Hermes::Ord> *e,
          Func<Ord> **ext) const;

        virtual MatrixFormVol<Scalar>* clone() const;

      protected:
        /// Norm used.
        ProjNormType projNormType;

        /// L2 error form.
        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain l2_error_form(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<SolFunctionDomain> *u,
          Func<SolFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext);

        /// H1 error form.
        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain h1_error_form(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<SolFunctionDomain> *u,
          Func<SolFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext);

        /// H1-seminorm error form.
        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain h1_error_semi_form(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<SolFunctionDomain> *u,
          Func<SolFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext);

        /// H-div error form.
        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain hdiv_error_form(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<SolFunctionDomain> *u,
          Func<SolFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext);

        /// H-curl error form.
        template<typename TestFunctionDomain, typename SolFunctionDomain>
        static SolFunctionDomain hcurl_error_form(int n, double *wt, Func<SolFunctionDomain> *u_ext[], Func<SolFunctionDomain> *u,
          Func<SolFunctionDomain> *v, Geom<TestFunctionDomain> *e, Func<SolFunctionDomain> **ext);
      };

      /// Sets user defined bilinear form which is used to calculate error.
      /** By default, all inherited class should set default bilinear forms for each element (i.e. i = j).
      *  This method can change it or set forms that can combine components (e.g. energy error).
      *  \param[in] i The first component index.
      *  \param[in] j The second component index.
      *  \param[in] bi_form A bilinear form which calculates value.
      *  \param[in] bi_ord A bilinear form which calculates order. */
      void set_error_form(int i, int j, MatrixFormVolError* form);
      void set_error_form(MatrixFormVolError* form);   ///< i = j = 0

      void set_norm_form(int i, int j, MatrixFormVolError* form);
      void set_norm_form(MatrixFormVolError* form);   ///< i = j = 0

      /// Type-safe version of calc_err_est() for one solution.
      /// @param[in] solutions_for_adapt - if sln and rsln are the solutions error of which is used in the function adapt().
      double calc_err_est(Solution<Scalar>*sln, Solution<Scalar>*rsln, bool solutions_for_adapt = true,
        unsigned int error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL);

      /// Calculates the error of the solution. 'n' must be the same
      /// as 'num' in the constructor. After that, n coarse solution
      /// pointers are passed, followed by n fine solution pointers.
      /// @param[in] solutions_for_adapt - if slns and rslns are the solutions error of which is used in the function adapt().
      double calc_err_est(Hermes::vector<Solution<Scalar>*> slns, Hermes::vector<Solution<Scalar>*> rslns,
        Hermes::vector<double>* component_errors = NULL, bool solutions_for_adapt = true,
        unsigned int error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL);

      /// Type-safe version of calc_err_exact() for one solution.
      /// @param[in] solutions_for_adapt - if sln and rsln are the solutions error of which is used in the function adapt().
      double calc_err_exact(Solution<Scalar>*sln, Solution<Scalar>*rsln, bool solutions_for_adapt = true,
        unsigned int error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL);

      /// Calculates the error of the solution.
      /// @param[in] solutions_for_adapt - if slns and rslns are the solutions error of which is used in the function adapt().
      double calc_err_exact(Hermes::vector<Solution<Scalar>*> slns, Hermes::vector<Solution<Scalar>*> rslns,
        Hermes::vector<double>* component_errors = NULL, bool solutions_for_adapt = true,
        unsigned int error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL);

      /// Refines elements based on results from calc_err_est().
      /**
      *  \param[in] refinement_selectors Vector of selectors.
      *  \param[in] thr A threshold. The meaning of the threshold is defined by the parameter strat.
      *  \param[in] strat A strategy. It specifies a stop condition which quits processing elements in the Adapt::regular_queue. Possible values are 0, 1, 2, and 3.
      *  \param[in] regularize Regularizing of a mesh.
      *  \param[in] same_order True if all element have to have same orders after all refinements are applied.
      *  \param[in] to_be_processed Error which has to be processed. Used in strategy number 3.
      *  \return True if no element was refined. In usual case, this indicates that adaptivity is not able to refine anything and the adaptivity loop should end. */
      bool adapt(Hermes::vector<RefinementSelectors::Selector<Scalar>*> refinement_selectors, double thr, int strat = 0,
        int regularize = -1, double to_be_processed = 0.0);

      /// Refines elements based on results from calc_err_est().
      /**
      *  \param[in] refinement_selector A pointer to a selector which will select a refinement.
      *  \param[in] thr A threshold. The meaning of the threshold is defined by the parameter strat.
      *  \param[in] strat A strategy. It specifies a stop condition which quits processing elements in the Adapt::regular_queue. Possible values are 0, 1, 2, and 3.
      *  \param[in] regularize Regularizing of a mesh.
      *  \param[in] same_order True if all element have to have same orders after all refinements are applied.
      *  \param[in] to_be_processed Error which has to be processed. Used in strategy number 3.
      *  \return True if no element was refined. In usual case, this indicates that adaptivity is not able to refine anything and the adaptivity loop should end. */
      bool adapt(RefinementSelectors::Selector<Scalar>* refinement_selector, double thr, int strat = 0,
        int regularize = -1, double to_be_processed = 0.0);

      /// Returns a squared error of an element.
      /** \param[in] A component index.
      *  \param[in] An element index.
      *  \return Squared error. Meaning of the error depends on parameters of the function calc_errors_internal(). */
      double get_element_error_squared(int component, int id) const;
    protected:

      Exceptions::Exception* caughtException;

      /// A reference to an element.
      struct ElementReference {
        int id; ///< An element ID. Invalid if below 0.
        int comp; ///< A component which this element belongs to. Invalid if below 0.
        ElementReference(int id = -1, int comp = -1) : id(id), comp(comp) {}; ///< Constructor. It creates an invalid element reference.
      };

      /// Returns regular queue of elements
      /** \return A regular queue. */
      const Hermes::vector<ElementReference>& get_regular_queue() const;

      /// Apply a single refinement.
      /** \param[in] A refinement to apply. */
      void apply_refinement(const ElementToRefine& elem_ref);

      /// Apply a vector of refinements.
      /** \param[in] A vector of refinements to apply. */
      virtual void apply_refinements(std::vector<ElementToRefine>& elems_to_refine);

      /// Returns a vector of refinements generated during the last execution of the method adapt().
      /** \return A vector of refinements generated during the last execution of the method adapt(). The returned vector might change or become invalid after the next execution of the method adadpt(). */
      const std::vector<ElementToRefine>& get_last_refinements() const; ///< Returns last refinements.

      std::queue<ElementReference> priority_queue; ///< A queue of priority elements. Elements in this queue are processed before the elements in the Adapt::regular_queue.
      Hermes::vector<ElementReference> regular_queue; ///< A queue of elements which should be processes. The queue had to be filled by the method fill_regular_queue().
      std::vector<ElementToRefine> last_refinements; ///< A vector of refinements generated during the last finished execution of the method adapt().

      /// Fixes refinements of a mesh which is shared among multiple components of a multimesh.
      /** If a mesh is shared among components, it has to be refined similarly in order to avoid inconsistency.
      *  \param[in] meshes An array of meshes of components.
      *  \param[in] elems_to_refine A vector of refinements.
      *  \param[in] idx A 2D array that translates a pair (a component index, an element id) to an index of a refinement in the vector of refinements. If the index is below zero, a given element was not refined.
      *  \param[in] refinement_selector A selected used by the adaptivity. The selector is used to correct orders of modified refinements using RefinementSelectors::Selector::update_shared_mesh_orders(). */
      void fix_shared_mesh_refinements(Mesh** meshes, std::vector<ElementToRefine>& elems_to_refine, int** idx,
        RefinementSelectors::Selector<Scalar>*** refinement_selectors);

      /// Enforces the same order to an element of a mesh which is shared among multiple components.
      /** \param[in] meshes An array of meshes of components. */
      void homogenize_shared_mesh_orders(Mesh** meshes);

      int num;                              ///< Number of solution components (as in wf->neq).
      Hermes::vector<Space<Scalar>*> spaces;        ///< Spaces.
      bool **own_forms;
      int num_act_elems;                    ///< A total number of active elements across all provided meshes.
      Solution<Scalar>* sln[H2D_MAX_COMPONENTS];    ///< Coarse solution.
      Solution<Scalar>* rsln[H2D_MAX_COMPONENTS];   ///< Reference solutions.
      bool have_errors;                     ///< True if errors of elements were calculated.
      bool have_coarse_solutions;           ///< True if the coarse solutions were set.
      bool have_reference_solutions;        ///< True if the reference solutions were set.

      double* errors[H2D_MAX_COMPONENTS];   ///< Errors of elements. Meaning of the error depeds on flags used when the
      ///< method calc_errors_internal() was calls. Initialized in the method calc_errors_internal().
      double  errors_squared_sum;           ///< Sum of errors in the array Adapt::errors_squared. Used by a method adapt() in some strategies.

      double error_time;                    ///< Time needed to calculate the error.

      static const unsigned char HERMES_TOTAL_ERROR_MASK = 0x0F;    ///< A mask which masks-out total error type. Used by Adapt::calc_err_internal(). \internal
      static const unsigned char HERMES_ELEMENT_ERROR_MASK = 0xF0;  ///< A mask which masks-out element error type. Used by Adapt::calc_err_internal(). \internal

      MatrixFormVolError* error_form[H2D_MAX_COMPONENTS][H2D_MAX_COMPONENTS]; ///< Bilinear forms to calculate error.
      MatrixFormVolError* norm_form[H2D_MAX_COMPONENTS][H2D_MAX_COMPONENTS];  ///< Bilinear forms to calculate norm (set to error_form by default).

      /// Calculates error between a coarse solution and a reference solution and sorts components according to the error.
      /** If overridden, this method has to initialize errors (Array::errors), sum of errors (Array::error_sum), norms of components (Array::norm), number of active elements (Array::num_act_elems). Also, it has to fill the regular queue through the method fill_regular_queue().
      *  \param[in] error_flags Flags which calculates the error. It can be a combination of ::HERMES_TOTAL_ERROR_REL, ::HERMES_TOTAL_ERROR_ABS, ::HERMES_ELEMENT_ERROR_REL, ::HERMES_ELEMENT_ERROR_ABS.
      *  \return The total error. Interpretation of the error is specified by the parameter error_flags. */
      virtual double calc_err_internal(Hermes::vector<Solution<Scalar>*> slns, Hermes::vector<Solution<Scalar>*> rslns,
        Hermes::vector<double>* component_errors, bool solutions_for_adapt, unsigned int error_flags);

      /// One Space version.
      virtual double calc_err_internal(Solution<Scalar>* sln, Solution<Scalar>* rsln,
        Hermes::vector<double>* component_errors, bool solutions_for_adapt, unsigned int error_flags);

      /// Evaluates a square of an absolute error of an active element among a given pair of components.
      /** The method uses a bilinear forms to calculate the error. This is done by supplying a differences (f1 - v1) and (f2 - v2) at integration points to the bilinear form,
      *  where f1 and f2 are values of (coarse) solutions of the first and the second component respectively,
      *  v1 and v2 are values of reference solutions of the first and the second component respectively.
      *  \param[in] form A bilinear form.
      *  \param[in] sln1 A (coarse) solution corresponding to the first component.
      *  \param[in] sln2 A (coarse) solution corresponding to the second component.
      *  \param[in] rsln1 A reference solution corresponding to the first component.
      *  \param[in] rsln2 A reference solution corresponding to the second component.
      *  \param[in] rv1 A reference map of a (coarse) solution sln1.
      *  \param[in] rv2 A reference map of a (coarse) solution sln2.
      *  \param[in] rrv1 A reference map of a reference solution rsln1.
      *  \param[in] rrv2 A reference map of a reference solution rsln2.
      *  \return A square of an absolute error. */
      virtual double eval_error(MatrixFormVolError* form,
        MeshFunction<Scalar>*sln1, MeshFunction<Scalar>*sln2, MeshFunction<Scalar>*rsln1, MeshFunction<Scalar>*rsln2);

      /// Evaluates the square of a norm of an active element in the reference solution among a given pair of components.
      /** The method uses a bilinear forms to calculate the norm. This is done by supplying a v1 and v2 at integration points to the bilinear form,
      *  where v1 and v2 are values of reference solutions of the first and the second component respectively.
      *  \param[in] form A bilinear form.
      *  \param[in] rsln1 A reference solution corresponding to the first component.
      *  \param[in] rsln2 A reference solution corresponding to the second component.
      *  \param[in] rrv1 A reference map of a reference solution rsln1.
      *  \param[in] rrv2 A reference map of a reference solution rsln2.
      *  \return A square of a norm. */
      virtual double eval_error_norm(MatrixFormVolError* form,
        MeshFunction<Scalar>*rsln1, MeshFunction<Scalar>*rsln2);

      /// Builds an ordered queue of elements that are be examined.
      /** The method fills Adapt::standard_queue by elements sorted accordin to their error descending.
      *  The method assumes that Adapt::errors_squared contains valid values.
      *  If a special order of elements is requested, this method has to be overridden.
      *  /param[in] meshes An array of pointers to meshes of a (coarse) solution. An index into the array is an index of a component.
      *  /param[in] meshes An array of pointers to meshes of a reference solution. An index into the array is an index of a component. */
      virtual void fill_regular_queue(const Mesh** meshes);

    private:
      /// A functor that compares elements accoring to their error. Used by std::sort().
      class CompareElements
      {
      private:
        double** errors; ///< A 2D array of squared errors: the first index is an index of component, the second index is an element ID.
      public:
        CompareElements(double** errors); ///< Constructor.
        /// Compares two elements.
        /** \param[in] e1 A reference to the first element.
        *  \param[in] e1 A reference to the second element.
        *  \return True if a squared error of the first element is greater than a squared error of the second element. */
        bool operator ()(const ElementReference& e1, const ElementReference& e2) const;
      };
    };
  }
}
#endif