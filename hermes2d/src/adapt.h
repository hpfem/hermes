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

#include "norm.h"
#include "forms.h"
#include "space.h"
#include "tuple.h"
#include "weakform.h"
#include "integrals_h1.h"
#include "integrals_hcurl.h"
#include "integrals_hdiv.h"
#include "ref_selectors/selector.h"

/** \defgroup g_adapt Adaptivity
 *  \brief Adaptivity provides framework for modyfying elements in order to decrease errors of the solution.
 *
 *  Adaptivity classes calculates error of every element.
 *  An error of an element is calculated by comparing an
 *  coarse solution with a reference solution. Errors of
 *  elements defines an order in which elements are examined.
 *  During examining an element, a refinement is proposed
 *  and the element is refined if applicable. The refinement
 *  is proposed through refinement selectors, see \ref g_selectors.
 *
 *  All adaptivity classes have to be derived from the class Adapt.
 *  Currently available classes are:
 *  - H1Adapt
 *  - L2Adapt
 *    \if H2D_COMPLEX # -HcurlAdapt \endif
 */

#define H2D_MAX_COMPONENTS 10 ///< A maximum number of components.

H2D_API_USED_TEMPLATE(Tuple<Space*>); ///< Instantiated template. It is used to create a clean Windows DLL interface.
H2D_API_USED_TEMPLATE(Tuple<Solution*>); ///< Instantiated template. It is used to create a clean Windows DLL interface.

// Constant used by Adapt::calc_eror().
#define H2D_TOTAL_ERROR_REL  0x00  ///< A flag which defines interpretation of the total error. \ingroup g_adapt
                                   ///  The total error is divided by the norm and therefore it should be in a range [0, 1].
                                   ///  \note Used by Adapt::calc_elem_errors().. This flag is mutually exclusive with ::H2D_TOTAL_ERROR_ABS.
#define H2D_TOTAL_ERROR_ABS  0x01  ///< A flag which defines interpretation of the total error. \ingroup g_adapt
                                   ///  The total error is absolute, i.e., it is an integral over squares of differencies.
                                   ///  \note Used by Adapt::calc_elem_errors(). This flag is mutually exclusive with ::H2D_TOTAL_ERROR_REL.
#define H2D_ELEMENT_ERROR_REL 0x00 ///< A flag which defines interpretation of an error of an element. \ingroup g_adapt
                                   ///  An error of an element is a square of an error divided by a square of a norm of a corresponding component.
                                   ///  When norms of 2 components are very different (e.g. microwave heating), it can help.
                                   ///  Navier-stokes on different meshes work only when absolute error (see ::H2D_ELEMENT_ERROR_ABS) is used.
                                   ///  \note Used by Adapt::calc_elem_errors(). This flag is mutually exclusive with ::H2D_ELEMENT_ERROR_ABS.
#define H2D_ELEMENT_ERROR_ABS 0x10 ///< A flag which defines interpretation of of an error of an element. \ingroup g_adapt
                                   ///  An error of an element is a square of an asolute error, i.e., it is an integral over squares of differencies.
                                   ///  \note Used by Adapt::calc_elem_errors(). This flag is mutually exclusive with ::H2D_ELEMENT_ERROR_REL.

// Matrix forms for error calculation.
  typedef scalar (*matrix_form_val_t) (int n, double *wt, Func<scalar> *u_ext[], 
                                       Func<scalar> *u, Func<scalar> *v, Geom<double> *e, 
                                       ExtData<scalar> *); ///< A bilinear form callback function.
  typedef Ord (*matrix_form_ord_t) (int n, double *wt, Func<Ord> *u_ext[], 
                                    Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, 
                                    ExtData<Ord> *); ///< A bilinear form to estimate order of a function.

///< Structure to hold adaptivity parameters together.
struct AdaptivityParamType {
  double err_stop; 
  int ndof_stop;
  double threshold; 
  int strategy;
  int mesh_regularity;
  double to_be_processed;
  int total_error_flag;
  int elem_error_flag;

  Tuple<int> error_form_i;
  Tuple<int> error_form_j;
  Tuple<matrix_form_val_t> error_form_val;
  Tuple<matrix_form_ord_t> error_form_ord;

  AdaptivityParamType(double err_stop = 1.0, int ndof_stop = 50000,
	  	      double threshold = 0.3, int strategy = 0, 
                      int mesh_regularity = -1, double to_be_processed = 0.0,
                      int total_error_flag = H2D_TOTAL_ERROR_REL,
                      int elem_error_flag = H2D_ELEMENT_ERROR_REL)
  {
    this->err_stop = err_stop;
    this->ndof_stop = ndof_stop;
    this->threshold = threshold;
    this->strategy = strategy;
    this->mesh_regularity = mesh_regularity;
    this->to_be_processed = to_be_processed;
    this->total_error_flag = total_error_flag;
    this->elem_error_flag = elem_error_flag;
    error_form_i = Tuple<int>();
    error_form_j = Tuple<int>();
    error_form_val = Tuple<matrix_form_val_t>();
    error_form_ord = Tuple<matrix_form_ord_t>();
  }; 
  
  void set_error_form(int i, int j, matrix_form_val_t mfv, matrix_form_ord_t mfo) 
  {
    if (error_form_val.size() > 100) error("too many error forms in AdaptivityParamType::add_error_form().");
    this->error_form_i.push_back(i);
    this->error_form_j.push_back(j);
    this->error_form_val.push_back(mfv);
    this->error_form_ord.push_back(mfo);
  }
};

/// Evaluation of an error between a (coarse) solution and a refernece solution and adaptivity. \ingroup g_adapt
/** The class provides basic functionality necessary to adaptively refine elements.
 *  Given a reference solution and a coarse solution, it calculates error estimates
 *  and it acts as a container for the calculated errors.
 *  The class has to be inherited in order to be used.
 */
class H2D_API Adapt
{
public:
  /// Constructor. Suitable for problems where various solution components belong to different spaces (L2, H1, Hcurl, 
  /// Hdiv). If proj_norms are not specified, they are expected to be set later by set_error_form.
  Adapt(Tuple<Space *> spaces_, Tuple<int> proj_norms = Tuple<int>()); 
  virtual ~Adapt();  ///< Destructor. Deallocates allocated private data.

  /// Sets user defined bilinear form which is used to calculate error.
  /** By default, all inherited class should set default bilinear forms for each element (i.e. i = j).
   *  This method can change it or set forms that can combine components (e.g. energy error).
   *  \param[in] i The first component index.
   *  \param[in] j The second component index.
   *  \param[in] bi_form A bilinear form which calculates value.
   *  \param[in] bi_ord A bilinear form which calculates order. */
  void set_error_form(int i, int j, matrix_form_val_t bi_form, matrix_form_ord_t bi_ord);
  void set_error_form(matrix_form_val_t bi_form, matrix_form_ord_t bi_ord);   // i = j = 0

  /// Sets solutions and reference solutions.
  /** \param[in] solutions Coarse solutions. The number of solutions has to match a number of components.
   *  \param[in] ref_solutions Reference solutions. The number of reference solutions has to match a number of components. */
  void set_solutions(Tuple<Solution*> solutions, Tuple<Solution*> ref_solutions);
  void set_solutions(Solution* solution, Solution* ref_solution) 
  {
    set_solutions(Tuple<Solution*>(solution), Tuple<Solution*>(ref_solution));
  }

  /// Calculates error between a coarse solution and a reference solution and sorts components according to the error.
  /** If overrided, this method has to initialize errors (Array::errors), sum of errors (Array::error_sum), norms of components (Array::norm), number of active elements (Array::num_act_elems). Also, it has to fill the regular queue through the method fill_regular_queue().
   *  \param[in] error_flags Flags which calculates the error. It can be a combination of ::H2D_TOTAL_ERROR_REL, ::H2D_TOTAL_ERROR_ABS, ::H2D_ELEMENT_ERROR_REL, ::H2D_ELEMENT_ERROR_ABS.
   *  \return The total error. Interpretation of the error is specified by the parameter error_flags. */
  virtual double calc_elem_errors(unsigned int error_flags = H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_ABS);

  /// Refines elements based on results from calc_elem_errors().
  /** The behavior of adaptivity can be controlled through methods should_ignore_element()
   *  and can_refine_element() which are inteteded to be overriden if neccessary.
   *  \param[in] refinement_selector A point to a selector which will select a refinement.
   *  \param[in] thr A threshold. The meaning of the threshold is defined by the parameter strat.
   *  \param[in] strat A strategy. It specifies a stop condition which quits processing elements in the Adapt::regular_queue. Possible values are 0, 1, 2, and 3.
   *  \param[in] regularize Regularizing of a mesh.
   *  \param[in] same_order True if all element have to have same orders after all refinements are applied.
   *  \param[in] to_be_processed Error which has to be processed. Used in strategy number 3.
   *  \return True if no element was refined. In usual case, this indicates that adaptivity is not able to refine anything and the adaptivity loop should end. */
  bool adapt(Tuple<RefinementSelectors::Selector *> refinement_selectors, double thr, int strat = 0,
             int regularize = -1, double to_be_processed = 0.0);

  /// Unrefines the elements with the smallest error.
  /** \note This method is provided just for backward compatibility reasons. Currently, it is not used by the library.
   *  \param[in] thr A stop condition relative error threshold. */
  void unrefine(double thr);

  /// A reference to an element.
  struct ElementReference {
    int id; ///< An element ID. Invalid if below 0.
    int comp; ///< A component which this element belongs to. Invalid if below 0.
    ElementReference(int id = -1, int comp = -1) : id(id), comp(comp) {}; ///< Constructor. It creates an invalid element reference.
  };

  /// Returns a squared error of an element.
  /** \param[in] A component index.
   *  \param[in] An element index.
   *  \return Squared error. Meaning of the error depends on parameters of the function calc_elem_errors(). */
  double get_element_error_squared(int component, int id) const { error_if(!have_errors, "Element errors have to be calculated first, call calc_elem_errors()."); return errors_squared[component][id]; };

  /// Returns regular queue of elements
  /** \return A regular queue. */
  const std::vector<ElementReference>& get_regular_queue() const { return regular_queue; };

  /// Returns a total number of active elements.
  /** \return A total number of active elements. If below 0, errors were not calculated yet, see set_solutions() */
  int get_total_active_elements() const { return num_act_elems; };

  /// Apply a single refienement.
  /** \param[in] A refinement to apply. */
  void apply_refinement(const ElementToRefine& elem_ref);

  /// Apply a vector of refinements.
  /** \param[in] A vector of refinements to apply. */
  virtual void apply_refinements(std::vector<ElementToRefine>& elems_to_refine);

  /// Returns a vector of refinements generated during the last execution of the method adapt().
  /** \return A vector of refinements generated during the last execution of the method adapt(). The returned vector might change or become invalid after the next execution of the method adadpt(). */
  const std::vector<ElementToRefine>& get_last_refinements() const; ///< Returns last refinements.

protected: //adaptivity
  int num_act_elems; ///< A total number of active elements across all provided meshes.
  std::queue<ElementReference> priority_queue; ///< A queue of priority elements. Elements in this queue are processed before the elements in the Adapt::regular_queue.
  std::vector<ElementReference> regular_queue; ///< A queue of elements which should be processes. The queue had to be filled by the method fill_regular_queue().
  std::vector<ElementToRefine> last_refinements; ///< A vector of refinements generated during the last finished execution of the method adapt().

  /// Returns true if a given element should be ignored and not processed through refinement selection.
  /** Overload this method to omit some elements from processing.
   *  \param[in] inx_element An index of an element in the regular queue. -1 if the element cames from the priority queue.
   *  \param[in] mesh A mesh that contains the element.
   *  \parar[in] element A pointer to the element.
   *  \return True if the element should be skipped. */
  virtual bool should_ignore_element(const int inx_element, const Mesh* mesh, const Element* element) { return false; };

  /// Returns true if a given element can be refined using proposed refinement.
  /** Overload this method to
   *  - avoid application of a refinement even thought a selector considered this refinement as the optimal one,
   *  - suggest a new refinement despite that the selector was not able to select a refinement.
   *  \param[in] mesh A mesh that contains the element.
   *  \param[in] e A point to the element.
   *  \param[in] refined True if a refinement of the element was found.
   *  \param[in,out] elem_ref The proposed refinement. Change a value of this parameter to select a different refinement.
   *  \return True if the element should not be refined using the refinement. */
  virtual bool can_refine_element(Mesh* mesh, Element* e, bool refined, ElementToRefine& elem_ref) { return refined; };

  /// Fixes refinements of a mesh which is shared among multiple components of a multimesh.
  /** If a mesh is shared among components, it has to be refined similarly in order to avoid incosistency.
   *  \param[in] meshes An array of meshes of components.
   *  \param[in] elems_to_refine A vector of refinements.
   *  \param[in] idx A 2D array that translates a pair (a component index, an element id) to an index of a refinement in the vector of refinements. If the index is below zero, a given element was not refined.
   *  \param[in] refinement_selector A selected used by the adaptivity. The selector is used to correct orders of modified refinements using RefinementSelectors::Selector::update_shared_mesh_orders(). */
  void fix_shared_mesh_refinements(Mesh** meshes, std::vector<ElementToRefine>& elems_to_refine, AutoLocalArray2<int>& idx, 
       Tuple<RefinementSelectors::Selector *> refinement_selectors);

  /// Enforces the same order to an element of a mesh which is shared among multiple compoenets.
  /** \param[in] meshes An arrat of meshes of components. */
  void homogenize_shared_mesh_orders(Mesh** meshes);

protected: //object state
  bool have_errors; ///< True if errors of elements were calculated.
  bool have_solutions; ///< True if solutions were set.

protected: // spaces & solutions
  int neq;                              ///< Number of solution components (as in wf->neq).
  Tuple<Space*> spaces;                 ///< Spaces. 
  Solution* sln[H2D_MAX_COMPONENTS];    ///< Coarse solution. 
  Solution* rsln[H2D_MAX_COMPONENTS];   ///< Reference solutions. 

protected: // element error arrays
  double* errors_squared[H2D_MAX_COMPONENTS]; ///< Errors of elements. Meaning of the error depeds on flags used when the method calc_elem_errors() was calls. Initialized in the method calc_elem_errors().
  double  errors_squared_sum; ///< Sum of errors in the array Adapt::errors_squared. Used by a method adapt() in some strategies.

protected: //forms and error evaluation
  matrix_form_val_t form[H2D_MAX_COMPONENTS][H2D_MAX_COMPONENTS]; ///< Bilinear forms to calculate error
  matrix_form_ord_t ord[H2D_MAX_COMPONENTS][H2D_MAX_COMPONENTS];  ///< Bilinear forms to calculate error

  /// Evaluates a square of an absolute error of an active element among a given pair of components.
  /** The method uses a bilinear forms to calculate the error. This is done by supplying a differences (f1 - v1) and (f2 - v2) at integration points to the bilinear form,
   *  where f1 and f2 are values of (coarse) solutions of the first and the second component respectively,
   *  v1 and v2 are values of reference solutions of the first and the second component respectively.
   *  \param[in] bi_fn A bilinear form.
   *  \param[in] bi_ord A bilinear form which is used to calculate the order of a quadrature.
   *  \param[in] sln1 A (coarse) solution corresponding to the first component.
   *  \param[in] sln2 A (coarse) solution corresponding to the second component.
   *  \param[in] rsln1 A reference solution corresponding to the first component.
   *  \param[in] rsln2 A reference solution corresponding to the second component.
   *  \param[in] rv1 A reference map of a (coarse) solution sln1.
   *  \param[in] rv2 A reference map of a (coarse) solution sln2.
   *  \param[in] rrv1 A reference map of a reference solution rsln1.
   *  \param[in] rrv2 A reference map of a reference solution rsln2.
   *  \return A square of an absolute error. */
  virtual double eval_elem_error_squared(matrix_form_val_t bi_fn, matrix_form_ord_t bi_ord,
                    MeshFunction *sln1, MeshFunction *sln2, MeshFunction *rsln1, MeshFunction *rsln2,
                    RefMap *rv1,        RefMap *rv2,        RefMap *rrv1,        RefMap *rrv2);

  /// Evaluates a square of a norm of an active element in the reference solution among a given pair of components.
  /** The method uses a bilinear forms to calculate the norm. This is done by supplying a v1 and v2 at integration points to the bilinear form,
   *  where v1 and v2 are values of reference solutions of the first and the second component respectively.
   *  \param[in] bi_fn A bilinear form.
   *  \param[in] bi_ord A bilinear form which is used to calculate the order of a quadrature.
   *  \param[in] rsln1 A reference solution corresponding to the first component.
   *  \param[in] rsln2 A reference solution corresponding to the second component.
   *  \param[in] rrv1 A reference map of a reference solution rsln1.
   *  \param[in] rrv2 A reference map of a reference solution rsln2.
   *  \return A square of a norm. */
  virtual double eval_elem_norm_squared(matrix_form_val_t bi_fn, matrix_form_ord_t bi_ord,
                   MeshFunction *rsln1, MeshFunction *rsln2, RefMap *rrv1, RefMap *rrv2);

  /// Builds an ordered queue of elements that are be examined.
  /** The method fills Adapt::standard_queue by elements sorted accordin to their error descending.
   *  The method assumes that Adapt::errors_squared contains valid values.
   *  If a special order of elements is requested, this method has to be overriden.
   *  /param[in] meshes An array of pointers to meshes of a (coarse) solution. An index into the array is an index of a component.
   *  /param[in] meshes An array of pointers to meshes of a reference solution. An index into the array is an index of a component. */
  virtual void fill_regular_queue(Mesh** meshes, Mesh** ref_meshes);

private: 
  /// A functor that compares elements accoring to their error. Used by std::sort().
  class CompareElements {
  private:
    double** errors_squared; ///< A 2D array of squared errors: the first index is an index of component, the second index is an element ID.
  public:
    CompareElements(double** errors_squared): errors_squared(errors_squared) {}; ///< Constructor.
    /// Compares two elements.
    /** \param[in] e1 A reference to the first element.
     *  \param[in] e1 A reference to the second element.
     *  \return True if a squared error of the first element is greater than a squared error of the second element. */
    bool operator ()(const ElementReference& e1,const ElementReference& e2) const {
      return errors_squared[e1.comp][e1.id] > errors_squared[e2.comp][e2.id];
    };
  };
};

// Mesh is adapted to represent a given function with given accuracy
// in a given projection norm.
void adapt_to_exact_function(Space *space, int proj_norm, ExactFunction exactfn,
                    RefinementSelectors::Selector* selector, double threshold, int strategy,
                    int mesh_regularity, double err_stop, int ndof_stop, bool verbose,
                    Solution* sln);

#endif
