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

#ifndef __H2D_REFINEMENT_PROJ_BASED_SELECTOR_H
#define __H2D_REFINEMENT_PROJ_BASED_SELECTOR_H

#include "optimum_selector.h"

using namespace Hermes::Algebra::DenseMatrixOperations;
namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors {
      /// Error of an element of a candidate for various permutations of orders. \ingroup g_selectors
      /** If not noted otherwise, the first index is the horizontal order, the second index is the vertical order.
      *  The maximum allowed order is ::H2DRS_MAX_ORDER + 1. */
      typedef double CandElemProjError[H2DRS_MAX_ORDER + 2][H2DRS_MAX_ORDER + 2];

      /// A general projection-based selector. \ingroup g_selectors
      /** Calculates an error of a candidate as a combination of errors of
      *  elements of a candidate. Each element of a candidate is calculated
      *  separatelly.
      *
      *  \section s_override Expanding
      *  In order to implement a support for a new space or a new approach to calculation of squared error,
      *  implement following methods:
      *  - precalc_shapes()[optional]
      *  - precalc_ortho_shapes()[optional]
      *  - precalc_ref_solution()
      *  - build_projection_matrix()
      *  - evaluate_rhs_subdomain()
      *  - evaluate_error_squared_subdomain()
      */
      template<typename Scalar>
      class HERMES_API ProjBasedSelector : public OptimumSelector<Scalar> {
      protected:
        class TrfShapeExp;

      public: //API
        /// Destructor
        virtual ~ProjBasedSelector();

        /// Sets error weights.
        /** An error weight is a multiplicative coefficient that modifies an error of candidate.
        *  Error weights can be used to proritize refinements.
        *  Error weights are applied in the method evaluate_cands_error().
        *  \param[in] weight_h An error weight of H-candidate. The default value is ::H2DRS_DEFAULT_ERR_WEIGHT_H.
        *  \param[in] weight_p An error weight of P-candidate. The default value is ::H2DRS_DEFAULT_ERR_WEIGHT_P.
        *  \param[in] weight_aniso An error weight of ANISO-candidate. The default value is ::H2DRS_DEFAULT_ERR_WEIGHT_ANISO. */
        void set_error_weights(double weight_h = H2DRS_DEFAULT_ERR_WEIGHT_H, double weight_p = H2DRS_DEFAULT_ERR_WEIGHT_P, double weight_aniso = H2DRS_DEFAULT_ERR_WEIGHT_ANISO);

        double get_error_weight_h() const;
        double get_error_weight_p() const;
        double get_error_weight_aniso() const;

        /// Evaluated shapes for all possible transformations for all points. The first index is a transformation, the second index is an index of a shape function.
        typedef Hermes::vector<TrfShapeExp> TrfShape[H2D_TRF_NUM];

        bool* cached_shape_vals_valid; ///< True if shape values were already initialized.
        TrfShape* cached_shape_ortho_vals; ///< Precalculated valus of orthogonalized shape functions.
        TrfShape* cached_shape_vals; ///< Precalculate values of shape functions.

      protected: //evaluated shape basis
        /// A transform shaped function expansions.
        /** The contents of the class can be accessed through an array index operator.
        *  The first index is an index of the function expansion, the second index is
        *  and index of an integration point.
        *
        *  \important Use strictly through reference.
        *  This allocates an array internally through the method allocate().
        *  If the use attempt to assign an instance which is not empty (i.e., the array is not allocate),
        *  an error will be issued. Assigning an empty instance will cause
        *  deallocation of the internal structures. */
        class TrfShapeExp {
        public:
          /// A default contructor. Creates an empty instance.
          TrfShapeExp();

          /// Desructor.
          virtual ~TrfShapeExp();

          /// Assignment operator. Prevent unauthorized copying of the pointer.
          const TrfShapeExp& operator = (const TrfShapeExp& other)
          {
            delete [] values; values = NULL;
            if(other.values == NULL)
              throw Exceptions::Exception("Unable to assign a non-empty values. Use references instead.");
            return *this;
          }
        private:
          int num_gip; ///< A number of integration points.
          int num_expansion; ///< A number of expansions.
          double** values; ///< Values. The first index is index of a functions expansion, the second index is an index of a an integration point.

          /// Allocates a space for function expansions.
          /** \param[in] num_expansion A number of expansions.
          *  \param[in] num_gip A number of itegration points. */
          void allocate(int num_expansion, int num_gip);

          /// Operator for accessing of contents.
          /** \param[in] inx_expansion An index of a function expansion.
          *  \return A pointer to an array of values of a function expansion at integration points. The returned pointer should not be stored outside the calling method or deallocated. */
          inline double* operator[](int inx_expansion);

          /// Returns true if the instance is empty, i.e., the method allocate() was not called yet.
          /** \return True if the instance is empty, i.e., the method allocate() was not called yet. */
          inline bool empty() const;

          template<typename T> friend class ProjBasedSelector;
          template<typename T> friend class L2ProjBasedSelector;
          template<typename T> friend class H1ProjBasedSelector;
          template<typename T> friend class HcurlProjBasedSelector;
          template<typename T> friend class Adapt;
        };

        /// Calculates values of shape function at GIP for all transformations.
        /** Override this method to supply a pre-calculated vales of shape function expansions
        *  at integration points. If override, the method has to supply precalculate expansions
        *  for all transformations plus an identity transformation. An index of the identity
        *  transformation is ::H2D_TRF_IDENTITY.
        *  \param[in] gip_points Integration points. The first index is an index of an integration point, the second index is an element of the enum ::GIP2DIndices.
        *  \param[in] num_gip_points A number of integration points.
        *  \param[in] trfs A transformations. The array has ::H2D_TRF_NUM elements. The index of the identity transformation is ::H2D_TRF_IDENTITY.
        *  \param[in] num_noni_trfs A number of transformations which are not identity. This number might be lower than than ::H2D_TRF_NUM.
        *  \param[in] shapes Shape functions.
        *  \param[in] max_shape_inx A maximum index of a shape function. This is used to resize a sub-array of the array \a svals.
        *  \param[out] svals A precalculated values of shape functions. The user has to resize the array of shape functions for every
        *             used transformation including the identity to the a size defined by \a max_shape_inx. The system will assume that shape functions
        *             are precalculated if the array corresponding to the identity function is not empty.
        */
        virtual void precalc_shapes(const double3* gip_points, const int num_gip_points, const Trf* trfs, const int num_noni_trfs, const Hermes::vector<typename OptimumSelector<Scalar>::ShapeInx>& shapes, const int max_shape_inx, TrfShape& svals, ElementMode2D mode) {};

        /// Calculates values of orthogonalized shape function at GIP for all transformations.
        /** Override this method to supply a pre-calculated vales of orthonormalized shape function expansions
        *  at integration points. If override, the method has to supply precalculate expansions
        *  for all transformations plus an identity transformation. An index of the identity
        *  transformation is ::H2D_TRF_IDENTITY.
        *
        *  If overridden and if this method orthonormalizes shape functions at the integration points, it is suggested
        *  to use the method precalc_shapes() in order to obtain initial values of shape functions at integration points.
        *  \param[in] gip_points Integration points. The first index is an index of an integration point, the second index is an element of the enum ::GIP2DIndices.
        *  \param[in] num_gip_points A number of integration points.
        *  \param[in] trfs A transformations. The array has ::H2D_TRF_NUM elements. The index of the identity transformation is ::H2D_TRF_IDENTITY.
        *  \param[in] num_noni_trfs A number of transformations which are not identity. This number might be lower than than ::H2D_TRF_NUM.
        *  \param[in] shapes Shape functions.
        *  \param[in] max_shape_inx A maximum index of a shape function. This is used to resize a sub-array of the array \a svals.
        *  \param[out] svals A precalculated values of shape functions. The user has to resize the array of shape functions for every
        *             used transformation including the identity to the a size defined by \a max_shape_inx. The system will assume that shape functions
        *             are precalculated if the array corresponding to the identity function is not empty.
        */
        virtual void precalc_ortho_shapes(const double3* gip_points, const int num_gip_points, const Trf* trfs, const int num_noni_trfs, const Hermes::vector<typename OptimumSelector<Scalar>::ShapeInx>& shapes, const int max_shape_inx, TrfShape& ortho_svals, ElementMode2D mode) {};

      protected:
        /// Constructor.
        /** Intializes attributes, projection matrix cache (ProjBasedSelector::proj_matrix_cache), and allocates rhs cache (ProjBasedSelector::rhs_cache).
        *  \param[in] cand_list A predefined list of candidates.
        *  \param[in] max_order A maximum order which considered. If ::H2DRS_DEFAULT_ORDER, a maximum order supported by the selector is used.
        *  \param[in] shapeset A shapeset. It cannot be NULL.
        *  \param[in] vertex_order A range of orders for vertex functions. Use an empty range (i.e. Range()) to skip vertex functions.
        *  \param[in] edge_bubble_order A range of orders for edge and bubble functions. Use an empty range (i.e. Range()) to skip edge and bubble functions. */
        ProjBasedSelector(CandList cand_list, int max_order, Shapeset* shapeset, const Range& vertex_order, const Range& edge_bubble_order);

      protected: //internal logic
        /// True if the selector has already warned about possible inefficiency.
        /** If OptimumSelector::cand_list does not generate candidates with elements of
        *  non-uniform orders and if the precalculate orthonormal base is available,
        *  the user should not use non-uniform order in the meshed in
        *  order to gain efficiency. */
        bool warn_uniform_orders;

      protected: //error evaluation
        static const int H2DRS_VALCACHE_INVALID = 0; ///< State of value cache: item contains undefined or invalid value. \ingroup g_selectors
        static const int H2DRS_VALCACHE_VALID = 1; ///< State of value cache: item contains a valid value. \ingroup g_selectors
        static const int H2DRS_VALCACHE_USER = 2; ///< State of value cache: the first state ID which can be used by the user. \ingroup g_selectors

        /// An item of a value cache.
        template<typename T>
        struct ValueCacheItem {
          /// Returns true if value is mared as valid.
          /** \return True if the value is marked as valid. */
          bool is_valid() const { return state != H2DRS_VALCACHE_INVALID; };

          /// Marks a value.
          /** \param[in] new_state A new state of the value. By default, it marks the value as valid. */
          void mark(int new_state = H2DRS_VALCACHE_VALID) { state = new_state; };

          /// Sets a value.
          /** \param[in] new_value A new value. It does not change state of the value. */
          void set(T new_value) { value = new_value; };

          /// Returns the value. It does check the state of the value.
          /** \return A current value. */
          T get() const { return value; };

          /// Constructor.
          /** By default, the item is set as invalid.
          *  \param value A starting value.
          *  \param state A state of the value. */
          ValueCacheItem(const T& value = 0, const int state = H2DRS_VALCACHE_INVALID) : value(value), state(state) {}; ///< Default constructor. By default, it creates an item that contains invalid value.
        private:
          T value; ///< A value stored in the item.
          int state; ///< A state of the image: ::H2DRS_VALCACHE_INVALID or ::H2DRS_VALCACHE_VALID or any other user-defined value. The first user defined state has to have number ::H2DRS_VALCACHE_USER.
        };
        /// A projection matrix cache type.
        /** Defines a cache of projection matrices for all possible permutations of orders. */
        typedef double** ProjMatrixCache[H2DRS_MAX_ORDER + 2][H2DRS_MAX_ORDER + 2];

        /// An array of projection matrices.
        /** The first index is the mode (see the enum ElementMode2D). The second and the third index
        *  is the horizontal and the vertical order respectively.
        *
        *  All matrices are square dense matrices and they have to be created through the function new_matrix().
        *  If record is NULL, the corresponding matrix has to be calculated. */
        ProjMatrixCache proj_matrix_cache[H2D_NUM_MODES];

        double error_weight_h; ///< A coefficient that multiplies error of H-candidate. The default value is ::H2DRS_DEFAULT_ERR_WEIGHT_H.
        double error_weight_p; ///< A coefficient that multiplies error of P-candidate. The default value is ::H2DRS_DEFAULT_ERR_WEIGHT_P.
        double error_weight_aniso; ///< A coefficient that multiplies error of ANISO-candidate. The default value is ::H2DRS_DEFAULT_ERR_WEIGHT_ANISO.

        /// Calculates error of candidates.
        /** Overriden function. For details, see OptimumSelector::evaluate_cands_error(). */
        virtual void evaluate_cands_error(Hermes::vector<Cand>& candidates, Element* e, MeshFunction<Scalar>* rsln);

        /// Calculates projection errors of an elements of candidates for all permutations of orders.
        /** Errors are not normalized and they are squared.
        *  The range of orders is defined through parameters \a info_h, \a info_h, and \a info_aniso.
        *
        *  If defining a new evalution (e.g. using a difference space) of errors,
        *  follows instructions in \ref s_override.
        *  \param[in] e An element that is being examined by the selector.
        *  \param[in] info_h Information about H-candidates: range of orders, etc.
        *  \param[in] info_p Information about P-candidates: range of orders, etc.
        *  \param[in] info_aniso Information about ANISO-candidates: range of orders, etc.
        *  \param[in] rsln A reference solution.
        *  \param[out] herr An error of elements of H-candidates of various permutation of orders.
        *  \param[out] perr An error of elements of P-candidates of various permutation of orders.
        *  \param[out] anisoerr An error of elements of ANISO-candidates of various permutation of orders. */
        virtual void calc_projection_errors(Element* e, const typename OptimumSelector<Scalar>::CandsInfo& info_h, const typename OptimumSelector<Scalar>::CandsInfo& info_p, const typename OptimumSelector<Scalar>::CandsInfo& info_aniso, MeshFunction<Scalar>* rsln, CandElemProjError herr[H2D_MAX_ELEMENT_SONS], CandElemProjError perr, CandElemProjError anisoerr[H2D_MAX_ELEMENT_SONS]);

        /// Calculate projection errors of an element of an candidate considering multiple orders.
        /** An element of a candidate may span over multiple sub-domains. All integration uses the reference domain.
        *
        *  Modify this method in order to add ortho-adaptivity.
        *  \param[in] mode A mode (enum ElementMode2D).
        *  \param[in] gip_points Integration points in the reference domain.
        *  \param[in] num_gip_points A number of integration points.
        *  \param[in] num_sub A number of subdomains.
        *  \param[in] sub_domains Subdomains (elements of a reference mesh) that occupy the element of a candidate. The first index is an index of the subdomain.
        *  \param[in] sub_trfs Transformation from a reference domain of a subdomain to a reference domain of the element of a candidate. The first index is an index of the subdomain.
        *  \param[in] sub_rvals Values at integration points for every subdomain. Contents of this array (the second index) is defined by the method precalc_ref_solution(). The first index is an index of the subdomain.
        *  \param[in] sub_nonortho_svals
        *  \param[in] sub_ortho_svals
        *  \param[in] info Information about candidates: range of orders, etc.
        *  \param[out] errors_squared Calculated squared errors for all orders specified through \a info. */
        void calc_error_cand_element(const ElementMode2D mode, double3* gip_points, int num_gip_points, const int num_sub, Element** sub_domains, Trf** sub_trfs, Scalar*** sub_rvals, Hermes::vector<TrfShapeExp>** sub_nonortho_svals, Hermes::vector<TrfShapeExp>** sub_ortho_svals, const typename OptimumSelector<Scalar>::CandsInfo& info, CandElemProjError errors_squared);

      protected: //projection
        /// Projection of an element of a candidate.
        struct ElemProj {
          int* shape_inxs; ///< Used shape indices
          int num_shapes; ///< A number of used shape indices.
          Hermes::vector<TrfShapeExp>& svals; ///< A precalculated shape-function values. Empty is not defined.
          Scalar* shape_coeffs; ///< Coefficients of shape indices of a projection.
          int max_quad_order; ///< An encoded maximum order of the projection. If triangle, the vertical order is equal to the horizontal order.
        };

        /// Integration points in the reference domain of an element of a candidate.
        /** The structure assumes Gauss Integration Points (GIP). */
        struct ElemGIP {
          double3* gip_points; ///< Integration points and weights. The first index is an index of an integration point, the second index is defined through the enum GIP2DIndices.
          int num_gip_points; ///< A number of integration points.
          Scalar** rvals; ///< Values of a reference solution at the integration points. The first index is an index of the function expansion (f, df/dx, ...), the second index is an index of the integration point. The meaning of the second index is defined through the method precalc_ref_solution().
        };

        /// A transformation from a reference domain of a subdomain to a reference domain of an element of a candidate.
        struct ElemSubTrf {
          Trf* trf; ///< A transformation.
          double coef_mx; ///< A coefficient that scales df/dx for each subdomain. A coefficient represents effects of a transformation \a trf on df/dx.
          double coef_my; ///< A coefficient that scales df/dy for each subdomain. A coefficient represents effects of a transformation \a trf on df/dy.
        };

        /// A shape function on subdomain of an element.
        struct ElemSubShapeFunc {
          int inx; ///< An index of a shape function.
          TrfShapeExp& svals; ///< Evaluate values of a shape function. If TrfShapeExp::empty(), no precalculated values are available.
        };

        /// Returns an array of values of the reference solution at integration points.
        /** The method have to set an active element and an quadrature on its own.
        *
        *  Override to provide all necessary values. Provided pointers should stay valid through an execution
        *  of the method calc_projection_errors(). Since an explicit deallocation of these pointers is not done,
        *  it is suggested to provide pointers to attributes of the class rather than to dynamically allocate
        *  an array.
        *  The method can assume that the an element is refined into ::H2D_MAX_ELEMENT_SONS elements (sons) in the reference mesh.
        *  \param[in] inx_son An index of a son of an element, i.e., an index of a subdomain. The index is in a range[0, H2D_MAX_ELEMENT_SONS - 1].
        *  \param[in] rsln A reference solution.
        *  \param[in] element An element of the coarse solution. An element of both the same geometry and the same ID have to be present in the mesh of the reference solution.
        *  \param[in] intr_gip_order An order of quadrature integration. The number of quadrature points should be retrieved through a quadrature stored in the paremeter \a rsln.
        *  \return A pointer to 2D array. The first index is an index of the function expansion (f, df/dx, ...), the second index is an index of the integration point. */
        virtual Scalar** precalc_ref_solution(int inx_son, MeshFunction<Scalar>* rsln, Element* element, int intr_gip_order) = 0;

        /// Builds projection matrix using a given set of shapes.
        /** Override to calculate a projection matrix.
        *  \param[in] gip_points Integration points. The first index is an index of an integration point, the second index is defined through the enum GIP2DIndices.
        *  \param[in] num_gip_points A number of integration points.
        *  \param[in] shape_inx An array of shape indices.
        *  \param[in] num_shapes A number of shape indices in the array.
        *  \return A projection matrix. The matrix has to be allocated trought new_matrix(). The size of the matrix has to be \a num_shapes x \a num_shapes. */
        virtual double** build_projection_matrix(double3* gip_points, int num_gip_points, const int* shape_inx, const int num_shapes, ElementMode2D) = 0;

        /// Evaluates a value of the right-hande side in a subdomain.
        /** Override to calculate a value of the right-hand side.
        *  \param[in] sub_elem An element of a reference mesh that corresponds to a subdomain.
        *  \param[in] sub_gip Integration points. Locations of integration points are defined in the reference domain. Use \a sub_trf to transform it to the reference domain of an element of a candidate.
        *  \param[in] sub_trf A transformation from a reference domain of a subdomain to the reference domain of an element of a candidate.
        *  \param[in] sub_shape Information about a shape function: shape index and calculated expansions at integration points, if any.
        *  \return A value of the righ-hand size of a given shape function. */
        virtual Scalar evaluate_rhs_subdomain(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemSubShapeFunc& sub_shape) = 0;

        /// Evaluates an squared error of a projection of an element of a candidate onto subdomains.
        /** Override to calculate an error using a provided projection and subdomains.
        *  \param[in] sub_elem An element of a reference mesh that corresponds to a subdomain.
        *  \param[in] sub_gip Integration points. Locations of integration points are defined in the reference domain. Use \a sub_trf to transform it to the reference domain of an element of a candidate.
        *  \param[in] sub_trf A transformation from a reference domain of a subdomain to the reference domain of an element of a candidate.
        *  \param[in] elem_proj A projection of an element of a candidate on subdomains.
        *  \return A squared error of an element of a candidate. */
        virtual double evaluate_error_squared_subdomain(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemProj& elem_proj) = 0;
      };
    }
  }
}
#endif