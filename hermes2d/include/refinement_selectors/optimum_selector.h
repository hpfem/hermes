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

#ifndef __H2D_REFINEMENT_OPTIMUM_SELECTOR_H
#define __H2D_REFINEMENT_OPTIMUM_SELECTOR_H

#include <ostream>
#include "selector.h"
#include "../shapeset/shapeset.h"
#include "candidates.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors
    {
      /// A selector that chooses an optimal candidates based on a score. \ingroup g_selectors
      /** This is a base class for all selectors that chooses an candidate based on some
      *  evaluated criteria. Currently, the criteria is based on an error change per DOF. */
      template<typename Scalar>
      class HERMES_API OptimumSelector : public Selector<Scalar>
      {
      public:
        /// Destructor.
        virtual ~OptimumSelector();

        /// Set the score DOF exponent.
        void set_dof_score_exponent(double exponent);

      protected:
        /// Constructor.
        /** \note Parameters \a vertex_order and \a edge_bubble_order fixes the fact that a shapeset returns a valid index even though a given shape is not invalid in the space.
        *  \param[in] cand_list A predefined list of candidates.
        *  \param[in] max_order A maximum order which considered. If ::H2DRS_DEFAULT_ORDER, a maximum order supported by the selector is used.
        *  \param[in] shapeset A shapeset. It cannot be NULL.
        *  \param[in] vertex_order A range of orders for vertex functions. Use an empty range (i.e. Range()) to skip vertex functions.
        *  \param[in] edge_bubble_order A range of orders for edge and bubble functions. Use an empty range (i.e. Range()) to skip edge and bubble functions. */
        OptimumSelector(CandList cand_list, int max_order, Shapeset* shapeset, const Range& vertex_order, const Range& edge_bubble_order);

        //candidates
        /// Information about candidates.
        struct CandsInfo
        {
          bool uniform_orders; ///< True if all elements of all examined candidates have uniform orders.
          int min_quad_order; ///< Minimum quad order of all elements of all examined candidates.
          int max_quad_order; ///< Maximum quad order of all elements of all examined candidates. If less than zero, no candidate is generated.

          /// Default constructor. Creates info that declares no candidates and uniform orders.
          CandsInfo() : uniform_orders(true), min_quad_order(-1), max_quad_order(-1) {};

          /// Returns true if there are no candidates.
          /** \return True if there are no candidates. */
          bool is_empty() const { return (min_quad_order < 0 || max_quad_order < 0); };
        };

        CandList cand_list; ///< Allowed candidate types.

        /// Updates information about candidates. Initial information is provided.
        /** \param[in,out] info_h Information about all H-candidates.
        *  \param[in,out] info_p Information about all P-candidates.
        *  \param[in,out] info_aniso Information about all ANISO-candidates. */
        void update_cands_info(Hermes::vector<Cand>& candidates, CandsInfo& info_h, CandsInfo& info_p, CandsInfo& info_aniso) const;

        /// Appends cancidates of a given refinement and a given range of orders.
        /** If either borders or a ranges is invalid (i.e. smaller than zero)
        * or if a upper boundary is below the lower boudary, no candidate is appended.
        *  \param[in] start_quad_order A lower boundary of a range in a form of an encoded order.
        *  \param[in] last_order The upper boundery of a range in a form of an encoded order.
        *  \param[in] split A refinement, see the enum RefinementTypes.
        *  \param[in] iso_p True if both orders (horizontal and vertical) should be modified uniformly. Used in a case of a triangle. */
        void append_candidates_split(Hermes::vector<Cand>& candidates, const int start_quad_order, const int last_order, const int split, bool iso_p);

        /// Fill a list of candidates.
        /** Override to generate or adjust generated candidates. The method has to initialize the array OptimumSelector::candidates.
        *  If triangle, all generated candidates have to have the vertical order equal to the horizontal order.
        *  An order of any element of any candidate has to fit into a range[OptimumSelector::current_min_order, OptimumSelector::current_max_order].
        *  \param[in] e An element that is being refined.
        *  \param[in] quad_order An encoded order of the element. If triangle, the vertical order is equal to the horizontal order.
        *  \param[in] max_ha_quad_order A maximum encoded order of an element of a H-candidate or an ANISO-candidate. In the case of ANIO-candidates, the maximum is applied only to modified orders.
        *  \param[in] max_p_quad_order A maximum encoded order of an element of a P-candidate.
        *  \return A vector of candidates. The first candidate has to be equal to the original element with a refinement ::H2D_REFINEMENT_P.
         */
        virtual Hermes::vector<Cand> create_candidates(Element* e, int quad_order);

        /// Calculates error, dofs, and score of candidates.
        /** \param[in] e An element that is being refined.
        *  \param[in] rsln A reference solution which is used to calculate the error.
        *  \param[out] avg_error An average of \f$\log_{10} e\f$ where \f$e\f$ is an error of a candidate. It cannot be NULL.
        *  \param[out] dev_error A deviation of \f$\log_{10} e\f$ where \f$e\f$ is an error of a candidate. It cannot be NULL. */
        void evaluate_candidates(Hermes::vector<Cand>& candidates, Element* e, MeshFunction<Scalar>* rsln);

        /// Sorts and selects the best candidate and the best H-candidate according to the score.
        /** Any two candidates with the same score are skipped since it is not possible to decide between them.
        *  The method assumes that the candidate at the index 0 is the original element therefore
        *  it skips this candidate. The selected H-candidate is used to handle a case when a mesh is shared
        *  amond multiple components.
        *
        *  Override to redefined the algoritm of the selecting. If overridden, the implementation has to select both
        *  the best candidate and the best H-candidate.
        *  \param[in] e An element that is being refined.
        *  \param[in] avg_error An average of \f$\log_{10} e\f$ where \f$e\f$ is an error of a candidate.
        *  \param[in] dev_error A deviation of \f$\log_{10} e\f$ where \f$e\f$ is an error of a candidate.
        *  \param[out] best_candidates Four best candidates \
        *     0 - overall
        *     1 - 4 : indexed by enum RefinementType.
        */
        virtual void select_best_candidate(Hermes::vector<Cand>& candidates, Element* e, Cand*& best_candidate, Cand* best_candidates_specific_type[4]);

        /// Calculates error of candidates.
        /** This method has to be implemented in inherited classes.
        *  \param[in] e An element that is being refined.
        *  \param[in] rsln A reference solution which is used to calculate the error.
        */
        virtual void evaluate_cands_error(Hermes::vector<Cand>& candidates, Element* e, MeshFunction<Scalar>* rsln) = 0;

        /// Calculates DOF of candidates.
        /** It uses a list of shape indices (OptimumSelector::shape_indices) to
        *  count a number of DOFs. No number of DOFs cannot be zero.
        *  \param[in] e An element that is being refined.
        *  \param[in] rsln A reference solution which is used to calculate the error. */
        virtual void evaluate_cands_dof(Hermes::vector<Cand>& candidates, Element* e, MeshFunction<Scalar>* rsln);

        /// Evalutes score of candidates.
        /** It calculates score \f$s\f$ of a candidate as \f[s = \frac{\log_{10} e_0 - \log_{10} e}{(d - d_0)^c},\f]
        *  where \f$e\f$ and \f$e_0\f$ are errors of a candidate and the original element (i.e. a candidate at the index 0)
        *  respectively, \f$d\f$ and \f$d_0\f$ are number of DOFs of a candidate and the original element (i.e. a candidate at the index 0)
        *  respectively.
        *  If the score is zero, the score is invalid.
        *
        *  If overridden, the higher score the better candidate.
        *  \param[in] e An element that is being refined. */
        virtual void evaluate_cands_score(Hermes::vector<Cand>& candidates, Element* e);

        /// Number of shape functions for
        /// - mode
        /// - horizontal order + 1 (any)
        /// - vertical order + 1 (any)
        /// - shape function type
        int ****num_shapes;
      
        /// Compares scores. Used to sort scores ascending.
        /** \param[in] a The first candidate.
        *  \param[in] b The second candidate.
        *  \return True if score of \a a is greater than the score of \a b. */
        static bool compare_cand_score(const Cand& a, const Cand& b);

        //orders and their range
        /// Sets OptimumSelector::current_max_order and OptimumSelector::current_min_order.
        /** This method has to be implemented by derived classes and it is mean to be
        *  space dependent, i.e., it should differ in a case of H1, L2, and Hcurl.
        *  \param[in] element An element that is being refined. */
        virtual void get_current_order_range(Element* element, int& min_order, int& max_order) = 0;

        //shape functions
        /// A shape function type.
        enum ShapeType {
          H2DST_VERTEX = 0x01, ///< Vertex function.
          H2DST_HORIZ_EDGE = 0x02, ///< Horizontal edge function.
          H2DST_VERT_EDGE = 0x04, ///< Verical edge function.
          H2DST_TRI_EDGE = 0x08, ///< Triangle edge.
          H2DST_BUBBLE = 0x10 ///< Bubble function.
        };

        enum ShapeTypeInt {
          H2DSI_VERTEX,
          H2DSI_HORIZ_EDGE,
          H2DSI_VERT_EDGE,
          H2DSI_TRI_EDGE,
          H2DSI_BUBBLE,
          H2DSI_ANY
        };

        /// A shape index.
        /** Any element order higher than both the vertical and the horizontal direction will use a given shape function. */
        struct ShapeInx {
          int order_h; ///< A minimal horizonal order of an element that can use this shape function.
          int order_v; ///< A minimal vertical order of an element that can use this shape function.
          int inx; ///< An index of the shape function.
          ShapeType type; ///< A type of the shape function. It is used to calculate DOF in Optimum::evaluate_cands_dof().

          /// Constructor.
          /** \param[in] order_h A minimal horizonal order of an element that can use this shape function.
          *  \param[in] order_v A minimal vertical order of an element that can use this shape function.
          *  \param[in] inx An index of the shape function.
          *  \param[in] type A type of the shape function. */
          ShapeInx(int order_h, int order_v, int inx, ShapeType type) : order_h(order_h), order_v(order_v), inx(inx), type(type) {};
        };

        Shapeset *shapeset; ///< A shapeset used to calculate error.

        Hermes::vector<ShapeInx> shape_indices[H2D_NUM_MODES]; ///< Shape indices. The first index is a mode (ElementMode2D).
        int max_shape_inx[H2D_NUM_MODES]; ///< A maximum index of a shape function. The first index is a mode (ElementMode2D).
        int next_order_shape[H2D_NUM_MODES][H2DRS_MAX_ORDER+1]; ///< An index to the array OptimumSelector::shape_indices of a shape function of the next uniform order. The first index is a mode (ElementMode2D), the second index is an order.
        bool has_vertex_shape[H2D_NUM_MODES]; ///< True if the shapeset OptimumSelector::shapeset contains vertex functions. The index is a mode (ElementMode2D).
        bool has_edge_shape[H2D_NUM_MODES]; ///< True if the shapeset OptimumSelector::shapeset contains edge functions. The index is a mode (ElementMode2D).
        bool has_bubble_shape[H2D_NUM_MODES]; ///< True if the shapeset OptimumSelector::shapeset contains bubble functions. The index is a mode (ElementMode2D).

        /// Adds an index (or indices) of a bubble function of a given order if the shape index was not used yet.
        /** This function adds indices of bubble functions that were not added yet on a quadrilateral.
        *  Since a shapeset allows only to retrieve a list of all possible bubbles based on an order of an element,
        *  this functions allows to create a back-mapping table which converts shape index to the smallest element order that uses the shape function.
        *  The function assumes that all shapes of a lower element order than a given combinations were already added.
        *  Used by build_shape_indices().
        *  \param[in] order_h A horizontal order of an element.
        *  \param[in] order_v A vertical order of an element.
        *  \param[in,out] used_shape_index A vector of used shape indices. If a shape index is present in the map, a shape was already added and it will not be added again.
        *  \param[in,out] indices A vector of shape indices. The vector is updated by the function. */
        void add_bubble_shape_index(int order_h, int order_v, std::map<int, bool>& used_shape_index, Hermes::vector<ShapeInx>& indices, ElementMode2D mode);

        /// Builds shape index table OptimumSelector::shape_indices.
        /** The method fills the array OptimumSelector::shape_indices for a given mode.
        *  The method is called by the constructor.
        *  \param[in] mode A mode (ElementMode2D).
        *  \param[in] vertex_order A range of orders in which to search for vertex functions.
        *  \param[in] edge_bubble_order A range of order in which to search for edge and bubble functions. */
        void build_shape_indices(const ElementMode2D mode, const Range& vertex_order, const Range& edge_bubble_order);

        /// Returns a number of shapes that may be contained in an element of a given order.
        /** \param[in] mode A mode (ElementMode2D).
        *  \param[in] order_h A horizontal order of the element. If ::H2DRS_ORDER_ANY, any order is allowed.
        *  \param[in] order_v A horizontal order of the element. If ::H2DRS_ORDER_ANY, any order is allowed.
        *  \param[in] allowed_type_mask A combination of flags specifying which orders are allowed. Flags are defined in the enum ShapeType.
        *  \return Returns a number of shape functions that satisfies given parameters. */
        int calc_num_shapes(int mode, int order_h, int order_v, int allowed_type_mask);

        /// Selects a refinement.
        /** Overriden function. For details, see Selector::select_refinement(). */
        virtual bool select_refinement(Element* element, int quad_order, MeshFunction<Scalar>* rsln, ElementToRefine& refinement); ///< Selects refinement.

        /// Score DOF exponent. Used in evaluate_cands_score.
        double dof_score_exponent;
      };
    }
  }
}

#endif