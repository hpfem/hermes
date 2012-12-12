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
#include "order_permutator.h"
#include "selector.h"
#include "../shapeset/shapeset.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors {
      /// Predefined list of candidates. \ingroup g_selectors
      enum CandList {
        H2D_NONE,  ///< No adaptivity. (Used only in modules.)
        H2D_P_ISO, ///< P-candidates only. Hermes::Orders are modified uniformly.
        H2D_P_ANISO, ///< P-candidates only. Hermes::Orders are modified non-uniformly.
        H2D_H_ISO, ///< H-candidates only. Hermes::Orders are not modified.
        H2D_H_ANISO, ///< H- and ANISO-candidates only. Hermes::Orders are not modified.
        H2D_HP_ISO, ///< H- and P-candidates only. Hermes::Orders are modified uniformly.
        H2D_HP_ANISO_H, ///< H-, ANISO- and P-candidates. Hermes::Orders are modified uniformly.
        H2D_HP_ANISO_P, ///< H- and P-candidates only. Hermes::Orders are modified non-uniformly.
        H2D_HP_ANISO ///< H-, ANISO- and P-candidates. Hermes::Orders are modified non-uniformly.
      };

      /// Options of the selector. \ingroup g_selectors
      enum SelOption {
        H2D_PREFER_SYMMETRIC_MESH, ///< Prefer symmetric mesh when selection of the best candidate is done. If two or more candiates has the same score, they are skipped. This option is set by default.
        H2D_APPLY_CONV_EXP_DOF ///< Use \f$d^c - d_0^c\f$, where \f$c\f$ is the convergence exponent, instead of \f$(d - d_0)^c\f$ to evaluate the score in the method OptimumSelector::evaluate_cands_score(). This option is not set by default.
      };

      /// Returns a string representation of a predefined candidate list. \ingroup g_selectors
      /** Used for debugging and output purposes.
      *  \param cand_list A predefined list of candidates.
      *  \return A string representation of the enum value. */
      extern HERMES_API const char* get_cand_list_str(const CandList cand_list);

      /// Returns true if a predefined candidate list may contain candidates that are HP. \ingroup g_selectors
      /** \param cand_list A predefined list of candidates.
      *  \return True if a predefined candidate list may contain candidates that are HP. */
      extern HERMES_API bool is_hp(const CandList cand_list);

      /// Returns true if a predefined candidate list may contain candidates with an anisotropic change of orders. \ingroup g_selectors
      /** \param cand_list A predefined list of candidates.
      *  \return True if a predefined candidate list may contain candidates with an anisotropic change of orders. */
      extern HERMES_API bool is_p_aniso(const CandList cand_list);

      /// A selector that chooses an optimal candidates based on a score. \ingroup g_selectors
      /** This is a base class for all selectors that chooses an candidate based on some
      *  evaluated criteria. Currently, the criteria is based on an error change per DOF. */
      template<typename Scalar>
      class HERMES_API OptimumSelector : public Selector<Scalar> {
      protected: //options
        bool opt_symmetric_mesh; ///< True if ::H2D_PREFER_SYMMETRIC_MESH is set. True by default.
        bool opt_apply_exp_dof; ///< True if ::H2D_APPLY_CONV_EXP_DOF is set. False by default.
      public:
        /// Enables or disables an option.
        /** If overridden, the implementation has to call a parent implementation.
        *  \param[in] option An option that is going to be enabled. For possible values, see SelOption.
        *  \param[in] enable True to enable, false to disable. */
        virtual void set_option(const SelOption option, bool enable);

        /// A candidate.
        struct Cand {
          double error; ///< An error of this candidate.
          int dofs;  ///< An estimated number of DOFs.
          int split; ///< A refinement, see the enum RefinementType.
          int p[4]; ///< Encoded orders of sons, see ::H2D_MAKE_QUAD_ORDER. In a case of a triangle, the vertical order is equal to the horizontal one.
          double score; ///< A score of a candidate: the higher the better. If zero, the score is not valid and a candidate should be ignored. Evaluated in OptimumSelector::select_best_candidate.

          /// Constructor.
          /** \param[in] split A refinement, see the enum RefinementTypes.
          *  \param[in] order_elems Encoded orders for all element of candidate. If triangle, a vertical order has to be equal to the horizontal one. Unused elements of the array can be ignored. */
          Cand(const int split, const int order_elems[4])
            : dofs(-1), split(split), score(0) {
              p[0] = order_elems[0];
              p[1] = order_elems[1];
              p[2] = order_elems[2];
              p[3] = order_elems[3];
          };

          /// Constructor.
          /** \param[in] split A refinement, see the enum RefinementTypes.
          *  \param[in] order_elem0 Encoded order of the first element of the candidate. If triangle, a vertical order has to be equal to the horizontal one.
          *  \param[in] order_elem1 Encoded order of the second element of the candidate, if any. If triangle, a vertical order has to be equal to the horizontal one.
          *  \param[in] order_elem2 Encoded order of the third element of the candidate, if any. If triangle, a vertical order has to be equal to the horizontal one.
          *  \param[in] order_elem3 Encoded order of the fourth element of the candidate, if any. If triangle, a vertical order has to be equal to the horizontal one. */
          Cand(const int split, const int order_elem0, const int order_elem1 = 0, const int order_elem2 = 0, const int order_elem3 = 0)
            : dofs(-1), split(split), score(0) {
              p[0] = order_elem0;
              p[1] = order_elem1;
              p[2] = order_elem2;
              p[3] = order_elem3;
          };

          /// Returns a number of elements of a candidate.
          /** \return A number of elements of a candidate. */
          int get_num_elems() const {
            switch (split) {
            case H2D_REFINEMENT_H: return 4;
            case H2D_REFINEMENT_P: return 1;
            case H2D_REFINEMENT_ANISO_H:
            case H2D_REFINEMENT_ANISO_V:
              return 2;
            default:
              throw Hermes::Exceptions::Exception("Invalid refinement type %d.", split);
              return -1;
              break;
            }
          }
        };

        /// Returns a vector of the last generated candidates.
        /** \return A vector of last generated candidates. The vector will change if a new list is generated. */
        const Hermes::vector<Cand>& get_candidates() const { return candidates; };

        /// Number of shape functions for
        /// - mode
        /// - horizontal order + 1 (any)
        /// - vertical order + 1 (any)
        /// - shape function type
        int ****num_shapes;

      protected: //candidates
        /// Information about candidates.
        struct CandsInfo {
          bool uniform_orders; ///< True if all elements of all examined candidates have uniform orders.
          int min_quad_order; ///< Minimum quad order of all elements of all examined candidates. If less than zero, no candidate is generated.
          int max_quad_order; ///< Maximum quad order of all elements of all examined candidates. If less than zero, no candidate is generated.

          /// Default constructor. Creates info that declares no candidates and uniform orders.
          CandsInfo() : uniform_orders(true), min_quad_order(-1), max_quad_order(-1) {};

          /// Returns true if there are no candidates.
          /** \return True if there are no candidates. */
          bool is_empty() const { return (min_quad_order < 0 || max_quad_order < 0); };
        };

        CandList cand_list; ///< Allowed candidate types.
        double conv_exp; ///< Convergence power. Modifies difference between DOFs before they are used to calculate the score.
        Hermes::vector<Cand> candidates; ///< A vector of candidates. The first candidate has to be equal to the original element with a refinement ::H2D_REFINEMENT_P.

        /// Updates information about candidates. Initial information is provided.
        /** \param[in,out] info_h Information about all H-candidates.
        *  \param[in,out] info_p Information about all P-candidates.
        *  \param[in,out] info_aniso Information about all ANISO-candidates. */
        void update_cands_info(CandsInfo& info_h, CandsInfo& info_p, CandsInfo& info_aniso) const;

        /// Appends cancidates of a given refinement and a given range of orders.
        /** If either borders or a ranges is invalid (i.e. smaller than zero)
        * or if a upper boundary is below the lower boudary, no candidate is appended.
        *  \param[in] start_quad_order A lower boundary of a range in a form of an encoded order.
        *  \param[in] last_order The upper boundery of a range in a form of an encoded order.
        *  \param[in] split A refinement, see the enum RefinementTypes.
        *  \param[in] iso_p True if both orders (horizontal and vertical) should be modified uniformly. Used in a case of a triangle. */
        void append_candidates_split(const int start_quad_order, const int last_order, const int split, bool iso_p);

        /// Fill a list of candidates.
        /** Override to generate or adjust generated candidates. The method has to initialize the array OptimumSelector::candidates.
        *  If triangle, all generated candidates have to have the vertical order equal to the horizontal order.
        *  An order of any element of any candidate has to fit into a range[OptimumSelector::current_min_order, OptimumSelector::current_max_order].
        *  \param[in] e An element that is being refined.
        *  \param[in] quad_order An encoded order of the element. If triangle, the vertical order is equal to the horizontal order.
        *  \param[in] max_ha_quad_order A maximum encoded order of an element of a H-candidate or an ANISO-candidate. In the case of ANIO-candidates, the maximum is applied only to modified orders.
        *  \param[in] max_p_quad_order A maximum encoded order of an element of a P-candidate. */
        virtual void create_candidates(Element* e, int quad_order, int max_ha_quad_order, int max_p_quad_order);

        /// Calculates error, dofs, and score of candidates.
        /** \param[in] e An element that is being refined.
        *  \param[in] rsln A reference solution which is used to calculate the error.
        *  \param[out] avg_error An average of \f$\log_{10} e\f$ where \f$e\f$ is an error of a candidate. It cannot be NULL.
        *  \param[out] dev_error A deviation of \f$\log_{10} e\f$ where \f$e\f$ is an error of a candidate. It cannot be NULL. */
        void evaluate_candidates(Element* e, Solution<Scalar>* rsln, double* avg_error, double* dev_error);

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
        *  \param[out] selected_cand A pointer to a selected index of the best candidate. If the index is 0, the algorithm was not able to decide.
        *  \param[out] selected_h_cand A pointer to a selected index of the best H-candidate. If the index is 0, the algorithm was not able to decide.
        */
        virtual void select_best_candidate(Element* e, const double avg_error, const double dev_error, int* selected_cand, int* selected_h_cand);

        /// Calculates error of candidates.
        /** This method has to be implemented in inherited classes.
        *  \param[in] e An element that is being refined.
        *  \param[in] rsln A reference solution which is used to calculate the error.
        *  \param[out] avg_error An average of \f$\log_{10} e\f$ where \f$e\f$ is an error of a candidate. It cannot be NULL.
        *  \param[out] dev_error A deviation of \f$\log_{10} e\f$ where \f$e\f$ is an error of a candidate. It cannot be NULL. */
        virtual void evaluate_cands_error(Element* e, Solution<Scalar>* rsln, double* avg_error, double* dev_error) = 0;

        /// Calculates DOF of candidates.
        /** It uses a list of shape indices (OptimumSelector::shape_indices) to
        *  count a number of DOFs. No number of DOFs cannot be zero.
        *  \param[in] e An element that is being refined.
        *  \param[in] rsln A reference solution which is used to calculate the error. */
        virtual void evaluate_cands_dof(Element* e, Solution<Scalar>* rsln);

        /// Evalutes score of candidates.
        /** It calculates score \f$s\f$ of a candidate as \f[s = \frac{\log_{10} e_0 - \log_{10} e}{(d - d_0)^c},\f]
        *  where \f$e\f$ and \f$e_0\f$ are errors of a candidate and the original element (i.e. a candidate at the index 0)
        *  respectively, \f$d\f$ and \f$d_0\f$ are number of DOFs of a candidate and the original element (i.e. a candidate at the index 0)
        *  respectively, and \f$c\f$ is a convergence exponent (OptimumSelector::conv_exp).
        *  If the score is zero, the score is invalid.
        *
        *  If overridden, the higher score the better candidate.
        *  \param[in] e An element that is being refined. */
        virtual void evaluate_cands_score(Element* e);

      private:
        /// Compares scores. Used to sort scores ascending.
        /** \param[in] a The first candidate.
        *  \param[in] b The second candidate.
        *  \return True if score of \a a is greater than the score of \a b. */
        static bool compare_cand_score(const Cand& a, const Cand& b);

      protected: //orders and their range
        int current_max_order; ///< Current maximum order.
        int current_min_order; ///< Current minimum order.

        /// Sets OptimumSelector::current_max_order and OptimumSelector::current_min_order.
        /** This method has to be implemented by derived classes and it is mean to be
        *  space dependent, i.e., it should differ in a case of H1, L2, and Hcurl.
        *  \param[in] element An element that is being refined. */
        virtual void set_current_order_range(Element* element) = 0;

      protected: //shape functions
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

        /// Range of values.
        class Range
        {
        protected:
          int lower_bound;    ///< Lower boundary.
          int upper_bound;    ///< Upper boundary.
          bool empty_range; ///< intrue if range is empty.
          Range();
          Range(const int& lower_bound, const int& upper_bound);
          bool empty() const;
          const int& lower() const;
          const int& upper() const;
          bool is_in_closed(const Range& range) const;
          bool is_in_closed(const int& value) const;
          bool is_in_open(const int& value) const;
          void enlarge_to_include(const int& value);

          static Range make_envelope(const Range& a, const Range& b);
          template<typename T> friend class OptimumSelector;
          template<typename T> friend class ProjBasedSelector;
          template<typename T> friend class H1ProjBasedSelector;
          template<typename T> friend class L2ProjBasedSelector;
          template<typename T> friend class HcurlProjBasedSelector;
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

        /// Constructor.
        /** \note Parameters \a vertex_order and \a edge_bubble_order fixes the fact that a shapeset returns a valid index even though a given shape is not invalid in the space.
        *  \param[in] cand_list A predefined list of candidates.
        *  \param[in] conv_exp A conversion exponent, see evaluate_cands_score().
        *  \param[in] max_order A maximum order which considered. If ::H2DRS_DEFAULT_ORDER, a maximum order supported by the selector is used.
        *  \param[in] shapeset A shapeset. It cannot be NULL.
        *  \param[in] vertex_order A range of orders for vertex functions. Use an empty range (i.e. Range()) to skip vertex functions.
        *  \param[in] edge_bubble_order A range of orders for edge and bubble functions. Use an empty range (i.e. Range()) to skip edge and bubble functions. */
        OptimumSelector(CandList cand_list, double conv_exp, int max_order, Shapeset* shapeset, const Range& vertex_order, const Range& edge_bubble_order);

      public:
        /// Destructor.
        virtual ~OptimumSelector();
      protected:
        /// Selects a refinement.
        /** Overriden function. For details, see Selector::select_refinement(). */
        virtual bool select_refinement(Element* element, int quad_order, Solution<Scalar>* rsln, ElementToRefine& refinement); ///< Selects refinement.

        /// Generates orders of elements which will be created due to a proposed refinement in another component that shares the same a mesh.
        /** Overriden function. For details, see Selector::generate_shared_mesh_orders(). */
        virtual void generate_shared_mesh_orders(const Element* element, const int orig_quad_order, const int refinement, int tgt_quad_orders[H2D_MAX_ELEMENT_SONS], const int* suggested_quad_orders); ///< Updates orders of a refinement in another multimesh component which shares a mesh.
      };
    }
  }
}

#endif