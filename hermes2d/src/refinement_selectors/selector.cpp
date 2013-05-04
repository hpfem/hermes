#include "global.h"
#include "selector.h"
#include "candidates.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors
    {
      template<typename Scalar>
      bool HOnlySelector<Scalar>::select_refinement(Element* element, int order, MeshFunction<Scalar>* rsln, ElementToRefine& refinement)
      {
        refinement.split = H2D_REFINEMENT_H;
        refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][0] = 
          refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][1] = 
          refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][2] = 
          refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H][3] = 
          order;
        ElementToRefine::copy_orders(refinement.refinement_polynomial_order, refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_H]);
        return true;
      }

      template<typename Scalar>
      POnlySelector<Scalar>::POnlySelector(int max_order, int order_h_inc, int order_v_inc)
        : Selector<Scalar>(max_order), order_h_inc(order_h_inc), order_v_inc(order_v_inc)
      {
        if(order_h_inc < 0)
          throw Hermes::Exceptions::ValueException("horizontal increase", order_h_inc, 0);
        if(order_v_inc < 0)
          throw Hermes::Exceptions::ValueException("vertical increase", order_v_inc, 0);
      }

      template<typename Scalar>
      bool POnlySelector<Scalar>::select_refinement(Element* element, int order, MeshFunction<Scalar>* rsln, ElementToRefine& refinement)
      {
        refinement.split = H2D_REFINEMENT_P;

        //determin max. order
        int max_allowed_order = this->max_order;
        if(this->max_order == H2DRS_DEFAULT_ORDER)
          max_allowed_order = H2DRS_MAX_ORDER;

        //calculate new order
        int order_h = H2D_GET_H_ORDER(order), order_v = H2D_GET_V_ORDER(order);
        int new_order_h = std::min(max_allowed_order, order_h + order_h_inc);
        int new_order_v = std::min(max_allowed_order, order_v + order_v_inc);
        if(element->is_triangle())
          refinement.refinement_polynomial_order[0] = refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_P][0] = new_order_h;
        else
          refinement.refinement_polynomial_order[0] = refinement.best_refinement_polynomial_order_type[H2D_REFINEMENT_P][0] = H2D_MAKE_QUAD_ORDER(new_order_h, new_order_v);

        //decide if successful
        if(new_order_h > order_h || new_order_v > order_v)
          return true;
        else
          return false;
      }

      template class HERMES_API Selector<double>;
      template class HERMES_API Selector<std::complex<double> >;
      template class HERMES_API HOnlySelector<double>;
      template class HERMES_API HOnlySelector<std::complex<double> >;
      template class HERMES_API POnlySelector<double>;
      template class HERMES_API POnlySelector<std::complex<double> >;
    }
  }
}