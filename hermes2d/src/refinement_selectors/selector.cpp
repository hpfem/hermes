#include "global.h"
#include "solution.h"
#include "element_to_refine.h"
#include "selector.h"
namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors
    {
      template<typename Scalar>
      Selector<Scalar>* HOnlySelector<Scalar>::clone()
      {
        return new HOnlySelector();
      }

      template<typename Scalar>
      bool HOnlySelector<Scalar>::select_refinement(Element* element, int quad_order, Solution<Scalar>* rsln, ElementToRefine& refinement)
      {
        refinement.split = H2D_REFINEMENT_H;
        refinement.p[0] = refinement.p[1] = refinement.p[2] = refinement.p[3] = quad_order;
        refinement.q[0] = refinement.q[1] = refinement.q[2] = refinement.q[3] = quad_order;
        return true;
      }

      template<typename Scalar>
      void HOnlySelector<Scalar>::generate_shared_mesh_orders(const Element* element, const int orig_quad_order, const int refinement, int tgt_quad_orders[H2D_MAX_ELEMENT_SONS], const int* suggested_quad_orders)
      {
        if(suggested_quad_orders != NULL)
          for(int i = 0; i < H2D_MAX_ELEMENT_SONS; i++)
            tgt_quad_orders[i] = suggested_quad_orders[i];
        else
          for(int i = 0; i < H2D_MAX_ELEMENT_SONS; i++)
            tgt_quad_orders[i] = orig_quad_order;
      }

      template<typename Scalar>
      Selector<Scalar>* POnlySelector<Scalar>::clone()
      {
        return new POnlySelector(this->max_order, this->order_h_inc, this->order_v_inc);
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
      bool POnlySelector<Scalar>::select_refinement(Element* element, int quad_order, Solution<Scalar>* rsln, ElementToRefine& refinement)
      {
        refinement.split = H2D_REFINEMENT_P;

        //determin max. order
        int max_allowed_order = this->max_order;
        if(this->max_order == H2DRS_DEFAULT_ORDER)
          max_allowed_order = H2DRS_MAX_ORDER;

        //calculate new order
        int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);
        int new_order_h = std::min(max_allowed_order, order_h + order_h_inc);
        int new_order_v = std::min(max_allowed_order, order_v + order_v_inc);
        if(element->is_triangle())
          refinement.p[0] = refinement.q[0] = new_order_h;
        else
          refinement.p[0] = refinement.q[0] = H2D_MAKE_QUAD_ORDER(new_order_h, new_order_v);

        //decide if successful
        if(new_order_h > order_h || new_order_v > order_v)
          return true;
        else
          return false;
      }

      template<typename Scalar>
      void POnlySelector<Scalar>::generate_shared_mesh_orders(const Element* element, const int orig_quad_order, const int refinement, int tgt_quad_orders[H2D_MAX_ELEMENT_SONS], const int* suggested_quad_orders)
      {
        if(suggested_quad_orders != NULL)
          tgt_quad_orders[0] = suggested_quad_orders[0];
        else
          tgt_quad_orders[0] = orig_quad_order;
#ifdef _DEBUG
        for(int i = 1; i < H2D_MAX_ELEMENT_SONS; i++)
          tgt_quad_orders[i] = 0;
#endif
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