#include "../h2d_common.h"
#include "order_permutator.h"

namespace RefinementSelectors {

  OrderPermutator::OrderPermutator(int start_quad_order, int end_quad_order, bool iso_p, int* tgt_quad_order)
    : start_order_h(H2D_GET_H_ORDER(start_quad_order)), start_order_v(H2D_GET_V_ORDER(start_quad_order))
    , end_order_h(H2D_GET_H_ORDER(end_quad_order)), end_order_v(H2D_GET_V_ORDER(end_quad_order))
    , iso_p(iso_p), tgt_quad_order(tgt_quad_order) {
    assert_msg(start_order_h <= end_order_h && start_order_v <= end_order_v, "End orders (H:%d, V:%d) are below start orders (H:%d, V:%d).", end_order_h, end_order_v, start_order_h, start_order_v);
    reset();
  }

  bool OrderPermutator::next() {
    if (iso_p) {
      if (order_h >= end_order_h || order_v >= end_order_v)
        return false;

      order_h++; order_v++;
    }
    else {
      if (order_h >= end_order_h && order_v >= end_order_v)
        return false;

      order_h++;
      if (order_h > end_order_h) {
        order_h = start_order_h;
        order_v++;
      }
    }

    if (tgt_quad_order != NULL)
      *tgt_quad_order = H2D_MAKE_QUAD_ORDER(order_h, order_v);
    return true;
  }

  void OrderPermutator::reset() {
    order_h = start_order_h;
    order_v = start_order_v;
    if (tgt_quad_order != NULL)
      *tgt_quad_order = H2D_MAKE_QUAD_ORDER(order_h, order_v);
  }

}
