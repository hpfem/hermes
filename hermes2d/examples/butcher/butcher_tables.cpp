#include "hermes2d.h"

void init_butcher_table(ButcherTableType butcher_table, ButcherTable* bt) 
{
  double gamma = 1./sqrt(2.);

  switch (butcher_table) {
  case Explicit_RK_1: // Explicit Euler.
    bt->alloc(1);
    bt->set_B(0, 1.);
  break;

  case Implicit_RK_1: // Implicit Euler.
    bt->alloc(1);
    bt->set_A(0, 0, 1.);
    bt->set_B(0, 1.);
    bt->set_C(0, 1.);
  break;

  case Explicit_RK_2: // Explicit RK-2.
    bt->alloc(2);
    bt->set_A(1, 0, 2./3.);
    bt->set_A(1, 1, 0.);
    bt->set_B(0, 1./4.);
    bt->set_B(1, 3./4.);
    bt->set_C(0, 2./3.);
  break;

  case Implicit_Crank_Nicolson_2: // Implicit Crank Nicolson.
    bt->alloc(2);
    bt->set_A(0, 0, 1./2.);
    bt->set_A(0, 1, 1./2.);
    bt->set_B(0, 1./2.);
    bt->set_B(1, 1./2.);
    bt->set_C(0, 1.);
  break;

  case Implicit_SDIRK_2: // Implicit SDIRK-2 (second-order).
    bt->alloc(2);
    bt->set_A(0, 0, 1. - gamma);
    bt->set_A(0, 1, 0.);
    bt->set_A(1, 0, gamma);
    bt->set_A(1, 1, 1. - gamma);
    bt->set_B(0, gamma);
    bt->set_B(1, 1. - gamma);
    bt->set_C(0, 1. - gamma);
    bt->set_C(1, 1.);  
  break;

  case Implicit_Lobatto_IIIA_2: // Implicit Lobatto IIIA (second-order).
    bt->alloc(2);
    bt->set_A(1, 0, 1./2.);
    bt->set_A(1, 1, 1./2.);
    bt->set_B(0, 1./2.);
    bt->set_B(1, 1./2.);
    bt->set_C(1, 1.);
  break;

  case Implicit_Lobatto_IIIB_2: // Implicit Lobatto IIIB (second-order).
    bt->alloc(2);
    bt->set_A(0, 0, 1./2.);
    bt->set_A(0, 1, 1./2.);
    bt->set_B(0, 1./2.);
    bt->set_B(1, 1./2.);
    bt->set_C(0, 1./2.);
    bt->set_C(1, 1./2.);
  break;

  case Implicit_Lobatto_IIIC_2: // Implicit Lobatto IIIC (second-order).
    bt->alloc(2);
    bt->set_A(0, 0, 1./2.);
    bt->set_A(0, 1, -1./2.);
    bt->set_A(1, 0, 1./2.);
    bt->set_A(1, 1, 1./2.);
    bt->set_B(0, 1./2.);
    bt->set_B(1, 1./2.);
    bt->set_C(0, 0.);
    bt->set_C(1, 1.);
  break;

  case Explicit_RK_3: // Explicit RK-3.
    bt->alloc(3);
    bt->set_A(1, 0, 1./2.);
    bt->set_A(2, 0, -1.);
    bt->set_A(2, 1, 2.);
    bt->set_B(0, 1./6.);
    bt->set_B(1, 2./3.);
    bt->set_B(2, 1./6.);
    bt->set_C(1, 1./2.);
    bt->set_C(2, 1.);
  break;

  case Explicit_RK_4: // Explicit RK-4.
    bt->alloc(4);
    bt->set_A(1, 0, 1./2.);
    bt->set_A(2, 1, 1./2.);
    bt->set_A(3, 2, 1.);
    bt->set_B(0, 1./6.);
    bt->set_B(1, 1./3.);
    bt->set_B(2, 1./3.);
    bt->set_B(3, 1./6.);
    bt->set_C(1, 1./2.);
    bt->set_C(2, 1./2.);
    bt->set_C(3, 1.);
  break;

  case Implicit_Lobatto_IIIA_4: // Implicit Lobatto IIIA (fourth-order).
    bt->alloc(3);
    bt->set_A(1, 0, 5./24.);
    bt->set_A(2, 0, 1./6.);
    bt->set_A(1, 1, 1./3.);
    bt->set_A(2, 1, 2./3.);
    bt->set_A(1, 2, -1./24.);
    bt->set_A(2, 2, 1./6.);
    bt->set_B(0, 1./6.);
    bt->set_B(1, 2./3.);
    bt->set_B(2, 1./6.);
    bt->set_C(1, 1./2.);
    bt->set_C(2, 1.);
  break;

  case Implicit_Lobatto_IIIB_4: // Implicit Lobatto IIIB (fourth-order).
    bt->alloc(3);
    bt->set_A(0, 0, 1./6.);
    bt->set_A(1, 0, 1./6.);
    bt->set_A(2, 0, 1./6.);
    bt->set_A(0, 1, -1./6.);
    bt->set_A(1, 1, 1./3.);
    bt->set_A(2, 1, 5./6.);
    bt->set_B(0, 1./6.);
    bt->set_B(1, 2./3.);
    bt->set_B(2, 1./6.);
    bt->set_C(1, 1./2.);
    bt->set_C(2, 1.);
  break;

  case Implicit_Lobatto_IIIC_4: // Implicit Lobatto IIIC (fourth-order).
    bt->alloc(3);
    bt->set_A(0, 0, 1./6.);
    bt->set_A(1, 0, 1./6.);
    bt->set_A(2, 0, 1./6.);
    bt->set_A(0, 1, -1./3.);
    bt->set_A(1, 1, 5./12.);
    bt->set_A(2, 1, 2./3.);
    bt->set_A(0, 2, 1./6.);
    bt->set_A(1, 2, -1./12.);
    bt->set_A(2, 2, 1./6.);
    bt->set_B(0, 1./6.);
    bt->set_B(1, 2./3.);
    bt->set_B(2, 1./6.);
    bt->set_C(1, 1./2.);
    bt->set_C(2, 1.);
  break;

  default: error("Unknown Butcher's table.");
  }
}
