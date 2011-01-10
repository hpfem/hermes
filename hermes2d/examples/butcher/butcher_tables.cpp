// First-order.
ButcherTable HBT_Explicit_RK_1;
ButcherTable HBT_Implicit_RK_1;

// Second-order.
ButcherTable HBT_Explicit_RK_2;
ButcherTable HBT_Implicit_Crank_Nicolson_2;
ButcherTable HBT_Implicit_SDIRK_2;
ButcherTable HBT_Implicit_Lobatto_IIIA_2;
ButcherTable HBT_Implicit_Lobatto_IIIB_2;
ButcherTable HBT_Implicit_Lobatto_IIIC_2;

// Third-order.
ButcherTable HBT_Explicit_RK_3;

// Fourth-order.
ButcherTable HBT_Explicit_RK_4;
ButcherTable HBT_Implicit_Lobatto_IIIA_4;
ButcherTable HBT_Implicit_Lobatto_IIIB_4;
ButcherTable HBT_Implicit_Lobatto_IIIC_4;

void init_butcher_tables() 
{
  // Explicit Euler.
  HBT_Explicit_RK_1.alloc(1);
  HBT_Explicit_RK_1.set_B(0, 1.);

  // Implicit Euler.
  HBT_Implicit_RK_1.alloc(1);
  HBT_Implicit_RK_1.set_A(0, 0, 1.);
  HBT_Implicit_RK_1.set_B(0, 1.);
  HBT_Implicit_RK_1.set_C(0, 1.);

  // Explicit RK-2.
  HBT_Explicit_RK_2.alloc(2);
  HBT_Explicit_RK_2.set_A(1, 0, 2./3.);
  HBT_Explicit_RK_2.set_A(1, 1, 0.);
  HBT_Explicit_RK_2.set_B(0, 1./4.);
  HBT_Explicit_RK_2.set_B(1, 3./4.);
  HBT_Explicit_RK_2.set_C(0, 2./3.);

  // Implicit Crank Nicolson.
  HBT_Implicit_Crank_Nicolson_2.alloc(2);
  HBT_Implicit_Crank_Nicolson_2.set_A(0, 0, 1./2.);
  HBT_Implicit_Crank_Nicolson_2.set_A(0, 1, 1./2.);
  HBT_Implicit_Crank_Nicolson_2.set_B(0, 1./2.);
  HBT_Implicit_Crank_Nicolson_2.set_B(1, 1./2.);
  HBT_Implicit_Crank_Nicolson_2.set_C(0, 1.);

  // Implicit Lobatto IIIA (second-order).
  HBT_Implicit_Lobatto_IIIA_2.alloc(2);
  HBT_Implicit_Lobatto_IIIA_2.set_A(1, 0, 1./2.);
  HBT_Implicit_Lobatto_IIIA_2.set_A(1, 1, 1./2.);
  HBT_Implicit_Lobatto_IIIA_2.set_B(0, 1./2.);
  HBT_Implicit_Lobatto_IIIA_2.set_B(1, 1./2.);
  HBT_Implicit_Lobatto_IIIA_2.set_C(1, 1.);

  // Implicit Lobatto IIIB (second-order).
  HBT_Implicit_Lobatto_IIIB_2.alloc(2);
  HBT_Implicit_Lobatto_IIIB_2.set_A(0, 0, 1./2.);
  HBT_Implicit_Lobatto_IIIB_2.set_A(0, 1, 1./2.);
  HBT_Implicit_Lobatto_IIIB_2.set_B(0, 1./2.);
  HBT_Implicit_Lobatto_IIIB_2.set_B(1, 1./2.);
  HBT_Implicit_Lobatto_IIIB_2.set_C(0, 1./2.);
  HBT_Implicit_Lobatto_IIIB_2.set_C(1, 1./2.);

  // Implicit Lobatto IIIC (second-order).
  HBT_Implicit_Lobatto_IIIC_2.alloc(2);
  HBT_Implicit_Lobatto_IIIC_2.set_A(0, 0, 1./2.);
  HBT_Implicit_Lobatto_IIIC_2.set_A(0, 1, -1./2.);
  HBT_Implicit_Lobatto_IIIC_2.set_A(1, 0, 1./2.);
  HBT_Implicit_Lobatto_IIIC_2.set_A(1, 1, 1./2.);
  HBT_Implicit_Lobatto_IIIC_2.set_B(0, 1./2.);
  HBT_Implicit_Lobatto_IIIC_2.set_B(1, 1./2.);
  HBT_Implicit_Lobatto_IIIC_2.set_C(0, 0.);
  HBT_Implicit_Lobatto_IIIC_2.set_C(1, 1.);

  // Implicit SDIRK-2 (second-order).
  HBT_Implicit_SDIRK_2.alloc(2);
  double gamma = 1./sqrt(2.);
  HBT_Implicit_SDIRK_2.set_A(0, 0, 1. - gamma);
  HBT_Implicit_SDIRK_2.set_A(0, 1, 0.);
  HBT_Implicit_SDIRK_2.set_A(1, 0, gamma);
  HBT_Implicit_SDIRK_2.set_A(1, 1, 1. - gamma);
  HBT_Implicit_SDIRK_2.set_B(0, gamma);
  HBT_Implicit_SDIRK_2.set_B(1, 1. - gamma);
  HBT_Implicit_SDIRK_2.set_C(0, 1. - gamma);
  HBT_Implicit_SDIRK_2.set_C(1, 1.);  

  // Explicit RK-3.
  HBT_Explicit_RK_3.alloc(3);
  HBT_Explicit_RK_3.set_A(1, 0, 1./2.);
  HBT_Explicit_RK_3.set_A(2, 0, -1.);
  HBT_Explicit_RK_3.set_A(2, 1, 2.);
  HBT_Explicit_RK_3.set_B(0, 1./6.);
  HBT_Explicit_RK_3.set_B(1, 2./3.);
  HBT_Explicit_RK_3.set_B(2, 1./6.);
  HBT_Explicit_RK_3.set_C(1, 1./2.);
  HBT_Explicit_RK_3.set_C(2, 1.);

  // Explicit RK-4.
  HBT_Explicit_RK_4.alloc(4);
  HBT_Explicit_RK_4.set_A(0, 1, 1./2.);
  HBT_Explicit_RK_4.set_A(1, 2, 1./2.);
  HBT_Explicit_RK_4.set_A(2, 3, 1.);
  HBT_Explicit_RK_4.set_B(0, 1./6.);
  HBT_Explicit_RK_4.set_B(1, 1./3.);
  HBT_Explicit_RK_4.set_B(2, 1./3.);
  HBT_Explicit_RK_4.set_B(3, 1./6.);
  HBT_Explicit_RK_4.set_C(1, 1./2.);
  HBT_Explicit_RK_4.set_C(2, 1./2.);
  HBT_Explicit_RK_4.set_C(3, 1.);

  // Implicit Lobatto IIIA (fourth-order).
  HBT_Implicit_Lobatto_IIIA_4.alloc(3);
  HBT_Implicit_Lobatto_IIIA_4.set_A(1, 0, 5./24.);
  HBT_Implicit_Lobatto_IIIA_4.set_A(2, 0, 1./6.);
  HBT_Implicit_Lobatto_IIIA_4.set_A(1, 1, 1./3.);
  HBT_Implicit_Lobatto_IIIA_4.set_A(2, 1, 2./3.);
  HBT_Implicit_Lobatto_IIIA_4.set_A(1, 2, -1./24.);
  HBT_Implicit_Lobatto_IIIA_4.set_A(2, 2, 1./6.);
  HBT_Implicit_Lobatto_IIIA_4.set_B(0, 1./6.);
  HBT_Implicit_Lobatto_IIIA_4.set_B(1, 2./3.);
  HBT_Implicit_Lobatto_IIIA_4.set_B(2, 1./6.);
  HBT_Implicit_Lobatto_IIIA_4.set_C(1, 1./2.);
  HBT_Implicit_Lobatto_IIIA_4.set_C(2, 1.);

  // Implicit Lobatto IIIB (fourth-order).
  HBT_Implicit_Lobatto_IIIB_4.alloc(3);
  HBT_Implicit_Lobatto_IIIB_4.set_A(0, 0, 1./6.);
  HBT_Implicit_Lobatto_IIIB_4.set_A(1, 0, 1./6.);
  HBT_Implicit_Lobatto_IIIB_4.set_A(2, 0, 1./6.);
  HBT_Implicit_Lobatto_IIIB_4.set_A(0, 1, -1./6.);
  HBT_Implicit_Lobatto_IIIB_4.set_A(1, 1, 1./3.);
  HBT_Implicit_Lobatto_IIIB_4.set_A(2, 1, 5./6.);
  HBT_Implicit_Lobatto_IIIB_4.set_B(0, 1./6.);
  HBT_Implicit_Lobatto_IIIB_4.set_B(1, 2./3.);
  HBT_Implicit_Lobatto_IIIB_4.set_B(2, 1./6.);
  HBT_Implicit_Lobatto_IIIB_4.set_C(1, 1./2.);
  HBT_Implicit_Lobatto_IIIB_4.set_C(2, 1.);

  // Implicit Lobatto IIIC (fourth-order).
  HBT_Implicit_Lobatto_IIIC_4.alloc(3);
  HBT_Implicit_Lobatto_IIIC_4.set_A(0, 0, 1./6.);
  HBT_Implicit_Lobatto_IIIC_4.set_A(1, 0, 1./6.);
  HBT_Implicit_Lobatto_IIIC_4.set_A(2, 0, 1./6.);
  HBT_Implicit_Lobatto_IIIC_4.set_A(0, 1, -1./3.);
  HBT_Implicit_Lobatto_IIIC_4.set_A(1, 1, 5./12.);
  HBT_Implicit_Lobatto_IIIC_4.set_A(2, 1, 2./3.);
  HBT_Implicit_Lobatto_IIIC_4.set_A(0, 2, 1./6.);
  HBT_Implicit_Lobatto_IIIC_4.set_A(1, 2, -1./12.);
  HBT_Implicit_Lobatto_IIIC_4.set_A(2, 2, 1./6.);
  HBT_Implicit_Lobatto_IIIC_4.set_B(0, 1./6.);
  HBT_Implicit_Lobatto_IIIC_4.set_B(1, 2./3.);
  HBT_Implicit_Lobatto_IIIC_4.set_B(2, 1./6.);
  HBT_Implicit_Lobatto_IIIC_4.set_C(1, 1./2.);
  HBT_Implicit_Lobatto_IIIC_4.set_C(2, 1.);
}

