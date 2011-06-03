#include "hermes2d.h"
#include "traverse.h"



scalar bilinear_form(RealFunction* fu, RealFunction* fv, RefMap* ru, RefMap* rv)
  { return int_grad_u_grad_v_2(fu, fv, ru, rv); }

scalar linear_form(RealFunction* fv, RefMap* refmap)
  { return int_v(fv, refmap); }


int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);

  H1ShapesetOrtho shapeset;
  PrecalcShapeset pss(&shapeset);

  ///////

  Mesh mesh1;
  mesh1.load("../data/square1.mesh");
  mesh1.refine_element_id(0, 2);
  mesh1.refine_element_id(2, 2);
  /*mesh1.refine_element_id(4, 2);
  mesh1.refine_element_id(6, 2);*/

  H1Space space1(&mesh1, &shapeset);
  space1.set_uniform_order(3);
  space1.assign_dofs();

  EllipticProblem ep;
  ep.set_num_equations(1);
  ep.set_spaces(1, &space1);
  ep.set_pss(1, &pss);
  ep.set_bilinear_form(0, 0, bilinear_form);
  ep.set_linear_form(0, linear_form);

  Solution sln1;
  ep.create_stiffness_matrix();
  ep.assemble_stiffness_matrix_and_rhs();
  ep.solve_system(1, &sln1);

  ScalarView view1("Solution 1");
  view1.show(&sln1);

  ////////

  Mesh mesh2;
  mesh2.load("../data/square1.mesh");
  //mesh2.refine_element_id(0, 1);

  H1Space space2(&mesh2, &shapeset);
  space2.set_uniform_order(3);
  space2.assign_dofs();

  //BaseView bv;
  //bv.show(&space2);

  ep.set_spaces(1, &space2);
  ep.set_external_fns(1, &sln1);

  Solution sln2;
  ep.create_stiffness_matrix();
  ep.assemble_stiffness_matrix_and_rhs();
  ep.save_stiffness_matrix("mat2.m", "mat2");
  ep.solve_system(1, &sln2);

  ScalarView view2("Solution 2");
  view2.show(&sln2);

  ////////

  /*RefMap refmap;
  Quad2DStd quad;

  sln1.set_quad_2d(&quad);
  sln2.set_quad_2d(&quad);
  sln3.set_quad_2d(&quad);
  refmap.set_quad_2d(&quad);

  Mesh* meshes[] = { &mesh1, &mesh2, &mesh3, &mesh1 };
  Transformable* fn[] = { &sln1, &sln2, &sln3, &refmap };

  Traverse trav;
  trav.begin(4, meshes, fn);

  double result = 0.0;
  int o = make_quad_order(10, 10);
  Element** e;
  while ((e = trav.get_next_state()) != NULL)
  {
    sln1.set_quad_order(o, H2D_FN_VAL);
    sln2.set_quad_order(o, H2D_FN_VAL);
    sln3.set_quad_order(o, H2D_FN_VAL);
    double* val1 = sln1.get_fn_values();
    double* val2 = sln2.get_fn_values();
    double* val3 = sln3.get_fn_values();

    double3* pt = quad.get_points(o);
    double2x2* m = refmap.get_inv_ref_map(o);
    double* jac = refmap.get_jacobian(o);
    for (int i = 0; i < quad.get_num_points(o); i++, m++)
      result += pt[i][2] * jac[i] * (val1[i] + val2[i] - val3[i]);
  }

  trav.finish();
  printf("%g\n", result);*/

  hermes2d_finalize();
  return 0;
}
