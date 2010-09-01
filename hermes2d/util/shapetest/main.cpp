#include "hermes2d.h"

#ifndef H2D_COMPLEX

//H1ShapesetOrtho shapeset;
H1ShapesetBeuchler shapeset;

#else

HcurlShapesetLegendre shapeset;
//HcurlShapesetGradLeg shapeset;

#endif


PrecalcShapeset precalc(&shapeset);
Quad2DStd quad;

int nc = shapeset.get_num_components();

const double eps = 1e-13;
inline bool eq(double a, double b) { return fabs(a - b) < eps; }
//inline bool eq(double a, double b) { return a == b; }


double2x2 rot[2][4] =
{
  { {{1,0},{0,1}}, {{-0.5,0.5},{-0.5,-0.5}}, {{0,-1},{1,0}} },
  { {{1,0},{0,1}}, {{0,1},{-1,0}}, {{-1,0},{0,-1}}, {{0,-1},{1,0}} }
};


void test_edge_rotation()
{
  info("Testing edge rotation...");
  int mode = shapeset.get_mode();
  int ne = mode ? 4 : 3;

  for (int ori = 0; ori <= 1; ori++)
  {
    for (int order = 0; order <= shapeset.get_max_order(); order++)
    {
      double *e01, *e02, *ee1, *ee2;
      precalc.set_active_shape(shapeset.get_edge_index(0, ori, order));
      precalc.set_quad_order(quad.get_edge_points(0));
      e01 = precalc.get_fn_values(0);
      if (nc > 1) e02 = precalc.get_fn_values(1);

      for (int e = 1; e < ne; e++)
      {
        precalc.set_active_shape(shapeset.get_edge_index(e, ori, order));
        precalc.set_quad_order(quad.get_edge_points(e));
        ee1 = precalc.get_fn_values(0);
        if (nc > 1) ee2 = precalc.get_fn_values(1);

        int np = quad.get_num_points(quad.get_edge_points(0));
        if (nc == 1)
        {
          for (int i = 0; i < np; i++)
            if (!eq(e01[i], ee1[i]))
            {
              info("order=%d, ori=%d, edge=%d -- not equal to edge 0", order, ori, e);
            }
        }
        else
        {
          for (int i = 0; i < np; i++)
          {
            double x = rot[mode][e][0][0] * ee1[i] + rot[mode][e][0][1] * ee2[i];
            double y = rot[mode][e][1][0] * ee1[i] + rot[mode][e][1][1] * ee2[i];
            if (!eq(e01[i], x) || !eq(e02[i], y))
            {
              info("order=%d, ori=%d, edge=%d -- not equal to edge 0", order, ori, e);
              printf("x comp: 0-ta %g, %d-ta %g\n", e01[i], e, x);
              printf("y comp: 0-ta %g, %d-ta %g\n\n", e02[i], e, y);
            }
          }
        }
      }
    }
  }
}




int main(int argc, char* argv[])
{
  hermes2d_initialize(&argc, argv);

  info("SHAPESET TESTER");
  info("num_components = %d", shapeset.get_num_components());
  info("max_order = %d", shapeset.get_max_order());

  precalc.set_quad_2d(&quad);
  for (int mode = 0; mode <= 1; mode++)
  {
    shapeset.set_mode(mode);
    quad.set_mode(mode);
    precalc.set_mode(mode);
    info(mode ? "\nTESTING QUADS\n" : "\nTESTING TRIANGLES\n");

    //test_orders(&shapeset);
    test_edge_rotation();
    //test_edge_orientation(&shapeset);
    //test_num_bubbles(&shapeset);
  }

  //info("\nALL OK!\n");

  printf("\n");
  hermes2d_finalize();
  return 0;
}
