#include "hermes2d.h"
#include "solver_pardiso.h"

const double Re = 300;
const double tau = 0.05;


int main(int argc, char* argv[])
{
  Mesh mesh;
  mesh.load("airfoil.mesh");
  mesh.refine_towards_boundary(3, 5);

  H1Space xvel(&mesh);
  xvel.set_simple_bc(11, BC_ESSENTIAL, 1.0);
  xvel.set_simple_bc(1, BC_NONE);
  xvel.set_uniform_order(2);

  H1Space yvel(&mesh);
  yvel.set_simple_bc(1, BC_NONE, 1.0);
  yvel.set_simple_bc(ANY_MARKER, BC_ESSENTIAL, 0.0);
  yvel.set_uniform_order(2);

  L2Space press(&mesh);
  press.set_simple_bc(ANY_MARKER, BC_NONE);
  press.set_uniform_order(1);

  Solution xsln, ysln, press;
  xsln.set_zero(&mesh);
  ysln.set_zero(&mesh);

  WeakForm wf(3);
  wf.set_const(2, "Re", Re, "tau", tau);
  wf.set_base(3, "u1", "u2", "p");
  wf.set_test(3, "v1", "v2", "q");
  wf.set_ext_fn(2, "X", &xprev, "Y", &yprev);
  wf.set_eqn(0, "(u1,v1)/tau + [u1,v1]/Re + (X*u1_x + Y*u1_y, v1) - (p,v1_x) = (X,v1)/tau");
  wf.set_eqn(1, "(u1,v1)/tau + [u2,v2]/Re + (X*u2_x + Y*u2_y, v2) - (p,v2_y) = (Y,v2)/tau");
  wf.set_eqn(2, "(u1_x,q) + (u2_y,q) = 0");

  PardisoSolver solver;
  LinSystem sys(&wf, &solver);
  sys.set_spaces(3, &xvel, &yvel, &press);

  VectorView vview("Velocity");
  ScalarView pview("Pressure");

  for (int i = 1; i <= 1000; i++)
  {
    sys.assemble();
    sys.solve(3, &xsln, &ysln, &press);

    vview.show(&xsln, &ysln);
    pview.show(&press);
  }

  return 0;
}
