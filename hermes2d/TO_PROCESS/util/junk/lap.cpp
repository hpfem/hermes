#include "hermes2d.h"

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  mesh.load("domain.mesh");
  mesh.refine_element_id(0);

  // create an H1 space
  H1Space space(&mesh);
  space.set_uniform_order(5);

  // set up the weak formulation
  WeakForm wf;
  wf.set_eqn("[u,v] = (2*v)");

  UmfpackSolver solver;
  Solution sln;

  // assemble and solve the linear system
  LinSystem sys(&wf, &solver);
  sys.set_spaces(1, &space);
  sys.assemble();
  sys.solve(1, &sln);

  // show the result
  ScalarView view("Solution");
  view.show(&sln);

  return 0;
}
