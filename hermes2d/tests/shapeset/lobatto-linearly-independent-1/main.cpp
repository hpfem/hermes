#include "hermes2d.h"
#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int P_INIT = 10;
double EPS = 10e-14;
// This test makes sure the Lobatto shape 
// functions are linearly independent.

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  // We load the mesh on a (-1, 1)^2 domain.
  mloader.load("ref_square.mesh", &mesh);            

  // Create an H1 space with default shapeset,
  // natural BC, and linear elements.
  H1Space space(&mesh, NULL, NULL, P_INIT);
  // The type of element, mesh_mode = 4 means a rectangle element.
  int mesh_mode = 4;
  int n = get_num_dofs(&space);

  printf("\n.........................\n");

  int *fn_idx = new int [n];
  int m = 0;
  int order = P_INIT;

  // Get the vertex fns index.
  info("Get the vertex fns index.");
  for (int i = 0; i < mesh_mode; i++, m++)
  {
    fn_idx[m] = space.get_shapeset()->get_vertex_index(i);
    printf("m = %d, get_vertex_index(m) = %d\n",
          m, space.get_shapeset()->get_vertex_index(m));
  }
  // Get the edge fns index.
  info("Get the edge fns index.");
  for (int edge_order = 2; edge_order <= order; edge_order++)
  {
    for (int j = 0; j < mesh_mode; j++, m++)
    {
      fn_idx[m] = space.get_shapeset()->get_edge_index(j, 0, edge_order);
      printf("m = %d, get_edge_index(m) = %d\n", m, fn_idx[m]);
    }
  }
  // Get the bubble fns index.
  info("Get the bubble fns index.");
  int number_bubble = space.get_shapeset()->get_num_bubbles(H2D_MAKE_QUAD_ORDER(order, order));
  int *bubble_idx = space.get_shapeset()->get_bubble_indices(H2D_MAKE_QUAD_ORDER(order, order));
  for (int i = 0; i < number_bubble; i++, m++ )
  {
    fn_idx[m] = bubble_idx[i];
    printf("m = %d, get_bubble_index(m) = %d\n", m, fn_idx[m]);
  }
  
  info("Assembling matrix ...");
  CooMatrix mat(n);
  Vector* rhs = new AVector(n);
  CommonSolverSciPyUmfpack solver;

  info("ndof = %d", n);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < order; j++)
    {
      double x1 = -1 + (2.0/order)*j;
      double y1 = -1;
      double value = space.get_shapeset()->get_fn_value(fn_idx[i], x1, y1, 0);
      mat.add(i, j, value);
      printf("get fn[%d] value = %f  ", i, value); 
      printf("x1 = %f, y1 = %f\n", x1, y1);
    }
    for (int j = 0; j < order; j++)
    {
      double y2 = -1 + (2.0/order)*j;;
      double x2 = 1;
      double value = space.get_shapeset()->get_fn_value(fn_idx[i], x2, y2, 0);
      mat.add(i, j+order, value);
      printf("get fn[%d] value = %f  ", i, value); 
      printf("x2 = %f, y2 = %f\n", x2, y2);
    }
    for (int j = 0; j < order; j++)
    {
      double x3 = 1 + (-1.0)*(2.0/order)*j;
      double y3 = 1;
      double value = space.get_shapeset()->get_fn_value(fn_idx[i], x3, y3, 0);
      mat.add(i, j+2*order, value);
      printf("get fn[%d] value = %f  ", i, value); 
      printf("x3 = %f, y3 = %f\n", x3, y3);
    }
    for (int j = 0; j < order; j++)
    {
      double x4 = -1;
      double y4 = 1 + (-1.0)*(2.0/order)*j;
      double value = space.get_shapeset()->get_fn_value(fn_idx[i], x4, y4, 0);
      mat.add(i, j+3*order, value);
      printf("get fn[%d] value = %f  ", i, value); 
      printf("x4 = %f, y4 = %f\n", x4, y4);
    }
  }

  int bubble = 0;
  for (int i = order*4; i < n; i++)
  {
      double x = 0.0 + (1.0/number_bubble)*bubble;
      double y = 0.0 + (1.0/number_bubble)*bubble;
      double value = space.get_shapeset()->get_fn_value(fn_idx[i], x, y, 0);
      mat.add(i, i, value);
      printf("get fn[%d] value = %f  ", i, value); 
      printf("x = %f, y = %f\n", x, y);
      bubble++;
  }

  info("ndof = %d", n);

  printf("Adding the rhs\n");
  for (int i = 0; i < n; i++)
  {
    rhs->add(i, 0.0);
  }

  info("Solution ...");
  solver.solve(&mat, rhs);

  for (int i = 0; i < n; i++)
  {
    // Get the value of the matrix solution by calling Vertex::get().
    if (rhs->get(i) >= EPS)
    {
      printf("Shape functions are not linearly independent\n");
      return ERROR_FAILURE; 
    }
  }
  printf("Success!\n");
 
  delete rhs;
  delete [] fn_idx;
  return ERROR_SUCCESS;
}

