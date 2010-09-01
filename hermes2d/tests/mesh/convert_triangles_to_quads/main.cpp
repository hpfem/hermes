#include "hermes2d.h"

// This example shows how to convert triangle elements into quadrangular elements.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    printf("please input as this format: convert_to_quads  meshfile.mesh \n");
    return ERROR_FAILURE;
  }

  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  Element* e;
  mloader.load(argv[1], &mesh);

  // Calculate the number of elements after refinement, starting from 0
  int element_num = 0;
  printf("The number of base elements is %d\n", mesh.get_max_element_id());
  for_all_elements(e, &mesh)
  {
    printf("e->id = %d  ", e->id);
    if (e->is_quad())
    {
      printf("type : quadrangle ");
      // one quad element refined into four elements
      element_num += 1;
    }
    else
    {
      printf("type : triangle   ");
      // one triangle element refined into three elements
      element_num += 3;
    }
    if (e->is_curved())
      printf("  curved\n");
    else
      printf("\n");
  }

  // convert the mesh
  mesh.convert_triangles_to_quads();
  if (element_num != mesh.get_max_element_id())
  {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }

  printf("The number of refined elements is %d\n", mesh.get_max_element_id());
  for_all_elements(e, &mesh)
  {
    printf("e->id = %d  ", e->id);
    if (e->is_quad())
    {
      printf("type : quadrangle ");
    }
    else
    {
      printf("type : triangle   ");
    }
    if (e->is_curved())
      printf("  curved\n");
    else
      printf("\n");
  }
  printf("Success!\n");
  return ERROR_SUCCESS;
}
