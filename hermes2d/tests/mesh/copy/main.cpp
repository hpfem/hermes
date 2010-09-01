#include "hermes2d.h"

// This test make sure that a copy of another mesh are correct,
// including the number of elements, the type of elements,
// and find out if the elements are curvilinear.

#undef ERROR_SUCCESS
#undef ERROR_FAILURE
#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    printf("please input as this format: copy  meshfile.mesh  number  \n");
    printf("number 1 : copy\n ");
    printf("number 2 : copy_base\n ");
    printf("number 3 : copy_refine\n ");
    return ERROR_FAILURE;
  }
  // load the mesh file
  Mesh mesh;
  Mesh dup;
  Element* e;
  H2DReader mloader;
  mloader.load(argv[1], &mesh);
  int mesh_base_element = mesh.get_num_base_elements();

  mesh.refine_all_elements();
  // Calculate the number of elements after refinement, starting from 0
  printf("Elements (count =  %d, the refined mesh elements)\n", mesh.get_max_element_id());

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

  switch (atoi(argv[2]))
  {
    // Creates a copy of another mesh.
    case 1:      dup.copy(&mesh);
    if (mesh.get_max_element_id() != dup.get_max_element_id())
    {
      printf("Failure!\n");
      return ERROR_FAILURE;
    }break;

    // Copies the coarsest elements of another mesh.
    case 2:      dup.copy_base(&mesh);
    if (mesh.get_num_base_elements() != dup.get_max_element_id())
    {
      printf("Failure!\n");
      return ERROR_FAILURE;
    }break;

    // Copies the refined active elements of another mesh.
    case 3:      dup.copy_converted(&mesh);
    if (mesh.get_max_element_id() != dup.get_max_element_id() + mesh_base_element)
    {
      printf("Failure!\n");
      return ERROR_FAILURE;
    }break;
    default: break;
  }

  printf("Elements (count =  %d, the copyed mesh elements)\n", dup.get_max_element_id());
  for_all_elements(e, &dup)
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

