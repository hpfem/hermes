#include "hermes2d.h"

// This test make sure that all elements are refined correct,
// including the number of elements, the type of elements,
// and find out if the elements are curvilinear.

#undef ERROR_SUCCESS
#undef ERROR_FAILURE
#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    printf("please input as this format: refinements meshfile.mesh \n");
    return ERROR_FAILURE;
  }

  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  Element* e;
  mloader.load(argv[1], &mesh);

  // default refine type '0': one quad to four quads
  //         refine type '1': one quad to two 'horizontal' quads
  //         refine type '2': one quad to two 'vertical' quads
  int refine_type = 0;
  printf("refine type %d\n", refine_type);

  // Calculate the number of elements after refinement, starting from 0
  int element_num = 0;
  printf("Elements (count =  %d, the coarse mesh elements)\n", mesh.get_max_element_id());

  for_all_elements(e, &mesh)
  {
    printf("e->id = %d  ", e->id);
    if (e->is_quad())
    {
      printf("type : quadrangle ");
      if (refine_type == 0)
      {
        // refine type '0'
        element_num += 4;
      }
      else
      {
        // refine type '1','2'
        element_num += 2;
      }
    }
    else
    {
      printf("type : triangle   ");
      // one triangle element refined into four elements
      element_num += 4;
    }
    if (e->is_curved())
      printf("  curved\n");
    else
      printf("\n");
  }

  // Refine all elements.
  mesh.refine_all_elements(refine_type);

  if (element_num != mesh.get_max_element_id() - mesh.get_num_base_elements())
  {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }

  printf("Elements (count =  %d, including the coarse mesh elements)\n", mesh.get_max_element_id());
  for_all_elements(e, &mesh)
  {
    if( e->id >= mesh.get_num_base_elements() )
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
  }
  printf("Success!\n");
  return ERROR_SUCCESS;
}

