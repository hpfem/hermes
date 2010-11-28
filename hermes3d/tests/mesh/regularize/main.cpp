#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// Helpers.

int parse_reft(char *str) {
  if (strcasecmp(str, "x") == 0) return H3D_REFT_HEX_X;
  else if (strcasecmp(str, "y") == 0) return H3D_REFT_HEX_Y;
  else if (strcasecmp(str, "z") == 0) return H3D_REFT_HEX_Z;
  else if (strcasecmp(str, "xy") == 0 || strcasecmp(str, "yx") == 0) return H3D_H3D_REFT_HEX_XY;
  else if (strcasecmp(str, "xz") == 0 || strcasecmp(str, "zx") == 0) return H3D_H3D_REFT_HEX_XZ;
  else if (strcasecmp(str, "yz") == 0 || strcasecmp(str, "zy") == 0) return H3D_H3D_REFT_HEX_YZ;
  else if (strcasecmp(str, "xyz") == 0) return H3D_H3D_H3D_REFT_HEX_XYZ;
  else return H3D_REFT_HEX_NONE;
}

int main(int argc, char *args[]) 
{
  // Test variable.
  int success_test = 1;

  if (argc < 1) error("Not enough parameters.");

  Mesh mesh;
  H3DReader mloader;
  if (!mloader.load(args[1], &mesh)) error("Loading mesh file '%s'.", args[1]);

  // Apply refinements.
  bool ok = true;
  for (int k = 2; k < argc && ok; k += 2) {
    int elem_id, reft_id;
    sscanf(args[k], "%d", &elem_id);
    reft_id = parse_reft(args[k + 1]);
    ok = mesh.refine_element(elem_id, reft_id);
  }

  // Testing regularization.
  mesh.regularize();

  mesh.dump();

  return ERR_SUCCESS;
}
