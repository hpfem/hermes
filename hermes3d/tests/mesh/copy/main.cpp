#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// This test makes sure that meshes are copied correctly.

enum ECopyType {
	CT_NONE,
	CT_COPY,
	CT_COPY_BASE
};

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

ECopyType parse_copy_type(char *str) {
	if (strcasecmp(str, "copy") == 0) return CT_COPY;
	else if (strcasecmp(str, "copy_base") == 0) return CT_COPY_BASE;
	else return CT_NONE;
}


int main(int argc, char **args) 
{
  // Test variable.
  int success_test = 1;

  if (argc < 2) error("Not enough parameters.");

  ECopyType copy_type = parse_copy_type(args[1]);
  if (copy_type == CT_NONE) error("Unknown copy type (%s).\n", args[1]);

  Mesh mesh;
  H3DReader mloader;
  if (!mloader.load(args[2], &mesh)) error("Loading mesh file '%s'.", args[2]);

  // Apply refinements.
  bool ok = true;
  for (int k = 3; k < argc && ok; k += 2) {
    int elem_id, reft_id;
    sscanf(args[k], "%d", &elem_id);
    reft_id = parse_reft(args[k + 1]);
    ok = mesh.refine_element(elem_id, reft_id);
  }

  if (ok) {
    Mesh dup;
    switch (copy_type) {
      case CT_COPY:      dup.copy(mesh); break;
      case CT_COPY_BASE: dup.copy_base(mesh); break;
      default: break;
    }
    dup.dump();
  }
  else {
    warning("Unable to refine a mesh.");
    success_test = 0;
  }

  if (success_test) {
    return ERR_SUCCESS;
  }
  else {
    return ERR_FAILURE;
  }
}
