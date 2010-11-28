#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

#define ERROR_SUCCESS                                                           0
#define ERROR_FAILURE                                                           -1

// TODO: Improve the way how to compare two meshes (compare outputs from dump() is not enough).

// Tests themselves.
int test_mesh3d_loader(char *file_name)
{
  _F_
  Mesh mesh;
  H3DReader mloader;
  if (mloader.load(file_name, &mesh)) {
    mesh.dump();
    return ERR_SUCCESS;
  }
  else {
  printf("failed\n");
  return ERR_FAILURE;
  }
}

int test_hdf5_loader(char *file_name)
{
  _F_
  Mesh mesh;
  HDF5Reader mloader;
  if (mloader.load(file_name, &mesh)) {
    mesh.dump();
    return ERR_SUCCESS;
  }
  else {
    printf("failed\n");
    return ERR_FAILURE;
  }
}

int test_exodusii_loader(char *file_name)
{
  _F_
  Mesh mesh;
  ExodusIIReader mloader;
  if (mloader.load(file_name, &mesh)) {
    mesh.dump();
    return ERR_SUCCESS;
  }
  else {
    printf("failed\n");
    return ERR_FAILURE;
  }
}

int main(int argc, char *args[])
{
  _F_
  int ret = ERR_SUCCESS;

  if (argc < 3) error("Not enough parameters.");

  if (strcmp(args[1], "m3d") == 0)
    ret = test_mesh3d_loader(args[2]);
  else if (strcmp(args[1], "hdf5") == 0)
    ret = test_hdf5_loader(args[2]);
  else if (strcmp(args[1], "exoii") == 0)
    ret = test_exodusii_loader(args[2]);
	
  // If the correct behavior is a failure.
  if(argc == 4 && strcmp(args[3], "supp_to_fail") == 0)
    ret = (ret == ERR_FAILURE ? ERR_SUCCESS : ERR_FAILURE);

  if (ret == ERR_SUCCESS) {
    return ERR_SUCCESS;
  }
  else {
    return ERR_FAILURE;
  }
}
