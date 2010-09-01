#include "hermes2d.h"


int main(int argc, char* argv[])
{
  for (int i = 1; i < argc; i++)
  {
    Mesh mesh;
    mesh.load_old(argv[i]);
    printf("Converting %s ...\n", argv[i]);
    mesh.save(argv[i]);
  }

  return 0;
}
