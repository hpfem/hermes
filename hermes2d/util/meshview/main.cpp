#include "hermes2d.h"


int main(int argc, char* argv[])
{
  if (argc != 2)
  {
     printf("Usage: meshview mesh-filename\n");
     exit(1);
  }

  Mesh mesh;
  mesh.load(argv[1], true);

  MeshView mv;
  mv.show(&mesh);

  View::wait();
  return 0;
}

