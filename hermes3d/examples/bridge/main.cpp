#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include <hermes3d.h>

// This example shows the Bridge model
// a hexahedral mesh in CTU format.

int main(int argc, char **args) 
{
    // Load the mesh. 
    Mesh mesh;
    //CTUReader mloader;
    H3DReader mloader;

    std::cout << "Loading mesh ...\n" ;
    mloader.load("./ctu2h3d/mesh.mesh3d", &mesh);

    std::cout << "Saving mesh...\n" ;
    mloader.save("new.mesh.mesh3d", &mesh);

    return 0;
}
