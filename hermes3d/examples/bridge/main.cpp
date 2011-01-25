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
    H3DReader mloader;
    mloader.load("mesh.mesh3d", &mesh);

    mloader.save("new.mesh.mesh3d", &mesh);
    return 0;
}
