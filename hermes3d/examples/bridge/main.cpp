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
    CTUReader mloader;
    mloader.load("most-sup.top", &mesh);

    mloader.save_as_h3d("mesh.mesh3d", &mesh);
    return 0;
}
