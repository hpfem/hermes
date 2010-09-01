#include <stdio.h>

#include "hermes2d/_hermes2d_api.h"

int main(int argc, char **argv)
{

    Py_Initialize();
    PySys_SetArgv(argc, argv);
    if (import_hermes2d___hermes2d())
        printf("bad\n");
    else
        printf("imported\n");
    set_trace(5);
    printf("ok\n");
}
