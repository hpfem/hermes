#include <fstream>

#include "cpp_api.h"

void CppCallback::event(const char *msg)
{
    printf("CppCallback event(%s) got called.\n", msg);
}
