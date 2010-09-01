#include <stdexcept>

#include "stdcython.h"

void throw_exception(char *text)
{
    throw std::runtime_error(text);
}
