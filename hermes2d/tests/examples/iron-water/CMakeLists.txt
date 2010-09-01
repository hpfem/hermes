if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)

if (WITH_EXODUSII)

project(example-iron-water)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-iron-water ${BIN})

endif (WITH_EXODUSII)
