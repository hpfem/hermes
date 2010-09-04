if (NOT H2D_REAL)
    return()
endif (NOT H2D_REAL)

if (WITH_TRILINOS)

project(tutorial-43-trilinos-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-43-trilinos-adapt ${BIN})

endif (WITH_TRILINOS)

