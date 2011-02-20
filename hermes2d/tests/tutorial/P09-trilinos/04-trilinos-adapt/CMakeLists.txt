if (NOT H2D_REAL)
    return()
endif (NOT H2D_REAL)

if (WITH_TRILINOS)

project(tutorial-P09-04-trilinos-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(tutorial-P09-04-trilinos-adapt ${BIN})
set_tests_properties(tutorial-P09-04-trilinos-adapt PROPERTIES LABELS slow)

endif (WITH_TRILINOS)

