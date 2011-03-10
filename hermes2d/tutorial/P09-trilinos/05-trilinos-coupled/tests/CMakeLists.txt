if (NOT H2D_REAL)
    return()
endif (NOT H2D_REAL)

if (WITH_TRILINOS)

project(test-tutorial-P09-05-trilinos-coupled)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-tutorial-P09-05-trilinos-coupled ${BIN})
set_tests_properties(test-tutorial-P09-05-trilinos-coupled PROPERTIES LABELS slow)

endif (WITH_TRILINOS)

