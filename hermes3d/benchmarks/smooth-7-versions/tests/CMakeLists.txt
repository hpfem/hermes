if(NOT H3D_REAL)
    return()
endif(NOT H3D_REAL)

project(test-benchmark-adapt-smooth)
add_executable(${PROJECT_NAME}	main.cpp)

include(${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME}-x ${BIN} x)
add_test(${PROJECT_NAME}-y ${BIN} y)
add_test(${PROJECT_NAME}-z ${BIN} z)
add_test(${PROJECT_NAME}-xy ${BIN} xy)
add_test(${PROJECT_NAME}-xz ${BIN} xz)
add_test(${PROJECT_NAME}-yz ${BIN} yz)
add_test(${PROJECT_NAME}-xyz ${BIN} xyz)

set_tests_properties(${PROJECT_NAME}-x PROPERTIES LABELS slow)
set_tests_properties(${PROJECT_NAME}-y PROPERTIES LABELS slow)
set_tests_properties(${PROJECT_NAME}-z PROPERTIES LABELS slow)
set_tests_properties(${PROJECT_NAME}-xy PROPERTIES LABELS slow)
set_tests_properties(${PROJECT_NAME}-xz PROPERTIES LABELS slow)
set_tests_properties(${PROJECT_NAME}-yz PROPERTIES LABELS slow)
set_tests_properties(${PROJECT_NAME}-xyz PROPERTIES LABELS slow)
