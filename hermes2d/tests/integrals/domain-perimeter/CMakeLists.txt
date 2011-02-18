if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)

project(integrals-domain-perimeter)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(integrals-domain-perimeter-1 "${BIN}" domain.mesh0 8.0)
add_test(integrals-domain-perimeter-2 "${BIN}" domain.mesh1 6.650281539872986)
add_test(integrals-domain-perimeter-3 "${BIN}" domain.mesh2 10.681415028389544)
