IF(H3D_COMPLEX)

project(maxwell-ring)

add_executable(${PROJECT_NAME} main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES3D_CPLX})

#if(WITH_TESTS)
#  add_subdirectory(tests)
#endif(WITH_TESTS)

ENDIF(H3D_COMPLEX)
