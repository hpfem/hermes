IF(H2D_COMPLEX)

project(screen)

add_executable(${PROJECT_NAME} main.cpp fresnl.c)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})
ENDIF(H2D_COMPLEX)
