IF(H3D_COMPLEX)

project(bessel)

add_executable(${PROJECT_NAME}
	main.cpp
	bessel.cpp
)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES3D_CPLX})

if(WITH_TESTS)
  add_subdirectory(tests)
endif(WITH_TESTS)

ENDIF(H3D_COMPLEX)
