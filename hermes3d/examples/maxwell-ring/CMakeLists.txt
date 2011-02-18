project(maxwell-ring)

add_executable(${PROJECT_NAME}
	main.cpp
)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES3D_CPLX})
