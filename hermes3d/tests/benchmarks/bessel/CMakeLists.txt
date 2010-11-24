project(benchmark-bessel)

add_executable(${PROJECT_NAME}
	main.cpp
	bessel.cpp
)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES_CPLX_BIN})

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME} ${BIN})