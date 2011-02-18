project(calc-hex-hcurl-cplx)
add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES3D_CPLX})

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME}-1 ${BIN} hex1.mesh3d 3)
foreach(i RANGE 0 7)
    add_test(${PROJECT_NAME}-2-${i} ${BIN} hex2-ori${i}.mesh3d 3)
endforeach(i RANGE 0 7)
add_test(${PROJECT_NAME}-4  ${BIN} hex4.mesh3d 3)
add_test(${PROJECT_NAME}-4x ${BIN} hex4x.mesh3d 3)
add_test(${PROJECT_NAME}-4y ${BIN} hex4y.mesh3d 3)
add_test(${PROJECT_NAME}-4z ${BIN} hex4z.mesh3d 3)
add_test(${PROJECT_NAME}-8  ${BIN} hex8.mesh3d 3)
