project(calc-hex-hcurl)

include(CMake.vars OPTIONAL)

add_executable(${PROJECT_NAME}
	hex-hcurl.cpp
	${HERMES_COMMON_DIR}/timer.cpp
	${HERMES_COMMON_DIR}/error.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
include_directories(${hermes3d_SOURCE_DIR})

target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES}) 
target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME}-1 ${BIN} ${MESHES_DIR}/mesh3d/hex1.mesh3d 3)
foreach(i RANGE 0 7)
    add_test(${PROJECT_NAME}-2-${i} ${BIN} ${MESHES_DIR}/mesh3d/hex2-ori${i}.mesh3d 3)
endforeach(i RANGE 0 7)
add_test(${PROJECT_NAME}-4  ${BIN} ${MESHES_DIR}/mesh3d/hex4.mesh3d 3)
add_test(${PROJECT_NAME}-4x ${BIN} ${MESHES_DIR}/mesh3d/hex4x.mesh3d 3)
add_test(${PROJECT_NAME}-4y ${BIN} ${MESHES_DIR}/mesh3d/hex4y.mesh3d 3)
add_test(${PROJECT_NAME}-4z ${BIN} ${MESHES_DIR}/mesh3d/hex4z.mesh3d 3)
add_test(${PROJECT_NAME}-8  ${BIN} ${MESHES_DIR}/mesh3d/hex8.mesh3d 3)

