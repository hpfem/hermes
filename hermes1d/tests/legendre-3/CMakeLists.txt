project(legendre-3)

add_executable(${PROJECT_NAME} main.cpp)
include (../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(legendre-3 ${BIN})

target_link_libraries(${PROJECT_NAME} ${HERMES_BIN})
