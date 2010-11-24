project(lobatto-1)

add_executable(${PROJECT_NAME} main.cpp)
include (../../examples/CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(lobatto-1 ${BIN})

target_link_libraries(${PROJECT_NAME} ${HERMES_BIN})
