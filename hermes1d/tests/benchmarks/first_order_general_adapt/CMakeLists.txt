project(benchmark-first_order_general_adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-first_order_general_adapt ${BIN})

target_link_libraries(${PROJECT_NAME} ${HERMES_BIN})
