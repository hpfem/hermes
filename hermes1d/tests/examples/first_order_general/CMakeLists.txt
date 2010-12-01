project(example-first_order_general)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-first_order_general ${BIN})

target_link_libraries(${PROJECT_NAME} ${HERMES_BIN})
