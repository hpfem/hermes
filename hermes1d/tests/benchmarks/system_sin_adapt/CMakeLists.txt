project(benchmark-system_sin_adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-system_sin_adapt ${BIN})

target_link_libraries(${PROJECT_NAME} ${HERMES_BIN})
