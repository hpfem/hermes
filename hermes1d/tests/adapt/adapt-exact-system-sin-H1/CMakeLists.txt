project(adapt-exact-system-sin-H1)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(adapt-exact-system-sin-H1 ${BIN})

target_link_libraries(${PROJECT_NAME} ${HERMES_BIN})
