project(example-system_neutronics_fixedsrc)

add_executable(${PROJECT_NAME} main.cpp)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(example-system_neutronics_fixedsrc ${BIN})

target_link_libraries(${PROJECT_NAME} ${HERMES_BIN})
