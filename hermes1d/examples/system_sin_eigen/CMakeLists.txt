project(system_sin_eigen)

add_executable(${PROJECT_NAME} main.cpp)
include(../CMake.common)

target_link_libraries(${PROJECT_NAME} ${HERMES_BIN})