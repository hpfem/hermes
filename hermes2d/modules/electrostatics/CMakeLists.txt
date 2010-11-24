project(module-electrostatics)

add_library(hermes_electrostatics SHARED electrostatics.cpp)
target_link_libraries(hermes_electrostatics ${HERMES_REAL_BIN})

add_executable(${PROJECT_NAME} main.cpp disc.cpp)
target_link_libraries(${PROJECT_NAME} hermes_electrostatics)
include (../CMake.common)
