CYTHON_ADD_MODULE_PYX(electrostatics)
include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_BINARY_DIR}/modules/electrostatics)
add_python_library(electrostatics electrostatics.cpp)
target_link_libraries(electrostatics hermes_electrostatics)

