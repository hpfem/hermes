CYTHON_ADD_MODULE_PYX(hermes2d)
CYTHON_ADD_MODULE_PYX(numpy_utils)

include_directories(${PROJECT_BINARY_DIR}/src)
add_python_library(numpy_utils numpy_utils.cpp)
add_python_library(hermes2d hermes2d.cpp)
target_link_libraries(hermes2d ${HERMES_REAL_BIN})

add_subdirectory(modules)
