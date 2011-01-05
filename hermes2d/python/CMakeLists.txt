include(UseCython)

# if the user specified PYTHON_INCLUDE_PATH, let's use that, otherwise let's
# try to guess using cmake/FindPythonLibs.cmake
if (NOT PYTHON_INCLUDE_PATH)
    find_package(PythonLibs REQUIRED)
endif (NOT PYTHON_INCLUDE_PATH)
find_package(NumPy REQUIRED)

include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
include_directories(${TRILINOS_INCLUDE_DIR})
include_directories(${hermes2d_SOURCE_DIR}/src/)
include_directories(${hermes2d_SOURCE_DIR}/hermes_common/)

add_subdirectory(hermes2d)

add_custom_target(test-python
    COMMAND ctest -R python-tests -V
    )
add_test(python-tests ./test)
