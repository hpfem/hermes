project(hermes1d)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(WITH_PYTHON yes)
set(DEBUG yes)
set(RELEASE no)
set(WITH_EXAMPLES yes)
set(WITH_TESTS yes)

# Doxygen related
set(DOXYGEN_BINARY doxygen)
set(DOXYGEN_CONFIG_FILE ${PROJECT_SOURCE_DIR}/doc/Doxyfile)

set(SPHINX_BUILD_BIN sphinx-build)

# allow to override the default values in CMake.vars
if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
    include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)


if(DEBUG)
        set(HERMES_BIN hermes1d-debug)
endif(DEBUG)

if(RELEASE)
        set(HERMES_BIN hermes1d)
endif(RELEASE)


add_subdirectory(hermes_common)
add_subdirectory(src)

find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

if(WITH_PYTHON)
    include(UseCython)

    # if the user specified PYTHON_INCLUDE_PATH, let's use that, otherwise let's
    # try to guess using cmake/FindPythonLibs.cmake
    if (NOT PYTHON_INCLUDE_PATH)
        find_package(PythonLibs REQUIRED)
    endif (NOT PYTHON_INCLUDE_PATH)
    find_package(NumPy REQUIRED)

    include_directories(${hermes1d_SOURCE_DIR}/hermes_common/)
    include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
    include_directories(src)
    add_subdirectory(hermes1d)
endif(WITH_PYTHON)
if(WITH_EXAMPLES)
    add_subdirectory(examples)
endif(WITH_EXAMPLES)

if(WITH_TESTS)
    enable_testing()
	add_subdirectory(tests)
endif(WITH_TESTS)

add_subdirectory(doc)

add_custom_target(doc DEPENDS ${PROJECT_BINARY_DIR}/doc/xml/index.xml)
add_custom_command(
	TARGET    doc
	COMMAND   ${SPHINX_BUILD_BIN} -N -b html ${PROJECT_SOURCE_DIR}/doc ${PROJECT_BINARY_DIR}/doc/html
)

add_custom_target(doc-tex DEPENDS ${PROJECT_BINARY_DIR}/doc/xml/index.xml)
add_custom_command(
	TARGET    doc-tex
	COMMAND   ${SPHINX_BUILD_BIN} -N -b latex ${PROJECT_SOURCE_DIR}/doc ${PROJECT_BINARY_DIR}/doc/tex
)

add_custom_command(
	OUTPUT    ${PROJECT_BINARY_DIR}/doc/xml/index.xml
	COMMAND   ${DOXYGEN_BINARY} ${DOXYGEN_CONFIG_FILE}
)


add_custom_target(test-quick
	COMMAND   /usr/bin/ctest -E \"hnnd|qorder|lobatto|hex-h1-1face\"
)

message("\nConfiguration results:")
message("----------------------")
message("Build with Python: ${WITH_PYTHON}")
message("Build with debug: ${DEBUG}")
message("Build with release: ${RELEASE}")
message("Build with tests: ${WITH_TESTS}")
message("\n")
