project(hermes1d)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../hermes_common/cmake)

IF(H1D_DEBUG)
	set(DEBUG YES)
ELSE(H1D_DEBUG)
	set(DEBUG NO)
ENDIF(H1D_DEBUG)

IF(H1D_RELEASE)
	set(RELEASE NO)
ELSE(H1D_RELEASE)
	set(RELEASE NO)
ENDIF(H1D_RELEASE)


IF(H1D_WITH_EXAMPLES)
	set(WITH_EXAMPLES YES)
ELSE(H1D_WITH_EXAMPLES)
	set(WITH_EXAMPLES NO)
ENDIF(H1D_WITH_EXAMPLES)

IF(H1D_WITH_BENCHMARKS)
	set(WITH_BENCHMARKS YES)
ELSE(H1D_WITH_BENCHMARKS)
	set(WITH_BENCHMARKS NO)
ENDIF(H1D_WITH_BENCHMARKS)

IF(H1D_WITH_TESTS)
	set(WITH_TESTS YES)
ELSE(H1D_WITH_TESTS)
	set(WITH_TESTS NO)
ENDIF(H1D_WITH_TESTS)

# Doxygen related
set(DOXYGEN_BINARY doxygen)
set(DOXYGEN_CONFIG_FILE ${PROJECT_SOURCE_DIR}/doc/Doxyfile)

set(SPHINX_BUILD_BIN sphinx-build)

# For use in the common code.
set(H1D_REAL YES)
add_definitions(-DH1D_REAL)

if(DEBUG)
        set(HERMES1D hermes1d-debug)
endif(DEBUG)

if(RELEASE)
        set(HERMES1D hermes1d)
endif(RELEASE)

# Python inclusion is different in H1D, H2D, H3D.
if(WITH_PYTHON)
  find_package(Python REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(Cython REQUIRED)
  include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
  set(PYTHON_INSTALL_PATH lib/python/site-packages)
  # This is needed to make hermes_common working properly:
  SET(CYTHON_INCLUDE_DIRECTORIES . -I${CMAKE_CURRENT_SOURCE_DIR}/hermes1d/hermes_common)
  add_subdirectory(hermes1d)
endif(WITH_PYTHON)

# Subdirectories.
add_subdirectory(src)

if(WITH_EXAMPLES)
    add_subdirectory(examples)
endif(WITH_EXAMPLES)

if(WITH_BENCHMARKS)
    add_subdirectory(benchmarks)
endif(WITH_BENCHMARKS)

if(WITH_TESTS)
    enable_testing()
	add_subdirectory(tests)
endif(WITH_TESTS)

add_subdirectory(doc.old)

add_custom_target(doc.old DEPENDS ${PROJECT_BINARY_DIR}/doc.old/xml/index.xml)
add_custom_command(
	TARGET    doc.old
	COMMAND   ${SPHINX_BUILD_BIN} -N -b html ${PROJECT_SOURCE_DIR}/doc.old ${PROJECT_BINARY_DIR}/doc.old/html
)

add_custom_target(doc-tex DEPENDS ${PROJECT_BINARY_DIR}/doc.old/xml/index.xml)
add_custom_command(
	TARGET    doc-tex
	COMMAND   ${SPHINX_BUILD_BIN} -N -b latex ${PROJECT_SOURCE_DIR}/doc.old ${PROJECT_BINARY_DIR}/doc.old/tex
)

add_custom_command(
	OUTPUT    ${PROJECT_BINARY_DIR}/doc.old/xml/index.xml
	COMMAND   ${DOXYGEN_BINARY} ${DOXYGEN_CONFIG_FILE}
)