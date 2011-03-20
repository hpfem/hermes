cmake_minimum_required(VERSION 2.6)

project(hermes1d)

IF(H1D_DEBUG)
  set(DEBUG YES)
ELSE(H1D_DEBUG)
  set(DEBUG NO)
ENDIF(H1D_DEBUG)

IF(H1D_RELEASE)
  set(RELEASE YES)
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


#
# Version configuration.
#

if (NOT DEBUG AND NOT RELEASE)
  message(FATAL_ERROR "WITH_H1D is set to YES, but no versions are specified. 
                       Re-run with DEBUG and/or RELEASE set to YES.")
endif (NOT DEBUG AND NOT RELEASE)

# Depending on the setting of H1D_DEBUG and H1D_RELEASE, the user will have 
# two variables at his disposal: HERMES1D_RELEASE and/or HERMES1D_DEBUG
# (corresponding to the two libraries built on his system). He can use these
# to selectively link his project to either the optimized version or
# one with debugging symbols by adding the appropriate second argument to 
# set_common_target_properties in his project's CMakeLists.txt. If he calles
# that function with just the name of the project, the Debug version will be 
# automatically selected if present (i.e. H1D_DEBUG was set to YES). Otherwise,
# the Release version will be selected. The appropriate version of hermes_common
# will be linked as well.
if(MSVC)
    set(HERMES1D_RELEASE "hermes1d")
    set(HERMES1D_DEBUG "hermes1d-debug")
    if(RELEASE)
		set(HERMES1D ${HERMES1D_RELEASE})
	endif(RELEASE)
	if(DEBUG)
		set(HERMES1D ${HERMES1D_DEBUG})
	endif(DEBUG)
else(MSVC)
  if(RELEASE)
    set(HERMES1D_RELEASE "hermes1d")
    set(HERMES1D ${HERMES1D_RELEASE})
  
    # set(CMAKE_BUILD_TYPE Release) # This does not work with multi-version configurations,
                                    # because it is linked to the project (i.e. we would have
                                    # to have separate projects for Debug and Release)
  endif(RELEASE)
  if(DEBUG)
    set(HERMES1D_DEBUG "hermes1d-debug")
    set(HERMES1D ${HERMES1D_DEBUG})
    #set(CMAKE_BUILD_TYPE Debug)
  endif(DEBUG)
endif(MSVC)

#
# Python
#

if(WITH_PYTHON)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/hermes1d/h1d_wrapper/)
  set(CYTHON_INCLUDE_DIRECTORIES -I.)
  set(CYTHON_INCLUDE_DIRECTORIES ${CYTHON_INCLUDE_DIRECTORIES} -I${CMAKE_SOURCE_DIR}/hermes_common)
  add_subdirectory(hermes1d)
endif(WITH_PYTHON)

#
# Subdirectories.
#

add_subdirectory(src)

if(WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif(WITH_TESTS)

if(WITH_EXAMPLES)
    add_subdirectory(examples)
endif(WITH_EXAMPLES)

if(WITH_BENCHMARKS)
    add_subdirectory(benchmarks)
endif(WITH_BENCHMARKS)

#
# Documentation.
#

# Doxygen related
set(DOXYGEN_BINARY doxygen)
set(DOXYGEN_CONFIG_FILE ${PROJECT_SOURCE_DIR}/doc/Doxyfile)

set(SPHINX_BUILD_BIN sphinx-build)

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
