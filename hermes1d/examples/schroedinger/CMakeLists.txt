project(cppschroedinger)

include_directories(${hermes1d_SOURCE_DIR}/src)
SET(CYTHON_INCLUDE_DIRECTORIES -I${CMAKE_SOURCE_DIR}/hermes_common)
SET(CYTHON_INCLUDE_DIRECTORIES ${CYTHON_INCLUDE_DIRECTORIES} -I${hermes1d_SOURCE_DIR})

cython_add_module(_forms forms.cpp)
target_link_libraries(_forms ${HERMES1D})

if(WITH_TESTS)
  if(WITH_PYTHON)
#    FIXME: not adapted to the latest API
#    add_subdirectory(tests)
  endif(WITH_PYTHON)
endif(WITH_TESTS)

