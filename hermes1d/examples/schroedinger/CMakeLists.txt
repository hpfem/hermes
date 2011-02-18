project(cppschroedinger)

include_directories(${hermes1d_SOURCE_DIR}/hermes1d/h1d_wrapper/)
include_directories(${hermes1d_SOURCE_DIR}/../hermes_common/)
include_directories(${hermes1d_SOURCE_DIR}/src)

SET(CYTHON_INCLUDE_DIRECTORIES ${CYTHON_INCLUDE_DIRECTORIES} -I../../)
cython_add_module(_forms forms.cpp)
target_link_libraries(_forms ${HERMES1D})
