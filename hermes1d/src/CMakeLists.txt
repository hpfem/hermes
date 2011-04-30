set(SRC
    iterator.cpp 
    lobatto.cpp legendre.cpp
    discrete_problem.cpp solution.cpp space.cpp
    ogprojection.cpp
    ogprojection2.cpp
    linearizer.cpp quad_std.cpp transforms.cpp
    adapt.cpp graph.cpp h1_polys.cpp
    weakform.cpp
  )


#
# Build and install the library.
#

include(BuildAndInstallScripts)
include(PickRealOrCplxLibs)

# Add a target for a specified version of the Hermes library (in case of H1D, just RELEASE/DEBUG),
# and set its build properties (compile and link flags, installation directories). Due to
# different optional libraries for H1D, H2D and H1D, this macro is dimension-specific.
macro(BUILD_LIB HERMES_LIB HERMES_COMMON_LIB BUILD_FLAGS)
  add_library(${HERMES_LIB} SHARED ${SRC})
  
  add_definitions(-DCOMPLEX=std::complex<double>)
  set_target_properties(${HERMES_LIB} PROPERTIES COMPILE_FLAGS ${BUILD_FLAGS})
  PICK_REAL_OR_CPLX_INCS(${HERMES_COMMON_LIB} ${HERMES_LIB})

  target_link_libraries(  ${HERMES_LIB} ${HERMES_COMMON_LIB}  )
endmacro(BUILD_LIB)

set(FLAGS "-DHERMES_COMMON_REAL -DH1D_REAL")

if(MSVC)
  BUILD_LIB(${HERMES1D} ${HERMES_COMMON_LIB_REAL} ${FLAGS})
  ADD_MSVC_BUILD_FLAGS(${HERMES1D} ${HERMES1D_DEBUG} ${HERMES1D_RELEASE})
  INSTALL_LIB(${HERMES1D})
else(MSVC)
  if(DEBUG)
    set(BUILD_FLAGS "${FLAGS} ${DEBUG_FLAGS}")
    BUILD_LIB(${HERMES1D_DEBUG} ${HERMES_COMMON_LIB_REAL_DEBUG} ${BUILD_FLAGS})
    INSTALL_LIB(${HERMES1D_DEBUG})
  endif(DEBUG)
  if(RELEASE)
    set(BUILD_FLAGS "${FLAGS} ${RELEASE_FLAGS}")
    BUILD_LIB(${HERMES1D_RELEASE} ${HERMES_COMMON_LIB_REAL_RELEASE} ${BUILD_FLAGS})
    INSTALL_LIB(${HERMES1D_RELEASE})
  endif(RELEASE)
endif(MSVC)

# Install header files.

file(GLOB header_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${header_files} DESTINATION include/hermes/hermes1d)

