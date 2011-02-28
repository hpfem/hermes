include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#
# Source files for the Hermes3D library.
#
set(SRC
  hermes3d.cpp
  filter.cpp
  weakform/forms.cpp
  function.cpp
  mesh.cpp
  discrete_problem.cpp
  ogprojection.cpp
  loader/ctuReader.cpp
  loader/exodusii.cpp
  loader/mesh3d.cpp
  loader/hdf5.cpp
  norm.cpp
  output/gmsh.cpp
  output/vtk.cpp
  output/graph.cpp
  quadcheb.cpp
  quadstd.cpp
  refdomain.cpp
  refmap.cpp
  shapefn.cpp
  shapeset/shapeset.cpp
  shapeset/lobatto.cpp
  shapeset/h1lobattotetra.cpp
  shapeset/h1lobattotetradx.cpp
  shapeset/h1lobattotetrady.cpp
  shapeset/h1lobattotetradz.cpp
  shapeset/h1lobattohex.cpp
  shapeset/hcurllobattohex.cpp
  shapeset/refmapss.cpp
  solution.cpp
  space/space.cpp
  space/h1.cpp
  space/hcurl.cpp
  transform.cpp
  traverse.cpp
  weakform/weakform.cpp
        adapt/adapt.cpp
  adapt/proj.cpp
  adapt/projipol.cpp
  adapt/h1proj.cpp
  adapt/hcurlproj.cpp
  adapt/h1projipol.cpp
)

#
# Build and install the library.
#

include(BuildAndInstallScripts)
include(PickRealOrCplxLibs)

# Add a target for a specified version of the Hermes library (REAL/COMPLEX, RELEASE/DEBUG),
# and set its build properties (compile and link flags, installation directories). Due to
# different optional libraries for H1D, H2D and H3D, this macro is dimension-specific.
macro(BUILD_LIB HERMES_LIB HERMES_COMMON_LIB BUILD_FLAGS)
  add_library(${HERMES_LIB} SHARED ${SRC})
 
  set_target_properties(${HERMES_LIB} PROPERTIES COMPILE_FLAGS ${BUILD_FLAGS})
  PICK_REAL_OR_CPLX_INCS(${HERMES_COMMON_LIB} ${HERMES_LIB})

  target_link_libraries(  ${HERMES_LIB} 
    ${HERMES_COMMON_LIB}
    ${LAPACK_LIBRARY}
    ${CLAPACK_LIBRARY} ${BLAS_LIBRARY}
    ${PYTHON_LIBRARY}
  )
endmacro(BUILD_LIB)

# Build and install real version:
if(H3D_REAL)
  set(FLAGS "-DHERMES_COMMON_REAL -DH3D_REAL")

  if(MSVC)
    BUILD_LIB(${HERMES3D_REAL} ${HERMES_COMMON_LIB_REAL} ${FLAGS})
    ADD_MSVC_BUILD_FLAGS(${HERMES3D_REAL} ${HERMES3D_REAL_DEBUG} ${HERMES3D_REAL_RELEASE})
    INSTALL_LIB(${HERMES3D_REAL})
  else(MSVC)
    if(DEBUG)
      set(BUILD_FLAGS "${FLAGS} ${DEBUG_FLAGS}")
      BUILD_LIB(${HERMES3D_REAL_DEBUG} ${HERMES_COMMON_LIB_REAL_DEBUG} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES3D_REAL_DEBUG})
    endif(DEBUG)
    if(RELEASE)
      set(BUILD_FLAGS "${FLAGS} ${RELEASE_FLAGS}")
      BUILD_LIB(${HERMES3D_REAL_RELEASE} ${HERMES_COMMON_LIB_REAL_RELEASE} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES3D_REAL_RELEASE})
    endif(RELEASE)
  endif(MSVC)
endif(H3D_REAL)

# Build and install complex version:
if(H3D_COMPLEX)
  set(FLAGS "-DHERMES_COMMON_COMPLEX -DH3D_COMPLEX")

  if(MSVC)
    BUILD_LIB(${HERMES3D_CPLX} ${HERMES_COMMON_LIB_CPLX} ${FLAGS})
    ADD_MSVC_BUILD_FLAGS(${HERMES3D_CPLX} ${HERMES3D_CPLX_DEBUG} ${HERMES3D_CPLX_RELEASE})
    INSTALL_LIB(${HERMES3D_CPLX})
  else(MSVC)
    if(DEBUG)
      set(BUILD_FLAGS "${FLAGS} ${DEBUG_FLAGS}")
      BUILD_LIB(${HERMES3D_CPLX_DEBUG} ${HERMES_COMMON_LIB_CPLX_DEBUG} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES3D_CPLX_DEBUG})
    endif(DEBUG)
    if(RELEASE)
      set(BUILD_FLAGS "${FLAGS} ${RELEASE_FLAGS}")
      BUILD_LIB(${HERMES3D_CPLX_RELEASE} ${HERMES_COMMON_LIB_CPLX_RELEASE} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES3D_CPLX_RELEASE})
    endif(RELEASE)
  endif(MSVC)
endif(H3D_COMPLEX)

# Install header files.

file(GLOB INC           "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
file(GLOB INC_SPACE     "${CMAKE_CURRENT_SOURCE_DIR}/space/*.h")
file(GLOB INC_SHAPESET  "${CMAKE_CURRENT_SOURCE_DIR}/shapeset/*.h")
file(GLOB INC_LOADER    "${CMAKE_CURRENT_SOURCE_DIR}/loader/*.h")
file(GLOB INC_OUTPUT    "${CMAKE_CURRENT_SOURCE_DIR}/output/*.h")
file(GLOB INC_INTEGRALS "${CMAKE_CURRENT_SOURCE_DIR}/integrals/*.h")
file(GLOB INC_ADAPT     "${CMAKE_CURRENT_SOURCE_DIR}/adapt/*.h")
file(GLOB INC_WEAKFORM  "${CMAKE_CURRENT_SOURCE_DIR}/weakform/*.h")

install(FILES ${INC}           DESTINATION include/hermes/hermes3d/src)
install(FILES ${INC_SPACE}     DESTINATION include/hermes/hermes3d/src/space)
install(FILES ${INC_SHAPESET}  DESTINATION include/hermes/hermes3d/src/shapeset)
install(FILES ${INC_LOADER}    DESTINATION include/hermes/hermes3d/src/loader)
install(FILES ${INC_OUTPUT}    DESTINATION include/hermes/hermes3d/src/output)
install(FILES ${INC_INTEGRALS} DESTINATION include/hermes/hermes3d/src/integrals)
install(FILES ${INC_ADAPT}     DESTINATION include/hermes/hermes3d/src/adapt)
install(FILES ${INC_WEAKFORM}  DESTINATION include/hermes/hermes3d/src/weakform)
