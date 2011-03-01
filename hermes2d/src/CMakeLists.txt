#
# Source files for the Hermes2d library.
#
set(SRC
       adapt/adapt.cpp
       adapt/kelly_type_adapt.cpp

       function/solution.cpp 
       function/filter.cpp
       function/norm.cpp
       function/forms.cpp
	   
       linearizer/linear1.cpp 
       linearizer/linear2.cpp 
       linearizer/linear3.cpp 

       mesh/refmap.cpp 
       mesh/curved.cpp
       mesh/refinement_type.cpp 
       mesh/element_to_refine.cpp
       mesh/exodusii.cpp 
       mesh/hash.cpp 
       mesh/h2d_reader.cpp
       mesh/mesh.cpp 
       mesh/regul.cpp 
       mesh/transform.cpp 
       mesh/traverse.cpp
       mesh/trans.cpp
       mesh/python_reader.cpp

       quadrature/limit_order.cpp
       quadrature/quad_std.cpp

       ref_selectors/selector.cpp 
       ref_selectors/order_permutator.cpp 
       ref_selectors/optimum_selector.cpp 
       ref_selectors/proj_based_selector.cpp 
       ref_selectors/l2_proj_based_selector.cpp 
       ref_selectors/h1_proj_based_selector.cpp 
       ref_selectors/hcurl_proj_based_selector.cpp

       shapeset/shapeset.cpp 
       shapeset/shapeset_h1_ortho.cpp 
       shapeset/shapeset_h1_jacobi.cpp 
       shapeset/shapeset_h1_quad.cpp
       shapeset/shapeset_hc_legendre.cpp 
       shapeset/shapeset_hc_gradleg.cpp
       shapeset/shapeset_hd_legendre.cpp
       shapeset/shapeset_l2_legendre.cpp
       shapeset/precalc.cpp 

       space/space.cpp 
       space/space_h1.cpp 
       space/space_hcurl.cpp 
       space/space_l2.cpp
       space/space_hdiv.cpp

       views/base_view.cpp 
       views/mesh_view.cpp 
       views/order_view.cpp 
       views/scalar_view.cpp 
       views/stream_view.cpp 
       views/vector_base_view.cpp 
       views/vector_view.cpp 
       views/view.cpp 
       views/view_data.cpp 
       views/view_support.cpp

       weakform/weakform.cpp 

       neighbor.cpp
       numerical_flux.cpp
       graph.cpp
       ogprojection.cpp
       h2d_common.cpp  
       discrete_problem.cpp
       runge_kutta.cpp
       boundaryconditions/boundaryconditions.cpp
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

  set(_FLAGS ${BUILD_FLAGS})

  if(REPORT_RUNTIME_CONTROL)
    set(_FLAGS "${_FLAGS} -DHERMES_REPORT_RUNTIME_CONTROL")
  endif(REPORT_RUNTIME_CONTROL)
  if(WITH_VIEWER_GUI)
    set(_FLAGS "${_FLAGS} -DENABLE_VIEWER_GUI")
  endif(WITH_VIEWER_GUI)

  set_target_properties(${HERMES_LIB} PROPERTIES COMPILE_FLAGS ${_FLAGS})
  PICK_REAL_OR_CPLX_INCS(${HERMES_COMMON_LIB} ${HERMES_LIB})

  target_link_libraries(  ${HERMES_LIB} 
    ${HERMES_COMMON_LIB}
    ${GLUT_LIBRARY} ${GLEW_LIBRARY}
    ${ANTTWEAKBAR_LIBRARY}
    ${LAPACK_LIBRARY}
    ${CLAPACK_LIBRARY} ${BLAS_LIBRARY}
    ${PYTHON_LIBRARY}
  )
endmacro(BUILD_LIB)

# Build and install real version:
if(H2D_REAL)
  set(FLAGS "-DHERMES_COMMON_REAL -DH2D_REAL")

  if(MSVC)
    BUILD_LIB(${HERMES2D_REAL} ${HERMES_COMMON_LIB_REAL} ${FLAGS})
    ADD_MSVC_BUILD_FLAGS(${HERMES2D_REAL} ${HERMES2D_REAL_DEBUG} ${HERMES2D_REAL_RELEASE})
    INSTALL_LIB(${HERMES2D_REAL})
  else(MSVC)
    if(DEBUG)
      set(BUILD_FLAGS "${FLAGS} ${DEBUG_FLAGS}")
      BUILD_LIB(${HERMES2D_REAL_DEBUG} ${HERMES_COMMON_LIB_REAL_DEBUG} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES2D_REAL_DEBUG})
    endif(DEBUG)
    if(RELEASE)
      set(BUILD_FLAGS "${FLAGS} ${RELEASE_FLAGS}")
      BUILD_LIB(${HERMES2D_REAL_RELEASE} ${HERMES_COMMON_LIB_REAL_RELEASE} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES2D_REAL_RELEASE})
    endif(RELEASE)
  endif(MSVC)
endif(H2D_REAL)

# Build and install complex version:
if(H2D_COMPLEX)
  set(FLAGS "-DHERMES_COMMON_COMPLEX -DH2D_COMPLEX")

  if(MSVC)
    BUILD_LIB(${HERMES2D_CPLX} ${HERMES_COMMON_LIB_CPLX} ${FLAGS})
    ADD_MSVC_BUILD_FLAGS(${HERMES2D_CPLX} ${HERMES2D_CPLX_DEBUG} ${HERMES2D_CPLX_RELEASE})
    INSTALL_LIB(${HERMES2D_CPLX})
  else(MSVC)
    if(DEBUG)
      set(BUILD_FLAGS "${FLAGS} ${DEBUG_FLAGS}")
      BUILD_LIB(${HERMES2D_CPLX_DEBUG} ${HERMES_COMMON_LIB_CPLX_DEBUG} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES2D_CPLX_DEBUG})
    endif(DEBUG)
    if(RELEASE)
      set(BUILD_FLAGS "${FLAGS} ${RELEASE_FLAGS}")
      BUILD_LIB(${HERMES2D_CPLX_RELEASE} ${HERMES_COMMON_LIB_CPLX_RELEASE} ${BUILD_FLAGS})
      INSTALL_LIB(${HERMES2D_CPLX_RELEASE})
    endif(RELEASE)
  endif(MSVC)
endif(H2D_COMPLEX)


# Install header files.

file(GLOB INC           "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
file(GLOB INC_SPACE     "${CMAKE_CURRENT_SOURCE_DIR}/space/*.h")
file(GLOB INC_MESH      "${CMAKE_CURRENT_SOURCE_DIR}/mesh/*.h")
file(GLOB INC_BC        "${CMAKE_CURRENT_SOURCE_DIR}/boundaryconditions/*.h")
file(GLOB INC_SHAPESET  "${CMAKE_CURRENT_SOURCE_DIR}/shapeset/*.h")
file(GLOB INC_FUNCTION  "${CMAKE_CURRENT_SOURCE_DIR}/function/*.h")
file(GLOB INC_INTEGRALS "${CMAKE_CURRENT_SOURCE_DIR}/integrals/*.h")
file(GLOB INC_LINEARIZER "${CMAKE_CURRENT_SOURCE_DIR}/linearizer/*.h")
file(GLOB INC_ADAPT     "${CMAKE_CURRENT_SOURCE_DIR}/adapt/*.h")
file(GLOB INC_QUAD      "${CMAKE_CURRENT_SOURCE_DIR}/quadrature/*.h")
file(GLOB INC_REF_SEL   "${CMAKE_CURRENT_SOURCE_DIR}/ref_selectors/*.h")
file(GLOB INC_VIEWS     "${CMAKE_CURRENT_SOURCE_DIR}/views/*.h")
file(GLOB INC_WEAKFORM  "${CMAKE_CURRENT_SOURCE_DIR}/weakform/*.h")

install(FILES ${INC}           DESTINATION include/hermes/hermes2d/src/)
install(FILES ${INC_SPACE}     DESTINATION include/hermes/hermes2d/src/space)
install(FILES ${INC_MESH}      DESTINATION include/hermes/hermes2d/src/mesh)
install(FILES ${INC_BC}        DESTINATION include/hermes/hermes2d/src/boundaryconditions)
install(FILES ${INC_SHAPESET}  DESTINATION include/hermes/hermes2d/src/shapeset)
install(FILES ${INC_FUNCTION}  DESTINATION include/hermes/hermes2d/src/function)
install(FILES ${INC_INTEGRALS} DESTINATION include/hermes/hermes2d/src/integrals)
install(FILES ${INC_LINEARIZER} DESTINATION include/hermes/hermes2d/src/linearizer)
install(FILES ${INC_ADAPT}     DESTINATION include/hermes/hermes2d/src/adapt)
install(FILES ${INC_QUAD}      DESTINATION include/hermes/hermes2d/src/quadrature)
install(FILES ${INC_REF_SEL}   DESTINATION include/hermes/hermes2d/src/ref_selectors)
install(FILES ${INC_VIEWS}     DESTINATION include/hermes/hermes2d/src/views)
install(FILES ${INC_WEAKFORM}  DESTINATION include/hermes/hermes2d/src/weakform)
