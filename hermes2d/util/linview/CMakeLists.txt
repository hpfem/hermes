if(NOT WITH_GLUT)
    message(FATAL_ERROR "util/linview must be compiled with GLUT")
endif(NOT WITH_GLUT)

if(NOT H2D_REAL)
    message(STATUS "skipping util/linview (real version is not being built)")
    return()
endif(NOT H2D_REAL)

macro(LINVIEW_SETTINGS VIEW)
    add_executable(${PROJECT_NAME} main.cpp)
    
    include_directories(${hermes2d_SOURCE_DIR}/src)
    include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
	  if(MSVC)
		  include_directories(${DEP_ROOT}/include)
	  endif(MSVC)
	  
	  PICK_REAL_OR_CPLX_LIBS(${HERMES_REAL_BIN} ${PROJECT_NAME})
    get_property(FLAGS TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS)

    # Should debugging information be included?
    if(DEBUG)
      set(FLAGS "${FLAGS} -g")
    endif(DEBUG)

    set(FLAGS "${FLAGS} -DWHATVIEW=${VIEW}")
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS ${FLAGS})
        
    # Do we need JUDY and GLUT libraries here?
    target_link_libraries(  ${PROJECT_NAME} 
      ${HERMES_REAL_BIN}
      ${PYTHON_LIBRARIES}
      ${PETSC_LIBRARIES}
      ${UMFPACK_LIBRARIES}
      ${PARDISO_LIBRARY}
      ${MUMPS_LIBRARIES}
      ${TRILINOS_LIBRARIES}
      ${LAPACK_LIBRARIES} 
      ${EXODUSII_LIBRARY} 
      ${HDF5_LIBRARY}
    )
    add_dependencies(${PROJECT_NAME} ${HERMES_REAL_BIN})
endmacro(LINVIEW_SETTINGS)

add_definitions(-DH2D_REAL)

include(PickRealOrCplxLibs)

project(linview)
LINVIEW_SETTINGS(ScalarView)
set_target_properties(linview PROPERTIES COMPILE_FLAGS "${FLAGS} -DCONTOURS")

project(vecview)
LINVIEW_SETTINGS(VectorView)

project(ordview)
LINVIEW_SETTINGS(OrderView)
