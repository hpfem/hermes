if(NOT H2D_WITH_GLUT)
    message(FATAL_ERROR "util/linview must be compiled with GLUT")
endif(NOT H2D_WITH_GLUT)

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
	  
	PICK_REAL_OR_CPLX_LIBS(${HERMES2D_REAL} ${PROJECT_NAME})
    get_property(FLAGS TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS)

    # Should debugging information be included?
    if(DEBUG)
      set(FLAGS "${FLAGS} -g")
    endif(DEBUG)

    set(FLAGS "${FLAGS} -DWHATVIEW=${VIEW}")
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS ${FLAGS})
        
    target_link_libraries(${PROJECT_NAME} ${HERMES2D_REAL} ${TRILINOS_LIBRARIES})
    add_dependencies(${PROJECT_NAME} ${HERMES2D_REAL})
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
