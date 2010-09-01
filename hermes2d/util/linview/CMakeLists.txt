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
    include_directories(${hermes2d_SOURCE_DIR}/hermes_common/)
    include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
    include_directories(${TRILINOS_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} hermes_common ${PYTHON_LIBRARIES} ${HERMES_REAL_BIN} ${LAPACK_LIBRARIES} ${TRILINOS_LIBRARIES})
    # if having Trilinos with UMFPACK support we need to link against UMFPACK
    # if we do not have UMFPACK support at all, these variables are empty and cause no harm
    target_link_libraries(${PROJECT_NAME} ${AMD_LIBRARY} ${UMFPACK_LIBRARY})
	if(MSVC)
		include_directories(${DEP_ROOT}/include)
	endif(MSVC)
    target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS
        -DWHATVIEW=${VIEW})
endmacro(LINVIEW_SETTINGS)

project(linview)
LINVIEW_SETTINGS(ScalarView)
set_target_properties(linview PROPERTIES
    COMPILE_FLAGS "-DWHATVIEW=ScalarView -DCONTOURS")

project(vecview)
LINVIEW_SETTINGS(VectorView)

project(ordview)
LINVIEW_SETTINGS(OrderView)
