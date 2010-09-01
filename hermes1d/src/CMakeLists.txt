set(SRC
    common.cpp iterator.cpp qsort.cpp
    lobatto.cpp legendre.cpp
    discrete.cpp solution.cpp mesh.cpp
    linearizer.cpp quad_std.cpp transforms.cpp
    adapt.cpp graph.cpp h1_polys.cpp
    )

add_definitions(-DCOMPLEX=std::complex<double>)

include_directories(${hermes1d_SOURCE_DIR}/hermes_common/)
include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})

add_library(${HERMES_BIN} SHARED ${SRC})
if(DEBUG)
    set_target_properties(hermes1d-debug PROPERTIES COMPILE_FLAGS
        ${CMAKE_CXX_FLAGS_DEBUG})
endif(DEBUG)
if(RELEASE)
    set_target_properties(hermes1d PROPERTIES COMPILE_FLAGS
        ${CMAKE_CXX_FLAGS_RELEASE})
endif(RELEASE)

target_link_libraries(${HERMES_BIN} hermes_common)

install(TARGETS ${HERMES_BIN}
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

file(GLOB header_files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${header_files} DESTINATION include/hermes1d)

