#
# TODO: 
# - COMPLEX version
# - install the library system-wide
#

set(PYTHON_LIB_NAME hermes3dpy)

include_directories(${PYTHON_INCLUDE_PATH})
include_directories(${hermes3d_SOURCE_DIR})

add_custom_target(hermes3d.cpp
	COMMAND ${CYTHON_BIN} ${CMAKE_CURRENT_SOURCE_DIR}/"hermes3d.pyx" -o "${CMAKE_CURRENT_BINARY_DIR}/hermes3d.cpp"
	COMMENT "Generating hermes3d.cpp with cython"
)

add_custom_command(
	OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/hermes3d.cpp
	COMMENT "---"
)

add_library(${PYTHON_LIB_NAME} SHARED hermes3d.cpp)
add_dependencies(${PYTHON_LIB_NAME} hermes3d.cpp)

target_link_libraries(${PYTHON_LIB_NAME} hermes3d-debug)
target_link_libraries(${PYTHON_LIB_NAME} Judy)

