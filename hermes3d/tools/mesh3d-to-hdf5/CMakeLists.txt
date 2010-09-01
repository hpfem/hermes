project(mesh3d-to-hdf5)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake)

include(CMake.vars)

if(DEBUG)
	set(CMAKE_BUILD_TYPE Debug)
else(DEBUG)
	set(CMAKE_BUILD_TYPE Release)
endif(DEBUG)


add_executable(${PROJECT_NAME}
	main.cpp
	mesh3dreader.cpp
	mesh.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../..)

# Judy
find_package(JUDY REQUIRED)
if(JUDY_FOUND)
	include_directories(${JUDY_INCLUDE_DIR})
	target_link_libraries(${PROJECT_NAME} ${JUDY_LIBRARY})
endif(JUDY_FOUND)

# HDF5
find_package(HDF5 REQUIRED)
if(HDF5_FOUND)
	include_directories(${HDF5_INCLUDE_DIR})
	target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARY})
endif(HDF5_FOUND)

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

