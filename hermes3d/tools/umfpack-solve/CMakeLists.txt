project(umfpack-solve)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake)

include(CMake.vars)

if(DEBUG)
	set(CMAKE_BUILD_TYPE Debug)
else(DEBUG)
	set(CMAKE_BUILD_TYPE Release)
endif(DEBUG)


add_executable(${PROJECT_NAME}
	main.cpp
)

find_package(UMFPACK REQUIRED)
if(UMFPACK_FOUND)
	include_directories(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})
	target_link_libraries(${PROJECT_NAME} ${UMFPACK_LIBRARY} ${AMD_LIBRARY})
endif(UMFPACK_FOUND)



