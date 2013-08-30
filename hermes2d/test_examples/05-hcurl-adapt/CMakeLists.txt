project(05-hcurl-adapt)

add_executable(${PROJECT_NAME} main.cpp)

if(NOT MSVC)
  set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${HERMES_FLAGS})
endif()

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

# This is a mystery, if this line is not here, OpenMP is not used.
if(MSVC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
endif(MSVC)
