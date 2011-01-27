project(module-schroedinger-oscillator)

include_directories(${hermes2d_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} hermes_schroedinger)
