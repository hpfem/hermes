project(eigenvalue-h2p)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)
target_link_libraries(${PROJECT_NAME} hermes_schroedinger)
