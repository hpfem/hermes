project(module-basic)

add_library(hermes_basic SHARED basic.cpp)

add_executable(${PROJECT_NAME} main.cpp disc.cpp)
target_link_libraries(${PROJECT_NAME} hermes_basic)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
install(TARGETS hermes_basic
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
