project(module-basic)

add_library(hermes_basic SHARED basic.cpp)
target_link_libraries(hermes_basic ${HERMES2D_REAL})

add_executable(${PROJECT_NAME} main.cpp disc.cpp)
target_link_libraries(${PROJECT_NAME} hermes_basic)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(hermes_basic)
set_common_target_properties(module-basic)
install(TARGETS hermes_basic
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
