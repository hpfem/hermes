project(module-compressible-euler)

add_library(hermes_compressible_euler SHARED compressible-euler.cpp)
target_link_libraries(hermes_compressible_euler ${HERMES_REAL_BIN})

add_executable(${PROJECT_NAME} main.cpp disc.cpp)
target_link_libraries(${PROJECT_NAME} hermes_compressible_euler)
include (../CMake.common)

install(TARGETS hermes_compressible_euler
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
