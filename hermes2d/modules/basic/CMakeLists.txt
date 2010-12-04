project(module-basic)

add_library(hermes_basic SHARED basic.cpp)
target_link_libraries(hermes_basic ${HERMES_REAL_BIN})

add_executable(${PROJECT_NAME} main.cpp disc.cpp)
target_link_libraries(${PROJECT_NAME} hermes_basic)
include (../CMake.common)

install(TARGETS hermes_basic
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)
