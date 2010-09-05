include_directories(${hermes_common_SOURCE_DIR})

project(leaks)
add_executable(${PROJECT_NAME} main.cpp cpp_api.cpp)
target_link_libraries(${PROJECT_NAME} ${HERMES_COMMON})
add_dependencies(${PROJECT_NAME} ${HERMES_COMMON})

cython_add_module(my_api2 cpp_api.cpp)

set_source_files_properties(main.cpp PROPERTIES
    OBJECT_DEPENDS ${leaks_SOURCE_DIR}/my_api2_api.h
    )
add_custom_command(
    OUTPUT ${leaks_SOURCE_DIR}/my_api2_api.h
    COMMAND cython my_api2.pyx
    )


# tests:
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(leaks ${BIN})
