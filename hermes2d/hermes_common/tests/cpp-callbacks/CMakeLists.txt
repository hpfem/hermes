include_directories(${hermes_common_SOURCE_DIR})

project(cpp-callbacks)
add_executable(${PROJECT_NAME} main.cpp cpp_api.cpp)
target_link_libraries(${PROJECT_NAME} ${HERMES_COMMON})
add_dependencies(${PROJECT_NAME} ${HERMES_COMMON})

cython_add_module(my_api cpp_api.cpp)

set_source_files_properties(main.cpp PROPERTIES
    OBJECT_DEPENDS ${cpp-callbacks_SOURCE_DIR}/my_api_api.h
    )
add_custom_command(
    OUTPUT my_api_api.h
    COMMAND cython my_api.pyx
    )

# tests:
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(cpp-callbacks ${BIN})
