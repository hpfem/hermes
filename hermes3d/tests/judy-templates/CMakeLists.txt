project(judy-templates)
add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} ${JUDY_LIBRARY}) 

# Tests
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(judy-templates-1 ${BIN})

