include_directories(${hermes_common_SOURCE_DIR})

project(matrix)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES} ${HERMES_COMMON})



# tests:
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(matrix ${BIN})
