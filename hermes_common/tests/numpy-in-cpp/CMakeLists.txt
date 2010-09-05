include_directories(${hermes_common_SOURCE_DIR})

project(numpy-in-cpp)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${HERMES_COMMON} ${PYTHON_LIBRARIES})
add_dependencies(${PROJECT_NAME} ${HERMES_COMMON})




# tests:
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(numpy-in-cpp ${BIN})
