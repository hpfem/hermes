include_directories(${hermes_common_SOURCE_DIR})

project(python-in-cpp)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${HERMES_COMMON})
add_dependencies(${PROJECT_NAME} ${HERMES_COMMON})




# tests:
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(python-in-cpp ${BIN})
