include_directories(${hermes_common_SOURCE_DIR})

project(vector)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES} ${HERMES_COMMON})
set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS "-g")


# tests:
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(vector ${BIN})
