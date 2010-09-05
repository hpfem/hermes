project(timer-1)

include_directories(${hermes_common_SOURCE_DIR})
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES} ${HERMES_COMMON})


set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(timer-1 ${BIN})

