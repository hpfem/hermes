project(test-solution-save_and_load)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-solution-save_and_load ${BIN})