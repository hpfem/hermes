project(test-tutorial-P04-03-complex)

add_executable(${PROJECT_NAME} main.cpp ../definitions.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-tutorial-P04-03-complex ${BIN})
set_tests_properties(test-tutorial-P04-03-complex PROPERTIES LABELS slow)

