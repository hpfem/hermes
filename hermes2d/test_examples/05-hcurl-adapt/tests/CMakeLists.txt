project(test-tutorial-P04-04-hcurl)

add_executable(${PROJECT_NAME} main.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-tutorial-P04-04-hcurl ${BIN})
set_tests_properties(test-tutorial-P04-04-hcurl PROPERTIES LABELS slow)

