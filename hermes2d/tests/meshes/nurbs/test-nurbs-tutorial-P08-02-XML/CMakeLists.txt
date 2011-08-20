project(test-nurbs-tutorial-P08-02-XML)

add_executable(${PROJECT_NAME} main.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})
target_link_libraries(${PROJECT_NAME} ${HERMES2D})

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-tutorial-P08-02-nurbs-1 ${BIN} 1)
add_test(test-tutorial-P08-02-nurbs-2 ${BIN} 2)
add_test(test-tutorial-P08-02-nurbs-3 ${BIN} 3)