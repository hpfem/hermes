project(test-time_period)

add_executable(${PROJECT_NAME} main.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES_COMMON_LIB} ${TRILINOS_LIBRARIES})

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})