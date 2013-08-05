project(08-nonlinearity)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${HERMES_FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})
