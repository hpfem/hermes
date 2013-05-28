project(00-quickShow)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp definitions.h)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})
