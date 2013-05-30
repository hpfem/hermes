project(13-FCT-solid_body_rotation)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp lumped_projection.cpp hp_adapt.cpp highOrder.cpp lowOrder.cpp fct.cpp reg_estimator.cpp prev_solution.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})
