project(HSS_vs_pMultigrid)

add_executable(${PROJECT_NAME} main.cpp definitions.cpp ../euler_util.cpp algorithms.cpp)

set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${HERMES_FLAGS})

target_link_libraries(${PROJECT_NAME} ${HERMES2D})
