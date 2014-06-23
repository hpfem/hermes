project(mhd-example)

add_executable(${PROJECT_NAME} main.cpp numerical_flux.cpp numerical_flux.h weak_formulation.cpp weak_formulation.h util.cpp util.h)

if(NOT MSVC)
  set_property(TARGET ${PROJECT_NAME} PROPERTY COMPILE_FLAGS ${HERMES_FLAGS})
endif()

target_link_libraries(${PROJECT_NAME} ${HERMES2D})
