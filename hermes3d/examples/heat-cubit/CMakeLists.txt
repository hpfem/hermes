if(WITH_EXODUSII)

project(heat-cubit)
add_executable(${PROJECT_NAME}	main.cpp)

include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

endif(WITH_EXODUSII)
