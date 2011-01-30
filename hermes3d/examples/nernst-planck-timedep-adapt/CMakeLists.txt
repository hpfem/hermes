if(WITH_EXODUSII)

project(newton-np-timedep-adapt-system-3d)
add_executable(${PROJECT_NAME}  main.cpp)
include (${hermes3d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})

endif(WITH_EXODUSII)

