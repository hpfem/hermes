if (WITH_EXODUSII)

project(iron-water)

add_executable(${PROJECT_NAME} main.cpp)
include (../CMake.common)

endif (WITH_EXODUSII)
