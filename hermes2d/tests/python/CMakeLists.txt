project(python-in-cpp)

add_executable(${PROJECT_NAME} main1.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(python-in-cpp ${BIN})

project(numpy-in-cpp)

add_executable(${PROJECT_NAME} main2.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME}) 
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(numpy-in-cpp ${BIN})
