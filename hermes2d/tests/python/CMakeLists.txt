project(python-in-cpp)

add_executable(${PROJECT_NAME} main1.cpp)
include (../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(python-in-cpp ${BIN})

project(numpy-in-cpp)

add_executable(${PROJECT_NAME} main2.cpp)
include (../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(numpy-in-cpp ${BIN})
