project(lobatto-zero-values-1)

add_executable(${PROJECT_NAME} 
        main.cpp
)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(lobatto-zero-values-1 ${BIN})
