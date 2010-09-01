find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})

# examples
add_subdirectory(lobatto-kernel-function-1)
add_subdirectory(lobatto-linearly-independent-1)
add_subdirectory(lobatto-zero-values-1)
add_subdirectory(lobatto-zero-values-2)
