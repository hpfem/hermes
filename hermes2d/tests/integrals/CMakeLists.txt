find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})

# examples
add_subdirectory(domain-perimeter)
