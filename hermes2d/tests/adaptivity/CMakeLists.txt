find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})

# adaptivity tests
add_subdirectory(cand_proj)
