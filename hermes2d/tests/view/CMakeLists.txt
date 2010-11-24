find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})

# views features
IF(NOT NOGLUT)
  add_subdirectory(zoom-to-fit)
ENDIF(NOT NOGLUT)

