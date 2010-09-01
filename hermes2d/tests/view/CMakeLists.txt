find_package(JUDY REQUIRED)
include_directories(${JUDY_INCLUDE_DIR})

# views features
IF(NOT NOGLUT)
    # FIXME: disable this for now, as it fails to compile if Trilinos is
    # enabled (http://github.com/hpfem/hermes/issues#issue/1):
    add_subdirectory(zoom-to-fit)
ENDIF(NOT NOGLUT)

