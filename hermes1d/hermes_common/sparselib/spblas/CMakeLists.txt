include_directories(..)

add_library(spblas STATIC
    spmm.cc
    spsm.cc
)

set_property(TARGET spblas PROPERTY COMPILE_FLAGS "-g -fPIC")
