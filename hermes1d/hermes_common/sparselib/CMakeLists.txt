cmake_minimum_required(VERSION 2.6)

include_directories(mv)
include_directories(spblas)
add_definitions(-DCOMPLEX=std::complex<double>)

add_library(sparse
	compcol_double.cc
	comprow_double.cc
	coord_double.cc
	iotext_double.cc
	diagpre_double.cc
	icpre_double.cc
	ilupre_double.cc
	qsort_double.cc
	qsort_int.cc
	iohb_double.cc
	iohb.c
)

set_property(TARGET sparse PROPERTY COMPILE_FLAGS "-g -fPIC")

add_subdirectory(spblas)
add_subdirectory(mv)
