IF(H2D_COMPLEX)

project(benchmark-screen)

# use the complex version of the library:
set(HERMES ${HERMES_CPLX_BIN})

add_executable(${PROJECT_NAME} main.cpp fresnl.c)
include (../../CMake.common)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(benchmark-screen ${BIN})
set_tests_properties(benchmark-screen PROPERTIES LABELS slow)

ENDIF(H2D_COMPLEX)
