if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(example-butcher)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

#Possible Butcher's tables types are: 
#   Explicit_RK_1,               // Explicit Runge-Kutta RK-1, or explicit Euler method.
#   Implicit_RK_1,               // Implicit Runge-Kutta RK-1, or implicit Euler method.
#   Explicit_RK_2,               // Explicit Runge-Kutta RK-2 method.
#   Implicit_Crank_Nicolson_2_2, // Implicit Crank_Nicolson method.
#   Implicit_SDIRK_2_2,          // Implicit SDIRK-2-2 method.
#   Implicit_Lobatto_IIIA_2_2,   // Implicit Lobatto IIIA-2 method.
#   Implicit_Lobatto_IIIB_2_2,   // Implicit Lobatto IIIB-2 method.
#   Implicit_Lobatto_IIIC_2_2,   // Implicit Lobatto IIIB-2 method.
#   Explicit_RK_3,               // Explicit Runge-Kutta RK-3 method.
#   Explicit_RK_4,               // Explicit Runge-Kutta RK-4 method.
#   Implicit_Lobatto_IIIA_3_4,   // Implicit Lobatto IIIA-4 method.
#   Implicit_Lobatto_IIIB_3_4,   // Implicit Lobatto IIIB-4 method.
#   Implicit_Lobatto_IIIC_3_4,   // Implicit Lobatto IIIB-4 method.
#   Implicit_Radau_IIA_3_5       // Implicit Radau IIA-5 method.
#   Implicit_SIRK_2_2,           // Implicit SIRK-2-2 method.
#   Implicit_ESIRK_2_2,          // Implicit ESIRK-2-2 method.
#   Implicit_SDIRK_4_5           // Implicit SDIRK-4-5 method.


add_test(example-butcher-Explicit_RK_1 ${BIN} 1)
add_test(example-butcher-Implicit_RK_1 ${BIN} 2)

add_test(example-butcher-Explicit_RK_2 ${BIN} 3)
add_test(example-butcher-Implicit_Crank_Nicolson_2_2 ${BIN} 4)
add_test(example-butcher-Implicit_SDIRK_2_2 ${BIN} 5)
add_test(example-butcher-Implicit_Lobatto_IIIA_2_2 ${BIN} 6)
add_test(example-butcher-Implicit_Lobatto_IIIB_2_2 ${BIN} 7)
add_test(example-butcher-Implicit_Lobatto_IIIC_2_2 ${BIN} 8)

add_test(example-butcher-Explicit_RK_3 ${BIN} 9)

add_test(example-butcher-Explicit_RK_4 ${BIN} 10)
add_test(example-butcher-Implicit_Lobatto_IIIA_3_4 ${BIN} 11)
add_test(example-butcher-Implicit_Lobatto_IIIB_3_4 ${BIN} 12)
add_test(example-butcher-Implicit_Lobatto_IIIC_3_4 ${BIN} 13)

add_test(example-butcher-Implicit_Radau_IIA_3_5 ${BIN} 14)

add_test(example-butcher-Implicit_SIRK_2_2 ${BIN} 15)
add_test(example-butcher-Implicit_ESIRK_2_2 ${BIN} 16)
add_test(example-butcher-Implicit_SDIRK_4_5 ${BIN} 17)

