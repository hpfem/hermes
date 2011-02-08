if(NOT H2D_REAL)
    return()
endif(NOT H2D_REAL)
project(example-butcher-adapt)

add_executable(${PROJECT_NAME} main.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})

#Possible Butcher's tables types are: 
# Embedded explicit methods:
#   Explicit_HEUN_EULER_2_12_embedded,          // Heun-Euler (orders 1 and 2),
#   Explicit_BOGACKI_SHAMPINE_4_23_embedded,    // Bogacki-Shampine (orders 2 and 3),
#   Explicit_FEHLBERG_6_45_embedded,            // Fehlberg (orders 4 and 5),
#   Explicit_CASH_KARP_6_45_embedded,           // Cash-Karp (orders 4 and 5),
#   Explicit_DORMAND_PRINCE_7_45_embedded,      // Dormand-Prince (orders 4 and 5),
#
# Embedded implicit methods:
#   Implicit_ESDIRK_TRBDF2_3_23_embedded, // From the paper Analysis and implementation of TR-BDF2 by Hosea and Shampine.
#   Implicit_ESDIRK_TRX2_3_23_embedded,   // From the paper Analysis and implementation of TR-BDF2 by Hosea and Shampine.
#   Implicit_SDIRK_CASH_3_23_embedded,    // From the paper Diagonally Implicit Runge-Kutta Formulae with Error Estimates by J.R. Cash
#   Implicit_SDIRK_CASH_5_24_embedded,    // From the paper Diagonally Implicit Runge-Kutta Formulae with Error Estimates by J.R. Cash
#   Implicit_SDIRK_CASH_5_34_embedded,    // From the paper Diagonally Implicit Runge-Kutta Formulae with Error Estimates by J.R. Cash
#   Implicit_DIRK_7_45_embedded           // Implicit embedded DIRK method pair of orders four in five (from the paper 
#                                         // Fudziah Ismail et all: Embedded Pair of Diagonally Implicit Runge-Kutta  
#                                         // Method for Solving Ordinary Differential Equations). The method has
#                                         // 7 stages but the first one is explicit. BEWARE - THE B2 ROW IS PROBABLY WRONG!

add_test(example-butcher-adapt-Explicit_HEUN_EULER_2_12_embedded ${BIN} 1)
add_test(example-butcher-adapt-Explicit_BOGACKI_SHAMPINE_4_23_embedded ${BIN} 2)
add_test(example-butcher-adapt-EExplicit_FEHLBERG_6_45_embedded ${BIN} 3)
add_test(example-butcher-adapt-IExplicit_CASH_KARP_6_45_embedded ${BIN} 4)
add_test(example-butcher-adapt-Explicit_DORMAND_PRINCE_7_45_embedded ${BIN} 5)

add_test(example-butcher-adapt-Implicit_ESDIRK_TRBDF2_3_23_embedded ${BIN} 6)
add_test(example-butcher-adapt-Implicit_ESDIRK_TRX2_3_23_embedded ${BIN} 7)
add_test(example-butcher-adapt-Implicit_SDIRK_CASH_3_23_embedded ${BIN} 8)
add_test(example-butcher-adapt-Implicit_SDIRK_CASH_5_24_embedded ${BIN} 9)
add_test(example-butcher-adapt-Implicit_SDIRK_CASH_5_34_embedded ${BIN} 10)
add_test(example-butcher-adapt-Implicit_DIRK_7_45_embedded ${BIN} 11)


