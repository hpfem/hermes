add_subdirectory("00-quickShow")

add_subdirectory("01-poisson")

add_subdirectory("02-poisson-newton")

add_subdirectory("03-navier-stokes")

add_subdirectory("04-complex-adapt")

add_subdirectory("05-hcurl-adapt")

add_subdirectory("06-system-adapt")

add_subdirectory("07-newton-heat-rk")

add_subdirectory("08-nonlinearity")

IF(WITH_TRILINOS)
	add_subdirectory("09-trilinos-nonlinear")
ENDIF(WITH_TRILINOS)

add_subdirectory("10-linear-advection-dg-adapt")

add_subdirectory("11-transient-adapt")

add_subdirectory("12-picard")

add_subdirectory("13-FCT")

add_subdirectory("14-error-calculation")