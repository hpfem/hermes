# MSVC (Win) helper macros

# Makes Win32 path from Unix-style patch which is used by CMAKE. Used when a path is provided to an OS utility.
macro(MAKE_PATH PATH_OUT PATH_IN)
  if(WIN32)
    string(REPLACE "/" "\\" ${PATH_OUT} ${PATH_IN})
  else(WIN32)
    set(${PATH_OUT} ${PATH_IN})
  endif(WIN32)
endmacro(MAKE_PATH)

# This ensures that a .dll library is built for both debug and release configurations under MSVC.
macro(ADD_MSVC_BUILD_FLAGS LIB LIB_DEBUG LIB_RELEASE)
  get_target_property(FLAGS ${LIB} COMPILE_FLAGS)  
  set_target_properties(${LIB} PROPERTIES COMPILE_FLAGS "-DEXPORT_HERMES_DLL ${FLAGS}")
  set_target_properties(${LIB} PROPERTIES DEBUG_OUTPUT_NAME "${LIB_DEBUG}")
  set_target_properties(${LIB} PROPERTIES RELEASE_OUTPUT_NAME ${LIB_RELEASE})
endmacro(ADD_MSVC_BUILD_FLAGS)

# Installs a library to directories relative to CMAKE_INSTALL_PREFIX.
macro(INSTALL_LIB LIB)
	install(TARGETS ${LIB} 
				RUNTIME DESTINATION ${TARGET_ROOT}/bin 
				LIBRARY DESTINATION ${TARGET_ROOT}/lib 
				ARCHIVE DESTINATION ${TARGET_ROOT}/lib)
	if (MSVC)
	  MAKE_PATH(TARGET_DIR "${TARGET_ROOT}/bin")
    get_target_property(SOURCE_DEBUG_FILE ${LIB} LOCATION_Debug)
    MAKE_PATH(SOURCE_DEBUG_FILE ${SOURCE_DEBUG_FILE})
    get_target_property(SOURCE_RELEASE_FILE ${LIB} LOCATION_Release)
    MAKE_PATH(SOURCE_RELEASE_FILE ${SOURCE_RELEASE_FILE})
    add_custom_command(TARGET ${LIB}
      POST_BUILD
			COMMAND if not exist ${TARGET_DIR} mkdir ${TARGET_DIR}
      COMMAND if exist ${SOURCE_DEBUG_FILE} copy /Y ${SOURCE_DEBUG_FILE} ${TARGET_DIR}
      COMMAND if exist ${SOURCE_RELEASE_FILE} copy /Y ${SOURCE_RELEASE_FILE} ${TARGET_DIR})
    unset(TARGET_DIR)
	endif(MSVC)
endmacro(INSTALL_LIB)
