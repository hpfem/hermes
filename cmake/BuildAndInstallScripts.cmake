# Macro for generating classes for XML mesh parsing according to XSD
macro(GENERATE_XSD_FILES PROJECT_NAME HEADER_XML_FILE_OUTPUT SOURCE_XML_FILE_OUTPUT XSD_FILE TARGET_DIR)
  add_custom_target(${PROJECT_NAME} ALL DEPENDS ${HEADER_XML_FILE_OUTPUT} "src/${SOURCE_XML_FILE_OUTPUT}")

IF(WIN32)
  MAKE_PATH(PATH_FOR_MOVE "${PROJECT_SOURCE_DIR}/include/${SOURCE_XML_FILE_OUTPUT}")
  ADD_CUSTOM_COMMAND(
        COMMAND   ${XSD_BIN} ARGS cxx-tree --generate-doxygen --generate-ostream --hxx-suffix .h --cxx-suffix .cpp --root-element-first --generate-polymorphic --generate-serialization --output-dir include/${TARGET_DIR} ${XSD_FILE}
        COMMAND   move ARGS "/Y" "${PATH_FOR_MOVE}" "${PROJECT_SOURCE_DIR}\\src\\${TARGET_DIR}"
        OUTPUT   ${HEADER_XML_FILE_OUTPUT} "src/${SOURCE_XML_FILE_OUTPUT}")
ELSE()
  ADD_CUSTOM_COMMAND(
        COMMAND   ${XSD_BIN} ARGS cxx-tree --generate-doxygen --generate-ostream --hxx-suffix .h --cxx-suffix .cpp --root-element-first --generate-polymorphic --generate-serialization --output-dir include/${TARGET_DIR} ${XSD_FILE}
        COMMAND   mv ARGS "-f" "${PROJECT_SOURCE_DIR}/include/${SOURCE_XML_FILE_OUTPUT}" "${PROJECT_SOURCE_DIR}/src/${TARGET_DIR}/"
        OUTPUT   ${HEADER_XML_FILE_OUTPUT} "src/${SOURCE_XML_FILE_OUTPUT}")
ENDIF()
ADD_CUSTOM_COMMAND(
      TARGET    ${PROJECT_NAME}
      DEPENDS   ${HEADER_XML_FILE_OUTPUT} "src/${SOURCE_XML_FILE_OUTPUT}")

endmacro(GENERATE_XSD_FILES)

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
  set_target_properties(${LIB} PROPERTIES COMPILE_FLAGS "-DEXPORT_HERMES_DLL ${HERMES_FLAGS}")
  set_target_properties(${LIB} PROPERTIES DEBUG_OUTPUT_NAME "${LIB_DEBUG}")
  set_target_properties(${LIB} PROPERTIES RELEASE_OUTPUT_NAME ${LIB_RELEASE})
endmacro(ADD_MSVC_BUILD_FLAGS)

# Installs a library to directories relative to CMAKE_INSTALL_PREFIX.
macro(INSTALL_LIB LIB)
	install(TARGETS ${LIB} 
                RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
                LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
                ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endmacro(INSTALL_LIB)
