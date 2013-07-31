SET(XSD_INCLUDE_SEARCH_PATH
	${XSD_ROOT}/include
	/usr/include
	/usr/local/include
)

FIND_PATH(XSD_INCLUDE_DIR xsd/cxx/pre.hxx xsd/cxx/xml/dom/parsing-source.hxx xsd/cxx/post.hxx xsd/cxx/xml/sax/std-input-source.hxx xsd/cxx/tree/error-handler.hxx ${XSD_INCLUDE_SEARCH_PATH})

FIND_PROGRAM(XSD_BIN NAMES xsd xsdcxx)

# Report the found libraries, quit with fatal error if any required library has not been found.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(XSD DEFAULT_MSG XSD_INCLUDE_DIR XSD_BIN)