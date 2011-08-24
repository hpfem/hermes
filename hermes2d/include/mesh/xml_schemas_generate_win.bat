xsd cxx-tree --generate-doxygen --generate-ostream --hxx-suffix .h --cxx-suffix .cpp --root-element-first --generate-polymorphic --generate-serialization mesh_h2d_xml.xsd
xsd cxx-tree --generate-doxygen --generate-ostream --hxx-suffix .h --cxx-suffix .cpp --generate-serialization subdomains_h2d_xml.xsd
move /Y "*.cpp" ../../src/mesh