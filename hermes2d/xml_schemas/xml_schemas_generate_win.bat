xsd cxx-tree --generate-doxygen --generate-ostream --hxx-suffix .h --cxx-suffix .cpp --root-element-first --generate-polymorphic --generate-serialization --output-dir ../include/mesh mesh_h2d_xml.xsd
xsd cxx-tree --generate-doxygen --generate-ostream --hxx-suffix .h --cxx-suffix .cpp --root-element-first --generate-polymorphic --generate-serialization --output-dir ../include/mesh mesh_h1d_xml.xsd
xsd cxx-tree --generate-doxygen --generate-ostream --hxx-suffix .h --cxx-suffix .cpp --root-element-first --generate-polymorphic --generate-serialization --output-dir ../include/mesh subdomains_h2d_xml.xsd
move /Y "..\include\mesh\*.cpp" ..\src\mesh
xsd cxx-tree --generate-doxygen --generate-ostream --hxx-suffix .h --cxx-suffix .cpp --generate-serialization --output-dir ../include/space space_h2d_xml.xsd
move /Y "..\include\space\*.cpp" ..\src\space
xsd cxx-tree --generate-doxygen --generate-ostream --hxx-suffix .h --cxx-suffix .cpp --generate-serialization --output-dir ../include/function solution_h2d_xml.xsd
move /Y "..\include\function\*.cpp" ..\src\function