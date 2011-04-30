// Windows DLL export/import definitions
#if defined(WIN32) || defined(_WINDOWS)
  // Visual Studio 2010.
  #if defined(EXPORT_HERMES_DLL)
    // when building DLL (target project defines this macro)
    #define TEUCHOS_LIB_DLL_EXPORT __declspec(dllexport)
  #else  
    // when using the DLL by a client project
    #define TEUCHOS_LIB_DLL_EXPORT __declspec(dllimport)
  #endif
#else 
  #define TEUCHOS_LIB_DLL_EXPORT
#endif