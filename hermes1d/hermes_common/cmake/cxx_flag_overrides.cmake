# This overrides CXX flags for MSVC
if(MSVC)
    set(MSVC_DEFINES "/DWIN32 /D_WINDOWS /Dpopen=_popen /Dpclose=_pclose /D__value=_value /Dfinite=_finite /Dhypot=_hypot /Disatty=_isatty /Dfileno=_fileno /D_CRT_SECURE_NO_WARNINGS /DYY_NO_UNISTD_H /D_USE_MATH_DEFINES /DIMPLELENT_C99 /wd4275 /wd4251")
    set(CMAKE_CXX_FLAGS_DEBUG_INIT          "/D_DEBUG /Od /Ob2 /MDd /Zi /RTC1 ${MSVC_DEFINES}")
    set(CMAKE_CXX_FLAGS_MINSIZEREL_INIT     "/O2 /Ob2 /MD ${MSVC_DEFINES}")
    set(CMAKE_CXX_FLAGS_RELEASE_INIT        "/O2 /Ob2 /MD ${MSVC_DEFINES}")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "/O2 /Ob2 /MD /Zi ${MSVC_DEFINES}")
else(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG_INIT          "${CMAKE_CXX_FLAGS_DEBUG_INIT} -D_DEBUG")
endif(MSVC)

