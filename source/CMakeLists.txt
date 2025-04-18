find_package(ROOT 6 CONFIG REQUIRED )

option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." TRUE)
if (${FORCE_COLORED_OUTPUT})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
       add_compile_options (-fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
       add_compile_options (-fcolor-diagnostics)
    endif ()
endif ()

if(ROOT_VERSION VERSION_LESS "6.24")
  message ( WARNING " Support for ROOT version " ${ROOT_VERSION} " will reach end soon!")
endif() 


if(ROOT_VERSION VERSION_GREATER_EQUAL "6.23.01")
    message ( "----> Use OSTAP/CMAKE setting for ROOT  >6.23")
    include(CMakeROOT_6_23.cmake)
elseif(ROOT_VERSION VERSION_GREATER_EQUAL "6.22")
    message ( FATAL_ERROR "FATAL ERROR: ----> ROOT version " ${ROOT_VERSION} " is not supported (new PyROOT)!" )
elseif(ROOT_VERSION VERSION_EQUAL "6.16.0")
    message ( "----> Use OSTAP/CMAKE setting for ROOT  6.16")
    include(CMakeROOT_6_16.cmake)
else()
    message ( "----> Use OSTAP/CMAKE setting for ROOT  6.14")
    include(CMakeROOT_6_14.cmake)
endif()


