find_package(ROOT 6 CONFIG REQUIRED )

message  ("ROOT 1 version:" ${ROOT_VERSION})
message  ("ROOT 2 version:" ${ROOT_FIND_VERSION})

option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." TRUE)
if (${FORCE_COLORED_OUTPUT})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
       add_compile_options (-fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
       add_compile_options (-fcolor-diagnostics)
    endif ()
endif ()

if(ROOT_VERSION VERSION_EQUAL "6.16.0")
    include(CMakeROOT_6_16.cmake)
else()
    include(CMakeROOT_6_14.cmake)
endif()
