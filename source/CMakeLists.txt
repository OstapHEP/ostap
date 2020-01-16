find_package(ROOT 6 CONFIG REQUIRED )

message  ("ROOT 1 version:" ${ROOT_VERSION})
message  ("ROOT 2 version:" ${ROOT_FIND_VERSION})

if( ROOT_VERSION VERSION_LESS "6.16.0")
    include(CMakeROOT_6_14.cmake)
else()
    include(CMakeROOT_6_16.cmake)
endif()
