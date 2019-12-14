find_package(ROOT 6 CONFIG REQUIRED )

message  ("ROOT version:" ${ROOT_VERSION})
if(ROOT_VERSION VERSION_EQUAL "6.16.0")
    include(CMakeROOT_6_16.cmake)
else()
    include(CMakeROOT_6_14.cmake)
endif()
