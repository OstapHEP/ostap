cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

## ## try to ccache when available 
## find_program(CCACHE_PROGRAM ccache)
## if(CCACHE_PROGRAM)
##     set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
## endif()

project(ostap LANGUAGES CXX)
include(CTest)

set(OSTAP_VERSION_MAJOR 1)
set(OSTAP_VERSION_MINOR 4)
set(OSTAP_VERSION_PATCH 8)
set(OSTAP_VERSION_TWEAK 7)

set(OSTAP_VERSION ${OSTAP_VERSION_MAJOR}.${OSTAP_VERSION_MINOR}.${OSTAP_VERSION_PATCH}.${OSTAP_VERSION_TWEAK})

set(OSTAP_PYTHON_SITE_DIR "")

# Use compiler from PATH, not the system one
# http://stackoverflow.com/a/29904501/504346
find_program(CMAKE_C_COMPILER   NAMES $ENV{CC}  gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
# - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
# - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  set ( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/install CACHE PATH "..." FORCE )
endif ( )

add_subdirectory(source)
add_subdirectory(ostap)
add_subdirectory(examples)
add_subdirectory(scripts)
add_subdirectory(docs)
configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR})

### where to put the libraries
set(LIBRARY_OUTPUT_PATH    "${CMAKE_BINARY_DIR}/lib")
### where to put the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set (CPACK_PACKAGE_VERSION_MAJOR "${OSTAP_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${OSTAP_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${OSTAP_VERSION_PATCH}")
set (CPACK_PACKAGE_VERSION_TWEAK "${OSTAP_VERSION_TWEAK}")
set (CPACK_PACKAGE_ICON           ${CMAKE_CURRENT_SOURCE_DIR}/docs/ostap.jpg )
include (CPack)

