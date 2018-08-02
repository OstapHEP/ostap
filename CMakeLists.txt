# Use compiler from PATH, not the system one
# http://stackoverflow.com/a/29904501/504346
find_program(CMAKE_C_COMPILER   NAMES $ENV{CC}  gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(Ostap)

set(OSTAP_VERSION_MAJOR 1)
set(OSTAP_VERSION_MINOR 2)
set(OSTAP_VERSION_PATCH 0)
set(OSTAP_VERSION ${OSTAP_VERSION_MAJOR}.${OSTAP_VERSION_MINOR}.${OSTAP_VERSION_PATCH})

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
# - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
# - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

### where to put the libraries
set(LIBRARY_OUTPUT_PATH     "${CMAKE_BINARY_DIR}/lib")
### where to put the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

add_subdirectory(ostap)
add_subdirectory(source)

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

find_program(DOXYGEN_EXECUTABLE NAMES doxygen)
if ( EXISTS "${DOXYGEN_EXECUTABLE}" ) 
 set(DOXYGEN_IN  ${CMAKE_CURRENT_SOURCE_DIR}/docs/OstapDoxy.cfg.in)
 set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/OstapDoxy.cfg)
 configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
 add_custom_target(doxy_doc  
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen (could be *very* time consuming)"
        VERBATIM )
endif() 
