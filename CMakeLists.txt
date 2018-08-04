cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(Ostap LANGUAGES CXX)

set(OSTAP_VERSION_MAJOR 1)
set(OSTAP_VERSION_MINOR 2)
set(OSTAP_VERSION_PATCH 0)
set(OSTAP_VERSION ${OSTAP_VERSION_MAJOR}.${OSTAP_VERSION_MINOR}.${OSTAP_VERSION_PATCH})


# Use compiler from PATH, not the system one
# http://stackoverflow.com/a/29904501/504346
find_program(CMAKE_C_COMPILER   NAMES $ENV{CC}  gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
# - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
# - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})


add_subdirectory(ostap)
add_subdirectory(source)

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

### where to put the libraries
set(LIBRARY_OUTPUT_PATH    "${CMAKE_BINARY_DIR}/lib")
### where to put the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

find_program(DOXYGEN_EXECUTABLE NAMES doxygen)
if ( EXISTS "${DOXYGEN_EXECUTABLE}" ) 
 find_package (ROOT) 
 find_program (DOT_EXECUTABLE dot)
 if ( EXISTS "${DOT_EXECUTABLE}" )
 set (HAVE_DOT YES)
 else() 
 set (HAVE_DOT NO)
 endif()
 string(APPEND ROOT_DOXYTAG ${ROOT_VERSION_MAJOR} ${ROOT_VERSION_MINOR})
 find_program (WGET_EXECUTABLE wget) 
 add_custom_target ( root_doxy_tag
                     COMMAND ${WGET_EXECUTABLE} -q http://root.cern.ch/doc/v${ROOT_DOXYTAG}/ROOT.tag 
                     COMMENT "(w)get ROOT.tag file" )
 set(DOXYGEN_IN  ${CMAKE_CURRENT_SOURCE_DIR}/docs/OstapDoxy.cfg.in)
 set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/OstapDoxy.cfg)
 configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)              
 add_custom_target(doxy_doc 
        DEPENDS root_doxy_tag
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT} 2> grep -e 'Fontconfig error: Cannot load default config file' -v 
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen (could be time consuming)"
        VERBATIM)
 ##install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doxygen_doc/html DESTINATION doc)
 install(FILES     ${CMAKE_CURRENT_BINARY_DIR}/ROOT.tag         DESTINATION doc)
 install(FILES     ${CMAKE_CURRENT_BINARY_DIR}/ostap_doxy.tag   DESTINATION doc)
endif() 
