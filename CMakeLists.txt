# Use compiler from PATH, not the system one
# http://stackoverflow.com/a/29904501/504346
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(ostap)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

### where to put the libraries
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
### where to put the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")


add_subdirectory(source)
