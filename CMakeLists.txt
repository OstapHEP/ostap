cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
message ( "----> CMAKE  version    : " ${CMAKE_VERSION} ) 

string(TIMESTAMP OSTAP_BUILD_TIMESTAMP   "%b %d, %Y, %H:%M [UTC]" UTC)
string(TIMESTAMP OSTAP_BUILD_UTCUNIX     "%s"                     UTC)

set(CMAKE_COLOR_DIAGNOSTIC ON)
set(CMAKE_COLOR_MAKEFILE   ON)

if ( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.15" )
   ## https://cmake.org/cmake/help/latest/policy/CMP0094.html
   cmake_policy ( SET CMP0094 NEW ) 
endif() 

if ( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.31" )
   ## https://cmake.org/cmake/help/latest/policy/CMP0177.html
   cmake_policy ( SET CMP0177 NEW ) 
endif() 

if ( ${CMAKE_CXX_COMPILER_LAUNCHER} )
else () 
  find_program(CCACHE_PROGRAM ccache)
  if(CCACHE_PROGRAM)
    # Support Unix Makefiles and Ninja
    ## set_property(GLOBAL PROPERTY CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM} )
  endif()
endif() 

project(ostap LANGUAGES CXX)
include(CTest)

## https://patorjk.com/software/taag/#p=display&f=Roman&t=%20Ostap3
     
message ( " ********************************************************************** " ) 
message ( " *                                                                    * " )
message ( " *      .oooooo.                .                          .oooo.     * " ) 
message ( " *     d8P'  `Y8b             .o8                        .dP\"\"Y88b    * " )
message ( " *    888      888  .oooo.o .o888oo  .oooo.   oo.ooooo.        ]8P'   * " )
message ( " *    888      888 d88(  \"8   888   `P  )88b   888' `88b     <88b.    * " )
message ( " *    888      888 `\"Y88b.    888    .oP\"888   888   888      `88b.   * " )  
message ( " *    `88b    d88' o.  )88b   888 . d8(  888   888   888 o.   .88P    * " ) 
message ( " *     `Y8bood8P'  8\"\"888P'   \"888\" `Y888\"\"8o  888bod8P' `8bd88P'     * " ) 
message ( " *                                             888                    * " ) 
message ( " *                                            o888o                   * " )
message ( " *                                                                    * " ) 
message ( " ********************************************************************** " )


set(OSTAP_VERSION_MAJOR  3)
set(OSTAP_VERSION_MINOR  2)
set(OSTAP_VERSION_PATCH  1)
set(OSTAP_VERSION_TWEAK  0)

set(OSTAP_VERSION ${OSTAP_VERSION_MAJOR}.${OSTAP_VERSION_MINOR}.${OSTAP_VERSION_PATCH}.${OSTAP_VERSION_TWEAK})

message ( "----> OSTAP           version    : " ${OSTAP_VERSION} )
message ( "----> OSTAP           timestamp  : " ${OSTAP_BUILD_TIMESTAMP} )

set(OSTAP_PYTHON_SITE_DIR "")

cmake_host_system_information(RESULT OSTAP_HOSTNAME QUERY HOSTNAME)
cmake_host_system_information(RESULT OSTAP_DOMAIN   QUERY FQDN    )

message ( "----> Host name                  : " ${OSTAP_HOSTNAME}  )
message ( "----> Domain                     : " ${OSTAP_DOMAIN}    )

cmake_host_system_information(RESULT OSTAP_PROCESSOR_NAME            QUERY PROCESSOR_NAME)
cmake_host_system_information(RESULT OSTAP_PROCESSOR_DESCRIPTION     QUERY PROCESSOR_DESCRIPTION)
cmake_host_system_information(RESULT OSTAP_TOTAL_PHYSICAL_MEMORY     QUERY TOTAL_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT OSTAP_AVAILABLE_PHYSICAL_MEMORY QUERY AVAILABLE_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT OSTAP_NUMBER_OF_LOGICAL_CORES   QUERY NUMBER_OF_LOGICAL_CORES )
cmake_host_system_information(RESULT OSTAP_NUMBER_OF_PHYSICAL_CORES  QUERY NUMBER_OF_PHYSICAL_CORES )
    
message ( "----> Processor name             : " ${OSTAP_PROCESSOR_NAME}            )
message ( "----> Processor description      : " ${OSTAP_PROCESSOR_DESCRIPTION}     )
message ( "----> Logical  cores             : " ${OSTAP_NUMBER_OF_LOGICAL_CORES}   )
message ( "----> Physical cores             : " ${OSTAP_NUMBER_OF_PHYSICAL_CORES}  )
message ( "----> Total     memory [MiB]     : " ${OSTAP_TOTAL_PHYSICAL_MEMORY}     )
message ( "----> Available memory [MiB]     : " ${OSTAP_AVAILABLE_PHYSICAL_MEMORY} )
        
cmake_host_system_information(RESULT OSTAP_IS_64BIT QUERY IS_64BIT)
message ( "----> 64-bit     ?               : " ${OSTAP_IS_64BIT} )

cmake_host_system_information(RESULT OSTAP_HAS_FPU QUERY HAS_FPU)
message ( "----> FPU        ?               : " ${OSTAP_HAS_FPU}  ) 
    
cmake_host_system_information(RESULT OSTAP_HAS_MMX QUERY HAS_MMX)
message ( "----> MMX        ?               : " ${OSTAP_HAS_MMX} )

cmake_host_system_information(RESULT OSTAP_HAS_MMX_PLUS QUERY HAS_MMX_PLUS)
message ( "----> MMX+       ?               : " ${OSTAP_HAS_MMX_PLUS} ) 

cmake_host_system_information(RESULT OSTAP_HAS_SSE QUERY HAS_SSE)
message ( "----> SSE        ?               : " ${OSTAP_HAS_SSE} )

cmake_host_system_information(RESULT OSTAP_HAS_SSE2 QUERY HAS_SSE2)
message ( "----> SSE2       ?               : " ${OSTAP_HAS_SSE2} )

cmake_host_system_information(RESULT OSTAP_HAS_SSE_FP QUERY HAS_SSE_FP)
message ( "----> SSE_FP     ?               : " ${OSTAP_HAS_SSE_FP} ) 

cmake_host_system_information(RESULT OSTAP_HAS_SSE_MMX QUERY HAS_SSE_MMX)
message ( "----> SSE_MMX    ?               : " ${OSTAP_HAS_SSE_MMX} )

cmake_host_system_information(RESULT OSTAP_AMD_3DNOW QUERY HAS_AMD_3DNOW)
message ( "----> AMD 3DNOW  ?               : " ${OSTAP_AMD_3DNOW} ) 

cmake_host_system_information(RESULT OSTAP_AMD_3DNOW_PLUS QUERY HAS_AMD_3DNOW_PLUS)
message ( "----> AMD 3DNOW+ ?               : " ${OSTAP_AMD_3DNOW_PLUS} )
            
cmake_host_system_information(RESULT OSTAP_OS_NAME     QUERY OS_NAME)
cmake_host_system_information(RESULT OSTAP_OS_PLATFORM QUERY OS_PLATFORM)
cmake_host_system_information(RESULT OSTAP_OS_VERSION  QUERY OS_VERSION)
cmake_host_system_information(RESULT OSTAP_OS_RELEASE  QUERY OS_RELEASE)

message ( "----> OS name                    : " ${OSTAP_OS_NAME}     )
message ( "----> OS platform                : " ${OSTAP_OS_PLATFORM} )
message ( "----> OS version                 : " ${OSTAP_OS_VERSION}  )
message ( "----> OS release                 : " ${OSTAP_OS_RELEASE}  )

    
# Use compiler from PATH, not the system one
# http://stackoverflow.com/a/29904501/504346
find_program(CMAKE_C_COMPILER   NAMES $ENV{CC}  gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

message ( "----> C-compiler                 : " $ENV{CC}  )
message ( "----> C++-compiler               : " $ENV{CXX} )
if ( CMAKE_C_COMPILER_ID ) 
  message ( "----> C-compiler      ID         : " ${CMAKE_C_COMPILER_ID}   )
endif() 
if ( CMAKE_CXX_COMPILER_ID ) 
  message ( "----> C++-compiler    ID         : " ${CMAKE_CXX_COMPILER_ID} )
endif() 
if ( CMAKE_C_COMPILER_VERSION ) 
  message ( "----> C-compiler      version    : " ${CMAKE_C_COMPILER_VERSION}   )
endif() 
if ( CMAKE_CXX_COMPILER_VERSION ) 
  message ( "----> C++-compiler    version    : " ${CMAKE_CXX_COMPILER_VERSION} ) 
endif()
if ( CMAKE_CXX_COMPILER_LAUNCHER )
  message ( "----> C++             launcher   : " ${CMAKE_CXX_COMPILER_LAUNCHER} ) 
endif() 
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
# - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
# - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
  set ( CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/install CACHE PATH "..." FORCE )
endif ( )

add_subdirectory(scripts)
add_subdirectory(source)
add_subdirectory(ostap)
## add_subdirectory(docs)

configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR})

### where to put the libraries
set(LIBRARY_OUTPUT_PATH    "${CMAKE_BINARY_DIR}/lib")
### where to put the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

# ===============================================================================
# build a CPack driven installer package
# ===============================================================================

include (InstallRequiredSystemLibraries)

set (CPACK_PACKAGE_VENDOR              "Ostap team" )
set (CPACK_PACKAGE_CONTACT             "Ivan.Belyaev@cern.ch"                  ) 
set (CPACK_OUTPUT_FILE_PREFIX          ${CMAKE_INSTALL_PREFIX}/packages        )
set (CPACK_RESOURCE_FILE_LICENSE       "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE"   )
set (CPACK_RESOURCE_FILE_README        "${CMAKE_CURRENT_SOURCE_DIR}/README.md" )
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Ostap: user-friendly (Py)ROOT-based analysis environment"  )
set (CPACK_PACKAGE_VERSION_MAJOR       "${OSTAP_VERSION_MAJOR}" )
set (CPACK_PACKAGE_VERSION_MINOR       "${OSTAP_VERSION_MINOR}" )
set (CPACK_PACKAGE_VERSION_PATCH       "${OSTAP_VERSION_PATCH}" )
set (CPACK_PACKAGE_VERSION_TWEAK       "${OSTAP_VERSION_TWEAK}" )
set (CPACK_PACKAGE_VERSION             "${CPACK_PACKAGE_VERSION_MAJOR}-${CPACK_PACKAGE_VERSION_MINOR}-${CPACK_PACKAGE_VERSION_PATCH}-${CPACK_PACKAGE_VERSION_TWEAK}" )

set (CPACK_PACKAGE_ICON             ${CMAKE_CURRENT_SOURCE_DIR}/docs/ostap.jpg )
if ( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.7" )
   set ( CPACK_PACKAGE_CHECKSUM  "SHA1" ) 
endif()

set (CPACK_COMPONENTS_ALL examples export headers libraries python scripts) 

set (CPACK_RPM_PACKAGE_LICENSE      "BSD-3-Clause" )
set (CPACK_RPM_PACKAGE_RELEASE_DIST  OK )

set (CPACK_RPM_FILE_NAME             "RPM-DEFAULT"  )
set (CPACK_DEBIAN_FILE_NAME          "DEB-DEFAULT"  )

## find_program( OSTAP_DPKH_EXECUTABLE NAMES dpkg     )        
## find_program( OSTAP_RPM_EXECUTABLE  NAMES rpmbuild )        

## if ( OSTAP_DPKG_EXECUTABLE AND OSTAP_RPM_EXECUTABLE )
##  set(CPACK_GENERATOR 7Z STGZ TBZ2 TGZ TZ TZST ZIP DEB RPM)
## elseif ( OSTAP_DPKG_EXECUTABLE )
##  set(CPACK_GENERATOR 7Z STGZ TBZ2 TGZ TZ TZST ZIP DEB )
## elseif ( OSTAP_RPM_EXECUTABLE  )
##  ## set(CPACK_GENERATOR 7Z STGZ TBZ2 TGZ TZ TZST ZIP RPM )
##  set(CPACK_GENERATOR STGZ RPM )
## else   ()
set(CPACK_GENERATOR 7Z STGZ TBZ2 TGZ TZ TZST ZIP )  
##endif ()

## RPM & DEB
set(CPACK_THREADS   0    )
 
include (CPack)
