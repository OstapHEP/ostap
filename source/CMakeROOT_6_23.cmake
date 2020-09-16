find_package(ROOT 6 CONFIG REQUIRED )

message ( "----> ROOT   version   : " ${ROOT_VERSION} )

# =============================================================================
## Locate proper Python/PythonLibs 
find_program(ROOT_CONFIG_EXECUTABLE NAMES root-config)


execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --python2-version
                 OUTPUT_VARIABLE PY2VERSION_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --python3-version
                 OUTPUT_VARIABLE PY3VERSION_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

## message ('ROOT libraries: ' ${ROOT_LIBRARIES} )
 
message ('-Python2 version: ' ${PY2VERSION_ROOT} ) 
message ('-Python3 version: ' ${PY3VERSION_ROOT} ) 

add_library ( root_pyroot INTERFACE IMPORTED) 
if     (PY2VERSION_ROOT)
  find_package(Python2 ${PY2VERSION_ROOT} COMPONENTS Interpreter Development NumPy)
  set(Python_VERSION ${Python2_VERSION}  )
  if  (Python2_FOUND)
    message ( "Python version: " ${Python3_VERSION} '/' ${Python_VERSION})
  else()
    message ( "Python2 is NOT FOUND!" ) 
  endif()
  set_target_properties( root_pyroot PROPERTIES INTERFACE_LINK_LIBRARIES "ROOT::PyROOT2;Python2::Python")
elseif (PY3VERSION_ROOT)
  find_package(Python3 ${PY3VERSION_ROOT} COMPONENTS Interpreter Development NumPy)
  set(Python_VERSION ${Python3_VERSION}  )
  if  (Python3_FOUND)
    message ( "Python version: " ${Python3_VERSION} '/' ${Python_VERSION})
  else()
    message ( "Python3 is NOT FOUND!" ) 
  endif()  
  set_target_properties( root_pyroot PROPERTIES INTERFACE_LINK_LIBRARIES "ROOT::PyROOT3;Python3::Python")
endif () 

set(PYTHON_VERSION ${Python_VERSION}  )

# =============================================================================
## Locate GSL 
## since we are using MathMore, GSL *must be* somewhere around
## try with gsl-config in the path 
if ( EXISTS "$ENV{GSL_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{GSL_ROOT_DIR}" GSL_ROOT_DIR )
  set( GSL_ROOT_DIR "${GSL_ROOT_DIR}" CACHE PATH "Prefix for GSL installation." )
endif()
if ( NOT EXISTS "${GSL_ROOT_DIR}" )
 find_program( GSL_CONFIG_EXECUTABLE NAMES gsl-config )        
 if( EXISTS "${GSL_CONFIG_EXECUTABLE}" )
  execute_process(
    COMMAND "${GSL_CONFIG_EXECUTABLE}" --prefix
    OUTPUT_VARIABLE GSL_ROOT_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE )
 endif()
endif()

find_package(GSL REQUIRED GSL_ROOT_DIR COMPONENTS gsl)
message ('GSL version:'     ${GSL_VERSION})
message ('GSL includes:'    ${GSL_INCLUDE_DIRS})
message ('GSL libraries:'   ${GSL_LIBRARIES})

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  #message(STATUS "YES" "${CMAKE_CXX_COMPILER_ID}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-register")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build.config.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/build.config"
)
                                                           
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/Ostap/Config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Ostap/Config.h"
  )

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})


function( MAKE_DICT name header selection )
   REFLEX_GENERATE_DICTIONARY( ${name} ${header} SELECTION ${selection})
   add_library( ${name}Dict MODULE ${name}.cxx)
   add_dependencies(${name}Dict ${name}-dictgen ostap ROOT::MathMore ROOT::GenVector root_pyroot)
   target_link_libraries   ( ${name}Dict ostap ROOT::MathMore ROOT::GenVector root_pyroot )
endfunction( MAKE_DICT )

## include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})


include(CMake_OSTAP.cmake) 


## set_target_properties(ostap
##     PROPERTIES
##     NO_SYSTEM_FROM_IMPORTED ON
##     )

execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --has-cxx17
                 OUTPUT_VARIABLE CXX17_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --has-cxx14
                 OUTPUT_VARIABLE CXX14_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --has-cxx11
                 OUTPUT_VARIABLE CXX11_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

if     ( ${CXX17_ROOT} STREQUAL "yes" ) 
target_compile_features (ostap PUBLIC cxx_std_17 )
message ( C++17 ) 
elseif ( ${CXX17_ROOT} STREQUALS "yes" ) 
target_compile_features (ostap PUBLIC cxx_std_14 )
message ( C++14 ) 
elseif ( ${CXX11_ROOT} STREQUALS "yes" ) 
target_compile_features (ostap PUBLIC cxx_std_11 )
message ( C++11 ) 
endif() 

target_link_libraries   (ostap ROOT::MathMore ROOT::ROOTVecOps ROOT::GenVector ROOT::ROOTTPython root_pyroot ROOT::RooFit ROOT::Hist ROOT::Tree ROOT::TreePlayer ROOT::RIO ROOT::TMVA ROOT::ROOTDataFrame ROOT::Core GSL::gsl )
###add_dependencies (ostap ROOT::ROOTVecOps ROOT::ROOTDataFrame ROOT::RooFit ROOT::MathMore ROOT::ROOTTPython root_pyroot ROOT::Core GSL::gsl )

target_include_directories (ostap
    PUBLIC 
        $<INSTALL_INTERFACE:include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src ${Python2_INCLUDE_DIRS} 
## ${GSL_INCLUDE_DIRS} 
## /cvmfs/sft-nightlies.cern.ch/lcg/nightlies/dev3python3/Fri/vdt/0.4.3/x86_64-centos7-gcc9-opt/include/
)

get_target_property(incdirs1 ROOT::MathMore INTERFACE_INCLUDE_DIRECTORIES)
message ( INCDIRS1 ${incdirs1} )
get_target_property(incdirs2 root_pyroot    INTERFACE_INCLUDE_DIRECTORIES)
message ( INCDIRS2 ${incdirs2} )
get_target_property(incdirs3 root_pyroot    INTERFACE_LINK_LIBRARIES)
message ( INCDIRS3 ${incdirs3} )

## include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${incdirs1} ${PYTHON_INCLUDE_DIRS} )

## ## make dictionaries 
MAKE_DICT (ostap src/dict/Ostap.hh src/dict/Ostap.xml )

## For clang we need to define __MATH_LONG_DOUBLE_CONSTANTS for
## defining M_El, M_PIl math constants.
## set(COMPILE_DEFINITIONS -D__MATH_LONG_DOUBLE_CONSTANTS)
## add_definitions(-D__MATH_LONG_DOUBLE_CONSTANTS)
## set_target_properties(ostap ostapDict PROPERTIES
## SUFFIX ".so"
## COMPILE_DEFINITIONS __MATH_LONG_DOUBLE_CONSTANTS)

install ( TARGETS ostap     EXPORT   ostap-export 
                            LIBRARY  DESTINATION lib 
                            INCLUDES DESTINATION include )

install ( TARGETS ostapDict LIBRARY  DESTINATION lib )

install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/Ostap     DESTINATION include       )
install ( FILES     ${CMAKE_CURRENT_BINARY_DIR}/Ostap/Config.h    DESTINATION include/Ostap )
install ( FILES     ${CMAKE_CURRENT_BINARY_DIR}/ostap_rdict.pcm   DESTINATION lib           )
install ( FILES     ${CMAKE_CURRENT_BINARY_DIR}/ostapDict.rootmap DESTINATION lib           )

install(EXPORT ostap-export
  FILE         OstapTargets.cmake
  NAMESPACE    ostap::
  DESTINATION  cmake/Ostap
)
