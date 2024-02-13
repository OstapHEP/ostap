find_package(ROOT 6 CONFIG REQUIRED )

message ( "----> ROOT   version   : " ${ROOT_VERSION} )

# =============================================================================
## Locate proper Python/PythonLibs 
find_program(ROOT_CONFIG_EXECUTABLE NAMES root-config)


if(ROOT_VERSION VERSION_LESS "6.31.01")

execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --python2-version
                 OUTPUT_VARIABLE PY2VERSION_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

endif() 

execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --python3-version
                 OUTPUT_VARIABLE PY3VERSION_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

## message ('ROOT libraries: ' ${ROOT_LIBRARIES} )
## message ('PY2VERSION_ROOT:' ${PY2VERSION_ROOT} )
## message ('PY3VERSION_ROOT:' ${PY3VERSION_ROOT} )

add_library ( root_pyroot INTERFACE IMPORTED) 
if     (PY2VERSION_ROOT)
  find_package(Python2 ${PY2VERSION_ROOT} COMPONENTS Interpreter Development NumPy)
  set(Python_VERSION ${Python2_VERSION}  )
  if  (Python2_FOUND)
    message ( "Python version: " ${Python2_VERSION} '/' ${Python_VERSION})
  else()
    message ( "Python2 is NOT FOUND!" ) 
  endif()
  message ( "----> Python version    : " ${Python2_VERSION}      )
  message ( "----> Python executable : " ${Python2_EXECUTABLE}   )
  message ( "----> Python include    : " ${Python2_INCLUDE_DIRS} )
  message ( "----> Python libraries  : " ${Python2_LIBRARIES}    )
  set( PYTHON_INCLUDE_DIRS ${Python2_INCLUDE_DIRS} )
  set( PYTHON_LIBRARIES    ${Python2_LIBRARIES}    )
  set( PYTHON_VERSION      ${Python2_VERSION}      PARENT_SCOPE)
  set( PYTHON_EXECUTABLE   ${Python2_EXECUTABLE}   PARENT_SCOPE)
  set_target_properties( root_pyroot PROPERTIES INTERFACE_LINK_LIBRARIES "ROOT::PyROOT2;ROOT::ROOTTPython;Python2::Python")
elseif (PY3VERSION_ROOT)
  find_package(Python3 ${PY3VERSION_ROOT} COMPONENTS Interpreter Development NumPy)
  ## find_package(Python3 ${PY3VERSION_ROOT} COMPONENTS Interpreter Development)
  set(Python_VERSION ${Python3_VERSION}  )
  if  (Python3_FOUND)
    message ( "Python version: " ${Python3_VERSION} '/' ${Python_VERSION})
  else()
    message ( "Python3 is NOT FOUND!" ) 
  endif()  
  message ( "----> Python version    : " ${Python3_VERSION}      )
  message ( "----> Python executable : " ${Python3_EXECUTABLE}   )
  message ( "----> Python include    : " ${Python3_INCLUDE_DIRS} )
  message ( "----> Python libraries  : " ${Python3_LIBRARIES}    )
  set( PYTHON_INCLUDE_DIRS ${Python3_INCLUDE_DIRS} )
  set( PYTHON_LIBRARIES    ${Python3_LIBRARIES}    )
  set( PYTHON_VERSION      ${Python3_VERSION}      PARENT_SCOPE)
  set( PYTHON_EXECUTABLE   ${Python3_EXECUTABLE}   PARENT_SCOPE)

  if(ROOT_VERSION VERSION_LESS "6.31.01")
    set_target_properties( root_pyroot PROPERTIES INTERFACE_LINK_LIBRARIES "ROOT::PyROOT3;ROOT::ROOTTPython;Python3::Python")
  else  ()
    set_target_properties( root_pyroot PROPERTIES INTERFACE_LINK_LIBRARIES "ROOT::ROOTTPython;Python3::Python")
  endif () 

endif () 


## message ( "----> ROOT   version   : " ${ROOT_VERSION} )
                   
## get_target_property( pyincdirc root_pyroot INTERFACE_INCLUDE_DIRECTORIES)
## get_target_property( pyincdirc   root_pyroot  INTERFACE_LINK_LIBRARIES)
## message ('ppincdirc'  ${pyincdirc} ) 
 

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
if (GSL_FOUND) 
message ( "----> GSL    version   : " ${GSL_VERSION} )
message ( "----> GSL    include   : " ${GSL_INCLUDE_DIRS} )
message ( "----> GSL    libraries : " ${GSL_LIBRARIES} )
else()
endif() 

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  #message(STATUS "YES" "${CMAKE_CXX_COMPILER_ID}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-register")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build.config.in"
  "${CMAKE_CURRENT_BINARY_DIR}/build.config"
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

if(ROOT_VERSION VERSION_LESS_EQUAL "6.31.00")

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
elseif ( ${CXX14_ROOT} STREQUAL "yes" ) 
target_compile_features (ostap PUBLIC cxx_std_14 )
elseif ( ${CXX11_ROOT} STREQUAL "yes" ) 
target_compile_features (ostap PUBLIC cxx_std_11 )
endif() 

else() 

execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --cxxstandard
                 OUTPUT_VARIABLE CXX_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

set_property(TARGET ostap PROPERTY CXX_STANDARD ${CXX_ROOT} )

if     ( ${CXX_ROOT} STREQUAL "26" ) 
target_compile_features (ostap PUBLIC cxx_std_26 )
                        message ( '26' ) 
elseif ( ${CXX_ROOT} STREQUAL "23" ) 
target_compile_features (ostap PUBLIC cxx_std_23 )
                        message ( '23' ) 
elseif ( ${CXX_ROOT} STREQUAL "20" ) 
target_compile_features (ostap PUBLIC cxx_std_20 )
                        message ( '20' ) 
elseif ( ${CXX_ROOT} STREQUAL "17" ) 
target_compile_features (ostap PUBLIC cxx_std_17 )
##                         message ( '17' ) 
else() 
##target_compile_features (ostap PUBLIC cxx_std_17 )
##                        message ( '17/0' ) 
endif() 

endif() 






target_link_libraries   ( ostap ROOT::MathMore ROOT::ROOTVecOps ROOT::GenVector root_pyroot ROOT::RooFit ROOT::Hist ROOT::Tree ROOT::TreePlayer ROOT::RIO ROOT::TMVA ROOT::ROOTDataFrame GSL::gsl )

target_include_directories (ostap
    PUBLIC 
        $<INSTALL_INTERFACE:include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src ${Python_INCLUDE_DIRS} 
)

## get_target_property(incdirs1 ROOT::MathMore INTERFACE_INCLUDE_DIRECTORIES)
## message ( INCDIRS1 ${incdirs1} )
## get_target_property(incdirs2 root_pyroot    INTERFACE_INCLUDE_DIRECTORIES)
## message ( INCDIRS2 ${incdirs2} )
## get_target_property(incdirs3 root_pyroot    INTERFACE_LINK_LIBRARIES)
## message ( INCDIRS3 ${incdirs3} )

## include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${incdirs1} ${PYTHON_INCLUDE_DIRS} )

## ## make dictionaries 
## MAKE_DICT (ostap src/dict/Ostap.hh src/dict/Ostap.xml )
## REFLEX_BUILD_DICTIONARY ( ostap src/dict/Ostap.hh SELECTION src/dict/Ostap.xml )
## ROOT_GENERATE_DICTIONARY( G__ostapDict ${CMAKE_CURRENT_SOURCE_DIR}/src/dict/Ostap.hh 
##                          LINKDEF      ${CMAKE_CURRENT_SOURCE_DIR}/src/dict/Ostap.xml
##                          MODULE       ostap )
## 

REFLEX_GENERATE_DICTIONARY( ostap     ${CMAKE_CURRENT_SOURCE_DIR}/src/dict/Ostap.hh 
                            SELECTION ${CMAKE_CURRENT_SOURCE_DIR}/src/dict/Ostap.xml )

add_library           ( ostapDict MODULE ostap.cxx)
add_dependencies      ( ostapDict ostap-dictgen ostap ROOT::MathMore ROOT::GenVector root_pyroot )
target_link_libraries ( ostapDict               ostap ROOT::MathMore ROOT::GenVector root_pyroot )


## REFLEX_GENERATE_DICTIONARY( ${name} ${header} SELECTION ${selection})
## add_library( ${name}Dict MODULE ${name}.cxx)
## add_dependencies(${name}Dict ${name}-dictgen ostap ROOT::MathMore ROOT::GenVector root_pyroot)
## target_link_libraries   ( ${name}Dict ostap ROOT::MathMore ROOT::GenVector root_pyroot )
## endfunction


## add_library      (ostapDict SHARED G__ostapDict.cxx )
## add_dependencies (ostapDict ostap  G__ostapDict ) 
## target_include_directories (G__ostapDict
##         PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
## )
## 
## For clang we need to define __MATH_LONG_DOUBLE_CONSTANTS for
## defining M_El, M_PIl math constants.
## set(COMPILE_DEFINITIONS -D__MATH_LONG_DOUBLE_CONSTANTS)
## add_definitions(-D__MATH_LONG_DOUBLE_CONSTANTS)
## set_target_properties(ostap ostapDict PROPERTIES
## SUFFIX ".so"
## COMPILE_DEFINITIONS __MATH_LONG_DOUBLE_CONSTANTS)
## 
install ( TARGETS ostap     EXPORT   ostap-export 
                            LIBRARY  DESTINATION lib 
                            INCLUDES DESTINATION include )

install ( TARGETS ostapDict LIBRARY  DESTINATION lib )

install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/Ostap     
                    DESTINATION include 
                    FILES_MATCHING
                    PATTERN "*.h"
                    PATTERN "*.hpp"
                    PATTERN "*.icpp"
                    PATTERN  "*#*" EXCLUDE )
install ( FILES     ${CMAKE_CURRENT_BINARY_DIR}/Ostap/Config.h    DESTINATION include/Ostap )
install ( FILES     ${CMAKE_CURRENT_BINARY_DIR}/ostap_rdict.pcm   DESTINATION lib           )
install ( FILES     ${CMAKE_CURRENT_BINARY_DIR}/ostapDict.rootmap DESTINATION lib           )
install ( FILES     ${CMAKE_CURRENT_BINARY_DIR}/build.config      DESTINATION lib           )

install(EXPORT ostap-export
  FILE         OstapTargets.cmake
  NAMESPACE    ostap::
  DESTINATION  cmake/Ostap
)
