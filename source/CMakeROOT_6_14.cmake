# root_generate_dictionary(ostap_dict ${CMAKE_CURRENT_SOURCE_DIR}/dict/Dict.h ${CMAKE_CURRENT_SOURCE_DIR}/dict/selections.xml)

## find_package(ROOT 6 CONFIG REQUIRED)
## find_package(ROOT 6 CONFIG REQUIRED COMPONENTS Smatrix Core MathCore MathMore Minuit2 GenVector Hist Matrix RIO TMVA Tree Thread TreePlayer RooFit RooFitCore PyROOT)
find_package(ROOT 6 CONFIG REQUIRED COMPONENTS Smatrix Core MathCore MathMore Minuit2 GenVector Hist Matrix RIO TMVA Tree Thread TreePlayer RooFit RooFitCore PyROOT)

if (ROOT_FOUND) 
message ( "----> ROOT   version   : " ${ROOT_VERSION} )
message ( "----> ROOT   include   : " ${ROOT_INCLUDE_DIRS} )
message ( "----> ROOT   libraries : " ${ROOT_LIBRARIES} )
endif()

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

find_package(GSL REQUIRED GSL_ROOT_DIR)
if (GSL_FOUND) 
message ( "----> GSL    version   : " ${GSL_VERSION} )
message ( "----> GSL    include   : " ${GSL_INCLUDE_DIRS} )
message ( "----> GSL    libraries : " ${GSL_LIBRARIES} )
else()
endif() 
# =============================================================================
## Locate Python/PythonLibs 
## since we are using MathMore, GSL *must be* somewhere around
## try with gsl-config in the path 
find_program(ROOT_CONFIG_EXECUTABLE NAMES root-config )        
execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --python-version
                 OUTPUT_VARIABLE PYVERSION_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

if ( ${CMAKE_VERSION}  VERSION_LESS "3.12") 
   find_package(Python ${PYVERSION_ROOT} COMPONENTS Interpreter Development )
   message ( "----> Python version    : " ${Python_VERSION}      )
   message ( "----> Python executable : " ${Python_EXECUTABLE}   )
   message ( "----> Python include    : " ${Python_INCLUDE_DIRS} )
   message ( "----> Python libraries  : " ${Python_LIBRARIES}    )
   set( PYTHON_INCLUDE_DIRS ${Python_INCLUDE_DIRS} )
   set( PYTHON_LIBRARIES    ${Python_LIBRARIES}    )
   set( PYTHON_VERSION      ${Python_VERSION}      PARENT_SCOPE)
   set( PYTHON_EXECUTABLE   ${Python_EXECUTABLE}   PARENT_SCOPE)
elseif ( ${PYVERSION_ROOT} VERSION_LESS "3.0" ) 
   find_package(Python2 ${PYVERSION_ROOT} COMPONENTS Interpreter Development NumPy)
   message ( "----> Python version    : " ${Python2_VERSION}      )
   message ( "----> Python executable : " ${Python2_EXECUTABLE}   )
   message ( "----> Python include    : " ${Python2_INCLUDE_DIRS} )
   message ( "----> Python libraries  : " ${Python2_LIBRARIES}    )
   set( PYTHON_INCLUDE_DIRS ${Python2_INCLUDE_DIRS} )
   set( PYTHON_LIBRARIES    ${Python2_LIBRARIES}    )
   set( PYTHON_VERSION      ${Python2_VERSION}      )
   set( PYTHON_VERSION      ${Python2_VERSION}      PARENT_SCOPE)
   set( PYTHON_EXECUTABLE   ${Python2_EXECUTABLE}   PARENT_SCOPE)
else() 
   find_package(Python3 ${PYVERSION_ROOT} COMPONENTS Interpreter Development NumPy)
   message ( "----> Python version    : " ${Python3_VERSION}      )
   message ( "----> Python executable : " ${Python3_EXECUTABLE}   )
   message ( "----> Python include    : " ${Python3_INCLUDE_DIRS} )
   message ( "----> Python libraries  : " ${Python3_LIBRARIES}    )
   set( PYTHON_INCLUDE_DIRS ${Python3_INCLUDE_DIRS} )
   set( PYTHON_LIBRARIES    ${Python3_LIBRARIES}    )
   set( PYTHON_VERSION      ${Python3_VERSION}      PARENT_SCOPE)
   set( PYTHON_EXECUTABLE   ${Python3_EXECUTABLE}   PARENT_SCOPE)
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
   REFLEX_GENERATE_DICTIONARY( ${name} ${header} SELECTION ${selection} OPTIONS -D__MATH_LONG_DOUBLE_CONSTANT -Wno-inconsistent-missing-override)
   add_library( ${name}Dict MODULE ${name}.cxx)
   add_dependencies(${name}Dict ${name}-dictgen ostap)
   ##target_compile_features ( ${name}Dict PUBLIC cxx_std_14 )
   target_link_libraries   ( ${name}Dict ostap )
endfunction( MAKE_DICT )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})


include(CMake_OSTAP.cmake) 

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
elseif ( ${CXX17_ROOT} STREQUALS "yes" ) 
target_compile_features (ostap PUBLIC cxx_std_14 )
elseif ( ${CXX11_ROOT} STREQUALS "yes" ) 
target_compile_features (ostap PUBLIC cxx_std_11 )
endif() 


##target_compile_features    (ostap PUBLIC cxx_std_14 )
target_link_libraries      (ostap ${ROOT_LIBRARIES} ${GSL_LIBRARIES} ${PYTHON_LIBRARIES})

target_include_directories (ostap
    PUBLIC 
        $<INSTALL_INTERFACE:include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR} ${ROOT_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS}
)

## make dictionaries 
MAKE_DICT (ostap src/dict/Ostap.hh src/dict/Ostap.xml )
 
target_include_directories (ostapDict
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS}
)

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
install ( FILES     ${CMAKE_CURRENT_BINARY_DIR}/build.config      DESTINATION lib           )

install(EXPORT ostap-export
  FILE         OstapTargets.cmake
  NAMESPACE    ostap::
  DESTINATION  cmake/Ostap
)
