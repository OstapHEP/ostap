find_package(ROOT 6 CONFIG REQUIRED )

if (ROOT_FOUND) 
message ( "----> ROOT   version   : " ${ROOT_VERSION} )
message ( "----> ROOT   include   : " ${ROOT_INCLUDE_DIRS} )
message ( "----> ROOT   libraries : " ${ROOT_LIBRARIES} )
endif()

# =============================================================================
## Locate Python/PythonLibs 
find_program(ROOT_CONFIG_EXECUTABLE NAMES root-config )        
execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --python-version
                 OUTPUT_VARIABLE PYVERSION_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

if ( ${CMAKE_VERSION}  VERSION_LESS "3.12") 
   find_package(PythonInterp ${PYVERSION_ROOT} REQUIRED )
   find_package(PythonLibs                     REQUIRED )
   message ( "----> Python version    : " ${PYTHON_VERSION_STRING}      )
   message ( "----> Python executable : " ${PYTHON_EXECUTABLE}   )
   message ( "----> Python include    : " ${PYTHON_INCLUDE_DIRS} )
   message ( "----> Python libraries  : " ${PYTHON_LIBRARIES}    )
   ## set( PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} )
   ## set( PYTHON_LIBRARIES    ${PYTHON_LIBRARIES}    )
   set( PYTHON_VERSION      ${PYTHON_VERSION_STRING}     PARENT_SCOPE)
   set( PYTHON_EXECUTABLE   ${PYTHON_EXECUTABLE}         PARENT_SCOPE)
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



## message ('ROOT   version:' ${ROOT_VERSION})
## message ('PYTHON version:' ${PYTHONLIBS_VERSION_STRING})

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  #message(STATUS "YES" "${CMAKE_CXX_COMPILER_ID}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-register")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/Ostap/Config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Ostap/Config.h"
  )

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})


function( MAKE_DICT name header selection )
   REFLEX_GENERATE_DICTIONARY( ${name} ${header} SELECTION ${selection} OPTIONS -D__MATH_LONG_DOUBLE_CONSTANT -Wno-inconsistent-missing-override)
   add_library( ${name}Dict MODULE ${name}.cxx)
   add_dependencies(${name}Dict ${name}-dictgen ostap ROOT::MathMore ROOT::GenVector ROOT::PyROOT)
   target_link_libraries   ( ${name}Dict ostap ROOT::MathMore ROOT::GenVector ROOT::PyROOT )
endfunction( MAKE_DICT )

## include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})


include(CMake_OSTAP.cmake)

set_target_properties(ostap
    PROPERTIES
    NO_SYSTEM_FROM_IMPORTED ON
    )
target_compile_features    (ostap PUBLIC cxx_std_17 )
target_link_libraries      (ostap ROOT::MathMore ROOT::ROOTVecOps ROOT::GenVector ROOT::PyROOT ROOT::RooFit ROOT::Hist ROOT::Tree ROOT::TreePlayer ROOT::TMVA ROOT::ROOTDataFrame ${PYTHON_LIBRARIES})

target_include_directories (ostap
    PUBLIC 
        $<INSTALL_INTERFACE:include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src ${PYTHON_INCLUDE_DIRS}
)

get_target_property(incdirs1 ROOT::MathMore INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${incdirs1} ${PYTHON_INCLUDE_DIRS} )

## make dictionaries 
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
