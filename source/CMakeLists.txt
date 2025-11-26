find_package(ROOT 6 CONFIG REQUIRED )

option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." TRUE)
if (${FORCE_COLORED_OUTPUT})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
       add_compile_options (-fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
       add_compile_options (-fcolor-diagnostics)
    endif ()
endif ()

if(ROOT_VERSION VERSION_LESS "6.26")
  message ( FATAL_ERROR " ROOT   version " ${ROOT_VERSION} "<6.26 is not supported anymore" )
endif() 
    
message ( "----> ROOT            version    : " ${ROOT_VERSION} )
       
find_program ( ROOT_CONFIG_EXECUTABLE NAMES root-config )
find_program ( ROOT_EXECUTABLE        NAMES root        )        
find_program ( ROOT_EXE_EXECUTABLE    NAMES root.exe    )        
message ( "----> ROOT            executable : " ${ROOT_EXECUTABLE} )
message ( "----> ROOT.EXE        executable : " ${ROOT_EXE_EXECUTABLE} )
message ( "----> ROOT config     executable : " ${ROOT_CONFIG_EXECUTABLE} )
message ( "----> ROOT C++        flags      : " ${ROOT_CXX_FLAGS} )
    
# =============================================================================
## Locate proper Python/PythonLibs 

execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --python3-version
                 OUTPUT_VARIABLE PY3VERSION_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

if (PY3VERSION_ROOT VERSION_LESS "3.8" )
  message ( FATAL_ERROR " Python version " ${PY3VERSION}   "<3.8 is no tsupported anymore")
endif() 

add_library ( root_pyroot INTERFACE IMPORTED)
    
find_package(Python3 ${PY3VERSION_ROOT} COMPONENTS Interpreter Development NumPy)
if  (Python3_FOUND)
else()
  message ( FATAL_ERROR "Python3 is NOT FOUND!" ) 
endif()  
if (Python3_VERSION  VERSION_LESS "3.8" )
  message ( FATAL_ERROR " Python version " ${Python3_VERSION} "<3.8 is not supported anymore")
endif() 

set(Python_VERSION ${Python3_VERSION}  )
message ( "----> PYTHON          version    : " ${Python3_VERSION}      )
message ( "----> PYTHON          executable : " ${Python3_EXECUTABLE}   )
message ( "----> PYTHON          include    : " ${Python3_INCLUDE_DIRS} )
message ( "----> PYTHON          libraries  : " ${Python3_LIBRARIES}    )
set( PYTHON_INCLUDE_DIRS ${Python3_INCLUDE_DIRS} )
set( PYTHON_LIBRARIES    ${Python3_LIBRARIES}    )
set( PYTHON_VERSION      ${Python3_VERSION}      PARENT_SCOPE)
set( PYTHON_EXECUTABLE   ${Python3_EXECUTABLE}   PARENT_SCOPE)
if(ROOT_VERSION VERSION_LESS "6.31.01")
  set_target_properties( root_pyroot PROPERTIES INTERFACE_LINK_LIBRARIES "ROOT::PyROOT3;ROOT::ROOTTPython;Python3::Python")
else  ()
  set_target_properties( root_pyroot PROPERTIES INTERFACE_LINK_LIBRARIES "ROOT::ROOTTPython;Python3::Python")
endif () 

# =============================================================================
## Locate GSL 
## since we are using MathMore, GSL *must be* somewhere around
## try with gsl-config in the path 
if ( EXISTS "$ENV{GSL_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{GSL_ROOT_DIR}" GSL_ROOT_DIR )
  set( GSL_ROOT_DIR "${GSL_ROOT_DIR}" CACHE PATH "Prefix for GSL installation." )
  message ( "----> GSL_ROOT_DIR               : " ${GSL_ROOT_DIR} )   
endif()
if ( NOT EXISTS "${GSL_ROOT_DIR}" )
 find_program( GSL_CONFIG_EXECUTABLE NAMES gsl-config )        
 if( EXISTS "${GSL_CONFIG_EXECUTABLE}" )
  message ( "----> GSL config      executable : " ${GSL_CONFIG_EXECUTABLE} )
  execute_process(
    COMMAND "${GSL_CONFIG_EXECUTABLE}" --prefix
    OUTPUT_VARIABLE GSL_ROOT_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE )
 endif()
endif()

find_package(GSL REQUIRED GSL_ROOT_DIR COMPONENTS gsl)
if (GSL_FOUND) 
  message ( "----> GSL             version    : " ${GSL_VERSION} )
  message ( "----> GSL             include    : " ${GSL_INCLUDE_DIRS} )
  message ( "----> GSL             libraries  : " ${GSL_LIBRARIES} )
else()
  message ( SEND_ERROR "----> No GSL libraries are found!" )
endif() 


set(CMAKE_CXX_FLAGS ${ROOT_CXX_FLAGS} )

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  #message(STATUS "YES" "${CMAKE_CXX_COMPILER_ID}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-register")
## set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

message ( "----> C++             flags      : " ${CMAKE_CXX_FLAGS} )

  
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build.config.in"
  "${CMAKE_CURRENT_BINARY_DIR}/.ostap.build"
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

if ( ${ROOT_CXX_STANDARD} )
  set ( CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD} )
  message ( "----> ROOT/C++        standard   : " ${ROOT_CXX_STANDARD} )
  message ( "----> C++             satndard   : " ${CMAKE_CXX_STANDARD} )
endif ()


include(CMake_OSTAP.cmake) 


if ( ${ROOT_CXX_STANDARD} )
  
  set_property ( TARGET ostap PROPERTY CXX_STANDARD ${CMAKE_CXX_STANDARD}  )
  
elseif(ROOT_VERSION VERSION_LESS_EQUAL "6.31.00")

  execute_process ( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --has-cxx26
                    OUTPUT_VARIABLE CXX26_ROOT
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
  execute_process ( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --has-cxx23
                    OUTPUT_VARIABLE CXX23_ROOT
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
  execute_process ( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --has-cxx20
                    OUTPUT_VARIABLE CXX20_ROOT
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
  execute_process ( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --has-cxx17
                    OUTPUT_VARIABLE CXX17_ROOT
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
  execute_process ( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --has-cxx14
                    OUTPUT_VARIABLE CXX14_ROOT
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
  execute_process ( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --has-cxx11
                  OUTPUT_VARIABLE CXX11_ROOT
                  OUTPUT_STRIP_TRAILING_WHITESPACE )

  if   ( ${CXX26_ROOT} STREQUAL "yes" )
     target_compile_features (ostap PUBLIC cxx_std_26 )
     message ( "----> C++             standard   : C++26 " )
     set_property(TARGET ostap PROPERTY CXX_STANDARD 26 )
  elseif ( ${CXX23_ROOT} STREQUAL "yes" ) 
     target_compile_features (ostap PUBLIC cxx_std_23 )
     message ( "----> C++             standard   : C++23 " )
     set_property(TARGET ostap PROPERTY CXX_STANDARD 23 )
  elseif ( ${CXX20_ROOT} STREQUAL "yes" ) 
     target_compile_features (ostap PUBLIC cxx_std_20 )
     message ( "----> C++             standard   : C++20 " )
     set_property(TARGET ostap PROPERTY CXX_STANDARD 20 )
  elseif ( ${CXX17_ROOT} STREQUAL "yes" ) 
     target_compile_features (ostap PUBLIC cxx_std_17 )
     message ( "----> C++             standard   : C++17 " )
     set_property(TARGET ostap PROPERTY CXX_STANDARD 17 )
  elseif ( ${CXX14_ROOT} STREQUAL "yes" ) 
     target_compile_features (ostap PUBLIC cxx_std_14 )
     message ( "----> C++             standard   : C++14 " )  
     set_property(TARGET ostap PROPERTY CXX_STANDARD 14 )
  elseif ( ${CXX11_ROOT} STREQUAL "yes" ) 
     target_compile_features (ostap PUBLIC cxx_std_11 )
     message ( "----> C++             standard   : C++11 " )
     set_property(TARGET ostap PROPERTY CXX_STANDARD 11 )
   else()
     target_compile_features (ostap PUBLIC cxx_std_17 )
     message ( "----> C++             standard   : C++17 " )    
     set_property(TARGET ostap PROPERTY CXX_STANDARD 17 )
  endif() 

else() 

   execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --cxxstandard
                    OUTPUT_VARIABLE CXX_ROOT
                    OUTPUT_STRIP_TRAILING_WHITESPACE )
    message ( "----> C++/ROOT        standard   : " ${CXX_ROOT} )
   if     ( ${CXX_ROOT} STREQUAL "26" ) 
      target_compile_features (ostap PUBLIC cxx_std_26 )
      message ( "----> C++             standard   : C++26 " )    
      set_property(TARGET ostap PROPERTY CXX_STANDARD 26 )
   elseif ( ${CXX_ROOT} STREQUAL "23" ) 
      target_compile_features (ostap PUBLIC cxx_std_23 )
      message ( "----> C++             standard   : C++23 " )    
      set_property(TARGET ostap PROPERTY CXX_STANDARD 23 )
   elseif ( ${CXX_ROOT} STREQUAL "20" ) 
      target_compile_features (ostap PUBLIC cxx_std_20 )  
      message ( "----> C++             standard   : C++20 " )    
      set_property(TARGET ostap PROPERTY CXX_STANDARD 20 )
   elseif ( ${CXX_ROOT} STREQUAL "17" ) 
      target_compile_features (ostap PUBLIC cxx_std_17 )
      message ( "----> C++             standard   : C++17 " )    
      set_property(TARGET ostap PROPERTY CXX_STANDARD 17 )
   elseif ( ${CXX_ROOT} STREQUAL "14" ) 
      target_compile_features (ostap PUBLIC cxx_std_14 )
      message ( "----> C++             standard   : C++14 " )    
      set_property(TARGET ostap PROPERTY CXX_STANDARD 14 )
   else() 
      ##target_compile_features (ostap PUBLIC cxx_std_17 )
      message ( "----> C++             standard   : C++17 " )    
      set_property(TARGET ostap PROPERTY CXX_STANDARD 17 )
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

REFLEX_GENERATE_DICTIONARY( ostap     ${CMAKE_CURRENT_SOURCE_DIR}/src/dict/Ostap.hh 
                            SELECTION ${CMAKE_CURRENT_SOURCE_DIR}/src/dict/Ostap.xml )

add_library           ( ostapDict MODULE ostap.cxx )
add_dependencies      ( ostapDict ostap-dictgen ostap ROOT::MathMore ROOT::GenVector root_pyroot )
target_link_libraries ( ostapDict               ostap ROOT::MathMore ROOT::GenVector root_pyroot )

install ( TARGETS   ostap
          EXPORT    ostap-export
          COMPONENT export 
          LIBRARY   DESTINATION lib 
          INCLUDES  DESTINATION include )

install ( TARGETS     ostapDict
          LIBRARY DESTINATION lib
          COMPONENT   libraries )

install ( DIRECTORY   ${CMAKE_CURRENT_SOURCE_DIR}/include/Ostap     
          DESTINATION include 
          FILES_MATCHING
          COMPONENT   headers 
          PATTERN     "*.h"
          PATTERN     "*.hpp"
          PATTERN     "*.icpp"
          PATTERN     "*#*"    EXCLUDE )
    
install ( FILES     ${CMAKE_CURRENT_BINARY_DIR}/Ostap/Config.h     DESTINATION include/Ostap COMPONENT headers    )
install ( FILES     ${CMAKE_CURRENT_BINARY_DIR}/ostap_rdict.pcm    DESTINATION lib           COMPONENT libraries  )  
install ( FILES     ${CMAKE_CURRENT_BINARY_DIR}/ostapDict.rootmap  DESTINATION lib           COMPONENT libraries  )
install ( FILES     ${CMAKE_CURRENT_BINARY_DIR}/.ostap.build       DESTINATION .             COMPONENT libraries  )

install(EXPORT ostap-export
  FILE         OstapTargets.cmake
  COMPONENT    export 
  NAMESPACE    ostap::
  DESTINATION  cmake/Ostap
)

enable_testing()

add_test ( NAME ostap.test_aclic_cpp1
           WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary2
	   COMMAND ${ROOT_EXE_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/src/test_ACLiC.cpp    -q -b -l -x ) 
add_test ( NAME ostap.test_aclic_cpp2
           WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary2
	   COMMAND ${ROOT_EXE_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/src/test_ACLiC.cpp+   -q -b -l -x ) 
add_test ( NAME ostap.test_aclic_cpp3
           WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary2
	   COMMAND ${ROOT_EXE_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/src/test_ACLiC.cpp++  -q -b -l -x ) 
