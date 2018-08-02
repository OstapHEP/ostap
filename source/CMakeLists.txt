# root_generate_dictionary(ostap_dict ${CMAKE_CURRENT_SOURCE_DIR}/dict/Dict.h ${CMAKE_CURRENT_SOURCE_DIR}/dict/selections.xml)

find_package(ROOT       REQUIRED COMPONENTS Smatrix Core MathCore MathMore Minuit2 GenVector Hist Matrix RIO TMVA Tree Thread TreePlayer RooFit RooFitCore PyROOT)

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
                       
# =============================================================================
## Locate Python/PythonLibs 
## since we are using MathMore, GSL *must be* somewhere around
## try with gsl-config in the path 
find_program(ROOT_CONFIG_EXECUTABLE NAMES root-config )        
execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --python-version
                 OUTPUT_VARIABLE PYVERSION_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

find_package(PythonInterp ${PYVERSION_ROOT} REQUIRED )
find_package(PythonLibs                     REQUIRED )


if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  #message(STATUS "YES" "${CMAKE_CXX_COMPILER_ID}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

if(PYTHONLIBS_FOUND)
  set(OSTAP_USE_PYTHONLIBS 1)
else() 
  set(OSTAP_USE_PYTHONLIBS 0)
endif()               

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/Ostap/Config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Ostap/Config.h"
  )

function( MAKE_DICT name header selection )
   REFLEX_GENERATE_DICTIONARY( ${name} ${header} SELECTION ${selection} OPTIONS -D__MATH_LONG_DOUBLE_CONSTANTS)
   add_library( ${name}Dict SHARED ${name}.cxx )
   add_dependencies(${name}Dict ${name}-dictgen)
   # set_target_properties( ${name} PROPERTIES COMPILE_FLAGS ${ROOT_CXX_FLAGS} )
   target_link_libraries( ${name}Dict ostap ${ROOT_LIBRARIES})
endfunction( MAKE_DICT )

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})


if(PYTHONLIBS_FOUND)
 include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})
else() 
 include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})
endif()               

# REFLEX_GENERATE_DICTIONARY(ostap_dict src/dict/Ostap.hh SELECTION src/dict/Ostap.xml)
MAKE_DICT (ostap src/dict/Ostap.hh src/dict/Ostap.xml )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ostap_rdict.pcm   DESTINATION ${LIBRARY_OUTPUT_PATH})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ostapDict.rootmap DESTINATION ${LIBRARY_OUTPUT_PATH})
#---Create a shared library 
add_library(ostap SHARED src/format.cpp
                         src/gauss.cpp
                         src/BSpline.cpp
                         src/Bernstein.cpp
                         src/Bernstein1D.cpp
                         src/Bernstein2D.cpp
                         src/Bernstein3D.cpp
                         src/Binomial.cpp
                         src/BreitWigner.cpp
                         src/Choose.cpp
                         src/Combine.cpp
                         src/Chi2Fit.cpp
                         src/EigenSystem.cpp   
                         src/Error2Exception.cpp   
                         src/Exception.cpp
                         src/Faddeeva.cpp 
                         src/Formula.cpp   
                         src/Fourier.cpp   
                         src/Funcs.cpp   
                         src/GSL_sentry.cpp 
                         src/GSL_utils.cpp 
                         src/Hesse.cpp
                         src/HistoDump.cpp
                         src/HistoInterpolation.cpp
                         src/HistoMake.cpp
                         src/HistoProject.cpp
                         src/HistoStat.cpp
                         src/IFuncs.cpp
                         src/Interpolation.cpp
                         src/Iterator.cpp
                         src/Lomont.cpp
                         src/LorentzVectorWithError.cpp
                         src/Math.cpp
                         src/Models.cpp
                         src/Models2D.cpp
                         src/MoreMath.cpp
                         src/Mute.cpp
                         src/NStatEntity.cpp
                         src/Notifier.cpp
                         src/OstapDataFrame.cpp
                         src/Peaks.cpp
                         src/PDFs.cpp
                         src/PDFs2D.cpp
                         src/PDFs3D.cpp
                         src/PhaseSpace.cpp
                         src/Point3DWithError.cpp
                         src/Polynomials.cpp   
                         src/Printable.cpp
                         src/PyFuncs.cpp ## temporarily removed 
                         src/PyPdf.cpp   ## temporarily removed 
                         src/PyIterator.cpp
                         src/PySelector.cpp
                         src/PySelectorWithCuts.cpp
                         src/Polarization.cpp
                         src/SFactor.cpp
                         src/StatEntity.cpp
                         src/StatVar.cpp
                         src/StatusCode.cpp
                         src/Tee.cpp
                         src/Tensors.cpp
                         src/Tmva.cpp
                         src/UStat.cpp
                         src/Valid.cpp
                         src/ValueWithError.cpp
                         src/Vector3DWithError.cpp
                         src/Workspace.cpp    
                         src/WStatEntity.cpp    
                         src/nSphere.cpp      
                         src/owens.cpp      
                        )    

target_link_libraries(ostap ${ROOT_LIBRARIES} ${GSL_LIBRARIES})

# For clang we need to define __MATH_LONG_DOUBLE_CONSTANTS for
# defining M_El, M_PIl math constants.
#set(COMPILE_DEFINITIONS -D__MATH_LONG_DOUBLE_CONSTANTS)
# add_definitions(-D__MATH_LONG_DOUBLE_CONSTANTS)
set_target_properties(ostap ostapDict PROPERTIES
  SUFFIX ".so"
  COMPILE_DEFINITIONS __MATH_LONG_DOUBLE_CONSTANTS)
