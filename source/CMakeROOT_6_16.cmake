find_package(ROOT 6 CONFIG REQUIRED )

# =============================================================================
## Locate Python/PythonLibs 
find_program(ROOT_CONFIG_EXECUTABLE NAMES root-config )        
execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --python-version
                 OUTPUT_VARIABLE PYVERSION_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

find_package(PythonInterp ${PYVERSION_ROOT} REQUIRED )
find_package(PythonLibs                     REQUIRED )

## message ('ROOT   version:' ${ROOT_VERSION})
## message ('PYTHON version:' ${PYTHONLIBS_VERSION_STRING})

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  #message(STATUS "YES" "${CMAKE_CXX_COMPILER_ID}")
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


#---Create a shared library 
add_library(ostap SHARED src/format.cpp
                         src/gauss.cpp
                         src/AddBranch.cpp
                         src/AddVars.cpp
                         src/BLOB.cpp
                         src/BSpline.cpp
                         src/Bernstein.cpp
                         src/Bernstein1D.cpp
                         src/Bernstein2D.cpp
                         src/Bernstein3D.cpp
                         src/Binomial.cpp
                         src/BreitWigner.cpp
                         src/ChebyshevApproximation.cpp
                         src/Choose.cpp
                         src/Combine.cpp
                         src/Chi2Fit.cpp
                         src/Dalitz.cpp
                         src/DalitzIntegrator.cpp
                         src/DataFrameUtils.cpp
                         src/EigenSystem.cpp   
                         src/Error2Exception.cpp   
                         src/Exception.cpp
                         src/Faddeeva.cpp 
                         src/Formula.cpp   
                         src/Fourier.cpp   
                         src/Funcs.cpp   
                         src/GetWeight.cpp 
                         src/GSL_sentry.cpp 
                         src/GSL_utils.cpp 
                         src/Hesse.cpp
                         src/HistoDump.cpp
                         src/HistoInterpolation.cpp
                         src/HistoMake.cpp
                         src/HistoProject.cpp
                         src/HistoStat.cpp
                         src/IFuncs.cpp
                         src/Integrator.cpp
                         src/Interpolation.cpp
                         src/Iterator.cpp
                         src/Kinematics.cpp
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
                         src/Parameterization.cpp
                         src/Params.cpp
                         src/Peaks.cpp
                         src/PDFs.cpp
                         src/PDFs2D.cpp
                         src/PDFs3D.cpp
                         src/PhaseSpace.cpp
                         src/Point3DWithError.cpp
                         src/Polynomials.cpp   
                         src/Printable.cpp
                         src/PyBLOB.cpp
                         src/PyCallable.cpp 
                         src/PyFuncs.cpp 
                         src/PyPdf.cpp   
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
                         src/Topics.cpp
                         src/Tmva.cpp
                         src/UStat.cpp
                         src/Valid.cpp
                         src/ValueWithError.cpp
                         src/Vector3DWithError.cpp
                         src/Workspace.cpp    
                         src/WStatEntity.cpp    
                         src/nSphere.cpp      
                         src/owens.cpp      
                         src/hcubature.cpp                         
                         src/pcubature.cpp
                        )

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
