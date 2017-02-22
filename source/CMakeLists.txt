# root_generate_dictionary(ostap_dict ${CMAKE_CURRENT_SOURCE_DIR}/dict/Dict.h ${CMAKE_CURRENT_SOURCE_DIR}/dict/selections.xml)

find_package(ROOT  REQUIRED COMPONENTS Smatrix Core MathCore MathMore GenVector Hist Matrix RIO Tree Thread TreePlayer RooFit RooFitCore PyROOT)
find_package(GSL   REQUIRED)


if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  #message(STATUS "YES" "${CMAKE_CXX_COMPILER_ID}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

function( MAKE_DICT name header selection )
   REFLEX_GENERATE_DICTIONARY( ${name} ${header} SELECTION ${selection} OPTIONS -D__MATH_LONG_DOUBLE_CONSTANTS)
   add_library( ${name}Dict SHARED ${name}.cxx )
   add_dependencies(${name}Dict ${name}-dictgen)
   #set_target_properties( ${name} PROPERTIES COMPILE_FLAGS ${ROOT_CXX_FLAGS} )
   target_link_libraries( ${name}Dict ostap ${ROOT_LIBRARIES})
endfunction( MAKE_DICT )

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS})
# REFLEX_GENERATE_DICTIONARY(ostap_dict src/dict/Ostap.hh SELECTION src/dict/Ostap.xml)
MAKE_DICT (ostap src/dict/Ostap.hh src/dict/Ostap.xml )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ostap_rdict.pcm   DESTINATION ${LIBRARY_OUTPUT_PATH})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ostapDict.rootmap DESTINATION ${LIBRARY_OUTPUT_PATH})
#---Create a shared library with geneated dictionary
add_library(ostap SHARED src/format.cpp
                         src/BSpline.cpp
                         src/Bernstein.cpp
                         src/Binomial.cpp
                         src/Choose.cpp
                         src/Combine.cpp
                         src/Chi2Fit.cpp
                         src/EigenSystem.cpp   
                         src/Error2Exception.cpp   
                         src/Exception.cpp
                         src/Faddeeva.cpp 
                         src/Formula.cpp   
                         src/GSL_sentry.cpp 
                         src/GSL_utils.cpp 
                         src/Hesse.cpp
                         src/HistoDump.cpp
                         src/HistoInterpolation.cpp
                         src/HistoMake.cpp
                         src/HistoProject.cpp
                         src/HistoStat.cpp
                         src/Interpolation.cpp
                         src/Iterator.cpp
                         src/Lomont.cpp
                         src/LorentzVectorWithError.cpp
                         src/Math.cpp
                         src/Models.cpp
                         src/MoreMath.cpp
                         src/Mute.cpp
                         src/NStatEntity.cpp
                         src/PDFs.cpp
                         src/Point3DWithError.cpp
                         src/Polynomials.cpp
                         src/Printable.cpp
                         src/PyIterator.cpp
                         src/PySelector.cpp
                         src/PySelectorWithCuts.cpp
                         src/SFactor.cpp
                         src/StatEntity.cpp
                         src/StatVar.cpp
                         src/StatusCode.cpp
                         src/Tee.cpp
                         src/UStat.cpp
                         src/ValueWithError.cpp
                         src/Vector3DWithError.cpp
                         src/WStatEntity.cpp    
                         src/nSphere.cpp      
                        )    

target_link_libraries(ostap ${ROOT_LIBRARIES} ${GSL_LIBRARIES})

# For clang we need to define __MATH_LONG_DOUBLE_CONSTANTS for
# defining M_El, M_PIl math constants.
#set(COMPILE_DEFINITIONS -D__MATH_LONG_DOUBLE_CONSTANTS)
# add_definitions(-D__MATH_LONG_DOUBLE_CONSTANTS)
set_target_properties(ostap ostapDict PROPERTIES
  SUFFIX ".so"
  COMPILE_DEFINITIONS __MATH_LONG_DOUBLE_CONSTANTS)
# #---Create  a main program using the library
# add_executable(Main MainEvent.cxx)
# target_link_libraries(Main Event)
