# root_generate_dictionary(ostap_dict ${CMAKE_CURRENT_SOURCE_DIR}/dict/Dict.h ${CMAKE_CURRENT_SOURCE_DIR}/dict/selections.xml)
find_package(ROOT  REQUIRED)
find_package(Boost REQUIRED)
find_package(GSL   REQUIRED)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS})
REFLEX_GENERATE_DICTIONARY(ostap_dict src/dict/Dict.h SELECTION src/dict/selections.xml)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ostap_dict_rdict.pcm   DESTINATION ${LIBRARY_OUTPUT_PATH})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ostap_dictDict.rootmap DESTINATION ${LIBRARY_OUTPUT_PATH})
#---Create a shared library with geneated dictionary
FILE(GLOB_RECURSE sourceVar src/*.cpp)
##add_library(ostap SHARED ${sourceVar} ostap_dict.cxx)
add_library(ostap SHARED 
                  src/ValueWithError.cpp 
                  src/StatusCode.cpp 
                  ostap_dict.cxx)
target_link_libraries(ostap ${ROOT_LIBRARIES} ${GSL_LIBRARIES})

# #---Create  a main program using the library
# add_executable(Main MainEvent.cxx)
# target_link_libraries(Main Event)
