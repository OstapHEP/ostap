# root_generate_dictionary(ostap_dict ${CMAKE_CURRENT_SOURCE_DIR}/dict/Dict.h ${CMAKE_CURRENT_SOURCE_DIR}/dict/selections.xml)
REFLEX_GENERATE_DICTIONARY(ostap_dict dict/Dict.h SELECTION dict/selections.xml)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ostap_dict_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_PATH})
#---Create a shared library with geneated dictionary
add_library(ostap SHARED math/ValueWithError.cpp math/Lomont.cpp math/Math.cpp ostap_dict.cxx)
target_link_libraries(ostap ${ROOT_LIBRARIES})

# #---Create  a main program using the library
# add_executable(Main MainEvent.cxx)
# target_link_libraries(Main Event)