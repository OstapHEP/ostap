root_generate_dictionary(ostap_dict ${CMAKE_CURRENT_SOURCE_DIR}/dict/Dict.h LINKDEF LinkDef.h)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ostap_dict_rdict.pcm DESTINATION ${LIBRARY_OUTPUT_PATH})
#---Create a shared library with geneated dictionary
add_library(ostap SHARED math/ValueWithError.cpp math/Lomont.cpp ostap_dict.cxx)
target_link_libraries(ostap ${ROOT_LIBRARIES})

# #---Create  a main program using the library
# add_executable(Main MainEvent.cxx)
# target_link_libraries(Main Event)