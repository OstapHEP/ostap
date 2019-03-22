# root_generate_dictionary(ostap_dict ${CMAKE_CURRENT_SOURCE_DIR}/dict/Dict.h ${CMAKE_CURRENT_SOURCE_DIR}/dict/selections.xml)

find_package(ROOT 6 CONFIG REQUIRED )

if(ROOT_FIND_VERSION VERSION_LESS "6.16.0")
    include(${CMAKE_CURRENT_SOURCE_DIR}/lcg_94_93/CMakeList.txt)
else()
    include(${CMAKE_CURRENT_SOURCE_DIR}/lcg_95/CMakeList.txt)
endif()
