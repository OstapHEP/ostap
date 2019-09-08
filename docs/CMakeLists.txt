# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation?" ON)

find_package(Doxygen REQUIRED dot) 
## find_package(Sphinx  REQUIRED)

if (BUILD_DOC AND DOXYGEN_FOUND)

   find_package (ROOT)
   string(APPEND ROOT_DOXYTAG ${ROOT_VERSION_MAJOR} ${ROOT_VERSION_MINOR})

   find_file (ROOT_TAG_FILE ROOT.tag ${CMAKE_CURRENT_BINARY_DIR} )      
   if( NOT EXISTS ${ROOT_TAG_FILE} )  
     find_program (WGET_EXECUTABLE wget)
     set( ROOT_TAG_URL  http://root.cern.ch/doc/v${ROOT_DOXYTAG}/ROOT.tag ) 
     add_custom_target ( root_doxy_tag
                         DEPENDS ROOT.tag 
                         COMMAND ${WGET_EXECUTABLE} -q ${ROOT_TAG_URL} 
                         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                         COMMENT "(w)get ROOT.tag file from ${ROOT_TAG_URL}" )
     find_file (ROOT_TAG_FILE ROOT.tag ${CMAKE_CURRENT_BINARY_DIR} )
    else()
     add_custom_target ( root_doxy_tag COMMENT ${ROOT_TAG_FILE} ) 
    endif() 

  set(DOXYGEN_OUTPUT_DIR  ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
  set(DOXYGEN_INDEX_FILE  ${DOXYGEN_OUTPUT_DIR}/html/index.html)
  set(DOXYGEN_TAGFILE     ${DOXYGEN_OUTPUT_DIR}/ostap.tag)
  set(DOXYGEN_LATEX       ${DOXIGEN_OUTPUT_DIR}/latex/refman.tex)

  file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us

  set(DOXYFILE_IN  ${CMAKE_CURRENT_SOURCE_DIR}/OstapDoxy.cfg.in)
  set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/OstapDoxy.cfg)
 
  configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)              

  add_custom_target(DoxyGen
        DEPENDS root_doxy_tag ${DOXYFILE_IN} ostap 
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT} 2> grep -e 'Fontconfig error: Cannot load default config file' -v 
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen (could be time consuming O(15-20'))"
        MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
        VERBATIM)

  add_custom_target(PdfDoc 
        DEPENDS ${DOXIGEN_LATEX}
        COMMAND make pdf 
        WORKING_DIRECTORY ${DOXIGEN_OUTPUT_DIR}/latex
        COMMENT "Building PDF-documentation"
        VERBATIM)

## include(GNUInstallDirs)

if(EXISTS ${DOXYGEN_OUTPUT_DIR}/html)
 install(DIRECTORY ${DOXYGEN_OUTPUT_DIR}/html             DESTINATION docs )
endif()
if(EXISTS ${DOXYGEN_OUTPUT_DIR}/man)
 install(DIRECTORY ${DOXYGEN_OUTPUT_DIR}/man              DESTINATION docs )
endif()
if(EXISTS ${DOXYGEN_OUTPUT_DIR}/latex)
 install(DIRECTORY ${DOXYGEN_OUTPUT_DIR}/latex            DESTINATION docs )
endif()
if(EXISTS ${DOXYGEN_OUTPUT_DIR}/latex/refman.pdf)
 install(DIRECTORY ${DOXYGEN_OUTPUT_DIR}/latex/refman.pdf DESTINATION docs )
endif()
if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/ROOT.tag)
 install(FILES     ${CMAKE_CURRENT_BINARY_DIR}/ROOT.tag   DESTINATION docs )
endif()
if(EXISTS ${DOXYGEN_OUTPUT_DIR}/ostap.tag)
 install(FILES     ${DOXYGEN_OUTPUT_DIR}/ostap.tag        DESTINATION docs )
endif()

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

find_program(SPHINX_EXECUTABLE
             NAMES sphinx-build
             DOC "Path to sphinx-build executable")
if (SPHINX_EXECUTABLE)

 set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
 set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
 set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

 # Only regenerate Sphinx when:
 # - Doxygen has rerun
 # - Our doc files have been updated
 # - The Sphinx config has been updated
add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                   COMMAND 
                   ${SPHINX_EXECUTABLE} -b html
                   # Tell Breathe where to find the Doxygen output
                   -Dbreathe_projects.CatCutifier=${DOXYGEN_OUTPUT_DIR}/xml
                   ${SPHINX_SOURCE} ${SPHINX_BUILD}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   DEPENDS
                   # Other docs files you want to track should go here (or in some variable)
                   ${CMAKE_CURRENT_SOURCE_DIR}/index.rst
                   ${DOXYGEN_INDEX_FILE}
                   MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
                   COMMENT "Generating documentation with Sphinx")
 
# Nice named target so we can run the job easily
add_custom_target(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})


if(EXISTS (${SPHINX_BUILD} ) )
install(DIRECTORY ${SPHINX_BUILD} DESTINATION docs)
endif() 
 
endif(SPHINX_EXECUTABLE) 
endif(BUILD_DOC AND DOXYGEN_FOUND)



