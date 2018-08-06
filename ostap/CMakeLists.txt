

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py) 
install(DIRECTORY   ${CMAKE_CURRENT_SOURCE_DIR}
        DESTINATION .
        FILES_MATCHING PATTERN "*.py" )

## senseless? 
set( $ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}:$ENV{PYTHONPATH} ) 

