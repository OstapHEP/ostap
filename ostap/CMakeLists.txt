configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py) 
install(DIRECTORY   ${CMAKE_CURRENT_SOURCE_DIR}
        DESTINATION .
        FILES_MATCHING PATTERN "*.py" )

find_program(ROOT_CONFIG_EXECUTABLE NAMES root-config )
execute_process( COMMAND "${ROOT_CONFIG_EXECUTABLE}" --python-version
                 OUTPUT_VARIABLE PYVERSION_ROOT
                 OUTPUT_STRIP_TRAILING_WHITESPACE )

## senseless? 
set( $ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}:$ENV{PYTHONPATH} ) 

# first we can indicate the documentation build as an option and set it to ON by default
option(RUN_TESTS "Run tests?" ON)

enable_testing() 
file(GLOB files "${CMAKE_INSTALL_PREFIX}/ostap/*/tests/test_*.py")
foreach(filename ${files})
   get_filename_component(fname ${filename} NAME_WE)
   find_program(NOSE_EXECUTABLE NAMES nosetests)
   if( ${NOSE_EXECUTABLE-FOUND} )
     add_test (NAME ${fname} COMMAND ${NOSE_EXECUTABLE} ${filename} --process-timeout=5000 )
   else  () 
     if(PYVERSION_ROOT  VERSION_LESS "3.6")
             add_test (NAME ${fname} COMMAND python             ${filename} )
     else()
        add_test (NAME ${fname} COMMAND python3             ${filename} )
     endif()
   endif () 
   set_tests_properties(${fname} PROPERTIES TIMEOUT 5000 ) 
endforeach() 
## nosetests -w ../ostap/tools --exe --processes=16  --process-timeout=2000 -I '.*chopping.*'
