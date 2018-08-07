

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py) 
install(DIRECTORY   ${CMAKE_CURRENT_SOURCE_DIR}
        DESTINATION .
        FILES_MATCHING PATTERN "*.py" )

## senseless? 
set( $ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}:$ENV{PYTHONPATH} ) 

# first we can indicate the documentation build as an option and set it to ON by default
option(RUN_TESTS "Run tests?" ON)

enable_testing() 
message ('nose:' ${NOSE_EXECUTABLE} ) 
##file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*/tests/test_*.py")
file(GLOB files "${CMAKE_INSTALL_PREFIX}/ostap/*/tests/test_*.py")
foreach(filename ${files})
   get_filename_component(fname ${filename} NAME_WE)
   message( 'fname' ${fname} )
   find_program(NOSE_EXECUTABLE NAMES nosetests)
   if( ${NOSE_EXECUTABLE-FOUND} )
     add_test (NAME ${fname} COMMAND ${NOSE_EXECUTABLE} ${filename} --process-timeout=3000 )
   else  () 
     add_test (NAME ${fname} COMMAND python             ${filename} )
   endif () 
   set_tests_properties(${fname} PROPERTIES TIMEOUT 5000 ) 
endforeach() 
## nosetests -w ../ostap/tools --exe --processes=16  --process-timeout=2000 -I '.*chopping.*'


