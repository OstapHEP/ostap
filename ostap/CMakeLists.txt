configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py) 

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
  from distutils import sysconfig as sc
  print(sc.get_python_lib(prefix='', plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE
  OUTPUT_STRIP_TRAILING_WHITESPACE)

install(DIRECTORY   ${CMAKE_CURRENT_SOURCE_DIR}
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${PYTHON_SITE}
        FILES_MATCHING PATTERN "*.py" )


## senseless? 
set( $ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}:$ENV{PYTHONPATH} ) 

# first we can indicate the documentation build as an option and set it to ON by default
option(RUN_TESTS "Run tests?" ON)

enable_testing() 
file(GLOB files "${CMAKE_INSTALL_PREFIX}/${PYTHON_SITE}/ostap/*/tests/test_*.py")
find_program(NOSE_EXECUTABLE NAMES nosetests)
foreach(filename ${files})
   get_filename_component(fname ${filename} NAME_WE)
   if( ${NOSE_EXECUTABLE-FOUND} )
     add_test (NAME ${fname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary COMMAND ${NOSE_EXECUTABLE} ${filename} --process-timeout=5000 )
   else  () 
     add_test (NAME ${fname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary COMMAND python             ${filename} )
   endif () 
   set_tests_properties(${fname} PROPERTIES TIMEOUT 5000 ) 
endforeach() 
file(GLOB               allfiles "${CMAKE_INSTALL_PREFIX}/${PYTHON_SITE}/ostap/*/*.py")
list(SORT               allfiles)
list(REMOVE_DUPLICATES  allfiles)
list(FILTER             allfiles EXCLUDE REGEX ".*\\/__init__\\.py$")
list(FILTER             allfiles EXCLUDE REGEX ".*\\/tests/.*\\.py$")

foreach(filename ${allfiles})

   string(REPLACE FAKEfakeFAKEfake ququququ   test_name ${filename} )

   string(REPLACE fitting/param  fitting_param  test_name ${test_name})
   string(REPLACE histos/param   histos_param   test_name ${test_name})
   string(REPLACE math/param     math_param     test_name ${test_name})
   
   string(REPLACE fitting/utils  fitting_utils  test_name ${test_name})
   string(REPLACE logger/utils   logger_utils   test_name ${test_name})
   string(REPLACE utils/utils    utils_utils    test_name ${test_name})
   string(REPLACE core/utils    utils_utils    test_name ${test_name})

   string(REPLACE fitting/basic  fitting_basic  test_name ${test_name})
   string(REPLACE utils/basic    utils_basic    test_name ${test_name})

   string(REPLACE fitting/funcs  fitting_funcs  test_name ${test_name})
   string(REPLACE trees/funcs    trees_funcs    test_name ${test_name})

   string(REPLACE fitting/models fitting_models test_name ${test_name})
   string(REPLACE math/models    math_models    test_name ${test_name})

   get_filename_component(tname ${test_name} NAME_WE)
   ##string(PREPEND tname test_runit_)
   string(CONCAT tname test_runit_ ${tname})

   if( ${NOSE_EXECUTABLE-FOUND} )
     add_test (NAME ${tname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary  COMMAND ${NOSE_EXECUTABLE} ${filename} --process-timeout=5000 )
   else  () 
     add_test (NAME ${tname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary  COMMAND python             ${filename} )
   endif () 
   set_tests_properties(${tname} PROPERTIES TIMEOUT 5000 ) 

endforeach() 


## nosetests -w ../ostap/tools --exe --processes=16  --process-timeout=2000 -I '.*chopping.*'
