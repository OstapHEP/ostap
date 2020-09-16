
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py) 

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
  from distutils import sysconfig as sc
  print(sc.get_python_lib(prefix='', plat_specific=True))"
  OUTPUT_VARIABLE OSTAP_PYTHON_SITE_LOCAL
  OUTPUT_STRIP_TRAILING_WHITESPACE)

install(DIRECTORY   ${CMAKE_CURRENT_SOURCE_DIR}
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${OSTAP_PYTHON_SITE_LOCAL}
        FILES_MATCHING PATTERN "*.py" )

set(OSTAP_PYTHON_SITE_DIR ${OSTAP_PYTHON_SITE_LOCAL}  PARENT_SCOPE)


## senseless? 
set( $ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}:$ENV{PYTHONPATH} ) 

# first we can indicate the documentation build as an option and set it to ON by default
option(RUN_TESTS "Run tests?" ON)

enable_testing() 
file(GLOB files "${CMAKE_INSTALL_PREFIX}/${OSTAP_PYTHON_SITE_LOCAL}/ostap/*/tests/test_*.py")
find_program(NOSE_EXECUTABLE NAMES nosetests)
foreach(filename ${files})
   get_filename_component(tname ${filename} NAME_WE)
   string(CONCAT tname ostap. ${tname})
   if( ${NOSE_EXECUTABLE} )
     add_test (NAME ${tname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary COMMAND ${NOSE_EXECUTABLE} ${filename} --process-timeout=5000 )
   else  () 
     add_test (NAME ${tname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary COMMAND python ${filename} )
   endif () 
   set_tests_properties(${tname} PROPERTIES TIMEOUT 5000 ) 
endforeach() 
file(GLOB               allfiles "${CMAKE_INSTALL_PREFIX}/${OSTAP_PYTHON_SITE_LOCAL}/ostap/*/*.py")
list(SORT               allfiles)
list(REMOVE_DUPLICATES  allfiles)
list(FILTER             allfiles EXCLUDE REGEX ".*\\/__init__\\.py$")
list(FILTER             allfiles EXCLUDE REGEX ".*\\/tests/.*\\.py$")
## list(FILTER             allfiles EXCLUDE REGEX ".*\\/lzshelve\\.py$")
list(FILTER             allfiles EXCLUDE REGEX ".*\\/sp_.*\\.py$")
list(FILTER             allfiles EXCLUDE REGEX ".*\\/.*#.*\\.py$")

if (Python_VERSION VERSION_LESS "3.3")
  list(FILTER             allfiles EXCLUDE REGEX ".*\\/lzshelve\\.py$")
endif() 

set (PREFIX ${CMAKE_INSTALL_PREFIX}/${OSTAP_PYTHON_SITE_LOCAL}/ )

foreach(filename ${allfiles})

    string(REPLACE FAKEfakeFAKEfake ququququ   testname  ${filename} )
    string(REPLACE ${PREFIX} ""  testname ${testname})
    string(REPLACE "\.py" ""     testname ${testname})
    string(REPLACE "/" "_"       testname ${testname})
    string(CONCAT  testname ostap.runit_  ${testname})

    add_test (NAME ${testname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary  COMMAND python ${filename} )

    set_tests_properties(${testname} PROPERTIES TIMEOUT 5000 ) 

endforeach() 

install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}  DESTINATION . )

## nosetests -w ../ostap/tools --exe --processes=16  --process-timeout=2000 -I '.*chopping.*'
