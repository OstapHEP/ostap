
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py) 

macro(find_python_module module)
	string(TOUPPER ${module} module_upper)
	if(NOT PY_${module_upper})
		if    (ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
			set(PY_${module}_FIND_REQUIRED TRUE)
		elseif(ARGC GREATER 2 AND ARGV2 STREQUAL "REQUIRED")
			set(PY_${module}_FIND_REQUIRED TRUE)
		endif()
		if    (ARGC GREATER 1 AND ARGV1 STREQUAL "QUIET")
			set(PY_${module}_FIND_QUIETLY TRUE)
		elseif(ARGC GREATER 2 AND ARGV2 STREQUAL "QUIET")
			set(PY_${module}_FIND_QUIETLY TRUE)
		endif()
		# A module's location is usually a directory, but for binary modules
		# it's a .so file.
		execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c 
			"import re, ${module}; print(re.compile('/__init__.py.*').sub('',${module}.__file__))"
			RESULT_VARIABLE _${module}_status 
			OUTPUT_VARIABLE _${module}_location
			ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
		if(NOT _${module}_status)
			set(PY_${module_upper} ${_${module}_location} CACHE STRING 
				"Location of Python module ${module}")
		endif(NOT _${module}_status)
	endif(NOT PY_${module_upper})
	find_package_handle_standard_args(PY_${module} DEFAULT_MSG PY_${module_upper} )
	## find_package_handle_standard_args(PY_${module}   REQUIRED_VARS PY_${module_upper} )
endmacro(find_python_module)

macro(find_python_module_version module)
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c 
    "import ${module}; print ( getattr ( ${module} , '__version__' , '?.?' ) ) "
    RESULT_VARIABLE _${module}_status 
    OUTPUT_VARIABLE _${module}_version
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT _${module}_status)
        set( PY_${module}_VERSION ${_${module}_version} CACHE STRING 
        "Version of Python module ${module}")
    endif(NOT _${module}_status)
endmacro(find_python_module_version)

if(DEFINED ENV{VIRTUAL_ENV} OR DEFINED ENV{CONDA_PREFIX})
##   set(PIP_ARGS)
   set(PIP_ARGS "--user")
else()
   set(PIP_ARGS "--user")
endif()

# ========================================================================================
## NumPy & scipy
# ======================================================================================
    
find_python_module ( numpy QUIET )
if ( NOT ${PY_numpy_FOUND} )
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install numpy ${PIP_ARGS})
endif()
if ( NOT ${PY_numpy_FOUND} )
  message ( WARNING "----> warning                      numpy is NOT found!" )
endif()
find_python_module ( scipy QUIET )
if ( NOT ${PY_scipy_FOUND} )
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install scipy ${PIP_ARGS})
endif()
find_python_module ( scipy QUIET )
if ( NOT ${PY_scipy_FOUND} )
  message ( WARNINg "----> warning                      scipy is NOT found!" )
endif()

find_python_module ( numba QUIET )
if ( NOT ${PY_numba_FOUND} )
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install numba ${PIP_ARGS})
endif()
find_python_module ( numba QUIET )
if ( NOT ${PY_numba_FOUND} )
  message ( "----> warning (optional)           numba is NOT found!" ) 
endif()

find_python_module ( cffi  QUIET )
if ( NOT ${PY_cffi_FOUND} )
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install cffi ${PIP_ARGS})
endif()
find_python_module ( cffi  QUIET )
if ( NOT ${PY_cffi_FOUND} )
  message ( "----> warning (optional)          cffi is NOT found!" ) 
endif()


# ========================================================================================
## Check for PATHOS & friends
# ========================================================================================
find_python_module ( dill          QUIET )
if ( NOT ${PY_dill_FOUND} )
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install dill[graph]   ${PIP_ARGS})
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install dill[profile] ${PIP_ARGS})
endif()
find_python_module ( dill          QUIET )
if ( NOT ${PY_dill_FOUND} )
  message ( "----> warning (optional)            dill is NOT found!" ) 
endif()
##
find_python_module ( ppft          QUIET )
if ( NOT ${PY_ppft_FOUND} )
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install ppft[dill] ${PIP_ARGS})
endif()
find_python_module ( ppft          QUIET )
if ( NOT ${PY_ppft_FOUND} )
  message ( "----> warning (optional)            ppft is NOT found!" )
endif()
##
find_python_module ( multiprocess  QUIET )
if ( NOT ${PY_multiprocess_FOUND} )
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install multiprocess[dill] ${PIP_ARGS})
endif()
if ( NOT ${PY_multiprocess_FOUND} )
  message ( "----> warning (optional)    multiprocess is NOT found!" )
endif()
##
find_python_module ( pathos        QUIET )
if ( NOT ${PY_pathos_FOUND} )
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install pathos ${PIP_ARGS})
endif()
if ( NOT ${PY_pathos_FOUND} )
  message ( "----> warning (optional)          pathos is NOT found!" )
endif()
##

# ========================================================================================
## ipyparallel 
# ========================================================================================
find_python_module ( ipyparallel          QUIET )
if ( NOT ${PY_ipyparallel_FOUND} )
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install ipyparallel ${PIP_ARGS})
endif()
find_python_module ( ipyparallel          QUIET )
if ( NOT ${PY_ipyparallel_FOUND} )
  message ( "----> warning (optional)     ipyparallel is NOT found!" )
endif()

# =======================================================================================
## joblib: run some numpy-related code in parallel
# =======================================================================================
find_python_module ( joblib          QUIET )
if ( NOT ${PY_joblib_FOUND} )
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install joblib ${PIP_ARGS})
endif()
find_python_module ( joblib          QUIET )
if ( NOT ${PY_joblib_FOUND} )
  message ( "----> warning (optional)          joblib is NOT found!" ) 
endif() 

# ========================================================================================
## zstandard
# ========================================================================================
find_python_module ( zstandard QUIET )
if ( NOT ${PY_zstandard_FOUND} )
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install zstandard ${PIP_ARGS})    
endif()
find_python_module ( zstandard QUIET )
if ( NOT ${PY_zstandard_FOUND} )
  message ( "----> warning (optional)       zstandard   is NOT found!" )
endif()    
# ========================================================================================
    
# ========================================================================================
## terminaltables3 & terminaltables 
# ========================================================================================
if ( ${PYTHON_VERSION} VERSION_GREATER_EQUAL "3.9" )
   find_python_module ( terminaltables3 QUIET )
   if ( NOT ${PY_terminaltables3_FOUND} )
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install terminaltables3 ${PIP_ARGS})
   endif ()
   find_python_module ( terminaltables3 QUIET )
   if ( NOT ${PY_terminaltables3_FOUND} )
     message ( "----> warning (optional) terminaltables3 is NOT found!" )    
   endif ()                        
else  ()
   find_python_module ( terminaltables QUIET )
   if ( NOT ${PY_terminaltables_FOUND} )
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install terminaltables ${PIP_ARGS})
   endif ()
   find_python_module ( terminaltables QUIET )
   if ( NOT ${PY_terminaltables_FOUND} )
     message ( "----> warning (optional)  terminaltables is NOT found!" )    
   endif ()            
endif ()

# ========================================================================================
## tabulate 
# ========================================================================================
##  
find_python_module ( tabulate QUIET )
if ( NOT ${PY_tabulate_FOUND} )
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install tabulate ${PIP_ARGS})
endif()
find_python_module ( tabulate QUIET )
if ( NOT ${PY_tabulate_FOUND} )
  message ( "----> warning (optional)        tabulate is NOT found!" )
endif()

# ========================================================================================
## more_itertools 
# ========================================================================================
##   
find_python_module ( more_itertools QUIET )
if ( NOT ${PY_more_itertools_FOUND} )
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install more_itertools ${PIP_ARGS})
endif()
find_python_module ( more_itertools QUIET )
if ( NOT ${PY_more_itertools_FOUND} )
  message ( "----> warning (optional)  more_itertools is NOT found!" ) 
endif()
    
# ========================================================================================
## psutil 
# ========================================================================================
## 
find_python_module ( psutil QUIET )
if ( NOT ${PY_psutil_FOUND} )
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install psutil ${PIP_ARGS})
endif()
find_python_module ( psutil QUIET )
if ( NOT ${PY_psutil_FOUND} )
  message ( "----> warning (optional)         psutils is NOT found!" ) 
endif()

# ========================================================================================
## BerkeleyDB
# ========================================================================================
##
find_python_module ( berkeleydb QUIET )
if ( NOT ${PY_berkeleydb_FOUND} )
   if   ( ${PYTHON_VERSION} VERSION_LESS "3.7" )
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install "berkeleydb==18.1.4"  ${PIP_ARGS})
   elif ( ${PYTHON_VERSION} VERSION_LESS "3.8" )
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install "berkeleydb==18.1.6"  ${PIP_ARGS})
   elif ( ${PYTHON_VERSION} VERSION_LESS "3.9" )
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install "berkeleydb==18.1.10" ${PIP_ARGS})
   elif ( ${PYTHON_VERSION} VERSION_LESS "3.10" )
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install "berkeleydb==18.1.14" ${PIP_ARGS})
   else ()   
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install  berkeleydb           ${PIP_ARGS})
   endif()
endif()
find_python_module ( berkeleydb QUIET )
if ( NOT ${PY_berkeleydb_FOUND} )
  message ( "----> warning (optional)      berkeleydb is NOT found!" ) 
endif()

# ========================================================================================
## BDDDB3
# ========================================================================================
if ( ${PYTHON_VERSION} VERSION_LESS       "3.10" )
   find_python_module ( bsddb3 QUIET )
   if ( NOT ${PY_bsddb3_FOUND} )
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install bsddb3 ${PIP_ARGS})
   endif()
   find_python_module ( bsddb3 QUIET )
   if ( NOT ${PY_bsddb3_FOUND} )
     message ( "----> warning (optional)          bsddb3 is NOT found!" ) 
   endif()
endif()
    
# ========================================================================================
## lmdb : check&install only if mdb_stat exists ...
# ========================================================================================  
find_program(MDB_STAT_EXECUTABLE NAMES mdb_stat)
if ( ${MDB_STAT_EXECUTABLE} )
   find_python_module ( lmdb        QUIET )
   if ( NOT ${PY_lmdb_FOUND} )
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install lmdb ${PIP_ARGS})
   endif()    
   find_python_module ( lmdb        QUIET )
   if ( NOT ${PY_lmdb_FOUND} )
      message ( "----> warning (optional)            lmdb is NOT found!" )
   endif() 
endif()

# ========================================================================================
## tkrzw 
# ========================================================================================
##
find_python_module ( tkrzw QUIET )
if ( NOT ${PY_tkrzw_FOUND} )
   message ( "----> warning (optional)           tkrzw is NOT found!" ) 
endif()

# ========================================================================================
## sortedcontainers
# ========================================================================================
## 
find_python_module ( sortedcontainers QUIET )
if ( NOT ${PY_sortedcontainers_FOUND} )
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install sortedcontainers ${PIP_ARGS})
endif()
find_python_module ( sortedcontainers QUIET )
if ( NOT ${PY_sortedcontainers_FOUND} )
  message ( "----> warning              sortedcontainer is NOT found!" ) 
endif()

# ========================================================================================
## pragmastat 
# ========================================================================================
## 
find_python_module ( pragmastat QUIET )
if ( NOT ${PY_pragmastat_FOUND} )
   execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install pragmastat ${PIP_ARGS})
endif()
find_python_module ( pragmastat  QUIET )
if ( NOT ${PY_pragmastat_FOUND} )
  message ( "----> warning                 pragmastat is NOT found!" ) 
endif()

# ========================================================================================
## hep_ml
# ========================================================================================
if ( ${PYTHON_VERSION} VERSION_LESS       "3.12" )
   find_python_module ( hep_ml QUIET )
   if ( NOT ${PY_hep_ml_FOUND} )
      execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install hep_ml ${PIP_ARGS})
   endif()
   find_python_module ( hep_ml QUIET )
   if ( NOT ${PY_hep_ml_FOUND} )
     message ( "----> warning (optional)          hep_ml is NOT found!" ) 
   endif()
endif()

if ( ${PY_pathos_FOUND} AND ${PY_dill_FOUND} AND ${PY_multiprocess_FOUND} AND ${PY_ppft_FOUND} ) 
   message ( "----> Full pathos (pathos+dill+multiprocess+ppft) is found" )
else()
   message ( "----> warning                     pathos is NOT found!" )  
endif()

if ( ${PY_numpy_FOUND} )
  find_python_module_version ( numpy ) 
  message ( "----> numpy            version    : " ${PY_numpy_VERSION} )  
endif()
##
if ( ${PY_scipy_FOUND} )
  find_python_module_version ( scipy ) 
  message ( "----> scipy            version    : " ${PY_scipy_VERSION} )  
endif()

if ( ${PY_numba_FOUND} )
  find_python_module_version ( numba ) 
  message ( "----> numba            version    : " ${PY_numba_VERSION} )  
endif()

if ( ${PY_cffi_FOUND} )
  find_python_module_version ( cffi ) 
  message ( "----> cffi             version    : " ${PY_cffi_VERSION} )  
endif()

##
if ( ${PY_dill_FOUND} )
  find_python_module_version ( dill ) 
  message ( "----> dill             version    : " ${PY_dill_VERSION} )  
endif()
##
if ( ${PY_ppft_FOUND} )
  find_python_module_version ( ppft ) 
  message ( "----> ppft             version    : " ${PY_ppft_VERSION} )    
endif()
##
if ( ${PY_multiprocess_FOUND} )
  find_python_module_version ( multiprocess ) 
  message ( "----> multiprocess     version    : " ${PY_multiprocess_VERSION} )      
endif()
##
if ( ${PY_pathos_FOUND} )
  find_python_module_version ( pathos ) 
  message ( "----> pathos           version    : " ${PY_pathos_VERSION} )        
endif()
##
if ( ${PY_ipyparallel_FOUND} )
  find_python_module_version ( ipyparallel ) 
  message ( "----> ipyparallel      version    : " ${PY_ipyparallel_VERSION} )        
endif()
##
if ( ${PY_joblib_FOUND} )
  find_python_module_version ( joblib ) 
  message ( "----> joblib           version    : " ${PY_joblib_VERSION} )             
endif() 
##
if ( ${PY_zstandard_FOUND} )
  find_python_module_version ( zstandard ) 
  message ( "----> zstandard        version    : " ${PY_zstandard_VERSION} )        
endif()    
##
if ( ${PY_terminaltables3_FOUND} )
  find_python_module_version ( terminaltables3 ) 
  message ( "----> terminaltables3  version    : " ${PY_terminaltables3_VERSION} )        
endif ()                        
##
if ( ${PY_terminaltables_FOUND} )
  find_python_module_version ( terminaltables ) 
  message ( "----> terminaltables   version    : " ${PY_terminaltables_VERSION} )        
endif ()
##
if ( ${PY_tabulate_FOUND} )
  find_python_module_version ( tabulate ) 
  message ( "----> tabulate         version    : " ${PY_tabulate_VERSION} )        
endif()
##
if ( ${PY_more_itertools_FOUND} )
  find_python_module_version ( more_itertools ) 
  message ( "----> more_itertools   version    : " ${PY_more_itertools_VERSION} )             
endif()
## 
if ( ${PY_psutil_FOUND} )
  find_python_module_version ( psutil ) 
  message ( "----> psutil           version    : " ${PY_psutil_VERSION} )             
endif()
##
if (${PY_berkeleydb_FOUND} )
  find_python_module_version ( berkeleydb ) 
  message ( "----> berkeleydb       version    : " ${PY_berkeleydb_VERSION} )             
endif()
##
if ( ${PY_bsddb3_FOUND} )
  find_python_module_version ( bdsdb3 ) 
  message ( "----> bsddb3           version    : " ${PY_bsddb3_VERSION} )             
endif()
##
if ( ${PY_lmdb_FOUND} )
  find_python_module_version ( lmdb ) 
  message ( "----> lmdb             version    : " ${PY_lmdb_VERSION} )             
endif() 
##
if ( ${PY_tkrzw_FOUND} )
  find_python_module_version ( tkrzw ) 
  message ( "----> tkrzw            version    : " ${PY_tkrzw_VERSION} )             
endif()
##
if ( ${PY_sortedcontainers_FOUND} )
  find_python_module_version ( sortedcontainers ) 
  message ( "----> sortedcontainers version    : " ${PY_sortedcontainers_VERSION} )             
  endif()
##
if ( ${PY_pragmastat_FOUND} )
  find_python_module_version ( pragmastat ) 
  message ( "----> pragmastat       version    : " ${PY_pragmastat_VERSION} )             
endif()
##
if ( ${PY_hep_ml_FOUND} )
  find_python_module_version ( hep_ml ) 
  message ( "----> hep_ml           version    : " ${PY_hep_ml_VERSION} )             
endif()

# =========================================================================================
## INSTALLATION
# =========================================================================================

if ( ${PYTHON_VERSION} VERSION_GREATER_EQUAL "3.10"  )    
   execute_process(
      COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
      import sysconfig as sc
      print(sc.get_path( 'purelib' , vars = { 'base' : '' , 'userbase' : '' , 'platbase' : '' } ))"
      OUTPUT_VARIABLE OSTAP_PYTHON_SITE_LOCAL
      OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
   execute_process(
      COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
      from distutils import sysconfig as sc
      print(sc.get_python_lib(prefix='', plat_specific=True))"
      OUTPUT_VARIABLE OSTAP_PYTHON_SITE_LOCAL
      OUTPUT_STRIP_TRAILING_WHITESPACE)
endif() 


set(OSTAP_PYTHON_SITE_DIR ${OSTAP_PYTHON_SITE_LOCAL}  PARENT_SCOPE)

## senseless? 
set( $ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}:$ENV{PYTHONPATH} ) 

# first we can indicate the documentation build as an option and set it to ON by default
option(RUN_TESTS "Run tests?" ON)

enable_testing() 

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -Werror -c "import math"
  OUTPUT_VARIABLE OSTAP_PYTHON_CHECK_WARNINGS
  RESULT_VARIABLE OSTAP_PYTHON_CHECK_WARNINGS_RESULT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)

if ( OSTAP_PYTHON_CHECK_WARNINGS_RESULT ) 
   set ( OSTAP_TEST_PYTHON_WARNINGS "" ) 
else  () 
   set ( OSTAP_TEST_PYTHON_WARNINGS "" ) ## -Werror" ) 
endif () 

if ( OSTAP_TEST_PYTHON_WARNINGS ) 
    MESSAGE ( "All tests will be executed with " ${OSTAP_TEST_PYTHON_WARNINGS} " flag"  ) 
endif() 

# ========================================================================================
## main tests
# ========================================================================================
    
file(GLOB tests "${CMAKE_INSTALL_PREFIX}${OSTAP_PYTHON_SITE_LOCAL}/ostap/*/tests/test_*.py")

## option(USE_NOSE "Use notetests to run main tests"  ON ) 

find_program(NOSE_EXECUTABLE NAMES nosetests)
if ( NOSE_EXECUTABLE )
  message ( "Nosetest will be used for the main tests:" ${NOSE_EXECUTABLE} )
else ()
  message ( " No NOSETESTS: " ${USE_NOSE} "  and " ${NOSE_EXECUTABLE} ) 
endif()
    
foreach(filename ${tests})
  get_filename_component(tname ${filename} NAME_WE)
  string(CONCAT tname ostap. ${tname})

  string(FIND   ${tname}  parallel HAS_PARALLEL)

  if ( NOSE_EXECUTABLE AND HAS_PARALLEL LESS 0 )
    add_test (NAME ${tname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary COMMAND ${NOSE_EXECUTABLE} ${filename} --process-timeout=10000 )
    set_tests_properties(${tname} PROPERTIES TIMEOUT 10000 ENVIRONMENT OSTAP_CMAKE_TEST=1 ) 
  else   ()
    add_test (NAME ${tname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary COMMAND python ${OSTAP_TEST_PYTHON_WARNINGS} ${filename} )
    set_tests_properties(${tname} PROPERTIES TIMEOUT 10000 ENVIRONMENT OSTAP_CMAKE_TEST=1 ) 
  endif  ()
  
endforeach()


# ========================================================================================
## additional tests 
# ========================================================================================
file(GLOB               allfiles "${CMAKE_INSTALL_PREFIX}${OSTAP_PYTHON_SITE_LOCAL}/ostap/*/*.py")
list(SORT               allfiles)
list(REMOVE_DUPLICATES  allfiles)
list(FILTER             allfiles EXCLUDE REGEX ".*\\/__init__\\.py$")
list(FILTER             allfiles EXCLUDE REGEX ".*\\/tests/.*\\.py$")
list(FILTER             allfiles EXCLUDE REGEX ".*\\/sp_.*\\.py$")
list(FILTER             allfiles EXCLUDE REGEX ".*\\/.*#.*\\.py$")
##
if ( NOT ${PY_zstandard_FOUND} )
   list(FILTER          allfiles EXCLUDE REGEX ".*\\/zstshelve\\.py$")
endif()

set (PREFIX ${CMAKE_INSTALL_PREFIX}${OSTAP_PYTHON_SITE_LOCAL}/ )

foreach(filename ${allfiles})

   string(REPLACE FAKEfakeFAKEfake ququququ   testname  ${filename} )
   string(REPLACE ${PREFIX} ""  testname ${testname})
   string(REPLACE "\.py" ""     testname ${testname})
   string(REPLACE "/" "_"       testname ${testname})
   string(CONCAT  testname ostap.runit_  ${testname})

   add_test (NAME ${testname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary COMMAND python ${OSTAP_TEST_PYTHON_WARNINGS} ${filename} )

   set_tests_properties(${testname} PROPERTIES TIMEOUT 5000 ENVIRONMENT OSTAP_CMAKE_TEST=1 ) 

endforeach() 

add_test (NAME ostap.ostap        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary  COMMAND ostap     -b )
add_test (NAME ostap.dependencies WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary  COMMAND check-dependencies )
add_test (NAME ostap.sync_dirs    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary  COMMAND sync-dirs -h )
 
install ( DIRECTORY       ${CMAKE_CURRENT_SOURCE_DIR}
          DESTINATION     ${CMAKE_INSTALL_PREFIX}${OSTAP_PYTHON_SITE_LOCAL}
          FILES_MATCHING  PATTERN "*.py" 
          PATTERN  "*#*"  EXCLUDE 
 	  PERMISSIONS     OWNER_READ GROUP_READ  WORLD_READ )

install ( DIRECTORY       ${CMAKE_CURRENT_SOURCE_DIR}  
          DESTINATION     . 
          COMPONENT       python 
          FILES_MATCHING  PATTERN "*.py"          
          PATTERN  "*.#*" EXCLUDE      
	  PERMISSIONS     OWNER_READ GROUP_READ  WORLD_READ )

## install ( DIRECTORY       ${CMAKE_CURRENT_SOURCE_DIR}
##          DESTINATION     .${OSTAP_PYTHON_SITE_LOCAL}
##          COMPONENT       python 
##          FILES_MATCHING  PATTERN "*.py" 
##          PATTERN  "*#*"  EXCLUDE 
##          PERMISSIONS     OWNER_READ GROUP_READ  WORLD_READ )
	
## nosetests -w ../ostap/tools --exe --processes=16  --process-timeout=2000 -I '.*chopping.*'
