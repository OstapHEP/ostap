
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py.in ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py) 

# ==============================================================================
# Macro: find_python_module
# ==============================================================================
macro(find_python_module module)
    string(TOUPPER ${module} module_upper)

    if(NOT PY_${module_upper})
        if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
            set(PY_${module}_FIND_REQUIRED TRUE)
        elseif(ARGC GREATER 2 AND ARGV2 STREQUAL "REQUIRED")
            set(PY_${module}_FIND_REQUIRED TRUE)
        endif()

        if(ARGC GREATER 1 AND ARGV1 STREQUAL "QUIET")
            set(PY_${module}_FIND_QUIETLY TRUE)
        elseif(ARGC GREATER 2 AND ARGV2 STREQUAL "QUIET")
            set(PY_${module}_FIND_QUIETLY TRUE)
        endif()

        execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" -c 
            "import re, ${module}; print(re.compile('/__init__.py.*').sub('',${module}.__file__))"
            RESULT_VARIABLE _${module}_status 
            OUTPUT_VARIABLE _${module}_location
            ERROR_QUIET 
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        if(NOT _${module}_status)
            set(PY_${module_upper} "${_${module}_location}" CACHE PATH 
                "Location of Python module ${module}")
        endif()
    endif()

    find_package_handle_standard_args(PY_${module} DEFAULT_MSG PY_${module_upper})
endmacro(find_python_module)


# ==============================================================================
# Macro: find_python_module_version
# ==============================================================================
macro(find_python_module_version module)
    string(TOUPPER ${module} module_upper)

    if(NOT DEFINED PY_${module}_VERSION OR PY_${module}_VERSION STREQUAL "")
        execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" -c 
            "import ${module}; print(getattr(${module}, '__version__', '?.?.?'))"
            RESULT_VARIABLE _${module}_status 
            OUTPUT_VARIABLE _${module}_version
            ERROR_QUIET 
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        if(NOT _${module}_status)
            set(PY_${module}_VERSION "${_${module}_version}" CACHE STRING 
                "Version of Python module ${module}")
        endif()
    endif()
endmacro(find_python_module_version)


# ==============================================================================
# Macro: ensure_python_module_ex (возвращено в макросы для правильной видимости _FOUND)
# ==============================================================================
macro(ensure_python_module_ex module_name pip_package_spec)
    find_python_module(${module_name} QUIET)
    
    string(TOLOWER "${module_name}" module_lower)
    
    if(NOT PY_${module_name}_FOUND AND NOT PY_${module_lower}_FOUND)
        message(STATUS "----> ${module_name} not found. Installing via pip...")
        
        execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" -m pip install ${pip_package_spec} ${PIP_ARGS}
            OUTPUT_VARIABLE PIP_OUTPUT
            ERROR_VARIABLE  PIP_ERROR
            RESULT_VARIABLE PIP_RESULT
        )
        
        if(NOT PIP_RESULT EQUAL 0)
            message(WARNING "----> pip install ${module_name} failed:\n${PIP_ERROR}")
        endif()
        
        string(TOUPPER "${module_name}" module_upper)
        unset(PY_${module_upper} CACHE)
        unset(PY_${module_name} CACHE)
        
        find_python_module(${module_name} QUIET)
    endif()
endmacro(ensure_python_module_ex)


if(DEFINED ENV{VIRTUAL_ENV} OR DEFINED ENV{CONDA_PREFIX})
  set(PIP_ARGS) 
else()
  set(PIP_ARGS "--user")
endif()


# ========================================================================================
## NumPy & scipy
# ======================================================================================
    
ensure_python_module_ex ( numpy numpy )
if ( NOT PY_numpy_FOUND )
  message ( WARNING "----> warning                      numpy is NOT found!" )
endif()

ensure_python_module_ex ( scipy scipy )
if ( NOT PY_scipy_FOUND )
  message ( WARNING "----> warning                      scipy is NOT found!" )
endif()

ensure_python_module_ex ( numba numba )
if ( NOT PY_numba_FOUND )
  message ( "----> warning (optional)           numba is NOT found!" ) 
endif()

ensure_python_module_ex ( cffi cffi )
if ( NOT PY_cffi_FOUND )
  message ( "----> warning (optional)          cffi is NOT found!" ) 
endif()

# ========================================================================================
## cpuinfo & cpufeature 
# ========================================================================================
ensure_python_module_ex ( cpuinfo py-cpuinfo )
if ( NOT PY_cpuinfo_FOUND ) 
  message ( "----> warning                 cpuinfo is NOT found!" ) 
endif()

find_python_module ( cpufeature QUIET )
if ( NOT PY_cpufeature_FOUND ) 
  message ( "----> warning                 cpufeature is NOT found!" ) 
endif()

# ========================================================================================
## Check for PATHOS & friends
# ========================================================================================
ensure_python_module_ex ( dill "dill[graph,profile]" )
if ( NOT PY_dill_FOUND )
  message ( "----> warning (optional)            dill is NOT found!" ) 
endif()

ensure_python_module_ex ( ppft "ppft[dill]" )
if ( NOT PY_ppft_FOUND )
  message ( "----> warning (optional)            ppft is NOT found!" )
endif()

ensure_python_module_ex ( multiprocess "multiprocess[dill]" )
if ( NOT PY_multiprocess_FOUND )
  message ( "----> warning (optional)    multiprocess is NOT found!" )
endif()

ensure_python_module_ex ( pathos pathos )
if ( NOT PY_pathos_FOUND )
  message ( "----> warning (optional)          pathos is NOT found!" )
endif()

# ========================================================================================
## ipyparallel 
# ========================================================================================
ensure_python_module_ex ( ipyparallel ipyparallel )
if ( NOT PY_ipyparallel_FOUND )
  message ( "----> warning (optional)     ipyparallel is NOT found!" )
endif()

# =======================================================================================
## joblib
# =======================================================================================
ensure_python_module_ex ( joblib joblib )
if ( NOT PY_joblib_FOUND )
  message ( "----> warning (optional)          joblib is NOT found!" ) 
endif() 

# =======================================================================================
## threadpoolctl 
# =======================================================================================
ensure_python_module_ex ( threadpoolctl threadpoolctl )
if ( NOT PY_threadpoolctl_FOUND )
  message ( "----> warning (optional)   threadpoolctl is NOT found!" ) 
endif() 

# =======================================================================================
## cloudpickle
# =======================================================================================
ensure_python_module_ex ( cloudpickle cloudpickle )
if ( NOT PY_cloudpickle_FOUND )
  message ( "----> warning (optional)          cloudpickle is NOT found!" ) 
endif() 

# ========================================================================================
## zstandard
# ========================================================================================
ensure_python_module_ex ( zstandard zstandard )
if ( NOT PY_zstandard_FOUND )
  message ( "----> warning (optional)       zstandard   is NOT found!" )
endif()    

# ========================================================================================
## wcwidth 
# ========================================================================================
ensure_python_module_ex ( wcwidth wcwidth )
if ( NOT PY_wcwidth_FOUND )
  message ( "----> warning (optional)         wcwidth is NOT found!" )
endif()

# ========================================================================================
## terminaltables3 & terminaltables 
# ========================================================================================
if ( PYTHON_VERSION VERSION_GREATER_EQUAL "3.9.0" )
   ensure_python_module_ex ( terminaltables3 terminaltables3 )
   if ( NOT PY_terminaltables3_FOUND )
     message ( "----> warning (optional) terminaltables3 is NOT found!" )    
   endif()                        
else()
   ensure_python_module_ex ( terminaltables terminaltables )
   if ( NOT PY_terminaltables_FOUND )
     message ( "----> warning (optional)  terminaltables is NOT found!" )    
   endif()            
endif()

# ========================================================================================
## tabulate 
# ========================================================================================
ensure_python_module_ex ( tabulate tabulate )
if ( NOT PY_tabulate_FOUND )
  message ( "----> warning (optional)        tabulate is NOT found!" )
endif()

# ========================================================================================
## pathvalidate
# ========================================================================================
ensure_python_module_ex ( pathvalidate pathvalidate )
if ( NOT PY_pathvalidate_FOUND )
  message ( "----> warning (optional)    pathvalidate is NOT found!" )
endif()

# ========================================================================================
## more_itertools 
# ========================================================================================
ensure_python_module_ex ( more_itertools more_itertools )
if ( NOT PY_more_itertools_FOUND )
  message ( "----> warning (optional)  more_itertools is NOT found!" ) 
endif()

# ========================================================================================
## packaging
# ========================================================================================
ensure_python_module_ex ( packaging packaging )
if ( NOT PY_packaging_FOUND ) 
  message ( "----> warning                  packaging is NOT found!" ) 
endif()

# ========================================================================================
## psutil 
# ========================================================================================
ensure_python_module_ex ( psutil psutil )
if ( NOT PY_psutil_FOUND )
  message ( "----> warning (optional)         psutil is NOT found!" ) 
endif()

# ========================================================================================
## BerkeleyDB
# ========================================================================================
find_python_module ( berkeleydb QUIET )
if ( NOT PY_berkeleydb_FOUND )
  if     (  "3.10.0" VERSION_LESS PYTHON_VERSION  )
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install     berkeleydb           ${PIP_ARGS})
  elseif (  "3.9.0"  VERSION_LESS PYTHON_VERSION  )
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install -I "berkeleydb==18.1.14" ${PIP_ARGS})
  elseif (  "3.8.0"  VERSION_LESS PYTHON_VERSION  )
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install -I "berkeleydb==18.1.10" ${PIP_ARGS})
  elseif (  "3.7.0"  VERSION_LESS PYTHON_VERSION  )
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install -I "berkeleydb==18.1.6"  ${PIP_ARGS})
  elseif (  "3.6.0"  VERSION_LESS PYTHON_VERSION  )
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install -I "berkeleydb==18.1.4"  ${PIP_ARGS})
  endif() 
  unset(PY_BERKELEYDB CACHE)
  find_python_module ( berkeleydb QUIET )
endif()
if ( NOT PY_berkeleydb_FOUND )
  message ( "----> warning (optional)      berkeleydb is NOT found!" ) 
endif()

# ========================================================================================
## BDDDB3
# ========================================================================================
if ( PYTHON_VERSION VERSION_LESS "3.10.0" )
   ensure_python_module_ex ( bsddb3 bsddb3 )
   if ( NOT PY_bsddb3_FOUND )
     message ( "----> warning (optional)          bsddb3 is NOT found!" ) 
   endif()
endif()
    
# ========================================================================================
## duckdb & fasteners
# ========================================================================================
if ( PYTHON_VERSION VERSION_GREATER_EQUAL "3.10.0" )
   ensure_python_module_ex ( duckdb "duckdb[all]" )
   if ( NOT PY_duckdb_FOUND )
     message ( "----> warning (optional)          duckdb is NOT found!" )    
   endif ()
   
   ensure_python_module_ex ( fasteners fasteners )
   if ( NOT PY_fasteners_FOUND )
     message ( "----> warning (optional)       fasteners is NOT found!" )    
   endif ()                       
endif () 

# ========================================================================================
## tkrzw 
# ========================================================================================
find_python_module ( tkrzw QUIET )
if ( NOT PY_tkrzw_FOUND )  
   message ( "----> warning (optional)           tkrzw is NOT found!" ) 
endif()

# ========================================================================================
## sortedcontainers
# ========================================================================================
ensure_python_module_ex ( sortedcontainers sortedcontainers )
if ( NOT PY_sortedcontainers_FOUND )
  message ( "----> warning              sortedcontainers is NOT found!" ) 
endif()

# ========================================================================================
## statsmodels 
# ========================================================================================
ensure_python_module_ex ( statsmodels statsmodels )
if ( NOT PY_statsmodels_FOUND )
  message ( "----> warning                statsmodels is NOT found!" ) 
endif()

# ========================================================================================
## pragmastat 
# ========================================================================================
ensure_python_module_ex ( pragmastat pragmastat )
if ( NOT PY_pragmastat_FOUND )
  message ( "----> warning                 pragmastat is NOT found!" ) 
endif()

# ========================================================================================
## Pandas & sklearn
# ========================================================================================
ensure_python_module_ex ( pandas pandas )
if ( NOT PY_pandas_FOUND )
  message ( "----> warning (optional)        pandas is NOT found!" ) 
endif()

ensure_python_module_ex ( sklearn scikit-learn )
if ( NOT PY_sklearn_FOUND )
  message ( "----> warning (optional)       sklearn is NOT found!" ) 
endif()

# ========================================================================================
## sweetviz
# ========================================================================================
ensure_python_module_ex ( sweetviz sweetviz )
if ( NOT PY_sweetviz_FOUND )
  message ( "----> warning                   sweetviz is NOT found!" ) 
endif()

# ========================================================================================
## lightgbm
# ========================================================================================
ensure_python_module_ex ( lightgbm "lightgbm[pandas,scikit-learn]" )
if ( NOT PY_lightgbm_FOUND )
  message ( "----> warning                   lightgbm is NOT found!" )
endif()

# ========================================================================================
## xgboost
# ========================================================================================
ensure_python_module_ex ( xgboost xgboost )
if ( NOT PY_xgboost_FOUND )
   message ( "----> warning                 xgboost is NOT found!" ) 
endif()        

# ========================================================================================
## catboost
# ========================================================================================
if ( OSTAP_HAS_AVX AND OSTAP_HAS_AVX2 )
  ensure_python_module_ex ( catboost catboost )
  if ( NOT PY_catboost_FOUND ) 
    message ( "----> warning                catboost is NOT found!" ) 
  endif()
endif()

# ========================================================================================
## pytorch, tensorflow & keras 
# ========================================================================================
if ( OSTAP_HAS_AVX AND OSTAP_HAS_AVX2 )
  ensure_python_module_ex ( torch torch )
  if ( NOT PY_torch_FOUND ) 
    message ( "----> warning               (py)torch is NOT found!" ) 
  endif()

  ensure_python_module_ex ( tensorflow tensorflow )
  if ( NOT PY_tensorflow_FOUND ) 
    message ( "----> warning              tensorflow is NOT found!" ) 
  endif()

  ensure_python_module_ex ( keras keras )
  if ( NOT PY_keras_FOUND ) 
    message ( "----> warning                   keras is NOT found!" ) 
  endif()  
endif() 

# ========================================================================================
## hep_ml
# ========================================================================================
ensure_python_module_ex ( hep_ml hep_ml )
if ( NOT PY_hep_ml_FOUND )
  message ( "----> warning (optional)          hep_ml is NOT found!" ) 
endif()

# ========================================================================================
## snakeviz
# ========================================================================================
ensure_python_module_ex ( snakeviz snakeviz )
if ( NOT PY_snakeviz_FOUND )
  message ( "----> warning                   snakeviz is NOT found!" ) 
endif()

if ( PY_pathos_FOUND AND PY_dill_FOUND AND PY_multiprocess_FOUND AND PY_ppft_FOUND ) 
   message ( "----> Full pathos (pathos+dill+multiprocess+ppft) is found" )
else()
   message ( "----> warning                     pathos is NOT found!" )  
endif()

# ========================================================================================
## report (Sorted alphabetically)
# ========================================================================================

if ( PY_berkeleydb_FOUND )
  find_python_module_version ( berkeleydb ) 
  message ( "----> berkeleydb       version    : " ${PY_berkeleydb_VERSION} " \t" ${PY_BERKELEYDB} )              
endif()

if ( PY_bsddb3_FOUND )
  find_python_module_version ( bsddb3 ) 
  message ( "----> bsddb3           version    : " ${PY_bsddb3_VERSION} " \t" ${PY_BSDDB3} )                          
endif()

if ( OSTAP_HAS_AVX2 AND PY_catboost_FOUND )
  find_python_module_version ( catboost ) 
  message ( "----> catboost         version    : " ${PY_catboost_VERSION} " \t" ${PY_CATBOOST} ) 
endif()

if ( PY_cffi_FOUND )
  find_python_module_version ( cffi ) 
  message ( "----> cffi             version    : " ${PY_cffi_VERSION} " \t" ${PY_CFFI}  )  
endif()

if ( PY_cloudpickle_FOUND )
  find_python_module_version ( cloudpickle ) 
  message ( "----> cloudpickle      version    : " ${PY_cloudpickle_VERSION} " \t" ${PY_CLOUDPICKLE} )            
endif() 

if ( PY_cpufeature_FOUND )
  find_python_module_version ( cpufeature ) 
  message ( "----> cpufeature       version    : " ${PY_cpufeature_VERSION} " \t" ${PY_CPUFEATURE}  )  
endif()

if ( PY_cpuinfo_FOUND )
  find_python_module_version ( cpuinfo ) 
  message ( "----> cpuinfo          version    : " ${PY_cpuinfo_VERSION} " \t" ${PY_CPUINFO}  )  
endif()

if ( PY_dill_FOUND )
  find_python_module_version ( dill ) 
  message ( "----> dill             version    : " ${PY_dill_VERSION} " \t" ${PY_DILL} )  
endif()

if ( PYTHON_VERSION VERSION_GREATER_EQUAL "3.10.0" )
  if ( PY_duckdb_FOUND )
    find_python_module_version ( duckdb ) 
    message ( "----> duckdb           version    : " ${PY_duckdb_VERSION} " \t" ${PY_DUCKDB} )                                       
  endif()
  if ( PY_fasteners_FOUND )
    find_python_module_version ( fasteners ) 
    message ( "----> fasteners        version    : " ${PY_fasteners_VERSION} " \t" ${PY_FASTENERS} )                                       
  endif()
endif()

if ( PY_hep_ml_FOUND )
  find_python_module_version ( hep_ml ) 
  message ( "----> hep_ml           version    : " ${PY_hep_ml_VERSION} " \t" ${PY_HEP_ML} ) 
endif()

if ( PY_ipyparallel_FOUND )
  find_python_module_version ( ipyparallel ) 
  message ( "----> ipyparallel      version    : " ${PY_ipyparallel_VERSION} " \t" ${PY_IPYPARALLEL} )        
endif()

if ( PY_joblib_FOUND )
  find_python_module_version ( joblib ) 
  message ( "----> joblib           version    : " ${PY_joblib_VERSION} " \t" ${PY_JOBLIB} )            
endif() 

if ( OSTAP_HAS_AVX2 AND PY_keras_FOUND )
  find_python_module_version ( keras ) 
  message ( "----> keras            version    : " ${PY_keras_VERSION} " \t" ${PY_KERAS} ) 
endif()

if ( PY_lightgbm_FOUND )
  find_python_module_version ( lightgbm ) 
  message ( "----> lightgbm         version    : " ${PY_lightgbm_VERSION} " \t" ${PY_LIGHTGBM} )
  if ( PY_lightgbm_VERSION VERSION_GREATER_EQUAL "4.7.0" )
     find_python_module ( narwhals QUIET )
     if ( PY_narwhals_FOUND )
       find_python_module_version ( narwhals ) 
       message ( "----> narwhals         version    : " ${PY_narwhals_VERSION} " \t" ${PY_NARWHALS} )
     endif ()     
  endif () 
endif()

if ( PY_more_itertools_FOUND )
  find_python_module_version ( more_itertools ) 
  message ( "----> more_itertools   version    : " ${PY_more_itertools_VERSION} " \t" ${PY_MORE_ITERTOOLS} )                           
endif()

if ( PY_multiprocess_FOUND )
  find_python_module_version ( multiprocess ) 
  message ( "----> multiprocess     version    : " ${PY_multiprocess_VERSION} " \t" ${PY_MULTIPROCESS} )
endif()

if ( PY_numba_FOUND )
  find_python_module_version ( numba ) 
  message ( "----> numba            version    : " ${PY_numba_VERSION} " \t" ${PY_NUMBA} )  
endif()

if ( PY_numpy_FOUND )
  find_python_module_version ( numpy ) 
  message ( "----> numpy            version    : " ${PY_numpy_VERSION}   " \t" ${PY_NUMPY} )  
endif()

if ( PY_packaging_FOUND )
  find_python_module_version ( packaging ) 
  message ( "----> packaging        version    : " ${PY_packaging_VERSION} " \t" ${PY_PACKAGING} )                           
endif()

if ( PY_pandas_FOUND )
  find_python_module_version ( pandas ) 
  message ( "----> pandas           version    : " ${PY_pandas_VERSION} " \t" ${PY_PANDAS}  )  
endif()

if ( PY_pathos_FOUND )
  find_python_module_version ( pathos ) 
  message ( "----> pathos           version    : " ${PY_pathos_VERSION} " \t" ${PY_PATHOS} )        
endif()

if ( PY_pathvalidate_FOUND )
  find_python_module_version ( pathvalidate )   
  message ( "----> pathvalidate     version    : " ${PY_pathvalidate_VERSION} " \t" ${PY_PATHVALIDATE} )           
endif()

if ( PY_ppft_FOUND )
  find_python_module_version ( ppft ) 
  message ( "----> ppft             version    : " ${PY_ppft_VERSION} " \t" ${PY_PPFT} )    
endif()

if ( PY_pragmastat_FOUND )
  find_python_module_version ( pragmastat ) 
  message ( "----> pragmastat       version    : " ${PY_pragmastat_VERSION} " \t" ${PY_PRAGMASTAT} ) 
endif()

if ( PY_psutil_FOUND )
  find_python_module_version ( psutil ) 
  message ( "----> psutil           version    : " ${PY_psutil_VERSION} " \t" ${PY_PSUTIL} )              
endif()

if ( PY_scipy_FOUND )
  find_python_module_version ( scipy ) 
  message ( "----> scipy            version    : " ${PY_scipy_VERSION}  " \t" ${PY_SCIPY} )  
endif()

if ( PY_sklearn_FOUND )
  find_python_module_version ( sklearn ) 
  message ( "----> sklearn          version    : " ${PY_sklearn_VERSION} " \t" ${PY_SKLEARN}  )  
endif()

if ( PY_snakeviz_FOUND )
  find_python_module_version ( snakeviz ) 
  message ( "----> snakeviz         version    : " ${PY_snakeviz_VERSION} " \t" ${PY_SNAKEVIZ} ) 
endif()

if ( PY_sortedcontainers_FOUND )
  find_python_module_version ( sortedcontainers ) 
  message ( "----> sortedcontainers version    : " ${PY_sortedcontainers_VERSION} " \t" ${PY_SORTEDCONTAINERS} ) 
endif()

if ( PY_statsmodels_FOUND )
  find_python_module_version ( statsmodels ) 
  message ( "----> statsmodels      version    : " ${PY_statsmodels_VERSION} " \t" ${PY_STATSMODELS} ) 
endif()

if ( PY_sweetviz_FOUND )
  find_python_module_version ( sweetviz ) 
  message ( "----> sweetviz         version    : " ${PY_sweetviz_VERSION} " \t" ${PY_SWEETVIZ} ) 
endif()

if ( PY_tabulate_FOUND )
  find_python_module_version ( tabulate ) 
  message ( "----> tabulate         version    : " ${PY_tabulate_VERSION} " \t" ${PY_TABULATE} )               
endif()

if ( PY_terminaltables3_FOUND )
  find_python_module_version ( terminaltables3 ) 
  message ( "----> terminaltables3  version    : " ${PY_terminaltables3_VERSION} " \t" ${PY_TERMINALTABLES3} )        
endif()                        

if ( PY_terminaltables_FOUND )
  find_python_module_version ( terminaltables ) 
  message ( "----> terminaltables   version    : " ${PY_terminaltables_VERSION} " \t" ${PY_TERMINALTABLES} )                
endif()

if ( OSTAP_HAS_AVX2 AND PY_tensorflow_FOUND )
  find_python_module_version ( tensorflow ) 
  message ( "----> tensorflow       version    : " ${PY_tensorflow_VERSION} " \t" ${PY_TENSORFLOW} ) 
endif()

if ( PY_threadpoolctl_FOUND )
  find_python_module_version ( threadpoolctl )  
  message ( "----> threadpoolctl    version    : " ${PY_threadpoolctl_VERSION} " \t" ${PY_THREADPOOLCTL} )          
endif() 

if ( PY_tkrzw_FOUND )
  find_python_module_version ( tkrzw ) 
  message ( "----> tkrzw            version    : " ${PY_tkrzw_VERSION} " \t" ${PY_TKRZW} )                                       
endif()

if ( OSTAP_HAS_AVX2 AND PY_torch_FOUND )
  find_python_module_version ( torch ) 
  message ( "----> torch            version    : " ${PY_torch_VERSION} " \t" ${PY_TORCH} ) 
endif()

if ( PY_wcwidth_FOUND )
  find_python_module_version ( wcwidth )         
  message ( "----> wcwidth          version    : " ${PY_wcwidth_VERSION} " \t" ${PY_WCWIDTH} )        
endif()                        

if ( PY_xgboost_FOUND )
  find_python_module_version ( xgboost ) 
  message ( "----> xgboost          version    : " ${PY_xgboost_VERSION} " \t" ${PY_XGBOOST} ) 
endif()

if ( PY_zstandard_FOUND )
  find_python_module_version ( zstandard ) 
  message ( "----> zstandard        version    : " ${PY_zstandard_VERSION} " \t" ${PY_ZSTANDARD} )        
endif()    


# =========================================================================================
## INSTALLATION
# =========================================================================================

if ( PYTHON_VERSION VERSION_GREATER_EQUAL "3.10.0"  )    
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

set( $ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}:$ENV{PYTHONPATH} ) 

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
   set ( OSTAP_TEST_PYTHON_WARNINGS "-Werror -Wignore::DeprecationWarning" ) 
endif() 

set ( OSTAP_TEST_PYWARNINGS "error,ignore::DeprecationWarning,once::UserWarning,ignore::FutureWarning" ) 

if ( OSTAP_TEST_PYWARNINGS )  
  message ( "----> All tests will be executed with PYTHONWARNINGS=" ${OSTAP_TEST_PYWARNINGS}  )
endif()

set ( TEST_ENV1                                    "PYWARNINGS=${OSTAP_TEST_PYWARNINGS}" )
set ( TEST_ENV2                 "PYTHONBUFFERED=1" "PYWARNINGS=${OSTAP_TEST_PYWARNINGS}" )
set ( TEST_ENV3 "OSTAP_NCPUS=3" "PYTHONBUFFERED=1" "PYWARNINGS=${OSTAP_TEST_PYWARNINGS}" )

# ========================================================================================
## main tests
# ========================================================================================
    
file(GLOB tests "${CMAKE_INSTALL_PREFIX}${OSTAP_PYTHON_SITE_LOCAL}/ostap/*/tests/test_*.py")

option(USE_NOSE "Use nosetests to run main tests"  OFF ) 

if ( USE_NOSE )
  find_program(NOSE_EXECUTABLE NAMES nosetests)
  if ( NOSE_EXECUTABLE )
    message ( "Nosetest will be used for the main tests:" ${NOSE_EXECUTABLE} )
  else ()
    message ( " No NOSETESTS: " ${USE_NOSE} "  and " ${NOSE_EXECUTABLE} ) 
  endif()
endif() 

foreach(filename ${tests})
  get_filename_component(tname ${filename} NAME_WE)
  string(CONCAT tname ostap. ${tname})

  string(FIND   ${tname}  parallel TEST_PARALLEL )
  string(FIND   ${tname}  _pp      TEST_PP       )
  string(FIND   ${tname}  _pathos  TEST_PATHOS   )
  string(FIND   ${tname}  _ipypa   TEST_IPYPA    )

  if ( NOSE_EXECUTABLE AND TEST_PARALLEL LESS 0 )
    add_test (NAME ${tname}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary
      COMMAND ${NOSE_EXECUTABLE} ${filename} --process-timeout=6000 )    
  else   ()
    add_test (NAME ${tname}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary
      COMMAND ${PYTHON_EXECUTABLE} ${filename} )
  endif  ()

  if     ( TEST_IPYPA     GREATER_EQUAL 0  )
    set_tests_properties(${tname} PROPERTIES TIMEOUT 6000 ENVIRONMENT "${TEST_ENV3}" ) 
    message ( "----> Test will be executed with OSTAP_NCPUS=2 PYTHONBUFFERED=1 : " ${tname})  
  elseif ( TEST_PP     GREATER_EQUAL 0  )
    set_tests_properties(${tname} PROPERTIES TIMEOUT 6000 ENVIRONMENT "${TEST_ENV2}" ) 
    message ( "----> Test will be executed with PYTHONBUFFERED=1         : " ${tname})    
  elseif ( TEST_PATHOS GREATER_EQUAL 0 )
    set_tests_properties(${tname} PROPERTIES TIMEOUT 6000 ENVIRONMENT "${TEST_ENV2}" ) 
    message ( "----> Test will be executed with PYTHONBUFFERED=1         : " ${tname})    
  else   () 
    set_tests_properties(${tname} PROPERTIES TIMEOUT 6000 ENVIRONMENT "${TEST_ENV1}") 
  endif() 

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

if ( NOT PY_zstandard_FOUND )
   list(FILTER          allfiles EXCLUDE REGEX ".*\\/zstshelve\\.py$")
endif()

set (PREFIX ${CMAKE_INSTALL_PREFIX}${OSTAP_PYTHON_SITE_LOCAL}/ )

foreach(filename ${allfiles})
   string(REPLACE FAKEfakeFAKEfake ququququ   testname  ${filename} )
   string(REPLACE ${PREFIX} ""  testname ${testname})
   string(REPLACE "\.py" ""     testname ${testname})
   string(REPLACE "/" "_"       testname ${testname})
   string(CONCAT  testname ostap.runit_  ${testname})

   add_test (NAME ${testname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary COMMAND python ${filename} )
   set_tests_properties(${testname} PROPERTIES TIMEOUT 6000 ENVIRONMENT PYTHONWARNINGS=${OSTAP_TEST_PYWARNINGS} )
endforeach() 

add_test (NAME ostap.ostap        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary  COMMAND ostap     -b )
add_test (NAME ostap.dependencies WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary  COMMAND check-dependencies )
add_test (NAME ostap.sync_dirs1   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary  COMMAND sync-dirs -h )
add_test (NAME ostap.sync_dirs2   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary  COMMAND sync-dirs --progress --verbose ${CMAKE_BINARY_DIR}/Testing ${CMAKE_BINARY_DIR}/Testing2 )
add_test (NAME ostap.pplaunch     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Testing/Temporary  COMMAND pplaunch  -h )

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
