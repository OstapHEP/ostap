#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file parallel.py 
#  Useful utilities for multiprocessing and parallel processing for Ostap
#  Actualy it is just a little bit upgraded version of original
#  GaudiMP.Parallel module developed by Pere MATO Pere.Mato@cern.ch
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-02-23
# =============================================================================
""" Useful utilities for multiprocessing and parallel processing for Ostap
Actualy it is just a little bit upgraded version of original
GaudiMP.Parallel module developed by Pere MATO Pere.Mato@cern.ch
"""
# =============================================================================
__version__ = '$Revision$'
__author__  = 'Vanya BELYAEV Ivan.Belyaev@itep.ru'
__date__    = '2016-02-23'
__all__     = (
    'Task'          , ## the base class for task
    'WorkManager'   , ## task manager
    'GenericTask'   , ## very generic "template"  tasl
    'Checker'       , ## check of the object can be pickled/unpickled 
    )
# =============================================================================
from   ostap.parallel.task import Task, GenericTask
from   ostap.utils.env     import has_env, get_env, OSTAP_PARALLEL  
import sys, os, warnings 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.parallel.parallel')
else                      : logger = getLogger ( __name__         ) 
# =============================================================================
## possible types of workers 
workers = ( 'PATHOS'           , 
            'IPYPARALLEL'      , 
            'GAUDIMP'          ,
            'GAUDI'            ,
            'MULTIPROCESS'     ,
            'MULTIPRCCESSIONG' )
# =============================================================================
worker  = ''  
if has_env ( OSTAP_PARALLEL ) :
    # =========================================================================
    worker = get_env ( OSTAP_PARALLEL, '' ) .upper()
    if not worker in workers : worker = ''
# =============================================================================
if not worker : # =============================================================
    # =========================================================================
    import ostap.core.config as _CONFIG
    if 'PARALLEL' in _CONFIG.general :
        worker = _CONFIG.general.get( 'PARALLEL' , fallback = '' ).upper() 
        if not worker in workers : worker = ''

# =============================================================================
## Dill ?
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import dill
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    dill = None    

# =============================================================================
DILL_PY3_issue = False
# =============================================================================
if ( 3 , 6 ) <= sys.version_info and dill : # =================================
    # =========================================================================
    dill_version   =  getattr ( dill , '__version__' , '' )
    if not dill_version :  dill_version =  getattr ( dill , 'version' , '' )
    DILL_PY3_issue = dill_version < '0.3'
    if not DILL_PY3_issue : # =================================================
        from ostap.core.meta_info import root_info
        ## DILL_PY3_issue = root_info < ( 6 , 23 )
        DILL_PY3_issue = root_info < ( 6 , 24 , 6 )
    # =========================================================================
    if DILL_PY3_issue : worker = 'GAUDIMP'
    
# =============================================================================
## Check for ipyparallel
# =============================================================================
if not worker or 'IPYPARALLEL' == worker : # ==================================
    # =========================================================================
    if ( 3 , 6 ) <= sys.version_info : # ======================================
        # =====================================================================
        try : # ===============================================================
            # =================================================================
            with warnings.catch_warnings() :            
                warnings.simplefilter ( "ignore" )                
                import ipyparallel as _ipp
                worker = 'IPYPARALLEL' if ( 8 , 0 ) <= _ipp.version_info else ''
            # =================================================================
        except ImportError : # ================================================
            # =================================================================
            worker = ''
        # =====================================================================
    else : # ==================================================================
        # =====================================================================
        worker = ''

# =============================================================================
## Check for pathos
# =============================================================================
if not worker  or 'PATHOS' == worker : # ======================================
    # =========================================================================
    if not DILL_PY3_issue : # =================================================
        # =====================================================================
        try : # ===============================================================
            # =================================================================
            import pathos as _pathos 
            worker = 'PATHOS'
            # =================================================================
        except ImportError : # ================================================
            # =================================================================
            worker = ''
        # =====================================================================
    else : # ==================================================================
        # =====================================================================
        worker = ''
        
# ===============================================================================
WorkManager = None

# ===============================================================================
## Use ipyparallel?
# ===============================================================================
if 'IPYPARALLEL' == worker and ( 3 , 6 ) <= sys.version_info :
    # ===========================================================================
    try : # =====================================================================
        # =======================================================================
        from ostap.parallel.parallel_ipyparallel import WorkManager, Checker
        logger.debug ('Use WorkManager from ostap.parallel.ipyparallel')
        # =======================================================================
    except ImportError : # ======================================================
        # =======================================================================
        WorkManage, Checker = None, None
        worker      = 'PATHOS'

# ===============================================================================
## Use PATHOS?
# ===============================================================================
if 'PATHOS' == worker and not WorkManager and not DILL_PY3_issue :
    # ===========================================================================
    try : # =====================================================================
        from ostap.parallel.parallel_pathos import WorkManager, Checker
        logger.debug ('Use WorkManager from ostap.parallel.pathos')
        # =======================================================================
    except ImportError : # ======================================================
        # =======================================================================
        WorkManager, Checker = None, None

# ===============================================================================
## Use multiprocess/multiprocessing
# ===============================================================================
if not WorkManager : # ==========================================================
    # ===========================================================================
    from ostap.parallel.parallel_gaudi  import WorkManager, Checker
    logger.debug ('Use WorkManager from GaudiMP.Parallel' )
    worker = 'GAUDI'


if not Checker : from ostap.io.pickling import Checker

# =============================================================================
if '__main__' == __name__ :
    
    from ostap import banner
    logger.info ( __file__ + '\n' + banner )
    logger.info ( 80*'*' )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 


# =============================================================================
##                                                                      The END 
# =============================================================================
