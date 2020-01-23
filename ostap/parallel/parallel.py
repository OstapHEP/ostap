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
    'Task'        , ## the base class for task
    'WorkManager' , ## task manager
    'GenericTask' , ## very generic "template"  tasl 
    )
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.parallel.parallel')
else                      : logger = getLogger ( __name__         ) 
# =============================================================================
from ostap.parallel.task import  Task, GenericTask

import os

workers = 'PATHOS' , 'GAUDIMP'

worker  = '' 

if 'OSTAP_PARALLEL' in os.environ :
    
    worker  = os.environ['OSTAP_PARALLEL'].upper()
    if not  worker in workers : worker = ''

if not worker :

    import ostap.core.config as _CONFIG
    if 'PARALLEL' in _CONFIG.general :
        worker = _CONFIG.general.get('PARALLEL', fallback = '' ).upper() 
        if not  worker in workers : worker = ''

# ===============================================================================
from sys import version_info  as python_version
if 3 == python_version.major and 6 >= python_version.minor :
    ## for python 3.6 dill fails to serialize ROOT objects
    worker = 'GAUDIMP'
    
# ===============================================================================

if  'GAUDIMP' != worker :
    
    try :
        from ostap.parallel.mp_pathos import WorkManager 
        logger.debug ('Use TaskManager from ostap.parallel.pathos')
    except ImportError :
        from ostap.parallel.mp_gaudi  import WorkManager 
        logger.info  ('Use TaskManager from GaudiMP.Parallel'     )

else :
    
    from ostap.parallel.mp_gaudi  import WorkManager 
    logger.debug ('Use TaskManager from GaudiMP.Parallel'         )

    
    
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
# The END 
# =============================================================================
