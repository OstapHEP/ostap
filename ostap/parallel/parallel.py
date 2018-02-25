#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file parallel.py 
#
#  Useful utilities for multiprocessing and parallel processing for Ostap
#  Actualy it is just a little bit upgraded version of original
#  GaudiMP.Parallel module developed by Pere MATO Pere.Mato@cern.ch
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-02-23
#
# =============================================================================
"""   Useful utilities for multiprocessing and parallel processing for Ostap
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
    )
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'pstap.parallel.parallel')
else                      : logger = getLogger ( __name__         ) 
# =============================================================================
## try:
         
##     from ostap.parallel.mp_pathos import Task, WorkManager 
##     logger.info  ('Use Task and TaskManager from ostap.parallel.pathos')
    
## except ImportErorr :
    
##     logger.error ("Can't import ostap.parallel.mp_pathos:" )  ## , exc_info = True )
##     from ostap.parallel.mp_gaudi import Task, WorkManager 
##     logger.info  ('Use Task and TaskManager from GaudiMP.Parallel'    )


from ostap.parallel.mp_gaudi import Task, WorkManager 
logger.info  ('Use Task and TaskManager from GaudiMP.Parallel'    )
    
    
    
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
