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
    'GenericTask' , ## very generic "template"  tasl 
    )
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'pstap.parallel.parallel')
else                      : logger = getLogger ( __name__         ) 
# =============================================================================
import operator

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
## @class GenericTask
#  Generic ``temlated'' task for Parallel processing  
#    One needs to  define three functions/functors:
#    - processor   :<code>        output = processor   ( item )               </code>
#    - merger      :<code>updated_output = merger ( old_output , new_output ) </code>
#    - initializer :<code>        output = initializer (      )               </code> 
class GenericTask(Task) :
    """Generic ``temlated'' task for parallel processing.
    One needs to  define three functions/functors:
    - processor   :         output = processor   ( item ) 
    - merger      : updated_output = merger ( old_output , new_output )
    - initializer :         output = initializer (      )  
    """
    # =========================================================================
    def __init__ ( self                       ,
                   processor                  ,
                   merger      = operator.add ,
                   initializer = tuple        ) :
        """Generic task for parallel processing. One needs to  define three functions/functors
        - processor   :         output = processor   ( item ) 
        - merger      : updated_output = merger ( old_output , new_output )
        - initializer :         output = initializer (      )  
        """        
        self.__processor   = processor
        self.__merger      = merger
        self.__initializer = initializer
        
    # =========================================================================
    ## local initialization (executed once in parent process)
    def initializeLocal   ( self ) :
        """Local initialization (executed once in parent process)"""
        self.output = self.initializer () if self.initializer else () 
        
    # =========================================================================
    ## the actual processing of the single item 
    def process  ( self , item ) :
        """The actual processing of the single item"""
        self.output = self.processor ( item )
        
    # =========================================================================
    ## merge results 
    def _mergeResults ( self , result ) :
        """Merge processing results"""
        self.output = self.merger ( self.output , result )

    # =========================================================================
    @property
    def processor  ( self ) :
        """``processor'' : the actual function for each subprocess
        - Signature: output = processor ( item ) 
        """
        return self.__processor
    @property
    def merger     ( self ) :
        """``merger'' : the actual fuction to merge results
        - Signature: updated_output = merger ( old_output , new_output )         
        """
        return self.__merger
    @property
    def initializer ( self ) :
        """``initializer'' : the actual fuction to initialize local output  
        - Signature: output = initializer() 
        """
        return self.__initializer
    
    
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
