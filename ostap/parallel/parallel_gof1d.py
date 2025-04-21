#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file ostap/parallel/parallel_gof1d.py
#  Make paralllel toys for Goodness-of-fit estimates 
#  @see ostap.stats.gof1d 
#  @date   2024-11-21
#  @author Vanya  BELYAEV Ivan.Belyaev@cern.ch
# =============================================================================
""" Make paralllel toys for Goodness-of-fit estimates 
- see ostap.stats.gof1d
"""
# =============================================================================
__author__  = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
__date__    = "2020-01-18"
__version__ = '$Revision$'
__all__     = (
    )
# =============================================================================
from   ostap.parallel.parallel import Task, WorkManager
from   ostap.core.ostap_types  import string_types, integer_types
from   ostap.utils.basic       import numcpu
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.gof1d' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## @class GoF1DTask
#  Simple task object 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-01-18 
class GoF1DTask (Task) :
    """ The GoF task object
    """
    ## 
    def __init__ ( self , gof  ) : 
        
        self.__gof        = gof 
        self.__the_output = None 
        
    @property
    def the_output ( self ) :
        return self.__the_output
    @the_output.setter 
    def the_output ( self , value ) :
        self.__the_output = value 
    
    ## initialize the local task, no action 
    def initialize_local   ( self ) : self.__the_output = None 

    ## initialize the remote task, treat the random numbers  
    def initialize_remote  ( self , jobid = -1 ) :
        """ Initialize the remote task, treta the random numbers  
        """
        import random, ROOT
        from ostap.parallel.utils import random_random
        random_random ( jobid )        
        return self.initialize_local() 
        
    ## get the results 
    def results ( self ) :
        return self.__the_output
    
    ## merge results of toys 
    def merge_results ( self , result , jobid = -1 ) :
        """ Merge results of toys
        """
        if not self.__the_output : self.__the_output  = result
        else                     : self.__the_output += result 

    # =========================================================================
    ## the actual processing of toys 
    def process ( self , jobid , nToys ) :
        """ The actual processing of toys 
        """        
        from   ostap.stats.gof1d import GoF1DToys
        
        toys = GoF1DToys ( gof = self.__gof )
        toys.run ( nToys    = nToys ,
                   parallel = False ,
                   silent   = True  )
        
        ## del self.__gof        
        #
        return toys 

# =============================================================================
## Run GoF1D toys in parallel 
def parallel_gof1dtoys ( gof             ,
                         nToys    = 1000 ,
                         nSplit   = 0    ,
                         silent   = True , 
                         progress = True , **kwargs ) :

    assert isinstance ( nToys  , integer_types ) and 0 < nToys  ,\
        'Invalid "nToys"  argument %s/%s' % ( nToys  , type ( nToys  ) )
    
    if not nSplit : nSplit = max ( 2 , 2 * numcpu() ) 
    
    assert isinstance ( nSplit , integer_types ) and 0 < nSplit ,\
        'Invalid "nSplit" argument %s/%s' % ( nSplit , type ( nSplit ) )

    if nSplit < 2 or numcpu ()  < 2 :
        from   ostap.stats.gof1d import GoF1DToys        
        toys = GoF1DToys ( gof = gof )
        toys.run ( nToys    = nToys                  ,
                   parallel = False                  ,
                   silent   = silent or not progress ) 
        return toys
        
    ## create work manager 
    wmgr  = WorkManager ( silent   = silent and not progress ,
                          progress = progress or not silent  , **kwargs )
    
    task  = GoF1DTask ( gof = gof )

    if nToys <= nSplit : params = nToys * [ 1 ]
    else : 
        size , rem = divmod ( nToys , nSplit )
        if rem : params = [ size + rem ] + ( nSplit - 1 ) * [ size ]
        else   : params =                    nSplit       * [ size ]

    ## start parallel processing! 
    wmgr.process ( task , params  )
    
    return task.results () 

    
# =============================================================================
##                                                                      The END 
# =============================================================================
