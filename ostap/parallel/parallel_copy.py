#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file ostap/parallel/parallel_copy.py
#  Copy files in parallel
#  @date   2020-01-18
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
# =============================================================================
""" Make fitting toys in parallel
- see ostap.fitting.toys 
"""
# =============================================================================
__author__  = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
__date__    = "2020-01-18"
__version__ = '$Revision$'
__all__     = (
    'copy_files' , ## copy files in parallel
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
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.copy' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
## The simple task object for parallel copy 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2022-02-15 
class  CopyTask(Task) :
    """The simple task object for parallel copy
    """
    ## 
    def __init__ ( self , copier  ) :        
        self.__copier  = copier 
        self.__the_output   = () 

    @property
    def the_output ( self ) :
        return self.__the_output
    @the_output.setter 
    def the_output ( self , value ) :
        self.__the_output = value 
    
    def initialize_local   ( self ) : self.__the_output = ()

    ## get the results 
    def results ( self ) :
        return self.__the_output
    
    ## the actual processing 
    def process ( self , jobid , file_pairs ) :
        
        results = [] 
        for source , destination in file_pairs :
            result = self.__copier ( source , destination , progress = False )
            result = source , result            
            results.append ( result )
            
        self.the_output = tuple ( results ) 
        
        return self.results() 

    ## merge results of toys 
    def merge_results ( self , result , jobid = -1 ) :
        """Merge results of toys
        """
        self.the_output = self.the_output + result 
        
        return self.results () 

# =============================================================================
## Copy files in parallel:
#  @return sequence of (input,output) pairs
def copy_files ( file_pairs , progress = True , maxfiles = 5 , copier = None , **kwargs ) :
    """Copy files in parallel: 
    - return sequence of (input,output) pairs
    """
    if not copier :
        from ostap.utils.utils import copy_file
        copier = copy_file

    pairs = tuple ( p for p in file_pairs )
    # ========================================================================
    ## sequential copy 
    # =======================================================================
    nfiles = len ( pairs ) 
    if nfiles <= 1 or numcpu () <=1  :
        from ostap.utils.progress_bar import progress_bar
        silent = nfiles <= 1 or not progress 
        copied = [] 
        for f, nf in progress_bar ( pairs , silent = silent ) :
            output = copier ( f , nf , progress = progress and nfiles <=1 )
            result = f , output 
            copied.append ( result ) 
        return tuple ( copied )
    
    # ========================================================================
    ## start parallel processing
    # =======================================================================
    from ostap.utils.utils import chunked
    data   = chunked ( pairs , max ( maxfiles , 1 ) )
    task   = CopyTask ( copier = copier )    
    wmgr   = WorkManager ( silent = not progress , progress = progress , **kwargs )    
    wmgr.process ( task , data )
    copied = task.results ()

    ## check the final results 
    for f , nf in pairs :
        if os.path.exists ( nf ) and os.path.isfile ( nf ) : pass 
        else :
            logger.warning ( "copy_files: no expected output '%s'" % nf ) 

    return tuple ( copied )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
