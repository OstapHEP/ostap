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
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.parallel.copy' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
import ROOT
from   ostap.parallel.parallel import Task, WorkManager
from   ostap.core.ostap_types  import string_types, integer_types  
# =============================================================================
## The simple task object for parallel copy 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2022-02-15 
class  CopyTask(Task) :
    """The simple task object for parallel copy
    """
    ## 
    def __init__ ( self    ,
                   copier  ) :
        
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
## Copy files in parallel 
def copy_files ( file_pairs , progress = True , maxfiles = 5 , copier = None , **kwargs ) :
    """Copy files in parallel 
    """

    if not copier :
        from ostap.utils.utils import copy_file
        copier = copy_file
        
    task = CopyTask ( copier = copier )
    
    wmgr  = WorkManager ( silent = not progress , **kwargs )
    
    if maxfiles < 1 : maxfiles = 1
    
    from ostap.utils.utils import chunked
    data = chunked ( file_pairs , maxfiles )
    
    wmgr.process( task , data )

    copied = task.results () 
    return copied 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
