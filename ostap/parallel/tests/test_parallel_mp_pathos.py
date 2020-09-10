#!/usr/bin/env python
# ============================================================================
## @file test_parallel_mp_pathos.py
# Oversimplified script for parallel execution using mp_pathos
# ============================================================================
""" Oversimplified script for parallel execution using mp_pathos
"""
from   __future__        import print_function
import ROOT, time, sys 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_parallel_mp_pathos' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
from ostap.parallel.task      import Task, GenericTask
try :
    from ostap.parallel.mp_pathos import WorkManager 
except ImportError :
    logger.error ("Cannot import mp_pathos")
    WorkManager = None 

# =============================================================================
import ostap.histos.histos
from   ostap.utils.progress_bar import progress_bar 
# =============================================================================
## simple    function that created and  fill a histogram
def make_histos ( jobid , n ) :
    """Simple    function that creates and  fills a histogram
    """
    import ROOT, random 
    h1 = ROOT.TH1F ( 'h%d' %  jobid , '' , 100 , 0 , 10 )
    for i in range ( n ) : h1.Fill ( random.gauss (  5 ,  1 ) )
    return h1 

# =============================================================================
## simple "merger" for historgams 
def merge_histos  ( h1 , h2 ) :
    """Simple ``merger'' for historgams"""
    if h1 :
        h1.Add (  h2 )
        return h1
    return h2 

# ==============================================================================
## simple task to create and fill historgam 
class HTask(Task) :
    """Simple task to create and fill historgam
    """
    def __init__ (  self )                  : self.__histo = None
    def initialize_local  ( self )          : self.__histo = None
    def process  ( self  , jobid , n ) :        
        import ROOT, random 
        h1 = make_histos ( jobid , n )
        self.__histo = h1 
        return self.__histo 
    def merge_results ( self , result ) : 
        if not self.__histo : self.__histo = result
        else                : self.__histo.Add ( result ) 
    def results ( self ) :
        return self.__histo

# =============================================================================
## test parallel processing with mp_pathos
def test_mp_pathos () :
    """Test parallel processnig with mp_pathos
    """

    if not WorkManager :
        logger.error ("test_mp_pathos: cannot import WorkManager")
        return 
        
    ## start 25 jobs, and for each job create the histogram with 1000 entries 
    inputs = 25 * [ 1000 ]

    task     = HTask() 
    manager  = WorkManager ( silent = False )

    manager.process ( task ,  inputs ) 

    result = task.results() 
    
    logger.info ( "Histogram is %s" % result )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 

    
# =============================================================================
## test parallel processing with mp_pathos using generic task 
def test_mp_pathos_generic () :
    """Test parallel processnig with mp_paths + GenericTask
    """

    if not WorkManager :
        logger.error ("test_mp_pathos: cannot import WorkManager")
        return 

    ## start 25 jobs, and for each job create the histogram with 1000 entries 
    inputs = 25 * [ 1000 ]
    
    task     = GenericTask( processor = make_histos  ,
                            merger    = merge_histos )
    
    manager  = WorkManager ( silent = False )

    manager.process ( task ,  inputs ) 

    result = task.results() 
    
    logger.info ( "Histogram is %s" % result )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 
    

# =============================================================================
if '__main__' == __name__ :

    test_mp_pathos         () 
    test_mp_pathos_generic () 

        
# =============================================================================
##                                                                      The END 
# =============================================================================
