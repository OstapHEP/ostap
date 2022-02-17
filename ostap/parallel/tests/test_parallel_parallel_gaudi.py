#!/usr/bin/env python
# ============================================================================
## @file test_parallel_parallel_gaudi.py
# Oversimplified script for parallel execution using parallel_gaudi
# ============================================================================
""" Oversimplified script for parallel execution using parallel_gaudi
"""
# ============================================================================
import ROOT, time, sys 
# =============================================================================
from   itertools                import count 
import ostap.histos.histos
from   ostap.parallel.task      import Task, GenericTask
from   ostap.parallel.utils     import pool_context 
from   ostap.utils.progress_bar import progress_bar 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_parallel_parallel_gaudi' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
try :
    from ostap.parallel.parallel_gaudi import WorkManager 
except ImportError :
    logger.error ("Cannot import WorkManager from parallel_gaudi")
    WorkManager = None 
# =============================================================================


# =============================================================================
## simple    function that created and  fill a histogram
def make_histo  ( jobid , n ) :
    """Simple    function that creates and  fills a histogram
    """
    import ROOT, random 
    h1 = ROOT.TH1F ( 'h%d' %  jobid , '' , 100 , 0 , 10 )
    for i in range ( n ) : h1.Fill ( random.gauss (  5 ,  1 ) )
    return h1 

# =============================================================================
## simple    function that created and  fill a histogram
def make_histo2  ( param ) :
    """Simple    function that creates and fills a histogram
    """
    jobid , n  = param
    return make_histo ( jobid , n ) 

# =============================================================================
## simple "merger" for histograms
def merge_histos  ( h1 , h2 ) :
    """Simple ``merger'' for historgams"""
    if not h1 : return h2 
    h1.Add (  h2 )
    return h1

# =============================================================================
## start 5 jobs, and for each job create the histogram with 100 entries 
inputs = 5 * [ 100 ] 

# ==============================================================================
## simple task to create and fill historgam 
class HTask(Task) :
    """Simple task to create and fill historgam
    """
    def __init__ (  self )                  : self.__histo = None
    def initialize_local  ( self )          : self.__histo = None
    def process  ( self  , jobid , n ) :        
        import ROOT, random 
        h1 = ROOT.TH1F ( 'h%d' %  jobid , '' , 100 , 0 , 10 )
        for i in range ( n ) : h1.Fill ( random.gauss (  5 ,  1 ) )
        self.__histo = h1 
        return self.__histo 
    def merge_results ( self , result , jobid = -1 ) :        
        if not self.__histo : self.__histo = result
        else                : self.__histo.Add ( result ) 
    def results ( self ) :
        return self.__histo

# =============================================================================
## test parallel processing with parallel_gaudi (bare interface) 
def test_parallel_gaudi_mp_bare ( ) :
    """Test parallel processnig with parallel_gaudi (bare interface) 
    """
    logger  = getLogger ("ostap.test_parallel_gaudi_mp_bare")
    if not WorkManager :
        logger.error ("Failure to import WorkManager")
        return
    
    logger.info ('Test job submission with %s' % WorkManager  ) 

    
    ## vi = sys.version_info
    ## if 3<= vi.major and 6 <= vi.minor :
    ##     vip = '%s.%s.%s' % ( vi.major , vi.minor , vi.micro ) 
    ##     logger.warning ("test is disabled for Python %s (dill/ROOT issue)" % vip )
    ##     return
    
    ## create the manager 
    manager = WorkManager ( silent = False  )

    ## initialize  result 
    result   = None
    
    ## use the bare interface 
    for res in manager.iexecute ( make_histo2 , zip ( count() , inputs ) , progress = True ) :
        if result is None  : result = res
        else               : result.Add ( res )  

    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result


# =============================================================================
## test parallel processing with parallel_gaudi (task interface) 
def test_parallel_gaudi_mp_task ( ) :
    """Test parallel processnig with parallel_gaudi (task interface) 
    """
    logger  = getLogger ("ostap.test_parallel_gaudi_mp_task")
    if not WorkManager :
        logger.error ("Failure to import WorkManager")
        return
    
    logger.info ('Test job submission with %s' % WorkManager  ) 

    
    ## vi = sys.version_info
    ## if 3<= vi.major and 6 <= vi.minor :
    ##     vip = '%s.%s.%s' % ( vi.major , vi.minor , vi.micro ) 
    ##     logger.warning ("test is disabled for Python %s (dill/ROOT issue)" % vip )
    ##     return
    
    ## create the manager 
    manager = WorkManager ( silent = False  )

    ## create the task 
    task     = HTask()

    ## process the task 
    result   = manager.process ( task ,  inputs ) 
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with wait ( 1 ) , use_canvas ( 'test_parallel_gaudi_mp_task' ) : 
        result.draw (   ) 
    
    return result

    
# =============================================================================
## test parallel processing with parallel_gaudi (func interface) 
def test_parallel_gaudi_mp_func ( ) :
    """Test parallel processnig with parallel_gaudi (func interface) 
    """
    logger  = getLogger ("ostap.test_parallel_gaudi_mp_task")
    if not WorkManager :
        logger.error ("Failure to import WorkManager")
        return
    
    logger.info ('Test job submission with %s' % WorkManager  ) 

    
    ## vi = sys.version_info
    ## if 3<= vi.major and 6 <= vi.minor :
    ##     vip = '%s.%s.%s' % ( vi.major , vi.minor , vi.micro ) 
    ##     logger.warning ("test is disabled for Python %s (dill/ROOT issue)" % vip )
    ##     return
    
    ## create the manager 
    manager = WorkManager ( silent = False  )

    ## process the function  
    result   = manager.process ( make_histo , inputs , merger = merge_histos )
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with wait ( 1 ) , use_canvas ( 'test_parallel_gaudi_mp_func' ) : 
        result.draw (   ) 

    return result

    
# =============================================================================
## test parallel processing with parallel_gaudi (use generic task)
def test_parallel_gaudi_mp_generic ( ) :
    """Test parallel processnig with parallel_gaudi (use generic task)
    """
    logger  = getLogger ("ostap.test_parallel_gaudi_mp_generic")
    if not WorkManager :
        logger.error ("Failure to import WorkManager")
        return
    
    logger.info ('Test job submission with %s' % WorkManager  ) 

    ## vi = sys.version_info
    ## if 3<= vi.major and 6 <= vi.minor :
    ##     vip = '%s.%s.%s' % ( vi.major , vi.minor , vi.micro ) 
    ##     logger.warning ("test is disabled for Python %s (dill/ROOT issue)" % vip )
    ##     return
    
    ## create the manager 
    manager = WorkManager ( silent = False  )

    ## create the task 
    task    = GenericTask ( processor = make_histo   ,
                            merger    = merge_histos )    
    
    ## process the task 
    result   = manager.process ( task ,  inputs ) 
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with wait ( 1 ) , use_canvas ( 'test_parallel_gaudi_mp_generic' ) : 
        result.draw (   ) 
        
    return result


# =============================================================================
if '__main__' == __name__ :

    ## bare interface  
    test_parallel_gaudi_mp_bare    ()
    
    ## task interface  
    test_parallel_gaudi_mp_task    ()

    ## function interface 
    test_parallel_gaudi_mp_func    ()
    
    ## use generic task 
    test_parallel_gaudi_mp_generic ()
    
        
# =============================================================================
##                                                                      The END 
# =============================================================================
