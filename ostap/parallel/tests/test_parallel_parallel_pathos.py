#!/usr/bin/env python
# ============================================================================
## @file test_parallel_parallel_pathos.py
# Oversimplified script for parallel execution using parallel_pathos
# ============================================================================
""" Oversimplified script for parallel execution using parallel_pathos
"""
from   __future__        import print_function
import ROOT, time, sys 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_parallel_parallel_pathos' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
from   itertools            import count 
from   ostap.parallel.task  import Task, GenericTask
from   ostap.parallel.utils import pool_context 
import ostap.histos.histos
from   ostap.utils.progress_bar import progress_bar 
# =============================================================================
try :
    from ostap.parallel.parallel_pathos import WorkManager 
except ImportError :
    logger.error ("Cannot import WorkManager from parallel_pathos")
    WorkManager = None 
# =============================================================================
try : 
    import dill 
except ImportError :
    logger.error('Can not import dill')
    dill = None    

DILL_PY3_issue = False 
if ( 3 , 6 ) <= sys.version_info and dill :
    dill_version =  getattr ( dill , '__version__' , '' )
    if not dill_version :  dill_version =  getattr ( dill , 'version' , '' )
    DILL_PY3_issue = dill_version < '0.3'
    if not DILL_PY3_issue :
        from ostap.core.meta_info import root_info
        DILL_PY3_issue = root_info < ( 6 , 23 )
        
if DILL_PY3_issue : logger.warning ( "There is an issue with DILL/ROOT/PYTHON")

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
## test parallel processing with parallel_pathos (bare interface) 
def test_parallel_pathos_mp_bare ( ) :
    """Test parallel processnig with parallel_pathos (bare interface) 
    """
    logger  = getLogger ("test_parallel_pathos_mp_bare")
    if not WorkManager :
        logger.error ("Failure to import WorkManager")
        return
    
    logger.info ('Test job submission with %s' % WorkManager  ) 
    
    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return
    
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
## test parallel processing with parallel_pathos (bare interface) 
def test_parallel_pathos_pp_bare ( ) :
    """Test parallel processnig with parallel_pathos (bare interface) 
    """
    logger  = getLogger ("test_parallel_pathos_mp_bare")
    if not WorkManager :
        logger.error ("Failure to import WorkManager")
        return 
    
    logger.info ('Test job submission with %s' % WorkManager  ) 

    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return

    ## create the manager 
    manager  = WorkManager ( silent = False  , pp = True )
        
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
## test parallel processing with parallel_pathos (task interface) 
def test_parallel_pathos_mp_task ( ) :
    """Test parallel processnig with parallel_pathos (task interface) 
    """
    logger  = getLogger ("test_parallel_pathos_mp_task")
    if not WorkManager :
        logger.error ("Failure to import WorkManager")
        return
    
    logger.info ('Test job submission with %s' % WorkManager  ) 

    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return

    
    ## create the manager 
    manager = WorkManager ( silent = False  )

    ## create the task 
    task     = HTask()

    ## process the task 
    result   = manager.process ( task ,  inputs ) 
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result


# =============================================================================
## test parallel processing with parallel_pathos (task interface) 
def test_parallel_pathos_pp_task ( ) :
    """Test parallel processnig with parallel_pathos (task interface) 
    """
    logger  = getLogger ("test_parallel_pathos_pp_task")
    if not WorkManager :
        logger.error ("Failure to import WorkManager")
        return
    
    logger.info ('Test job submission with %s' % WorkManager  ) 

    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return

    ## create the manager 
    manager = WorkManager ( silent = False , pp = True  )

    ## create the task 
    task     = HTask()

    ## process the task 
    result   = manager.process ( task ,  inputs ) 
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result

    
# =============================================================================
## test parallel processing with parallel_pathos (func interface) 
def test_parallel_pathos_mp_func ( ) :
    """Test parallel processnig with parallel_pathos (func interface) 
    """
    logger  = getLogger ("test_parallel_pathos_mp_task")
    if not WorkManager :
        logger.error ("Failure to import WorkManager")
        return
    
    logger.info ('Test job submission with %s' % WorkManager  ) 
    
    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return
    
    ## create the manager 
    manager = WorkManager ( silent = False  )

    ## process the function  
    result   = manager.process ( make_histo , inputs , merger = merge_histos )
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result

    
# =============================================================================
## test parallel processing with parallel_pathos (func interface) 
def test_parallel_pathos_pp_func ( ) :
    """Test parallel processnig with parallel_pathos (func interface) 
    """
    logger  = getLogger ("test_parallel_pathos_pp_task")
    if not WorkManager :
        logger.error ("Failure to import WorkManager")
        return
    
    logger.info ('Test job submission with %s' % WorkManager  ) 

    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return

    ## create the manager 
    manager = WorkManager ( silent = False , pp = True )

    ## process the function  
    result   = manager.process ( make_histo , inputs , merger = merge_histos )
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result


# =============================================================================
## test parallel processing with parallel_pathos (use generic task)
def test_parallel_pathos_mp_generic ( ) :
    """Test parallel processnig with parallel_pathos (use generic task)
    """
    logger  = getLogger ("test_parallel_pathos_mp_generic")
    if not WorkManager :
        logger.error ("Failure to import WorkManager")
        return
    
    logger.info ('Test job submission with %s' % WorkManager  ) 

    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return
    
    ## create the manager 
    manager = WorkManager ( silent = False  )

    ## create the task 
    task    = GenericTask ( processor = make_histo   ,
                            merger    = merge_histos )    
    
    ## process the task 
    result   = manager.process ( task ,  inputs ) 
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result

# =============================================================================
## test parallel processing with parallel_pathos (use generic task)
def test_parallel_pathos_pp_generic ( ) :
    """Test parallel processnig with parallel_pathos (use generic task)
    """
    logger  = getLogger ("test_parallel_pathos_mp_generic")
    if not WorkManager :
        logger.error ("Failure to import WorkManager")
        return
    
    logger.info ('Test job submission with %s' % WorkManager  ) 

    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return

    ## create the manager 
    manager = WorkManager ( silent = False  , pp = True )

    ## create the task 
    task     = GenericTask ( processor = make_histo   ,
                             merger    = merge_histos )    

    ## process the task 
    result   = manager.process ( task ,  inputs ) 
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result




# =============================================================================
if '__main__' == __name__ :

    ## bare interface  
    test_parallel_pathos_mp_bare    ()
    test_parallel_pathos_pp_bare    ()
    
    ## task interface  
    test_parallel_pathos_mp_task    ()
    test_parallel_pathos_pp_task    ()

    ## function interface 
    test_parallel_pathos_mp_func    ()
    test_parallel_pathos_pp_func    ()
    
    ## use generic task 
    test_parallel_pathos_mp_generic ()
    test_parallel_pathos_pp_generic ()
    
        
# =============================================================================
##                                                                      The END 
# =============================================================================
