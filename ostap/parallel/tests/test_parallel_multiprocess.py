#!/usr/bin/env python
# ============================================================================
## @file test_parallel_multiprocess.py
# Oversimplified script for parallel execution using multiproces
# @see https://github.com/uqfoundation/multiprocess
# @see https://github.com/uqfoundation/dill
# ============================================================================
""" Oversimplified script for parallel execution using multiprocessing
"""
from   __future__        import print_function
import ROOT, time, sys 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_parallel_multiprocess' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
from   itertools                import count   
# =============================================================================
try : 
    import dill 
except ImportError :
    logger.error('Can not import dill')
    dill = None    
# =============================================================================
try : 
    import multiprocess
except ImportError :
    logger.error('Can not import multiprocess')
    multiprocess = None

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
import ostap.histos.histos
from   ostap.utils.progress_bar import progress_bar 
# =============================================================================
## simple    function that created and  fill a histogram
def make_histos ( item ) :
    """Simple    function that creates and  fills a histogram
    """
    i, n = item 
    import ROOT, random 
    h1 = ROOT.TH1F ( 'h%d' %  i , '' , 100 , 0 , 10 )
    for i in range ( n ) : h1.Fill ( random.gauss (  5 ,  1 ) )
    return h1 


# ===================================================================================
## @class MakeHisto
#  helper class to create a fill histograms
class MakeHisto(object) :
    """Helper class to create a fill histograms
    """
    def process  ( self , item ) :
        i, n = item 
        import ROOT, random 
        h1 = ROOT.TH1F ( 'h%d' %  i , '' , 100 , 0 , 10 )
        for i in range ( n ) : h1.Fill ( random.gauss (  5 ,  1 ) )
        return h1
    def __call__ ( self ,  item ) :
        return self.process ( item )

mh  = MakeHisto  ()

## start 5 jobs, and for each job create the histogram with 100 entries 
inputs = 5 * [ 100 ]

# =============================================================================
## test parallel processing with multiprocess
def test_multiprocess_function () :
    """Test parallel processnig with multiprocess
    """
    logger =    getLogger ("ostap.test_multiprocess_function")
    logger.info ('Test job submission with %s' %  multiprocess ) 
    
    if not dill :
        logger.error ( "dill is not available" )
        return
        
    if not multiprocess :
        logger.error ( "multiprocess is not available" )
        return 

    if DILL_PY3_issue : 
         logger.warning ("test is disabled for Python %s (DILL/ROOT/PY3 issue)" )
         return
    
    ncpus = multiprocess.cpu_count() 
    
    from multiprocess import Pool
    
    pool = Pool  ( ncpus ) 
    
    jobs = pool.imap_unordered ( make_histos , zip ( count() ,  inputs ) )
    
    result = None 
    for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
        if not result  : result = h
        else           : result.Add ( h )

    pool.close ()
    pool.join  ()
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 ) )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 


# =============================================================================
## test parallel processing with multiprocess
def test_multiprocess_callable  () :
    """Test parallel processnig with multiprocess
    """
    logger =    getLogger ("ostap.test_multiprocess_callable")
    logger.info ('Test job submission with %s' %  multiprocess ) 
    
    if not dill :
        logger.error ( "dill is not available" )
        return
        
    if not multiprocess :
        logger.error ( "multiprocess is not available" )
        return 
        
    if DILL_PY3_issue : 
         logger.warning ("test is disabled for Python %s (DILL/ROOT/PY3 issue)" )
         return
    
    ncpus = multiprocess.cpu_count() 
    
    from multiprocess import Pool
    
    pool = Pool  ( ncpus ) 
    
    jobs = pool.imap_unordered ( mh , zip ( count() ,  inputs ) )
    
    result = None 
    for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
        if not result  : result = h
        else           : result.Add ( h )

    pool.close ()
    pool.join  ()
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 ) )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 

    
# =============================================================================
if '__main__' == __name__ :

    test_multiprocess_function () 
    test_multiprocess_callable () 

        
# =============================================================================
##                                                                      The END 
# =============================================================================
