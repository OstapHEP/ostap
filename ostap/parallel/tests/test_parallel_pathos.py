#!/usr/bin/env python
# ============================================================================
## @file test_parallel_pathos.py
# Oversimplified script for parallel execution using Pathos 
# @see https://github.com/uqfoundation/pathos
# @see https://github.com/uqfoundation/dill
# ============================================================================
""" Oversimplified script for parallel execution using Pathos
"""
from   __future__        import print_function
import ROOT, time, sys 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_parallel_pathos' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
try : 
    import dill 
except ImportError :
    logger.error('Can not import dill')
    dill = None    
# =============================================================================
try : 
    import pathos
except ImportError :
    logger.error('Can not import pathos')
    pathos = None
    
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


# =============================================================================
## test parallel processing with pathos
def test_pathos () :
    """Test parallel processnig with pathos 
    """
    if not dill :
        logger.error ( "test_pathos: dill is not available" )
        return
        
    if not pathos :
        logger.error ( "test_pathos: pathos is not available" )
        return 
        
    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        logger.warning ("test_pathos is disabled for Python %s" % vi )
        return

    from pathos.helpers import cpu_count
    ncpus = cpu_count  ()
    
    from pathos.pools import ParallelPool
    pool = ParallelPool ( ncpus ) 
    pool.restart ( True ) 

    ## start 25 jobs, and for each job create the histogram with 1000 entries 
    inputs = 25 * [ 1000 ]

    # jobs = pool. imap ( make_histos ,  [  ( i , n )  for  ( i , n ) in enumerate ( inputs ) ] )
    jobs = pool.uimap ( make_histos ,  [  ( i , n )  for  ( i , n ) in enumerate ( inputs ) ] )
    
    result = None 
    for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
        if not result  : result = h
        else           : result.Add ( h )

    pool.close ()
    pool.join  ()
    pool.clear ()
    
    logger.info ( "Histogram is %s" % result )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 
    

# =============================================================================
if '__main__' == __name__ :

    ## pass
    test_pathos  () 

        
# =============================================================================
##                                                                      The END 
# =============================================================================
