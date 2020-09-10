#!/usr/bin/env python
# ============================================================================
## @file test_parallel_ppft.py
# Oversimplified script for parallel execution using Parallel Python
# @see https://github.com/uqfoundation/pathos
# @see https://github.com/uqfoundation/ppft
# ============================================================================
""" Oversimplified script for parallel execution using Parallel Python
"""
from   __future__        import print_function
import ROOT, time, sys 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_parallel_ppft' )
else : 
    logger = getLogger ( __name__ )
try : 
    import dill, ppft
except ImportError :
    logger.error('Can not import dill')
    dill = None    
# =============================================================================
try : 
    import ppft
except ImportError :
    logger.error('Can not import ppft')
    ppft = None
    
# =============================================================================
import ostap.histos.histos
from   ostap.utils.progress_bar import progress_bar 
# =============================================================================
## simple    function that created and  fill a histogram 
def make_histo  ( i , n ) :
    """Simple    function that creates and  fills a histogram
    """
    import ROOT
    import random
    h1 = ROOT.TH1F ( 'h%d' %  i , '' , 100 , 0 , 10 )
    for i in range ( n ) : h1.Fill ( random.gauss (  4 ,  1 ) )
    return h1 

# =============================================================================
## test parallel python 
def test_ppft() :
    """Test parallel python
    """

    if not dill :
        logger.error ( "test_ppft: dill is not available" )
        return 
    if not ppft :
        logger.error ( "test_ppft: ppft is not available" )
        return 
        
    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        logger.warning ("test_ppft is disabled for Python %s" % vi )
        return
    
    job_server = ppft.Server()
    
    ## start 100 jobs, and for each job create the histogram with 1000 entries 
    inputs = 100 * [ 1000 ]
    
    ## submit the jobs 
    jobs = [ ( i , job_server.submit ( make_histo , ( i , n ) , () , () ) ) for ( i , n ) in enumerate  ( inputs ) ]

    result = None 
    for input, job in progress_bar ( jobs ) :
        histo = job()
        if not result : result = histo
        else          :
            result.Add ( histo ) 
            del histo 

    logger.info ( "Histogram is %s" % result )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 
    

# =============================================================================
if '__main__' == __name__ :

    test_ppft () 

        
# =============================================================================
##                                                                      The END 
# =============================================================================
