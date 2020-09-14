#!/usr/bin/env python
# ============================================================================
## @file test_parallel_multiprocessing.py
# Oversimplified script for parallel execution using multiprocessing
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
    logger = getLogger ( 'test_parallel_multiprocessing' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
import multiprocessing
    
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

## start 10 jobs, and for each job create the histogram with 100 entries 
inputs = 10 * [ 100 ]

# =============================================================================
## test parallel processing with pathos
def test_multiprocessing () :
    """Test parallel processnig with multiprocessing
    """
    
    ncpus = multiprocessing.cpu_count() 
    
    from multiprocessing import Pool
    
    pool = Pool  ( ncpus ) 

    jobs = pool.imap_unordered ( make_histos ,  [  ( i , n )  for  ( i , n ) in enumerate ( inputs ) ] )
    
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

    ## pass
    test_multiprocessing  () 

        
# =============================================================================
##                                                                      The END 
# =============================================================================
