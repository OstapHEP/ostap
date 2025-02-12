#!/usr/bin/env python
# ============================================================================
## @file test_parallel_multiprocessing.py
# Oversimplified script for parallel execution using multiprocessing
# ============================================================================
""" Oversimplified script for parallel execution using multiprocessing
"""
# =============================================================================
from   itertools                import count    
from   ostap.utils.progress_bar import progress_bar 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait, batch_env 
import ostap.histos.histos
import multiprocessing
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
batch_env ( logger ) 
# =============================================================================

# =============================================================================
## simple    function that created and  fill a histogram
def make_histos ( item ) :
    """ Simple    function that creates and  fills a histogram
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
    """ Helper class to create a fill histoghrams
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

## start 50 jobs, and for each job create the histogram with 100 entries 
inputs = 50 * [ 100 ]

# =============================================================================
## test parallel processing with pathos
def test_multiprocessing_function () :
    """Test parallel processnig with multiprocessing
    """

    logger = getLogger ("test_multiprocessing_function")    
    logger.info ('Test job submission with module %s' %  multiprocessing )
    
    ncpus = multiprocessing.cpu_count() 
    
    from multiprocessing import Pool
    
    pool = Pool  ( ncpus ) 

    
    jobs = pool.imap_unordered ( make_histos , zip  ( count () , inputs ) )
    
    result = None 
    for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
        if not result  : result = h
        else           : result.Add ( h )

    pool.close ()
    pool.join  ()
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 ) )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 

    
    with wait ( 5 ) , use_canvas ( 'test_multiprocessing_function' ) : 
        result.draw (   ) 

    return result 


# =============================================================================
## test parallel processing
def test_multiprocessing_callable  () :
    """Test parallel processnig with multiprocessing
    """

    logger = getLogger ("test_multiprocessing_callable")
    
    logger.info ('Test job submission with module %s' %  multiprocessing )
    
    ncpus = multiprocessing.cpu_count() 
    
    from multiprocessing import Pool
    
    pool = Pool  ( ncpus ) 
    
    jobs = pool.imap_unordered ( mh , zip  ( count () , inputs ) )
    
    result = None 
    for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
        if not result  : result = h
        else           : result.Add ( h )

    pool.close ()
    pool.join  ()
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 ) )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with wait ( 5 ) , use_canvas ( 'test_multiprocessing_callable' ) : 
        result.draw (   ) 

    return result 

    
# =============================================================================
if '__main__' == __name__ :

    test_multiprocessing_function () 
    test_multiprocessing_callable () 

        
# =============================================================================
##                                                                      The END 
# =============================================================================
