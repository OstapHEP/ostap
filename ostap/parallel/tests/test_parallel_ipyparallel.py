#!/usr/bin/env python
# ============================================================================
## @file test_parallel_ipyparallel.py
# Oversimplified script for parallel execution using multiprocessing
# ============================================================================
""" Oversimplified script for parallel execution using multiprocessing
"""
# =============================================================================
from   itertools                import count    
from   ostap.utils.progress_bar import progress_bar 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
import ostap.histos.histos
import ROOT, time, sys, warnings  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_parallel_ipyparallel' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================

ipp = None 
if ( 3 , 6 )<= sys.version_info : 
    try :
        with warnings.catch_warnings() :
            warnings.simplefilter('ignore')
            import ipyparallel as ipp
    except ImportError : 
        ipp = None 

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
    """Helper class to create a fill histoghrams
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

## start 50 jobs, and for each job create the histogram with 1000 entries 
inputs = 50 * [ 1000 ]

# =============================================================================
## test parallel processing with ipyparallel
def test_ipyparallel_function () :
    """Test parallel processnig with ipyparallel
    """

    logger = getLogger ( "ostap.test_ipyparallel_function")    
    logger.info ('Test job submission with ipyparallel')
    
    if not (3,6)<= sys.version_info :
        logger.error ( "python3.6 is required for the test!")
        
    if not ipp :
        logger.error ( "ipyparallel module is not available")
        return
    
    if not (8,0) <= ipp.version_info :
        logger.error ( "ipyparallel module is too old %s" % str ( ipp.version_info ) ) 
        return

    result = None 
    with ipp.Cluster( silent = True ) as cluster :

        view    = cluster.load_balanced_view()
        
        results = view.map_async ( make_histos , zip  ( count () , inputs ) )
        
        for r in progress_bar ( results ) :
            if not result  : result = r
            else           : result.Add ( r )
            
    with wait ( 3 ) , use_canvas ( 'test_ipyparallel_function' ) : 
        logger.info ( "Histogram is %s" % result.dump ( 80 , 20 ) )
        logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
        result.draw (   ) 

    return result

    
# =============================================================================
if '__main__' == __name__ :

    test_ipyparallel_function () 
        
# =============================================================================
##                                                                      The END 
# =============================================================================
