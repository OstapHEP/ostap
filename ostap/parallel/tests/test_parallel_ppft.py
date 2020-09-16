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
# =============================================================================
try : 
    import dill
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

# ===================================================================================
## @class MakeHisto
#  helper class to create a fill histoghrams
class MakeHisto(object) :
    """Helper class to create a fill histoghrams
    """
    def process ( self , i , n ) :
        import ROOT
        import random
        h1 = ROOT.TH1F ( 'h%d' %  i , '' , 100 , 0 , 10 )
        for i in range ( n ) : h1.Fill ( random.gauss (  4 ,  1 ) )
        return h1 
    def __call__ ( self , i , n ) :
        return self.process ( i , n ) 

mh = MakeHisto()

## start 15 jobs, and for each job create the histogram with 100 entries 
inputs = 15 * [ 100 ]


# =============================================================================
## test parallel python with function 
def test_ppft_1 () :
    """Test parallel python (function)
    """

    if not dill :
        logger.error ( "test_ppft_1: dill is not available" )
        return 
    if not ppft :
        logger.error ( "test_ppft_1: ppft is not available" )
        return 
        
    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor : 
        vip = '%s.%s.%s' % ( vi.major , vi.minor , vi.micro ) 
        logger.warning ("test_ppft_1 is disabled for Python %s" % vip )
        return
    
    job_server = ppft.Server()
    
    jobs = [ ( i , job_server.submit ( make_histo , ( i , n ) , () , () ) ) for ( i , n ) in enumerate  ( inputs ) ]

    result = None 
    for input, job in progress_bar ( jobs ) :
        histo = job()
        if not result : result = histo
        else          :
            result.Add ( histo ) 
            del histo 

    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 


# =============================================================================
## test parallel python with functor
def test_ppft_2() :
    """Test parallel python  (functor) 
    """
    if not dill :
        logger.error ( "test_ppft_2: dill is not available" )
        return 
    if not ppft :
        logger.error ( "test_ppft_2: ppft is not available" )
        return 
        
    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        vip = '%s.%s.%s' % ( vi.major , vi.minor , vi.micro ) 
        logger.warning ("test_ppft_2 is disabled for Python %s" % vip )
        return
    
    job_server = ppft.Server()
    
    mh = MakeHisto() 
    jobs = [ ( i , job_server.submit ( mh.process , ( i , n ) , () , () ) ) for ( i , n ) in enumerate  ( inputs ) ]

    result = None 
    for input, job in progress_bar ( jobs ) :
        histo = job()
        if not result : result = histo
        else          :
            result.Add ( histo ) 
            del histo 

    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 
    
# =============================================================================
if '__main__' == __name__ :

    test_ppft_1 () 
    test_ppft_2 () 
        
# =============================================================================
##                                                                      The END 
# =============================================================================
