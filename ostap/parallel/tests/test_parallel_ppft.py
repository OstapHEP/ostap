#!/usr/bin/env python
# ============================================================================
## @file test_parallel_ppft.py
# Oversimplified script for parallel execution using Parallel Python
# @see https://github.com/uqfoundation/pathos
# @see https://github.com/uqfoundation/ppft
# ============================================================================
""" Oversimplified script for parallel execution using Parallel Python
"""
# =============================================================================
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.root_utils   import batch_env  
from   ostap.utils.progress_bar import progress_bar
from   ostap.parallel.utils     import uimap, fix_ppsrv 

import ostap.histos.histos
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
batch_env ( logger ) 
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import dill
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    logger.error('Can not import dill')
    dill = None
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import ppft
    fix_ppsrv ( ppft.Server )
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    logger.error('Can not import ppft')
    ppft = None

# =============================================================================
## simple  function that creates and  fills a histogram 
def make_histo  ( i , n ) :
    """ Simple    function that creates and fills a histogram
    """
    import ROOT
    import random
    h1 = ROOT.TH1F ( 'h%d' %  i , '' , 100 , 0 , 10 )
    for i in range ( n ) : h1.Fill ( random.gauss (  4 ,  1 ) )
    return h1 

# ===================================================================================
## @class MakeHisto
#  helper class to create and fill histograms
class MakeHisto(object) :
    """ Helper class to create and fill histoghrams
    """
    def process  ( self ,   *params ) :
        return make_histo ( *params )
    def __call__ ( self ,   *params ) :
        return make_histo ( *params )
    
mh  = MakeHisto  ()

## start NN jobs, and for each job create the histogram with 100 entries
NN     = 10  
inputs = NN * [ 100 ]

# =============================================================================
## test parallel python with with plain function 
def test_ppft_function () :
    """ Test parallel python with plain function
    """
    logger =    getLogger ("test_ppft_function")
    logger.info ('Test job submission with %s' %  ppft ) 
                  
    if not ppft :
        logger.error ( "ppdf is not available" )
        return 
    
    job_server = ppft.Server()
    
    jobs = [ ( i , job_server.submit ( make_histo , ( i , n ) ) ) for ( i , n ) in enumerate  ( inputs ) ]
    
    result = None 
    for input, job in progress_bar ( uimap ( jobs ) , max_value = len ( jobs ) ) :
        histo = job()
        if not result : result = histo
        else          :
            result.Add ( histo ) 
            del histo 

    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    job_server.print_stats()
    
    with use_canvas ( 'test_ppft_function' , wait = 1 ) : 
        result.draw (   ) 

    return result 

# =============================================================================
## test parallel python with object method  
def test_ppft_method() :
    """Test parallel python with object method  
    """
    logger =    getLogger ("test_ppft_method")
    logger.info ('Test job submission with %s' %  ppft ) 
    
    if not ppft :
        logger.error ( "ppft is not available" )
        return 
        
    job_server = ppft.Server()
    
    jobs = [ ( i , job_server.submit ( mh.process , ( i , n ) ) ) for ( i , n ) in enumerate  ( inputs ) ]

    result = None 
    for input, job in progress_bar ( uimap ( jobs ) , max_value = len ( jobs ) ) :
        histo = job()
        if not result : result = histo
        else          :
            result.Add ( histo ) 
            del histo 

    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    job_server.print_stats()

    with use_canvas ( 'test_ppft_method' , wait = 1 ) : 
        result.draw (   ) 

    return result 


# =============================================================================
## test parallel python with callable 
def test_ppft_callable1 () :
    """Test parallel python with callable  
    """
    logger = getLogger ("test_ppft_callable1")
    logger.info ('Test job submission with %s' %  ppft ) 
    
    if not ppft :
        logger.error ( "ppft is not available" )
        return 
        
    job_server = ppft.Server()
    
    jobs = [ ( i , job_server.submit ( mh.__call__  , ( i , n ) ) ) for ( i , n ) in enumerate  ( inputs ) ]

    result = None 
    for input, job in progress_bar ( jobs ) :
        histo = job()
        if not result : result = histo
        else          :
            result.Add ( histo ) 
            del histo 

    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with use_canvas ( 'test_ppft_callable1' , wait = 1 ) : 
        result.draw (   ) 

    return result 

# =============================================================================
## test parallel python with callable 
def test_ppft_callable2 () :
    """ Test parallel python with callable  
    """
    logger = getLogger ("test_ppft_callable2")
    logger.info ('Test job submission with %s' %  ppft ) 
    
    if not ppft :
        logger.error ( "ppft is not available" )
        return 
        
    job_server = ppft.Server()
    
    jobs = [ ( i , job_server.submit ( mh , ( i , n ) ) ) for ( i , n ) in enumerate  ( inputs ) ]

    result = None 
    for input, job in progress_bar ( jobs ) :
        histo = job()
        if not result : result = histo
        else          :
            result.Add ( histo ) 
            del histo 

    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with use_canvas ( 'test_ppft_callable2' , wait = 1 ) : 
        result.draw (   ) 

    return result 

# =============================================================================
if '__main__' == __name__ :

    import warnings 
    with warnings.catch_warnings( category = ResourceWarning ):
        warnings.simplefilter ( "ignore" )
        
        test_ppft_function  () 
        test_ppft_method    () 
        test_ppft_callable1 () 
        test_ppft_callable2 () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
