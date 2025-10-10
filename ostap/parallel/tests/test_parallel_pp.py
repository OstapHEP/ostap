#!/usr/bin/env python
# ============================================================================
## @file test_parallel_pp.py
# Oversimplified script for parallel execution using Parallel Python
# @see https://www.parallelpython.com/examples.php#CALLBACK
# ============================================================================
""" Oversimplified script for parallel execution using Parallel Python
- see https://www.parallelpython.com/examples.php#CALLBACK
"""
# =============================================================================
from   ostap.utils.progress_bar import progress_bar 
from   ostap.plotting.canvas    import use_canvas
from   ostap.parallel.utils     import uimap, fix_ppsrv 
from   ostap.utils.root_utils   import batch_env 
import ostap.histos.histos
import ROOT, random, time, sys, pp  
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_parallel_pp' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
batch_env ( logger ) 
# =============================================================================
fix_ppsrv ( pp.Server )
    
## simple    function that created and  fill a histogram 
def make_histo  ( i , n ) :
    """ Simple    function that creates and  fills a histogram
    """
    import ROOT
    import random
    h1 = ROOT.TH1F ( 'h%d' %  i , '' , 100 , 0 , 10 )
    for i in range ( n ) : h1.Fill ( random.gauss (  4 ,  1 ) )
    return h1 

# ===================================================================================
## @class MakeHisto
#  helper class to create a fill histograms
class MakeHisto(object) :
    """ Helper class to create a fill histoghrams
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
def test_pp_function () :
    """ Test parallel python with plain function
    """
    logger =    getLogger ("test_pp_function")
    logger.info ('Test job submission with %s' %  pp ) 

    job_server = pp.Server()
    
    jobs = [ ( i , job_server.submit ( make_histo , ( i , n ) ) ) for ( i , n ) in enumerate  ( inputs ) ]
    
    result = None 
    for input, job in progress_bar ( uimap ( jobs ) , max_value = len ( jobs ) ) :
        histo = job()
        if not result : result = histo
        else          :
            result.Add ( histo ) 
            del histo 

    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    job_server.print_stats()
    
    with use_canvas ( 'test_pp_function' , wait = 1 ) : 
        result.draw (   ) 

    return result 

# =============================================================================
## test parallel python with object method  
def test_pp_method() :
    """ Test parallel python with object method  
    """
    logger =    getLogger ("test_pp_method")
    logger.info ('Test job submission with %s' %  pp ) 
 
    job_server = pp.Server()    
    jobs = [ ( i , job_server.submit ( mh.process , ( i , n ) ) ) for ( i , n ) in enumerate  ( inputs ) ]

    result = None 
    for input, job in progress_bar ( uimap ( jobs ) , max_value = len ( jobs ) ) :
        histo = job()
        if not result : result = histo
        else          :
            result.Add ( histo ) 
            del histo 

    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    job_server.print_stats()

    with use_canvas ( 'test_pp_method' , wait = 1 ) : 
        result.draw (   ) 

    return result 

# =============================================================================
## test parallel python with callable 
def test_pp_callable1 () :
    """ Test parallel python with callable  
    """
    logger = getLogger ("test_pp_callable1")
    logger.info ('Test job submission with %s' %  pp ) 
    
    ## logger.warning ("test is disabled for UNKNOWN REASON")
    ## return

    job_server = pp.Server()
    
    jobs = [ ( i , job_server.submit ( mh.__call__ , ( i , n ) ) ) for ( i , n ) in enumerate  ( inputs ) ]

    result = None 
    for input, job in progress_bar ( uimap ( jobs ) , max_value = len ( jobs ) ) :
        histo = job()
        if not result : result = histo
        else          :
            result.Add ( histo ) 
            del histo 

    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with use_canvas ( 'test_pp_callable1' , wait = 1 ) : 
        result.draw (   ) 

    return result 

# =============================================================================
## test parallel python with callable2 
def test_pp_callable2 () :
    """ Test parallel python with callable2  
    """
    logger = getLogger ("test_pp_callable2")
    logger.info ('Test job submission with %s' %  pp ) 
    
    ## logger.warning ("test is disabled for UNKNOWN REASON")
    ## return

    job_server = pp.Server()
    
    jobs = [ ( i , job_server.submit ( mh , ( i , n ) ) ) for ( i , n ) in enumerate  ( inputs ) ]

    result = None 
    for input, job in progress_bar ( uimap ( jobs ) , max_value = len ( jobs ) ) :
        histo = job()
        if not result : result = histo
        else          :
            result.Add ( histo ) 
            del histo 

    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with use_canvas ( 'test_pp_callable2' , wait = 1 ) : 
        result.draw (   ) 

    return result 

# =============================================================================
if '__main__' == __name__ :

    test_pp_function  () 
    test_pp_method    () 
    test_pp_callable1 () 
    test_pp_callable2 () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
