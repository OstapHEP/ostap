#!/usr/bin/env python
# ============================================================================
## @file test_parallel_joblib.py
# Oversimplified script for parallel execution using joblib  
# ============================================================================
""" Oversimplified script for parallel execution using joblib
"""
# =============================================================================
from   ostap.utils.progress_bar import progress_bar 
from   ostap.plotting.canvas    import use_canvas
from   ostap.parallel.utils     import uimap, fix_ppsrv 
from   ostap.utils.root_utils   import batch_env 
import ostap.histos.histos
import joblib 
import ROOT, random, time, sys
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__  or '__builtin__' == __name__ : 
    logger = getLogger ( 'test_parallel_joblib' )
else : 
    logger = getLogger ( __name__ )
# =============================================================================
batch_env ( logger )
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import joblib # ===========================================================
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    joblib = None # ===========================================================
    logger.warning ( "No joblib module is available!" ) 

# =============================================================================    
## simple    function that creates and fills a histogram 
def make_histo  ( item ) :
    """ Simple    function that creates and  fills a histogram
    """
    i , n = item 
    import ROOT
    import random
    from ostap.core.core import hID 
    h1 = ROOT.TH1F ( hID() , '' , 100 , 0 , 10 )
    for i in range ( n ) : h1.Fill ( random.gauss (  4 ,  1 ) )
    return h1 

# ===================================================================================
## @class MakeHisto
#  helper class to create a fill histograms
class MakeHisto(object) :
    """ Helper class to create a fill histoghrams
    """
    def process  ( self ,   *args ) :
        return make_histo ( *args )
    def __call__ ( self ,   *args ) :
        return make_histo ( *args )

mh  = MakeHisto  ()

## start NN jobs, and for each job create the histogram with 100 entries
NN     = 10 
inputs = NN * [ 100 ]

def test_func ( *args , **kwargs ) :
    print ( 'ARGS' , args , kwargs )

# =============================================================================
config  = {}
if   joblib and '1.4.0' <= joblib.__version__ : config [ 'return_as' ] = 'generator_unordered'
elif joblib and '1.3.0' <= joblib.__version__ : config [ 'return_as' ] = 'generator'

# =============================================================================
## test joblib python with with plain function 
def test_joblib_function () :
    """ Test parallel python with plain function
    """
    logger =    getLogger ("test_joblib_function")
    logger.info ('Test job submission with %s' % joblib  ) 

    if not joblib :
        logger.error ( "No joblib module is available" )
        return 
        
    
    result = None

    with joblib.Parallel ( **config )  as executor:
        
        result = None 
        for histo in executor ( joblib.delayed ( make_histo )( ( i , a ) ) for i,a in enumerate ( inputs ) ) : 
            if result is None : result  = histo 
            else              : result += histo 
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with use_canvas ( 'test_joblib_function' , wait = 1 ) : 
        result.draw (   ) 

    return result 
                  
# =============================================================================
## test parallel python with object method  
def test_joblib_method() :
    """ Test parallel python with object method  
    """
    logger =    getLogger ("test_joblib_method")
    logger.info ('Test job submission with %s' % joblib ) 
    if not joblib :
        logger.error ( "No joblib module is available" )
        return 

    
    with joblib.Parallel ( **config )  as executor:
        
        result = None 
        for histo in executor ( joblib.delayed ( mh.process )( ( i , a ) ) for i,a in enumerate ( inputs ) ) : 
            if result is None : result  = histo 
            else              : result += histo 
    

    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with use_canvas ( 'test_joblib_method' , wait = 1 ) : 
        result.draw (   ) 

    return result 
                      
# =============================================================================
## test parallel python with callable object
def test_joblib_callable1 () :
    """ Test parallel python with callable object 
    """
    logger =    getLogger ("test_jobilb_callable1")
    logger.info ('Test job submission with %s' %  joblib  ) 
    if not joblib :
        logger.error ( "No joblib module is available" )
        return 
    
    with joblib.Parallel ( **config )  as executor:
        
        result = None 
        for histo in executor ( joblib.delayed ( mh )( ( i , a ) ) for i,a in enumerate ( inputs ) ) : 
            if result is None : result  = histo 
            else              : result += histo 
    
    with use_canvas ( 'test_joblib_callable1' , wait = 1 ) : 
        result.draw (   ) 

    return result 

# =============================================================================
## test parallel python with callable object
def test_joblib_callable2 () :
    """ Test parallel python with callable object 
    """
    logger =    getLogger ("test_joblib_callable2")
    logger.info ('Test job submission with %s' % joblib  ) 
    if not joblib :
        logger.error ( "No joblib module is available" )
        return 
    
    with joblib.Parallel ( **config )  as executor:
        
        result = None 
        for histo in executor ( joblib.delayed ( mh.__call__ )( ( i , a ) ) for i,a in enumerate ( inputs ) ) : 
            if result is None : result  = histo 
            else              : result += histo 

    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with use_canvas ( 'test_joblib_callable2' , wait = 1 ) : 
        result.draw (   ) 

    return result 

# =============================================================================
if '__main__' == __name__ :
    
    test_joblib_function  () 
    test_joblib_method    () 
    test_joblib_callable1 () 
    test_joblib_callable2 () 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
