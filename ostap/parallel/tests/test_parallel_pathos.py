#!/usr/bin/env python
# ============================================================================
## @file test_parallel_pathos.py
# Oversimplified script for parallel execution using Pathos 
# @see https://github.com/uqfoundation/pathos
# @see https://github.com/uqfoundation/dill
# ============================================================================
""" Oversimplified script for parallel execution using Pathos
"""
# ============================================================================
from   itertools                import count   
import ostap.histos.histos
from   ostap.utils.progress_bar import progress_bar 
from   ostap.parallel.utils     import pool_context
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.utils        import wait 
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
    logger.error ('Can not import pathos')
    pathos = None
    

DILL_PY3_issue = False 
if ( 3 , 6 ) <= sys.version_info and dill :
    dill_version =  getattr ( dill , '__version__' , '' )
    if not dill_version :  dill_version =  getattr ( dill , 'version' , '' )
    DILL_PY3_issue = dill_version < '0.3'
    if not DILL_PY3_issue :
        from ostap.core.meta_info import root_info
        ## DILL_PY3_issue = root_info < ( 6 , 23 )
        DILL_PY3_issue = root_info < ( 6 , 24 , 6  )

if DILL_PY3_issue : logger.warning ( "There is an issue with DILL/ROOT/PYTHON")
    
# =============================================================================
## simple function that creates and  fills a histogram
def make_histo ( item ) :
    """Simple    function that creates and  fills a histogram
    """
    i, n = item 
    import ROOT, random 
    h1 = ROOT.TH1F ( 'h%d' %  i , '' , 100 , 0 , 10 )
    for i in range ( n ) : h1.Fill ( random.gauss (  5 ,  1 ) )
    return h1 

# =============================================================================
## @class MakeHisto
#  class that creates and fill a histogram 
class MakeHisto(object)  :
    """Class that creates and fill a histogram"""
    def process  ( self , *args ) :
        return make_histo ( *args )
    def __call__ ( self , *args ) :
        return self.process ( *args ) 
    
## start 10 jobs, and for each job create the histogram with 100 entries 
inputs = 10 * [ 100 ]

# =============================================================================
## test parallel processing with pathos: ProcessPool
def test_pathos_mp_function () :
    """Test parallel processnig with pathos: ProcessPool
    """
    logger = getLogger("ostap.test_pathos_mp_function")
    if not pathos :
        logger.error ( "pathos is not available" )
        return 
    
    logger.info ('Test job submission with %s' %  pathos ) 
    
    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return
    
    from pathos.helpers import cpu_count
    ncpus = cpu_count  ()
    
    from pathos.pools import ProcessPool as Pool

    pool = Pool ( ncpus )
    logger.info ( "Pool is %s" % ( type ( pool ).__name__ ) )

    with pool_context   ( pool ) : 
        
        jobs = pool.uimap ( make_histo ,  zip ( count() , inputs ) )
        
        result = None 
        for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
            if not result  : result = h
            else           : result.Add ( h )
                
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with wait ( 1 ) , use_canvas ( 'test_pathos_mp_function' ) : 
        result.draw (   ) 

    return result 

# =============================================================================
## test parallel processing with pathos: ProcessPool
def test_pathos_mp_callable  () :
    """Test parallel processnig with pathos: ProcessPool
    """
    logger = getLogger("ostap.test_pathos_mp_callable")             
    if not pathos :
        logger.error ( "pathos is not available" )
        return 
        
    logger.info ('Test job submission with %s' %  pathos ) 

    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return
    
    from pathos.helpers import cpu_count
    ncpus = cpu_count  ()
    
    from pathos.pools import ProcessPool as Pool

    pool = Pool ( ncpus )
    logger.info ( "Pool is %s" % ( type ( pool ).__name__ ) )
    
    pool.restart ( True ) 
    
    jobs = pool.uimap ( make_histo ,  zip ( count() , inputs ) )
    
    result = None 
    for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
        if not result  : result = h
        else           : result.Add ( h )

    pool.close ()
    pool.join  ()
    pool.clear ()
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with wait ( 1 ) , use_canvas ( 'test_pathos_mp_callable' ) : 
        result.draw (   ) 

    return result 


# =============================================================================
## test parallel processing with pathos : ParallelPool 
def test_pathos_pp_function () :
    """Test parallel processnig with pathos: ParallelPool  
    """
    logger = getLogger("ostap.test_pathos_pp_function") 
    if not pathos :
        logger.error ( "pathos is not available" )
        return 
        
    logger.info ('Test job submission with %s' %  pathos ) 

    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return

    from pathos.helpers import cpu_count
    ncpus = cpu_count  ()
    
    from pathos.pools import ParallelPool as Pool 

    pool = Pool ( ncpus ,
                  secret='xxOGew', ppservers=[22602]
                  )
    logger.info ( "Pool is %s" %  ( type ( pool ).__name__ ) )

    pool.restart ( True ) 


    mh   = MakeHisto() 
    jobs = pool.uimap ( make_histo,  [  ( i , n )  for  ( i , n ) in enumerate ( inputs ) ] )
    
    result = None 
    for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
        if not result  : result = h
        else           : result.Add ( h )

    pool.close ()
    pool.join  ()
    pool.clear ()
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with wait ( 1 ) , use_canvas ( 'test_pathos_pp_function' ) : 
        result.draw (   ) 

    return result 


# =============================================================================
## test parallel processing with pathos : ParallelPool 
def test_pathos_pp_method () :
    """Test parallel processnig with pathos: ParallelPool  
    """
    logger = getLogger("ostap.test_pathos_pp_method ")
    if not pathos :
        logger.error ( "pathos is not available" )
        return 
    
    logger.info ('Test job submission with %s' %  pathos ) 
    
    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return

    from pathos.helpers import cpu_count
    ncpus = cpu_count  ()
    
    from pathos.pools import ParallelPool as Pool 

    pool = Pool ( ncpus )   
    logger.info ( "Pool is %s" %  ( type ( pool ).__name__ ) )

    pool.restart ( True ) 

    mh   = MakeHisto() 
    jobs = pool.uimap ( mh.process , [  ( i , n )  for  ( i , n ) in enumerate ( inputs ) ] )
    
    result = None 
    for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
        if not result  : result = h
        else           : result.Add ( h )

    pool.close ()
    pool.join  ()
    pool.clear ()
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with wait ( 1 ) , use_canvas ( 'test_pathos_pp_method' ) : 
        result.draw (   ) 

    return result 

# =============================================================================
## test parallel processing with pathos : ParallelPool 
def test_pathos_pp_callable () :
    """Test parallel processnig with pathos: ParallelPool  
    """
    logger = getLogger("ostap.test_pathos_pp_callable")         
    if not pathos :
        logger.error ( "pathos is not available" )
        return
    
    logger.info ('Test job submission with %s' %  pathos ) 
    
    if DILL_PY3_issue : 
        logger.warning ("test is disabled (DILL/ROOT/PY3 issue)" )
        return

    ## logger.warning ("test is disabled for UNKNOWN REASON")
    ## return

    from pathos.helpers import cpu_count
    ncpus = cpu_count  ()
    
    from pathos.pools import ParallelPool as Pool 

    pool = Pool ( ncpus )   
    logger.info ( "Pool is %s" %  ( type ( pool ).__name__ ) )

    pool.restart ( True ) 


    mh   = MakeHisto() 
    jobs = pool.uimap ( mh.process ,  [  ( i , n )  for  ( i , n ) in enumerate ( inputs ) ] )
    
    result = None 
    for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
        if not result  : result = h
        else           : result.Add ( h )

    pool.close ()
    pool.join  ()
    pool.clear ()
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    with wait ( 1 ) , use_canvas ( 'test_pathos_pp_callable' ) : 
        result.draw (   ) 

    return result 


# =============================================================================
if '__main__' == __name__ :

    test_pathos_mp_function ()
    test_pathos_mp_callable ()
    test_pathos_pp_function ()
    test_pathos_pp_method   ()
    test_pathos_pp_callable ()


# =============================================================================
##                                                                      The END 
# =============================================================================
