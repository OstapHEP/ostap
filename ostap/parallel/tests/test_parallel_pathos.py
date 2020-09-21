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
from   itertools                import count   
import ostap.histos.histos
from   ostap.utils.progress_bar import progress_bar 
from   ostap.parallel.utils     import pool_context

# =============================================================================
try : 
    import pathos
except ImportError :
    logger.error ('Can not import pathos')
    pathos = None
    
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
    
## start 10 jobs, and for each job create the histogram with 1000 entries 
inputs = 5 * [ 100 ]

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
            
    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        vip = '%s.%s.%s' % ( vi.major , vi.minor , vi.micro ) 
        logger.warning ("test is disabled for Python %s (dill/ROOT issue)" % vip )
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
    
    result.Draw (   ) 
    time.sleep  ( 2 )

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

    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        vip = '%s.%s.%s' % ( vi.major , vi.minor , vi.micro ) 
        logger.warning ("test is disabled for Python %s (dill/ROOT issue)" % vip )
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
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 


# =============================================================================
## test parallel processing with pathos : ParallelPool 
def test_pathos_pp_function () :
    """Test parallel processnig with pathos: ParallelPool  
    """
    logger = getLogger("test_pathos_pp_function") 
    if not pathos :
        logger.error ( "pathos is not available" )
        return 
        
    logger.info ('Test job submission with %s' %  pathos ) 

    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        vip = '%s.%s.%s' % ( vi.major , vi.minor , vi.micro ) 
        logger.warning ("test is disabled for Python %s (dill/ROOT issue)" % vip )
        return

    from pathos.helpers import cpu_count
    ncpus = cpu_count  ()
    
    from pathos.pools import ParallelPool as Pool 

    pool = Pool ( ncpus )   
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
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 


# =============================================================================
## test parallel processing with pathos : ParallelPool 
def test_pathos_pp_callable () :
    """Test parallel processnig with pathos: ParallelPool  
    """
    logger = getLogger("test_pathos_pp_callable") 
    logger.info ('Test job submission with %s' %  pathos ) 

    if not pathos :
        logger.error ( "pathos is not available" )
        return 
        
    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        vip = '%s.%s.%s' % ( vi.major , vi.minor , vi.micro ) 
        logger.warning ("test is disabled for Python %s (dill/ROOT issue)" % vip )
        return

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
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 

# =============================================================================
## test parallel processing with pathos : ParallelPool 
def test_pathos_pp_method () :
    """Test parallel processnig with pathos: ParallelPool  
    """
    logger = getLogger("test_pathos_pp_method ")
    if not pathos :
        logger.error ( "pathos is not available" )
        return 
    
    logger.info ('Test job submission with %s' %  pathos ) 
    
    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        vip = '%s.%s.%s' % ( vi.major , vi.minor , vi.micro ) 
        logger.warning ("test is disabled for Python %s (dill/ROOT issue)" % vip )
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
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 

# =============================================================================
## test parallel processing with pathos : ParallelPool 
def test_pathos_pp_callable () :
    """Test parallel processnig with pathos: ParallelPool  
    """
    logger = getLogger("test_pathos_pp_callable")         
    if not pathos :
        logger.error ( "pathos is not available" )
        return
    
    logger.info ('Test job submission with %s' %  pathos ) 
    
    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        vip = '%s.%s.%s' % ( vi.major , vi.minor , vi.micro ) 
        logger.warning ("test is disabled for Python %s (dill/ROOT issue)" % vip )
        return

    logger.warning ("test is disabled for UNKNOWN REASON")
    return

    from pathos.helpers import cpu_count
    ncpus = cpu_count  ()
    
    from pathos.pools import ParallelPool as Pool 

    pool = Pool ( ncpus )   
    logger.info ( "Pool is %s" %  ( type ( pool ).__name__ ) )

    pool.restart ( True ) 


    mh   = MakeHisto() 
    jobs = pool.uimap ( mh ,  [  ( i , n )  for  ( i , n ) in enumerate ( inputs ) ] )
    
    result = None 
    for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
        if not result  : result = h
        else           : result.Add ( h )

    pool.close ()
    pool.join  ()
    pool.clear ()
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 10 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

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
