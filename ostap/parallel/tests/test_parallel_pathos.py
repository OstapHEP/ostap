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
import ostap.histos.histos
from   ostap.utils.progress_bar import progress_bar 
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
    logger.error('Can not import pathos')
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
    def process (  self , *args ) :
        return make_histo ( *args )

## start 10 jobs, and for each job create the histogram with 1000 entries 
inputs = 10 * [ 100 ]

# =============================================================================
## test parallel processing with pathos: ProcessPool
def test_pathos () :
    """Test parallel processnig with pathos: ProcessPool
    """
    if not dill :
        logger.error ( "test_pathos: dill is not available" )
        return
        
    if not pathos :
        logger.error ( "test_pathos: pathos is not available" )
        return 
        
    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        logger.warning ("test_pathos is disabled for Python %s" % vi )
        return

    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        logger.warning ("test_pathos is disabled for Python %s" % vi )
        return
    
    from ostap.core.meta_info import root_version_int  as rv 
    if rv <= 62200 :
        logger.warning ("test_pathos is disabled for ROOT %s" % rv )
        return
        
    from pathos.helpers import cpu_count
    ncpus = cpu_count  ()
    
    from pathos.pools import ProcessPool as Pool

    pool = Pool ( ncpus )
    logger.info ( "test_pathos   : Pool is %s" % ( type ( pool ).__name__ ) )
    
    pool.restart ( True ) 

    # jobs = pool. imap ( make_histo ,  [  ( i , n )  for  ( i , n ) in enumerate ( inputs ) ] )
    # jobs = pool.uimap ( make_histo ,  [  ( i , n )  for  ( i , n ) in enumerate ( inputs ) ] )
    mh   = MakeHisto() 
    jobs = pool.uimap ( mh.process ,  [  ( i , n )  for  ( i , n ) in enumerate ( inputs ) ] )
    
    result = None 
    for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
        if not result  : result = h
        else           : result.Add ( h )

    pool.close ()
    pool.join  ()
    pool.clear ()
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 

# =============================================================================
## test parallel processing with pathos : ParallelPool 
def test_pathos_pp_1 () :
    """Test parallel processnig with pathos: ParallelPool  
    """
    if not dill :
        logger.error ( "test_pathos_pp_1: dill is not available" )
        return
        
    if not pathos :
        logger.error ( "test_pathos_pp_1: pathos is not available" )
        return 
        
    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        logger.warning ("test_pathos_pp_1 is disabled for Python %s" % vi )
        return

    from pathos.helpers import cpu_count
    ncpus = cpu_count  ()
    
    from pathos.pools import ParallelPool as Pool 

    pool = Pool ( ncpus )   
    logger.info ( "test_pathos_pp_1: Pool is %s" %  ( type ( pool ).__name__ ) )

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
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 

# =============================================================================
## test parallel processing with pathos : ParallelPool 
def test_pathos_pp_2 () :
    """Test parallel processnig with pathos: ParallelPool  
    """
    if not dill :
        logger.error ( "test_pathos_pp_2: dill is not available" )
        return
        
    if not pathos :
        logger.error ( "test_pathos_pp_2: pathos is not available" )
        return 
        
    vi = sys.version_info
    if 3<= vi.major and 6 <= vi.minor :
        logger.warning ("test_pathos_pp_2 is disabled for Python %s" % vi )
        return

    from pathos.helpers import cpu_count
    ncpus = cpu_count  ()
    
    from pathos.pools import ParallelPool as Pool 

    pool = Pool ( ncpus )   
    logger.info ( "test_pathos_pp_1: Pool is %s" %  ( type ( pool ).__name__ ) )

    pool.restart ( True ) 


    mh   = MakeHisto() 
    jobs = pool.uimap ( mh.process  ,  [  ( i , n )  for  ( i , n ) in enumerate ( inputs ) ] )
    
    result = None 
    for h in progress_bar ( jobs , max_value = len ( inputs ) ) :
        if not result  : result = h
        else           : result.Add ( h )

    pool.close ()
    pool.join  ()
    pool.clear ()
    
    logger.info ( "Histogram is %s" % result.dump ( 80 , 20 )  )
    logger.info ( "Entries  %s/%s" % ( result.GetEntries() , sum ( inputs ) ) ) 
    
    result.Draw (   ) 
    time.sleep  ( 2 )

    return result 

# =============================================================================
if '__main__' == __name__ :

    test_pathos      () 
    test_pathos_pp_1 () 
    test_pathos_pp_2 () 

        
# =============================================================================
##                                                                      The END 
# =============================================================================
