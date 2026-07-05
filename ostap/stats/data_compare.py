#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file data_compare.py 
#  Function to compare datasets
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2026-07-05  
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2026-07-05"
__all__     = (
    'numpy_compare' , ## compare two numpy arrays 
    'data_compare'  , ## compare two data sets 
)
# =============================================================================
from   ostap.math.math_base import FIRST_ENTRY , LAST_ENTRY 
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.data_compare' )
else                       : logger = getLogger ( __name__                   )
# =============================================================================

# =============================================================================
## Compare two numpy arrays
#  @code
#  nd1 = ...
#  nd2 = ...
#  comparator = ...
#  ## get t-value and p-value:
#  t_value , p_value = numpy_compare ( comparator , nd1 , nd2 ) 
#  @endcode 
def numpy_compare ( comparator ,
                    data1      ,
                    data2      , * , weight1 = None , weigth2 = None ) :
    """ Compare two numpy arrays
    >>> nd1 = ...
    >>> nd2 = ...
    >>> t_value , p_value = numpy_compare ( nd1 , nd2 )
    """
    ## retrun t and p-values  
    return comparator.pvalue ( data1   ,
                               data2   ,
                               weight1 = weight1 ,
                               weight2 = weight2 ) 

# =============================================================================
## Compare two datasets:
#  Each dataset is converted to numpy and then `numpy_compare` is invoked 
#  @see numpy_compare 
#  @code
#  data1 = ...
#  data2 = ...
#  ## get t-value and p-value:
#  t_value , p_value = data_compare ( nd1 , nd2 , 'X,y,z' , 'pt>1' )
#  @endcode 
def data_compare ( data         ,
                   data2        ,
                   comparator   , 
                   expressions  ,     ## variables in data1 
                   cuts         = ''   , * , ## cust for data1
                   expressions2 = None ,
                   cuts2        = None ,                            
                   first        = FIRST_ENTRY ,
                   last         =  LAST_ENTRY ,                                                                      
                   first2       = FIRST_ENTRY ,
                   last2        =  LAST_ENTRY ,                                                                      
                   cut_range    = ''          ,
                   cut_range2   = None        , 
                   use_frame    = False       , 
                   use_frame2   = None        , 
                   parallel     = False       , 
                   parallel2    = None        , 
                   progress     = False       ,
                   progress2    = None        , **config ) :
                
    if expressions2 is None : expressions2 = expressions
    if cuts2        is None : cuts2        = cuts
    
    if first2       is None : first2       = first
    if last2        is None : last2        = last
    
    if cut_range2   is None : cut_range2   = cut_range
    if use_frame2   is None : use_frame2   = use_frame
    if parallel2    is None : parallel2    = parallel
    if progress2    is None : progress2    = progress
    
    from   ostap.trees.cuts     import vars_and_cuts
    from   ostap.stats.statvars import data_slice
    
    varlst1 , cuts  , _  = vars_and_cuts ( expressions  , cuts )
    varlst1 = tuple ( sorted ( varlst1 ) )
    
    varlst2 , cuts2 , _  = vars_and_cuts ( expressions2 , cuts2 )
    varlst2 = tuple ( sorted ( varlst2 ) )
    
    assert varlst1 and len ( varlst1 ) == len ( varlst2 ) , "Different lengths for variable lists!"
    
    ## (1) create numpy datasets 
    nd1 , weight1 = data_slice ( data                   ,
                                 varlst1                , 
                                 cuts       = cuts      ,
                                 first      = first     ,
                                 last       = last      ,
                                 cut_range  = cut_range ,
                                 structured = True      , ## ATTENTION!  
                                 transpose  = True      ,
                                 progress   = progress  ,
                                 use_frame  = use_frame , 
                                 parallel   = parallel  )
    
    nd2 , weight2 = data_slice ( data2                   ,
                                 varlst2                 , 
                                 cuts       = cuts2      ,
                                 first      = first2     ,
                                 last       = last2      ,
                                 cut_range  = cut_range2 ,
                                 structured = True       , ## ATTENTION!  
                                 transpose  = True       ,
                                 progress   = progress2  ,
                                 use_frame  = use_frame2 , 
                                 parallel   = parallel2  )

    return numpy_compare ( comparator   ,
                           nd1          ,
                           nd2          ,
                           weigth1 = w1 ,
                           weigth2 = w2 )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
  
