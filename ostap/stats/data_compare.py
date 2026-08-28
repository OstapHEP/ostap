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
    ## 
    'numpy_compare'         , ## compare two (numpy) arrays 
    'data_compare'          , ## compare two datasets
    'ecdf_compare'          , ## visual (&chi2) comparison of two (W)ECDFs 
    'compare_variable'      , ## detailed compariso of two variables 
    ## 
    'ComparisonResult'      , ## helper class to report the comparison resuts 
    'HistoComparisonResult' , ## helper class to report the histogram comparison resuts 
)
# =============================================================================
from   collections            import namedtuple 
from   ostap.core.ostap_types import sequence_types, sized_types 
from   ostap.utils.core       import typename 
from   ostap.math.math_base   import FIRST_ENTRY , LAST_ENTRY
from   ostap.stats.utils      import weight_trivial, check_all 
from   ostap.stats.gof_utils  import clip_pvalue, the_header 
from   ostap.math.math_ve     import significance
from   ostap.stats.gof        import AGoFnp
from   ostap.stats.counters   import ECDF, WECDF
from   ostap.math.math_ve     import chi2_prob, significance
from   ostap.logger.symbols   import ( chi2ndf      as chi2ndf_symbol      , 
                                       infinity_pos as infinity_pos_symbol )
import ostap.logger.table     as     T 
import ostap.histos.histos 
import ostap.histos.graphs 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.data_compare' )
else                       : logger = getLogger ( __name__                   )
# =============================================================================
## comparison result  
ComparisonResult     = namedtuple ( 'ComparisonResult'  ,
                                    ( 'method'     ,
                                      'tvalue'     ,
                                      'pvalue'     ,
                                      'nsigma'     ,
                                      'importance' ,
                                      'table_row'  ) ,
                                    defaults = ( '<UNKNOWN>' , -1e+9 , -1 , -1 , None , '' ) )
# =============================================================================
## 1D-histogram comparison result  
HistoComparisonResult = namedtuple ( 'HistoComparisonResult'  ,
                                     ( 'histo1'    ,
                                       'histo2'    ,
                                       'chi2'      ,
                                       'nDoF'      , 
                                       'pvalue'    ,
                                       'nsigma'    ) )
# =============================================================================
## ECDF types 
ecdf_types = ECDF , WECDF
# =============================================================================
## Compare two numpy arrays
#  @code
#  nd1 = ...
#  nd2 = ...
#  comparator = ...
#  result     = numpy_compare ( comparator , nd1 , nd2 , importance = True ) 
#  @endcode 
def numpy_compare ( comparator ,
                    data1      ,
                    data2      , *     ,
                    weight1    = None  ,
                    weight2    = None  ,
                    importance = False ) :
    """ Compare two numpy arrays
    >>> data1      = ...
    >>> data2      = ...
    >>> comparator = ...
    >>> result     = numpy_compare ( comparator , data1 , data2 )
    """
    ## check input data 
    check_all ( data1 , data2 , weight1 , weight2 , "numpy_compare" ) 

    ## the sequence of comparators is provided 
    if isinstance ( comparator , sequence_types ) and all ( isinstance ( c , AGoFnp ) and c.two_samples for c in comparator ) :
        return tuple ( numpy_compare ( cmp        ,
                                       data1      = data1      ,
                                       data2      = data2      ,
                                       weight1    = weight1    ,
                                       weight2    = weight2    ,
                                       importance = importance ) for cmp in comparator )
    
    ## single comparator
    if not ( isinstance ( comparator , AGoFnp ) and comparator.two_samples ) :
        raise TypeError ( "Invalid comparator type: %s" % typename ( comparator )  ) 
    
    w1_trivial = weight_trivial ( weight1 )
    w2_trivial = weight_trivial ( weight2 )
    
    if   comparator.weights_supported : pass
    elif w1_trivial and w2_trivial    : pass
    else : raise TypeError ( "Comparator %s does not support weights!" % typename ( comparator ) )
    
    importance_features = None
    tvalue              = None 
    if importance and hasattr ( comparator , 'importance_features' ) :
        tvalue = comparator.tvalue ( data1      ,
                                     data2      ,
                                     weight1    = weight1    ,
                                     weight2    = weight2    ,
                                     importance = importance )        
        importance_features = {} 
        importance_features.update ( comparator.importance_features )
        
    ## get t&p-values 
    tvalue , pvalue  = comparator.pvalue ( data1   ,
                                           data2   ,
                                           tvalue  = tvalue  ,
                                           weight1 = weight1 ,
                                           weight2 = weight2 )
    
    pv      = clip_pvalue  ( pvalue ) 
    nsigma  = significance ( pv     ) ## convert  it to significance

    return ComparisonResult ( method     = comparator.method           ,
                              tvalue     = tvalue                      ,
                              pvalue     = pvalue                      ,
                              nsigma     = nsigma                      ,
                              importance = importance_features         , 
                              table_row  = comparator.the_row () [ 1 ] ) 

# =============================================================================
## Compare two datasets:
#  Each dataset is converted to numpy and then `numpy_compare` is invoked 
#  @see numpy_compare 
#  @code
#  data1      = ...
#  data2      = ...
#  comparator = ... 
#  results , table = data_compare ( comparator , nd1 , nd2 , 'X,y,z' , 'pt>1' , importance = True )
#  @endcode 
def data_compare ( comparator   ,  
                   data         ,
                   data2        ,
                   expressions  ,                   ## variables in data1 
                   cuts         = ''          , * , ## cuts for data1
                   expressions2 = None        ,
                   cuts2        = ''          ,                            
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
                   progress2    = None        ,
                   importance   = False       ,
                   ## 
                   title        = ''          , 
                   silent       = False       , **config  ) :
    """ Compare two datasets:
    Each dataset is converted to numpy and then `numpy_compare` is invoked 
    - see numpy_compare 
    >>> data1      = ...
    >>> data2      = ...
    >>> comparator = ... 
    >>> results , table = data_compare ( comparator , data1 , data2 , 'X,y,z' , 'pt>1' )
    """
    
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
    
    assert varlst1 and len ( varlst1 ) == len ( varlst2 ) , \
        "Different lenghts for variable lists!"

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

    ## use comparator(s)  
    results = numpy_compare ( comparator , 
                              nd1        ,
                              nd2        ,
                              weight1    = weight1 ,
                              weight2    = weight2 ,
                              importance = importance )

    ##
    
    if isinstance ( comparator , sequence_types ) and isinstance ( results , sequence_types ) and \
       isinstance ( comparator , sized_types    ) and isinstance ( results , sized_types    ) and \
       all ( isinstance ( r , ComparisonResult  ) for r in results ) :

        rr = results
        cc = comparator
        
    elif isinstance ( comparator , AGoFnp ) and isinstance ( results  , ComparisonResult ) :
        
        rr = results    , 
        cc = comparator ,
        
    else : raise TypeError ( "Inconsistent type for `results` and `comparator`"  )
    
    ## summary  table 
    if not silent :
        
        rows = [] 
        for r , c in zip ( rr , cc ) :
            
            _   , row = c.the_row ()
            row       = ( c.method , ) + row
            rows.append ( row )
            
        header = ( 'Method' , ) + the_header
        rows   = [ header ] + rows            
        title  = title if title else "Two-Sample comparison"
        rows   = T.remove_empty_columns ( rows )
        table  = T.table ( rows , title = title , prefix = "# " , alignment = 'lccccccccccc' )
        logger.info ( "%s:\n%s" % ( title , table ) )
        
    return results 

# =============================================================================
## Visual&chi2 comparison of (W)ECDFs
#  - create pooled sample
#  - find equidistant quantiles
#  - create two histograms accorfing to quantiels of pooled sample
#  - project (W)ECDFs into histogram
#  - normalize them as "density" (optional)
#  - draw them superimposed (optional)
#  - decorate & return the histograms
#  - calculate chi2-probability and convert it to #sigmas 
#  @code
#  ecdf1  = ...
#  ecdf2  = ...
#  result = ecdf_compare ( ecdf1 , ecdf2 , N = 50 ) 
#  @encode
def ecdf_compare ( ecdf1   ,
                   ecdf2   , *    ,
                   N       = 40   ,
                   density = True ,
                   fill    = True , 
                   draw    = True ) :
    """ Visual compare (W)ECDFs
    - create pooled sample
    - find equaidistant quantiles
    - create two histograms accorfing to quantiels of pooled sample
    - project (W)ECDFs into histogram
    - normalize them as "density" (optional)
    - draw them superimposed      (optional)
    - decorate & return the histograms
    - calculate chi2-probability and convert it to #sigmas 

    >>> ecdf1  = ...
    >>> ecdf2  = ...
    >>> result = ecdf_compare ( ecdf1 , ecdf2 , N = 50 ) 
    """
    if not isinstance ( ecdf1 , ecdf_types ) : raise TypeError  ( "Invalid `ecdf1` type: %s" % typename ( ecdf1 ) )
    if not isinstance ( ecdf2 , ecdf_types ) : raise TypeError  ( "Invalid `ecdf2` type: %s" % typename ( ecdf2 ) )
    if not ecdf1 . ok ()                     : raise ValueError ( "Invalid `ecdf1`!" )
    if not ecdf2 . ok ()                     : raise ValueError ( "Invalid `ecdf2`!" )
    
    if   isinstance   ( ecdf1 , WECDF ) : pooled = ecdf1.merge ( ecdf2 )
    else                                : pooled = ecdf2.merge ( ecdf1 )

    quantiles = pooled.quantiles_[N-1] ()
    
    import ostap.histos.axes      as     AXES
    axis      = AXES.axis_from_edges ( quantiles )
    
    total     = pooled.high_edge() - pooled.low_edge () 
    if total :
        ## avoid too wide bins 
        if 4 <= N : axis = axis.split_while ( total /  min ( N , 4 ) )
        ## avoid too narrow bins 
        axis             = axis.merge_while ( total / ( N * 10 ) )
        
    h1        = AXES.h1_axis ( axis , title = 'First dataset'  )
    h2        = AXES.h1_axis ( axis , title = 'Second dataset' )

    ecdf1.project ( h1 )
    ecdf2.project ( h2 )

    if density :
        
        h1 = h1.density ( silent = True )
        h2 = h2.density ( silent = True )
        
    if fill :
        
        h1.red  ( fill = True , opacity = 0.20 ) 
        h2.blue ( fill = True , opacity = 0.30 ) 

    if draw :

        if h2.GetMaximum() < h1.GetMaximum() : 
            
            h1.draw () 
            h2.draw ('same')
            
        else  :
             
            h2.draw ()
            h1.draw ( 'same' )
        
        if fill :
            
            h1.draw ( 'same hist' )
            h2.draw ( 'same hist' )


    chi2 = 0
    n    = 0 
    for i , _  , y1 in h1.items ()  :
        y2 = h2 [ i ]
        cov = y1.cov2() + y2.cov2()
        if 0 < cov :
            dy    = y1.value () - y2.value ()
            chi2 += dy * dy / cov
            n    += 1
            
    ## # of degrees of freedom 
    nDoF = ( n - 1 ) if density else n
    
    if 0 <= 0 and 1 <= nDoF :        
        pvalue = chi2_prob    ( chi2 , nDoF )
        nsigma = significance ( pvalue )
    else :
        pvalue , nsigma = float ( 'NaN' ) , float ( 'NaN' )
    
    return HistoComparisonResult ( histo1  = h1      ,
                                   histo2  = h2      ,
                                   chi2    = chi2    ,
                                   nDoF    = nDoF    ,   
                                   pvalue  = pvalue  , 
                                   nsigma  = nsigma  ) 


# ============================================================================
## Compare two variables from (presumably two) sources
#  @param what        (INPUT) the first variable name/expression
#  @param data        (INPUT) the first source:  TTree, TChain, RooDataSet,...
#  @param cuts        (INPUT) selection criteria 
#  @param cut_range   (INPUT) cut-range (for RooDataSet)
#  @param first       (INPUT) the first event in data to process 
#  @param last        (INPUT) the last  event in data to process 
#  @param var2        (INPUT) the second variable name/expression
#  @param data2       (INPUT) the second source:  TTree, TChain, RooDataSet,...
#  @param cuts2       (INPUT) selection criteria 
#  @param cut_range2  (INPUT) cut-range (for RooDataSet)
#  @param first2      (INPUT) the first event in data2 to process 
#  @param last2       (INPUT) the last  event in data2 to process 
#  @param N           (INPUT) number of bins in the histogram
#  @param draw        (INPUT) draw the comarion results ?
#  @param fill        (INPUT) fill the histogram area ?
#  @param density     (INPUT) convert historrams into density representaton?
#  @param comparators (INPUT) the list of Two-Sample/GoF estimators
#  @param progress    (INPUT) show the progress bar(s) ?
#  @param parallel    (INPUT) use parallelization if/when possible?
def compare_variable ( what        ,
                       data        , *           , 
                       cuts        = ''          , 
                       cut_range   = ''          ,
                       first       = FIRST_ENTRY ,
                       last        = LAST_ENTRY  ,
                       ## 
                       what2       = None        , 
                       data2       = None        ,
                       cuts2       = ''          ,
                       cut_range2  = ''          ,
                       first2      = FIRST_ENTRY ,
                       last2       = LAST_ENTRY  ,
                       ## 
                       comparators = ()          , ## COMPARATORS
                       ## (W)ECDF comparison 
                       N           = 30          , ## number of bins
                       draw        = True        , ## visualze results ?
                       fill        = True        ,
                       density     = True        ,
                       ##
                       silent      = False       ,
                       title       = ''          ,
                       progress    = False       ,
                       parallel    = False       ) :
    
    if what2      is None and \
       data2      is None and \
       cuts2      is None and \
       cut_range2 is None and \
       first2 == first    and \
       last2  == last        :  logger.error ( "Nothing to compare!" )

    if what2      is None : what2      = what 
    if data2      is None : data2      = data  
    if cuts2      is None : cuts2      = cuts 
    if cut_range2 is None : cut_range2 = cut_range

    from   ostap.stats.statvars import data_ECDF
    
    what , what2 = what.strip() , what2.strip () 
          
    ecdf1 = data_ECDF      ( data       = data       ,
                             expression = what       ,
                             cuts       = cuts       ,
                             cut_range  = cut_range  , 
                             first      = first      ,
                             last       = last       ,
                             progress   = progress   ,
                             parallel   = parallel   ) 
    
    ecdf2 = data_ECDF      ( data       = data2      ,
                             expression = what2      ,
                             cuts       = cuts2      ,
                             cut_range  = cut_range2 , 
                             first      = first2     ,
                             last       = last2      ,
                             progress   = progress   ,
                             parallel   = parallel   )

    if hasattr ( data  , 'soft_copy' ) : data  = data .soft_copy() 
    if hasattr ( data2 , 'soft_copy' ) : data2 = data2.soft_copy() 

    ## (1) vizual comparison of distributions 
    histos = ecdf_compare ( ecdf1   = ecdf1   ,
                            ecdf2   = ecdf2   ,
                            N       = N       , 
                            density = density ,
                            fill    = fill    ,
                            draw    = draw    )

    
    ## unpack (W)ECDFs 
    raw_data1 , weight1 = ecdf1.raw_data ()
    raw_data2 , weight2 = ecdf2.raw_data ()
    
    ## run comparators
    with_weight  = not weight_trivial ( weight1 ) or not weight_trivial ( weight2 )
    
    if isinstance ( comparators , AGoFnp ) : comparators = ( comparators , )
    comparators  = tuple (  c for c in comparators if isinstance ( c , AGoFnp ) and c.two_samples )
    if with_weight :
        comparators =  tuple ( c for c in comparators if c.weights_supported )

    results = []     
    rows    = []
    
    ## explicit loop over all comparators 
    for cmp in comparators :
        
        r = numpy_compare ( comparator = cmp       ,  
                            data1      = raw_data1 ,
                            data2      = raw_data2 ,
                            weight1    = weight1   ,
                            weight2    = weight2   )
        
        results.append ( r )

        if not silent : 
            _ , row = cmp.the_row ()
            row    = ( cmp.method , ) + tuple ( row )
            rows.append ( row ) 

    results = tuple ( results ) 

    ## summary table 
    if not silent :
        
        ## chi2-row 
        chi2ndf = '%.2f/%s'  % ( histos.chi2 , histos.nDoF) 
        pvalue  = '%6.2f'    % ( histos.pvalue * 100 )
        
        nsigma  = histos.nsigma 
        if  99 < float ( nsigma ) : nsigma = infinity_pos_symbol
        else                      : nsigma = '%5.2f' % nsigma
    
        row   = chi2ndf_symbol , chi2ndf   ,  '' , '' , '' , '' , '' , pvalue , nsigma
        
        header = ( 'Method' , ) + the_header 
        rows   = [ header ] + [ row ] + rows
        
        if   what2 == what : title = title if title else "Two-Sample test: %s"           %   what 
        else               : title = title if title else "Two-Sample test: (%s) vs (%s)" % ( what , what2 ) 
        
        rows   = T.remove_empty_columns ( rows )
        table  = T.table ( rows , title = title , prefix = "# " , alignment = 'lccccccccccc' )
        logger.info ( "%s:\n%s" % ( title , table ) )

    return histos , results 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
  
