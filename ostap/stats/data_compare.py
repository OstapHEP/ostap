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
    'numpy_compare'         , ## Compare two (numpy) arrays using Two-Samples comparators 
    'data_compare'          , ## Compare two datasets using Two-Samples comparators 
    'ecdf_compare'          , ## Cisual (&chi2) comparison of two (W)ECDFs 
    'compare_variable'      , ## detailed comparision of two variables
    'compare_projections'   , ## Use Random Axed/Cramer-Wold Two-Sample Test to compare datasets 
    ## 
    'ComparisonResult'      , ## helper class to report the comparison resuts 
    'HistoComparisonResult' , ## helper class to report the histogram comparison resuts 
)
# =============================================================================
from   ostap.core.ostap_types import sequence_types, sized_types, dictlike_types, integer_types   
from   ostap.utils.core       import typename 
from   ostap.math.math_base   import FIRST_ENTRY , LAST_ENTRY, axis_range 
from   ostap.stats.utils      import ( TwoDataSets            ,
                                       ComparisonResult       ,
                                       HistoComparisonResult  ,                                       
                                       prepare_data           , 
                                       weight_trivial         ,
                                       check_all              , 
                                       full_header            , 
                                       format_pvalue          ,
                                       format_nsigma          )  
from   ostap.stats.gof_utils  import clip_pvalue 
from   ostap.stats.gof        import AGoFnp
from   ostap.stats.counters   import ECDF, WECDF
from   ostap.math.math_ve     import chi2_prob, significance
from   ostap.logger.symbols   import chi2ndf  as chi2ndf_symbol  
import ostap.logger.table     as     T 
import ostap.histos.histos 
import ostap.histos.graphs 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.data_compare' )
else                       : logger = getLogger ( __name__                   )
# =============================================================================##
## ECDF types 
ecdf_types = ECDF , WECDF
# =============================================================================
## Compare two datasets:
#  Each dataset is converted to numpy and then `numpy_compare` is invoked 
#  @see numpy_compare 
#  @code
#  data1       = ...
#  data2       = ...
#  comparatorsа = ... 
#  results , table = data_compare ( comparators , nd1 , nd2 , 'X,y,z' , 'pt>1' , importance = True )
#  @endcode 
def data_compare ( comparators  ,  
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
                   prescale     = None        , 
                   prescale2    = None        , 
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
    >>> results , table = data_compare( comparator , data1 , data2 , 'X,y,z' , 'pt>1' )
    """
    ## convert data into numpy arrays 
    datasets = prepare_data ( data        = data          ,
                              data2       = data2         ,
                              ## 
                              expressions  = expressions  ,
                              cuts         = cuts         , 
                              first        = first        ,
                              last         = last         ,
                              cut_range    = cut_range    ,
                              prescale     = prescale     ,
                              use_frame    = use_frame    ,
                              parallel     = parallel     ,
                              progress     = progress     ,
                              ## 
                              expressions2 = expressions2 ,
                              cuts2        = cuts2        , 
                              first2       = first2       ,
                              last2        = last2        ,
                              cut_range2   = cut_range2   ,
                              prescale2    = prescale2    ,
                              use_frame2   = use_frame2   ,
                              parallel2    = parallel2    ,
                              progress2    = progress     )
    
    
    ## use comparator(s)  
    results = numpy_compare ( comparators                   , 
                              data1      = datasets.data1   ,  
                              data2      = datasets.data2   , 
                              weight1    = datasets.weight1 ,
                              weight2    = datasets.weight2 ,
                              importance = importance       )

    ##    
    if isinstance ( comparators , sequence_types ) and isinstance ( results , sequence_types ) and \
       isinstance ( comparators , sized_types    ) and isinstance ( results , sized_types    ) and \
       all ( isinstance ( r , ComparisonResult  ) for r in results ) :

        rr = results
        cc = comparators
        
    elif isinstance ( comparators , AGoFnp ) and isinstance ( results  , ComparisonResult ) :
        
        rr = results     , 
        cc = comparators ,
        
    else : raise TypeError ( "Inconsistent type for `results` and `comparators`"  )


    ## from   ostap.trees.cuts     import vars_and_cuts
    ## varlst1 , cuts  , _  = vars_and_cuts ( expressions  , cuts  ) 
    ## varlst2 , cuts2 , _  = vars_and_cuts ( expressions2 , cuts2 )
    ## varlst1 = tuple ( sorted ( varlst1 ) )    
    ## varlst2 = tuple ( sorted ( varlst2 ) )
    
    # =====================================================================================
    if False :
    ## for r in rr : # =======================================================================
        # =================================================================================
        if r.importance and isinstance ( r.importance , dictlike_types ) :            
            importance     = r.importance
            new_importance = {}
            for key, value in importance.items () :
                ## add it into new dictionary: 
                new_importance [ key ] = float ( value )
                if key.startswith ( 'f' ) :
                    # =====================================================================
                    try : # ===============================================================
                        # =================================================================
                        index = int ( key [ 1: ] )
                        if 0 <= index < len ( varlst1 ) :
                            name = varlst1 [ index ]
                            ## add proper name into dictionary 
                            new_importance [ name ] = float ( value )
                            ## remove old name 
                            del new_importance [ key ]
                        # =================================================================
                    except : # ============================================================
                        # =================================================================
                        pass
            ## replace 
            r.importance.clear()
            r.importance.update ( new_importance )
                            
    ## summary  table 
    if not silent :
        
        rows = [] 
        for r , c in zip ( rr , cc ) :
            
            _   , row = c.the_row ()
            ## print ( 'HERE:' , c.symbol, c.method , typename ( c ) , row , r ) 

            row       = ( c.symbol , c.method , ) + row
            rows.append ( row )
        header = ( '' , 'Method' , ) + full_header
        rows   = [ header ] + sorted ( rows ) 
        title  = title if title else "Two-Sample comparison"
        rows   = T.remove_empty_columns ( rows )
        table  = T.table ( rows , title = title , prefix = "# " , alignment = 'llccccccccccc' )
        logger.info ( "%s:\n%s" % ( title , table ) )
        
    return results 

# =============================================================================
## Compare two numpy datasets 
#  @code
#  data1 = ...
#  data2 = ...
#  comparator = ...
#  result     = numpy_compare ( comparator , data1 , data2 , importance = True ) 
#  @endcode 
def numpy_compare ( comparators ,
                    data1       ,
                    data2       , *     ,
                    weight1     = None  ,
                    weight2     = None  ,
                    importance  = False ) :
    """ Compare two numpy arrays
    >>> data1      = ...
    >>> data2      = ...
    >>> comparator = ...
    >>> result     = numpy_compare ( comparator , data1 , data2 )
    """
    ## check input data 
    check_all ( data1 , data2 , weight1 , weight2 , "numpy_compare" ) 

    ## the sequence of comparators is provided 
    if isinstance ( comparators , sequence_types ) and all ( isinstance ( c , AGoFnp ) and c.two_samples for c in comparators ) :
        return tuple ( numpy_compare ( cmp        ,
                                       data1      = data1      ,
                                       data2      = data2      ,
                                       weight1    = weight1    ,
                                       weight2    = weight2    ,
                                       importance = importance ) for cmp in comparators )
    ## single comparator
    comparator = comparators 
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

    return ComparisonResult ( symbol     = comparator.symbol           ,
                              method     = comparator.method           ,
                              tvalue     = tvalue                      ,
                              pvalue     = pvalue                      ,
                              nsigma     = nsigma                      ,
                              importance = importance_features         , 
                              table_row  = comparator.the_row () [ 1 ] )
                              
# =============================================================================
## Visual&chi2 comparison of (W)ECDFs
#  - create pooled sample  (with equal treatment of input ECDFs)
#  - find equidistant quantiles
#  - create two histograms according to quantiles of pooled sample
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
    - create pooled sample (with equal treatment of input ECDFs)
    - find equaidistant quantiles
    - create two histograms according to quantiles of pooled sample
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
    if not isinstance ( N , integer_types  ) : raise TypeError  ( "Invalid `N' type: %s"     % typename ( N     ) )
    if not 2 <= N                            : raise ValueError ( "Invalid `N'!"     )
    
    N = int ( N )
    
    pooled = WECDF ( False )
    sumw1  = ecdf1.sumw ()
    sumw2  = ecdf2.sumw ()
    
    pooled.add ( ecdf1 )
    pooled.add ( ecdf2 , sumw1 / sumw2 ) ## normalize to the same sum 
    
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
        
    if not fill :
        
        h1.red  () 
        h2.blue ()
         
    else :
         
        h1.red  ( fill = True , opacity = 0.20 ) 
        h2.blue ( fill = True , opacity = 0.30 ) 
        
    if draw :

        min1 , max1 = h1.minmax ( errors = True )
        min2 , max2 = h2.minmax ( errors = True ) 
        minv , maxv = axis_range ( min ( min1, min2 ) , max ( max1, max2 ) , delta = 0.1 )
        
        for h in ( h1, h2 ) : 
            h.SetMinimum ( minv )
            h.SetMaximum ( maxv )
        
        h1.draw () 
        h2.draw ('same')
            
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
    
    if 0 <= chi2 and 1 <= nDoF :        
        pvalue = chi2_prob    ( chi2   , nDoF )
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

    
    from   ostap.trees.cuts     import vars_and_cuts
    varlst , cuts  , _  = vars_and_cuts ( what , cuts )
    if 1 != len ( varlst ) : raise ValueError  ( "Invalid `what`: %s" % what )
    what = varlst [ 0 ]
    
    if not what2 and \
       data2      is None and \
       cuts2      is None and \
       cut_range2 is None and \
       first2 == first    and \
       last2  == last        :  logger.error ( "Nothing to compare!" )

    if not what2          : what2      = what 
    if data2      is None : data2      = data  
    if cuts2      is None : cuts2      = cuts 
    if cut_range2 is None : cut_range2 = cut_range
                    
    varlst , cuts2 , _  = vars_and_cuts ( what2 , cuts2 )
    if 1 != len ( varlst ) : raise ValueError  ( "Invalid `what2`: %s" % what2 )
    what2 = varlst [ 0 ]

    from   ostap.stats.statvars import data_ECDF
          
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
    if with_weight : comparators =  tuple ( c for c in comparators if c.weights_supported )

    results = []     
    rows    = []
    
    ## explicit loop over all comparators 
    for cmp in comparators :
        
        r = numpy_compare ( comparators = cmp       ,  
                            data1       = raw_data1 ,
                            data2       = raw_data2 ,
                            weight1     = weight1   ,
                            weight2     = weight2   )
        
        results.append ( r )

        if not silent : 
            _ , row = cmp.the_row ()
            row    = ( cmp.symbol , cmp.method , ) + tuple ( row )
            rows.append ( row ) 

    rows    = sorted ( rows    ) 
    results = tuple  ( results ) 

    ## summary table 
    if not silent :
        
        ## chi2-row 
        chi2ndf = '%.2f/%s'    % ( histos.chi2   , histos.nDoF) 
        
        pvalue  = format_pvalue  ( histos.pvalue )
        nsigma  = format_nsigma  ( histos.nsigma )
        
        row     = chi2ndf_symbol , chi2ndf_symbol , chi2ndf  ,  '' , '' , '' , '' , '' , pvalue , nsigma
        
        header  = ( '', 'Method' , ) + full_header 
        rows    = [ header ] + [ row ] + rows 
        
        if   what2 == what : title = title if title else "Two-Sample Test: %s"               %   what 
        else               : title = title if title else "Two-Sample Test: [ %s ] vs [ %s ]" % ( what , what2 ) 
        
        rows   = T.remove_empty_columns ( rows )
        table  = T.table ( rows , title = title , prefix = "# " , alignment = 'llcccccccccccc' )
        logger.info ( "%s:\n%s" % ( title , table ) )

        return histos , results

# =============================================================================
## Use Random Axes/Cramer-Wold Two-Sample Test to compare two datasets
## Compare two datasets:
#  Each dataset is converted to numpy and then `numpy_compare` is invoked 
#  @see numpy_compare 
#  @code
#  data1       = ...
#  data2       = ...
#  comparators = ... 
#  results , table = data_compare ( comparators , nd1 , nd2 , 'X,y,z' , 'pt>1' , importance = True )
#  @endcode 
def compare_projections ( comparator   ,  
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
                          prescale     = None        , 
                          prescale2    = None        , 
                          use_frame    = False       , 
                          use_frame2   = None        , 
                          parallel     = False       , 
                          parallel2    = None        , 
                          progress     = False       ,
                          progress2    = None        ,
                          ## 
                          title        = ''          , 
                          silent       = False       , **config  ) :
    """ Compare two datasets:
    Each dataset is converted to numpy and then `numpy_compare` is invoked 
    - see numpy_compare 
    >>> data1      = ...
    >>> data2      = ...
    >>> comparator = ... 
    >>> results , table = data_compare( comparator , data1 , data2 , 'X,y,z' , 'pt>1' )
    """
    import ostap.stats.random_axes as RA
    if not isinstance ( comparator , RA.RandomAxes ) :
        raise TypeError ( "Invalid comparator type : %s" % typename ( comparator ) )
    
    ## convert data into numpy arrays 
    datasets = prepare_data ( data        = data          ,
                              data2       = data2         ,
                              ## 
                              expressions  = expressions  ,
                              cuts         = cuts         , 
                              first        = first        ,
                              last         = last         ,
                              cut_range    = cut_range    ,
                              prescale     = prescale     ,
                              use_frame    = use_frame    ,
                              parallel     = parallel     ,
                              progress     = progress     ,
                              ## 
                              expressions2 = expressions2 ,
                              cuts2        = cuts2        , 
                              first2       = first2       ,
                              last2        = last2        ,
                              cut_range2   = cut_range2   ,
                              prescale2    = prescale2    ,
                              use_frame2   = use_frame2   ,
                              parallel2    = parallel2    ,
                              progress2    = progress     )
    
    ## use comparator
    _ , pvalue = comparator.pvalue ( data1   = datasets.data1   ,
                                     data2   = datasets.data2   ,
                                     weight1 = datasets.weight1 ,
                                     weight2 = datasets.weight2 )
    
    if not silent :
        title = title if title else '%s Random Axes/%s ' % ( comparator.symbol, comparator.method ) 
        logger.info ( '%s:\n%s' % ( title , comparator.report ( pvalue = pvalue , title = title , prefix = '# ' ) ) ) 

    ## 
    pv     = clip_pvalue  ( pvalue )
    nsigma = significance ( pv     )
    ## 
    return ComparisonResult ( symbol    = comparator.symbol   ,
                              method    = comparator.method    ,
                              tvalue    = None                 ,
                              pvalue    = pvalue               ,
                              nsigma    = nsigma               ,
                              importace = {}                   ,
                              table_row = comparator.table_row ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
  
