#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/tools/data_compare.py
#  Compare 1D (weighted) disributions 
#  - convert into (W)ECDF
#  - compare distributino using (unbinned) Two-Samples (weighted) tests
#  - make comparison plots 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2021-09-22
# =============================================================================
""" Compare 1D (weighted) disributions 
- convert into (W)ECDF
- compare distributino using (unbinned) Two-Samples (weighted) tests
- make comparison plots 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-06-07"
__all__     = (
) 
# =============================================================================
from   collections              import namedtuple 
from   ostap.utils.core         import typename
from   ostap.math.math_base     import FIRST_ENTRY , LAST_ENTRY
from   ostap.stats.utils        import weight_trivial
from   ostap.stats.counters     import ECDF, WECDF
from   ostap.stats.gof          import AGoFnp
from   ostap.math.math_ve       import chi2_prob, significance
from   ostap.logger.symbols     import ( greek_lower_sigma , script_t , script_p ,   
                                         chi2ndf      as chi2ndf_symbol      , 
                                         infinity_pos as infinity_pos_symbol ) 
import ostap.histos.axes        as     AXES
import ostap.logger.table       as     T 
import ostap.histos.histos 
import ostap.histos.graphs 
import numpy
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger( 'ostap.tools.data_compare' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## histogram comparison result  
HistoComparisonResult = namedtuple ( 'HistoComparisonResult'  ,
                                     ( 'histo1'  ,
                                       'histo2'  ,
                                       'chi2'    ,
                                       'nDoF'    , 
                                       'pvalue'  ,
                                       'nsigma'  ) )

# =============================================================================
ecdf_types = ECDF , WECDF
# =============================================================================
## Visual&chi2 comparison of (W)ECDFs
#  - create pooled sample
#  - find equidistant quantiles
#  - create two histograms accorfing to quantiels of pooled sample
#  - project (W)ECDFs into histogram
#  - normalize them as "density" (optional)
#  - draw them superimposed (optional)
#  - decorate & return the histograms
#  - calculate chi2-probability an dfocnvert it to #sigmas 
#  @code
#  ecdf1  = ...
#  ecdf2  = ...
#  result = compare_ecdfs ( ecdf1 , ecdf2 , N = 50 ) 
#  @encode
def compare_ecdfs ( ecdf1   ,
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
    - calculate chi2-probability an dfocnvert it to #sigmas 

    >>> ecdf1  = ...
    >>> ecdf2  = ...
    >>> result = compare_ecdfs ( ecdf1 , ecdf2 , N = 50 ) 
    """
    if not isinstance ( ecdf1 , ecdf_types ) : raise TypeError  ( "Invalid `ecdf1` type: %s" % typename ( ecdf1 ) )
    if not isinstance ( ecdf2 , ecdf_types ) : raise TypeError  ( "Invalid `ecdf2` type: %s" % typename ( ecdf2 ) )
    if not ecdf1 . ok ()                     : raise ValueError ( "Invalid `ecdf1`!" )
    if not ecdf2 . ok ()                     : raise ValueError ( "Invalid `ecdf2`!" )
    
    if   isinstance   ( ecdf1 , WECDF ) : pooled = ecdf1.merge ( ecdf2 )
    else                                : pooled = ecdf2.merge ( ecdf1 )

    quantiles = pooled.quantiles_[N-1] ()
    
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
        
        h1 = h1.density()
        h2 = h2.density()
        
    if fill :
        
        h1.red  ( fill = True , opacity = 0.25 ) 
        h2.blue ( fill = True , opacity = 0.25 ) 

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
        pvalue , nsigma = float( 'NaN' ) , float( 'NaN')
    
    return HistoComparisonResult ( histo1  = h1      ,
                                   histo2  = h2      ,
                                   chi2    = chi2    ,
                                   nDoF    = nDoF    ,   
                                   pvalue  = pvalue  , 
                                   nsigma  = nsigma  ) 

# ============================================================================
## Compare two variables from (presumably two) sources
#  @param var         (INPUT) the first variable name/expression
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
def compare_variables ( var         ,
                        data        , *           , 
                        cuts        = ''          , 
                        cut_range   = ''          ,
                        var2        = None        , 
                        data2       = None        ,
                        cuts2       = None        ,
                        cut_range2  = None        , 
                        first       = FIRST_ENTRY ,
                        last        = LAST_ENTRY  ,
                        first2      = FIRST_ENTRY ,
                        last2       = LAST_ENTRY  ,
                        comparators = ()          , ## COMPARATORS                         
                        N           = 20          , ## number of bins
                        draw        = True        , ## visualze results ?
                        fill        = True        ,
                        density     = True        ,
                        progress    = False       ,
                        parallel    = False       ) :
    
    if var2       is None and \
       data2      is None and \
       cuts2      is None and \
       cut_range2 is None and \
       first2 == first    and \
       last2  == last        :  logger.error ( "Nothing to compare!" )

    if var2       is None : var2       = var
    if data2      is None : data2      = data  
    if cuts2      is None : cuts2      = cuts 
    if cut_range2 is None : cut_range2 = cut_range

    from   ostap.stats.statvars import data_ECDF
    
    var , var2 = var.strip() , var2.strip() 
          
    ecdf1 = data_ECDF      ( data       = data       ,
                             expression = var        ,
                             cuts       = cuts       ,
                             cut_range  = cut_range  , 
                             first      = first      ,
                             last       = last       ,
                             progress   = progress   ,
                             parallel   = parallel   ) 
    
    ecdf2 = data_ECDF      ( data       = data2      ,
                             expression = var2       ,
                             cuts       = cuts2      ,
                             cut_range  = cut_range2 , 
                             first      = first2     ,
                             last       = last2      ,
                             progress   = progress   ,
                             parallel   = parallel   )
     
    ## (1) vizual comparison of distributions 
    histos = compare_ecdfs ( ecdf1   = ecdf1   ,
                             ecdf2   = ecdf2   ,
                             N       = N       , 
                             density = density ,
                             fill    = fill    ,
                             draw    = draw    )
    
    
    ## unpack (W)ECDFs 
    data1 , weight1 = ecdf1.raw_data()
    data2 , weight2 = ecdf2.raw_data()
    
    ## run comparators
    with_weight  = not weight_trivial ( weight1 ) or not weight_trivial ( weight2 )
    
    if isinstance ( comparators , AGoFnp ) : comparatros = ( comparators , )
    comparators  = tuple (  c for c in comparators if isinstance ( c , AGoFnp ) and c.two_samples )
    if with_weight : comparators =  tuple ( c for c in comparators if c.weights_supported )

    results = [] 
    header  = ( 'Method' , '%s-value' % script_t , '' , '' , '' , '' , '' , '%s-value [%%]' % script_p , '#%s' % greek_lower_sigma ) 
    
    chi2ndf = '%.2f/%s'  % ( histos.chi2 , histos.nDoF) 
    pvalue  = '%6.2f'    % ( histos.pvalue * 100 )  
    nsigma  = ( '%6.2f'  %   histos.nsigma ) if histos.nsigma < 100 else infinity_pos_symbol 
    row     = chi2ndf_symbol , chi2ndf   ,  '' , '' , '' , '' , '' , pvalue , nsigma  
            
    rows = [ row ] 
    
    from   ostap.stats.data_compare import numpy_compare        
    for c in comparators :
            
        r = numpy_compare ( c       ,
                            data1   = data1   ,
                            data2   = data2   ,
                            weight1 = weight1 ,
                            weight2 = weight2 )
        results.append ( r )
            
        header  , row = c.the_row ()
        header = ( 'Method' , ) + header 
        row    = ( c.method , ) + tuple ( row )

        rows.append ( row ) 

    rows  =  [ header ] + rows
     
    if   var2 == var : title = "Two-Sample test: %s"           %   var 
    else             : title = "Two-Sample test: (%s) vs (%s)" % ( var , var2 ) 
    
    rows  = T.remove_empty_columns ( rows ) 
    table = T.table ( rows , title = title , prefix = '# ' , alignment = 'lccccccccccc' )
    logger.info ( "%s:\n%s" % ( title , table ) )
        
            
    return histos , tuple ( results ) 

# ============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
        


# =============================================================================
##                                                                     The END 
# =============================================================================
