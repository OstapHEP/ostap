#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file sweet_viz.py 
#  Functino to analysse and compare datasets using sweetviz 1-liners  
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2026-07-05  
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2026-07-05"
__all__     = (
    'data_sweetviz_analyze'  , ## Analyze single dataset 
    'data_sweetviz_compare'  , ## Compare two datasets
    )
# =============================================================================
from   ostap.math.math_base import FIRST_ENTRY , LAST_ENTRY 
from   ostap.utils.cleanup  import CleanUp
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.stats.sweet_viz' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
## import sweetvis & pandas  
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import numpy    as np 
    import pandas   as pd
    import sweetviz as sv
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    np , pd , sv = None, None

# =============================================================================
## Use sweetviz to analyse single dataset 
#  @code
#  data = ....
#  data_sweetviz_analyse ( data , 'x,y,z' , 'pt>1' )
#  @endcode 
#  @see sweetviz
#  @see sweetviz.analyze
#  @see https://github.com/fbdesignpro/sweetviz
#  @see https://pypi.org/project/sweetviz/
def data_sweetviz_analyze ( data        ,
                            expressions ,
                            cuts        = ''   , *    ,
                            first       = FIRST_ENTRY ,
                            last        =  LAST_ENTRY ,                                                                      
                            cut_range   = ''          ,
                            use_frame   = False       , 
                            parallel    = False       , 
                            progress    = False       ,
                            anaconf     = {}          ,    ## parameters for sweetviz.analyze 
                            showconf    = {}          ) :  ## parameters for sweetviz.show_html
    
    """ Use sweetviz to analyse some statistics of data
    >>> data = ....
    >>> data_sweetviz_analyse ( data , 'x,y,z' , 'pt>1' )
    - see sweetviz
    - see https://github.com/fbdesignpro/sweetviz
    - see https://pypi.org/project/sweetviz/
    """
    
    ## check pre-requisites 
    if not np or not pd or not sv :
        logger.error ( "Numpy,pandas or sweetviz are not available, no action" )
        return

    if isinstance ( data , ROOT.RooAbsData ) and data.isWeighted() :
        import ostap.fitting.dataset 
        logger.warning ( 'The `weight`:%s will be ignored' % data.wname ()  )
        

    ## 
    from   ostap.trees.cuts     import vars_and_cuts
    from   ostap.stats.statvars import data_slice
    varlst , cuts, _  = vars_and_cuts ( expressions , cuts )
    varlst = tuple ( sorted ( varlst ) )
        
    assert varlst , "Empty variable list!"

    ## (1) create numpy dataset 
    ds , _ = data_slice ( data ,
                          varlst                 , 
                          cuts       = cuts      ,
                          first      = first     ,
                          last       = last      ,
                          cut_range  = cut_range ,
                          structured = True      , ## ATTENTION!  
                          transpose  = True      ,
                          progress   = progress  ,
                          use_frame  = use_frame , 
                          parallel   = parallel  )

    ## (2) convert it to Panda frame 
    df = pd.DataFrame ( ds )

    ## (3) make sweetviz report
    report   = sv.analyze ( df , **anaconf )

    if not 'filepath' in showconf : 
        tmp_path = CleanUp.tempfile ( suffix = '.html' , prefix = 'ostap-sweetviz-analyze-' )
        showconf [ 'filepath' ] = tmp_path
    
    ## (4) show report as HTML
    report.show_html ( **showconf )

    return showconf.get ( 'filepath' , '' )

# =============================================================================
## Use sweetviz to compare two datasets 
#  @code
#  data = ....
#  data_sweetviz_compare ( data1 , data2 , 'x,y,z' , 'pt>1' )
#  @endcode 
#  @see sweetviz
#  @see sweetviz.analyze
#  @see https://github.com/fbdesignpro/sweetviz
#  @see https://pypi.org/project/sweetviz/
def data_sweetviz_compare ( data         ,
                            data2        ,
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
                            progress2    = None        ,
                            name1        = 'First dataset'  ,
                            name2        = 'Second dataset' ,                            
                            cmpconf      = {}          ,    ## parameters for sweetviz.compare
                            showconf     = {}          ) :  ## parameters for sweetviz.show_html
    
    """ Use sweetviz to analyse some statistics of data
    >>> data = ....
    >>> data_sweetviz_compare( data , 'x,y,z' , 'pt>1' )
    - see sweetviz
    - see https://github.com/fbdesignpro/sweetviz
    - see https://pypi.org/project/sweetviz/
    """
    
    ## check pre-requisites 
    if not np or not pd or not sv :
        logger.error ( "Numpy,pandas or sweetviz are not available, no action" )
        return

    if isinstance ( data , ROOT.RooAbsData ) and data.isWeighted() :
        import ostap.fitting.dataset 
        logger.warning ( 'The `weight`:%s will be ignored' % data.wname ()  )
        

    if expressions2 is None : expressions2 = expressions
    if cuts2        is None : cuts2        = cuts

    if first2       is None : first2       = first
    if last2        is None : last2        = last
    
    if cut_range2   is None : cut_range2   = cut_range
    if use_frame2   is None : use_frame2   = use_frame
    if parallel2    is None : parallel2    = parallel
    if progress2    is None : progress2    = progress
    
    ## 
    from   ostap.trees.cuts     import vars_and_cuts
    from   ostap.stats.statvars import data_slice
    
    varlst1 , cuts  , _  = vars_and_cuts ( expressions  , cuts )
    varlst1 = tuple ( sorted ( varlst1 ) )
    
    varlst2 , cuts2 , _  = vars_and_cuts ( expressions2 , cuts2 )
    varlst2 = tuple ( sorted ( varlst2 ) )

    assert varlst1 and len ( varlst1 ) == len ( varlst2 ) , "Different lengths for variable lists!"
     
    ## (1) create numpy datasets 
    ds1 , _ = data_slice ( data                   ,
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
    
    ds2 , _ = data_slice ( data2                   ,
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
    
    ## (2) convert it to Panda frame 
    df1 = pd.DataFrame ( ds1 ) , name1 
    df2 = pd.DataFrame ( ds2 ) , name2 

    ## (3) make sweetviz report
    report   = sv.compare ( df1 , df2 , **cmpconf )

    if not 'filepath' in showconf : 
        tmp_path = CleanUp.tempfile ( suffix = '.html' , prefix = 'ostap-sweetviz-compare-' )
        showconf [ 'filepath' ] = tmp_path
    
    ## (4) show report as HTML
    report.show_html ( **showconf )

    return showconf.get ( 'filepath' , '' )


for t in ( ROOT.TTree      ,
           ROOT.RooDataSet ) : 
    t.sweetviz_analyze = data_sweetviz_analyze
    t.sweetviz_compare = data_sweetviz_compare


# ===============================================================================
_new_methods_ = (
    data_sweetviz_analyze , 
    data_sweetviz_compare ,
    ##
    ROOT.TTree      .sweetviz_analyze , 
    ROOT.RooDataSet .sweetviz_analyze ,
    ##
    ROOT.TTree      .sweetviz_compare , 
    ROOT.RooDataSet .sweetviz_compare
)         
_decorated_classes_ = (
    ROOT.TTree     ,
    ROOT.RooDataSet    
)

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    if not np : logger.error ( 'numpy    is not available!' )
    if not pd : logger.error ( 'pandas   is not available!' )
    if not sv : logger.error ( 'sweetviz is not available!' )

# =============================================================================
##                                                                      The END 
# =============================================================================
  
