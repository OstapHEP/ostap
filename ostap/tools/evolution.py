#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ==========================================================================================
## @file ostap/tools/evolution.py
#  Useful utility that allows to inspect non-parametric
#  evolution of the generic shape of distribution inm bins of othe variables
#  E.g. one can check if the shape of the resolution fuction depends on eg. pt
#  @date   2025-06-06
#  @author Vanya  BELYAEV Ivan.Belyaev@itep.ru
# =============================================================================
""" Useful utility that allows to inspect non-parametric
evolution of the generic shape of distribution inm bins of othe variables
E.g. one can check if the shape of the resolution fuction depends on eg. pt
"""
# =============================================================================
__author__  = 'Vanya BELYAEV  Ivan.Belyaev@itep.ru'
__date__    = "2025-06-06"
__version__ = '$Revision$'
# =============================================================================
__all__     = (
    'evlution' , ## the only one useful function
    )
# =============================================================================
from   ostap.core.ostap_types   import integer_types, sized_types  
from   ostap.core.pyrouts       import Ostap, hID
from   ostap.frames.frames      import frame_project, Frames_OK 
from   ostap.utils.valerrors    import ValWithErrors 
from   ostap.plotting.canvas    import use_canvas
from   ostap.utils.progress_bar import progress_bar
from   ostap.utils.timing       import timing 
import ostap.trees.trees
import ostap.histos.histos 
import ostap.histos.graphs 
import ostap.stats.moments
import ostap.stats.moment
import ROOT 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
# =============================================================================
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.tools.evolution' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
if Frames_OK : # ==============================================================
    # =========================================================================
    from ostap.frames.frames import frame_the_moment, frame_the_statVar
    def _X_stat ( tree , xvar, cuts ) :
        return frame_the_statVar ( tree , xvar , cuts , progress = False , report = False , lazy = False )
    ## get the momnet anc x
    def _X_and_Ymoment_ ( tree ,  N , yvar , xvar , cuts ) :
        mm = frame_the_moment  ( tree , N , yvar , cuts , progress = False , report = False , lazy = True )
        xs = frame_the_statVar ( tree ,     xvar , cuts , progress = False , report = False , lazy = True )
        mm = mm.GetValue()
        xs = xs.GetValue()
        return xs, mm 
    # ========================================================================= 
else : # ======================================================================
    #  ========================================================================
    from ostap.frames.frames import frame_moment, frame_statVar
    def _X_stat ( tree , xvar, cuts ) :
        return frame_statVar ( tree , xvar , cuts )
    ## get X and Y-moment the momnet anc x
    def _X_and_Ymoment_ ( tree ,  N , yvar , xvar , cuts ) :
        mm = frame_moment  ( tree , N , yvar , cuts )
        xs = frame_statVar ( tree ,     xvar , cuts )
        return xs, mm 
    
# ===========================================================================
## Helper function to study the non-parametric evolution of the distribution
#  shape of one variabl eas function of anmothe variable,
#  e.f. shape of resolution function
#  Generic features are studied
#  - mean value
#  - RMS 
#  - skewness 
#  - kurtosis
#  - 5th and 6th standartized moments
#  @code
#  tree = ...
#  result = evolution ( tree , 'm - m_true' , 'pt' , bins = 10 , draw = True , silent = False  )
#  
#  result['MEAN'    ][0].draw('ape1')
#  result['MEAN'    ][1].draw('pe1') 
#  result['RMS'     ][0].draw('ape1')
#  result['RMS'     ][1].draw('pe1') 
#  result['SKEWNESS'][0].draw('ape1')
#  result['SKEWNESS'][1].draw('pe1') 
#  result['KURTOSIS'][0].draw('ape1')
#  result['KURTOSIS'][1].draw('pe1') 
#  result['STD5'    ][0].draw('ape1')
#  result['STD5'    ][1].draw('pe1') 
#  result['STD6'    ][0].draw('ape1')
#  result['STD6'    ][1].draw('pe1') 
#  @endcode
#  @param tree (INPUT) input TTree
#  @param yvar (INOUT) variable/expression to be invstigated
#  @param xvar (INOUT) variable/expression 
#  @param bins (INPUT) how many bins in x ?
#  @param draw (INOUT) visualize the results in process?
#  @param silent (INPUT) silent processoing?
def evolution ( tree , yvar , xvar , cuts , bins = 10 , draw = True , silent = False ) :
    """ Helper function to study the non-parametric evolution of the distribution
     shape of one variabl eas function of anmothe variable,
    e.f. shape of resolution function
    Generic features are studied
    - mean value
    - RMS 
    - skewness 
    - kurtosis
    - 5th and 6th standartized moments
    >>> tree = ...
    >>> result = evolution ( tree , 'm - m_true' , 'pt' , bins = 10 , draw = True , silent = False  )
    
    >>> result['MEAN'    ][0].draw('ape1')
    >>> result['MEAN'    ][1].draw('pe1') 

    >>> result['RMS'     ][0].draw('ape1')
    >>> result['RMS'     ][1].draw('pe1') 

    >>> result['SKEWNESS'][0].draw('ape1')
    >>> result['SKEWNESS'][1].draw('pe1') 

    >>> result['KURTOSIS'][0].draw('ape1')
    >>> result['KURTOSIS'][1].draw('pe1') 

    >>> result['STD5'    ][0].draw('ape1')
    >>> result['STD5'    ][1].draw('pe1') 

    >>> result['STD6'    ][0].draw('ape1')
    >>> result['STD6'    ][1].draw('pe1') 
    """
    
    xvar = str ( ROOT.TCut ( xvar ) )
    yvar = str ( ROOT.TCut ( yvar ) )
    
    if isinstance ( bins , integer_types ) and 0 < bins :
        #
        with timing ( 'Choose proper binning scheme' , logger = logger ) : 
            xlims  = _X_stat ( tree , xvar , cuts )
            xmnmx  = xlims.minmax() 
            hh     = ROOT.TH1F ( hID() , '' , 200 * bins , *xmnmx )
            tree.fproject ( hh , xvar , cuts ) 
            from ostap.histos.histos import h1_axis 
            edges = hh.equal_edges ( bins )
            del hh 
            if not silent :
                vv = ', '.join ( '%.4g' % e for e in edges ) 
                logger.info ( 'Binning scheme: %s' % vv )
            hh    = h1_axis ( edges ) 
            bins  = tuple ( item for item in hh.bin_edges() )            
            del hh
            
    elif isinstance ( bins , ROOT.TAxis ) :
        bins   = tuple ( item  for item in bins.bin_edges() )             
    elif isinstance ( bins , ROOT.TH1  )  and 1 == bins.GetDimension () :
        bins   = tuple ( item  for item in bins.bin_edges() )             
    else :
        assert all ( isinstance ( item , sized_types ) and 2 == len ( item ) for item in bins ) , \
            "Invalid bins: %s" % str ( bins ) 

    ## max moment is 6, therefore 2*6=12 here 
    N      = 12 

    def calculate ( the_cuts ) :
        
        xs, mm  = _X_and_Ymoment_ ( tree , N , yvar, xvar , the_cuts )
        
        mean     = mm.mean           ()
        rms      = mm.rms            () 
        skewness = mm.skewness       () 
        kurtosis = mm.kurtosis       () 
        std5     = mm.std_moment_[5] () 
        std6     = mm.std_moment_[6] ()
        
        return mm , xs.mean() , mean, rms , skewness , kurtosis, std5 , std6 

    TGE  = ROOT.TGraphErrors
    TGAE = ROOT.TGraphAsymmErrors

    NB = len ( bins ) 
    g1_MEAN , g2_MEAN = TGAE ( NB ) , TGE ( 1 )
    g1_RMS  , g2_RMS  = TGAE ( NB ) , TGE ( 1 )
    g1_SKEW , g2_SKEW = TGAE ( NB ) , TGE ( 1 )
    g1_KURT , g2_KURT = TGAE ( NB ) , TGE ( 1 )
    g1_STD5 , g2_STD5 = TGAE ( NB ) , TGE ( 1 )
    g1_STD6 , g2_STD6 = TGAE ( NB ) , TGE ( 1 )

    for g in ( g1_MEAN, g1_RMS , g1_SKEW, g1_KURT , g1_STD5 , g1_STD6 ) : g.blue()
    for g in ( g2_MEAN, g2_RMS , g2_SKEW, g2_KURT , g2_STD5 , g2_STD6 ) : g.red ()

    #
    ## start with a single global bin
    with timing ( "Process the single global bin" , logger = logger ) :
        mm , xmean , mean, rms , skewness , kurtosis, std5 , std6 = calculate ( cuts )
        if not silent : 
            title = '%s in global single bin' % yvar 
            logger.info ( '%s:\n%s' % ( title , mm.table ( title = title , prefix = '# ' ) ) )

        g2_MEAN [ 0 ] = xmean , mean    
        g2_RMS  [ 0 ] = xmean , rms      
        g2_SKEW [ 0 ] = xmean , skewness 
        g2_KURT [ 0 ] = xmean , kurtosis 
        g2_STD5 [ 0 ] = xmean , std5     
        g2_STD6 [ 0 ] = xmean , std6     

    with timing ( "Process individual bins" , logger = logger ) :
        
        ibin  = 0 
        for xlow, xhigh in progress_bar ( bins , silent = silent ) : 
            
            cut_low  = ROOT.TCut ( '%.10g<=%s' % ( xlow , xvar  ) )
            cut_high = ROOT.TCut ( '%s<%.10g'  % ( xvar , xhigh ) )
            the_cut  = cut_low & cut_high & cuts 
            
            _ , xmean , mean, rms , skewness , kurtosis, std5 , std6 = calculate ( the_cut ) 

            xm = float ( xmean ) 
            xm = ValWithErrors ( xm , xlow - xm , xhigh - xm )
            
            g1_MEAN [ ibin ] = xm , mean     
            g1_RMS  [ ibin ] = xm , rms      
            g1_SKEW [ ibin ] = xm , skewness 
            g1_KURT [ ibin ] = xm , kurtosis 
            g1_STD5 [ ibin ] = xm , std5     
            g1_STD6 [ ibin ] = xm , std6     

            ibin += 1
            
            if draw and 1 == ibin :
                with use_canvas  ( "Mean     %s vs %s " % ( yvar, xvar ) ) :
                    g1_MEAN.draw ( 'ap' )
                    g2_MEAN.draw ( 'p'  )
                with use_canvas  ( "Rms      %s vs %s " % ( yvar, xvar ) ) :
                    g1_RMS .draw ( 'ap' )
                    g2_RMS .draw ( 'p'  )            
                with use_canvas  ( "Skewness %s vs %s " % ( yvar, xvar ) ) :
                    g1_SKEW.draw ( 'ap' )
                    g2_SKEW.draw ( 'p'  )                
                with use_canvas  ( "Kurtosis %s vs %s " % ( yvar, xvar ) ) :
                    g1_KURT.draw ( 'ap' )
                    g2_KURT.draw ( 'p'  )
                with use_canvas  ( "Std-5    %s vs %s " % ( yvar, xvar ) ) :
                    g1_STD5.draw ( 'ap' )
                    g2_STD5.draw ( 'p'  )
                with use_canvas  ( "Std-6    %s vs %s " % ( yvar, xvar ) ) :
                    g1_STD6.draw ( 'ap' )
                    g2_STD6.draw ( 'p'  )


    from ostap.utils.cidict import cidict, cidict_fun 
    graphs = cidict( transform = cidict_fun ) 
    
    graphs [ 'MEAN'     ] = g1_MEAN , g2_MEAN
    graphs [ 'RMS'      ] = g1_RMS  , g2_RMS 
    graphs [ 'SKEWNESS' ] = g1_SKEW , g2_SKEW
    graphs [ 'KURTOSIS' ] = g1_KURT , g2_KURT
    graphs [ 'STD5'     ] = g1_STD5 , g2_STD5
    graphs [ 'STD6'     ] = g1_STD6 , g2_STD6
    
    return graphs 


# =============================================================================
##                                                                      The END 
# =============================================================================
