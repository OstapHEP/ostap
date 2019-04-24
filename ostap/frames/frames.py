#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/frames/frames.py
#  Module with decoration of TDataFrame objects for efficient use in python
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-06-16
# =============================================================================
"""Decoration of Tree/Chain objects for efficient use in python"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'DataFrame'   , ## RDataFrame object 
    ) 
# =============================================================================
import ROOT
from ostap.core.core import cpp, Ostap 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.frames.frames' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some useful decorations for ROOT::RDataFrame objects')
# =============================================================================
try : 
    DataFrame = ROOT.ROOT.RDataFrame
except AttributeError :
    DataFrame = ROOT.ROOT.Experimental.TDataFrame 

Ostap.DataFrame    = DataFrame 

DataFrame.columns  = lambda s : tuple( s.GetColumnNames() ) 
DataFrame.branches = DataFrame.columns 

# =============================================================================
## Get the length/size of the data frame
#  @code
#  frame = ...
#  print len(frame)
#  @endcode 
def _fr_len_ ( f ) :
    """Get the length/size of the data frame
    >>> frame = ...
    >>> print len(frame)
    """
    return f.Count().GetValue() 
    

# =============================================================================
## Get the effective entries in data frame
#  @code
#  data = ...
#  neff = data.nEff('b1*b1')
#  @endcode
def _fr_nEff_  ( self , cuts = '' ) :
    """Get the effective entries in data frame 
    >>> data = ...
    >>> neff = data.nEff('b1*b1')
    """
    return Ostap.StatVar.nEff ( self , cuts )

# =============================================================================
## Get statistics for the  given expression in data frame
#  @code
#  data = ...
#  c1 = data.statVar( 'S_sw' , 'pt>10' ) 
#  c2 = data.statVar( 'S_sw' , 'pt>0'  )
#  @endcode
def _fr_statVar_ ( self  , expression ,  cuts = '' ) :
    """Get statistics for the  given expression in data frame
    >>> data = ...
    >>> c1 = data.statVar( 'S_sw' , 'pt>10' ) 
    >>> c2 = data.statVar( 'S_sw' )
    """
    return Ostap.StatVar.statVar( self , expression , cuts , )

# =============================================================================
## get the statistic for pair of expressions in DataFrame
#  @code
#  frame  = ...
#  stat1 , stat2 , cov2 , len = frame.statCov( 'x' , 'y' )
#  # apply some cuts 
#  stat1 , stat2 , cov2 , len = frame.statCov( 'x' , 'y' , 'z>0' )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-06-18
def _fr_statCov_ ( frame       ,
                   expression1 ,
                   expression2 ,
                   cuts = ''   ) :
    """Get the statistic for pair of expressions in DataFrame
    
    >>>  frame  = ...
    >>>  stat1 , stat2 , cov2 , len = framw.statCov( 'x' , 'y' )
    
    Apply some cuts:
    >>> stat1 , stat2 , cov2 , len = frame.statCov( 'x' , 'y' , 'z>0' )
    
    """
    import ostap.math.linalg 
    stat1  = Ostap.WStatEntity       ()
    stat2  = Ostap.WStatEntity       ()
    cov2   = Ostap.Math.SymMatrix2x2 ()

    length = Ostap.StatVar.statCov ( frame       ,
                                     expression1 ,
                                     expression2 ,
                                     cuts        ,    
                                     stat1       ,
                                     stat2       ,
                                     cov2        ) 
        
    return stat1 , stat2 , cov2 , length


# =============================================================================
## Simplified print out for the  frame 
#  @code 
#  frame = ...
#  print frame
#  @endcode 
def _fr_print_ ( t ) :
    """Simplified print out for the  frame 
    
    >>> frame = ...
    >>> print frame
    """
    ##
    res = "Frame Enries/#%d" %  len ( t )  
    ##
    _b          = t.columns ()
    res        +=        "\nBranches: %s" % list( _b )
    return res

# ==============================================================================
# decorate 
# ==============================================================================
if not hasattr ( DataFrame , '__len__') : DataFrame.__len__ = _fr_len_
DataFrame .__str__    = _fr_print_ 
DataFrame .__repr__   = _fr_print_

DataFrame .nEff       = _fr_nEff_
DataFrame .statVar    = _fr_statVar_
DataFrame .statCov    = _fr_statCov_


from ostap.stats.statvars import  data_decorate 
data_decorate ( DataFrame )
del data_decorate

from ostap.fitting.dataset import ds_draw , ds_project
DataFrame.draw    = ds_draw
DataFrame.project = ds_project

## # =============================================================================
## def _fproject_ ( tree , histo ,  what , cuts ) :

##     assert isinstance ( what , ( str , ROOT.TCut ) ) , "Invalid ``what'':%s" % what
##     assert isinstance ( cuts , ( str , ROOT.TCut ) ) , "Invalid ``cuts'':%s" % cuts
    
##     what = str ( what )
##     cuts = str ( cuts )

##     frame = DataFrame ( tree ).Filter ( cuts ) 

##     if isinstance ( histo , ROOT.TProfile2D ) :
##         histo  = ROOT.ROOT.RDF.TProfile2DModel ( histo ) 
##         return frame.Profile2D

    
##         return _fproject_ ( tree , ROOT.ROOT.RDF.TProfile2DModel ( histo ) , what , cuts )                             
##     elif isinstance ( histo , ROOT.TProfile   ) :
##         return _fproject_ ( tree , ROOT.ROOT.RDF.TProfile1DModel ( histo ) , what , cuts )
##     elif isinstance ( histo , ROOT.TH3 ) :
##         if isinstance ( histo,  ROOT.TH3D ) :
##             return _fproject_ ( tree , ROOT.ROOT.RDF.TH3DModel  ( histo ) , what , cuts )
##         return _fproject_ ( tree , ROOT.TH3D ( hd ) , what , cuts )
##     elif isinstance ( histo , ROOT.TH2 ) :
##         if isinstance ( histo,  ROOT.TH2D ) :
##             return _fproject_ ( tree , ROOT.ROOT.RDF.TH2DModel  ( histo ) , what , cuts )
##         return _fproject_ ( tree , ROOT.TH2D ( hd ) , what , cuts )
##     elif isinstance ( histo , ROOT.TH1 ) :
##         if isinstance ( histo,  ROOT.TH1D ) :
##             return _fproject_ ( tree , ROOT.ROOT.RDF.TH1DModel  ( histo ) , what , cuts )
##         return _fproject_ ( tree , ROOT.TH1D ( hd ) , what , cuts )

##     if not ROOT.ROOT.IsImplicitMTEnabled() :
##         ROOT.ROOT.EnableImplicitMT() :

##     frame = DataFrame ( tree )
    
##     h1    = frame.Filter ( cuts ).Histo1D
    


_decorated_classes_ = (
    DataFrame   ,
    )

_new_methods_       = (
    DataFrame.__len__          ,
    DataFrame.__repr__         ,
    DataFrame.__str__          ,
    DataFrame.columns          , 
    DataFrame.branches         ,
    #
    DataFrame.draw             , 
    DataFrame.project          ,
    #
    DataFrame.nEff             ,
    DataFrame.statVar          ,
    DataFrame.statCov          ,
    DataFrame.nEff             ,
    #
    DataFrame.get_moment       , 
    DataFrame.central_moment   , 
    DataFrame.mean             ,
    DataFrame.rms              ,
    DataFrame.skewness         ,
    DataFrame.kurtosis         ,
    DataFrame.quantile         ,
    DataFrame.median           ,
    DataFrame.quantiles        ,
    DataFrame.interval         ,
    DataFrame.terciles         ,
    DataFrame.quartiles        ,
    DataFrame.quintiles        ,
    DataFrame.deciles          ,
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
