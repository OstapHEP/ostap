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
    'report_prnt' , ## print the report 
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

# ==============================================================================
## modify constuctor for RDataFrame to enable/disable Implicit multithreading
#  @code
#  f1 = DataFrame ( ... , enable = True  ) ## default
#  f2 = DataFrame ( ... , enable = False ) ## default
#  @endcode
#  @see ROOT::EnableImplicitMT 
#  @see ROOT::DisableImplicitMT
#  @see ROOT::IsImplicitMTEnabled
def _fr_new_init_ ( self , *args , **kwargs ) :
    """Modify the DataFrame constuctor to allow (semi)automatic
    manipulations wth ROOT.ROOT.EnableImplicitMT/DisableImplicitMT
    - see ROOT.ROOT.EnableImplicitMT 
    - see ROOT.ROOT.DisableImplicitMT
    - see ROOT.ROOT.IsImplicitMTEnabled
    >>> f = DataFrame ( .... , enable = True  ) ## default 
    >>> f = DataFrame ( .... , enable = False )
    """
    
    mt = kwargs.pop ( 'enable' , True )
    
    if       mt and not ROOT.ROOT.IsImplicitMTEnabled() : ROOT.ROOT.EnableImplicitMT  ()
    elif not mt and     ROOT.ROOT.IsImplicitMTEnabled() : ROOT.ROOT.DisableImplicitMT ()
    
    return  self._fr_old_init_ ( *args , **kwargs )

if not hasattr ( DataFrame , '_fr_old_init_' ) :
    DataFrame._fr_old_init_ = DataFrame.__init__
    DataFrame.__init__      = _fr_new_init_
    
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

# ===============================================================================
## Print the frame report
def report_prnt (  report , prefix = '' ) :
    """Print a frame report
    """
    from ostap.core.core import binomEff
    table  = []
    lmax   = 5
    n0     = -1 
    for c in report :
        if  n0 <= 0 : n0 = c.GetAll () 
        name    = c.GetName ()
        passed  = c.GetPass ()
        all     = c.GetAll  ()
        eff1    = binomEff  ( passed , all ) * 100 
        eff2    = binomEff  ( passed ,  n0 ) * 100 
        table.append (  ( name , passed , all , eff1 , eff2 )  )
        lmax    = max ( len ( name ) , lmax , len ( 'Filter' ) )
    header   = '|    #input  |  #passed   |     efficiency [%]   | Cumulated efficiency [%] | ' 
    row      = '| %10d | %-10d | %8.3g +- %-8.3g | %8.3g +- %-8.3g     |'
    front    = '| %%-%ds ' % max ( lmax + 2 , len ( 'Selection' ) + 2 )
    prefix   = front % 'Filter'
    the_line = '\n# ' + '+' + ((len(prefix)-1)*'-') + '+' + (12*'-') + '+' + (12*'-') + '+' + (22*'-') + '+' + (26*'-') + '+'
    
    if prefix : text = str ( prefix ) + the_line 
    else      : text = the_line[1:]
    
    text    += '\n# ' + prefix  + header
    text    += the_line 
    for entry in table :
        n, p, a , e1 , e2 = entry
        line = row % ( a , p , e1.value() , e1.error() , e2.value() , e2.error() ) 
        text += '\n# ' + ( front % n ) + line
        
    return text + the_line 
        



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
