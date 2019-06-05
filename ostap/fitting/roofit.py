#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/roofit.py
#  Module with decoration of some RooFit objects for efficient use in python
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Decoration of some RooFit objects for efficient use in python
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'setStorage'    , ## define the default storage for  RooDataStore 
    'useStorage'    , ## define (as context) the default storage for  RooDataStore 
    'PDF_fun'       , ## wrapper of PDF to ``simple'' function 
    'SETVAR'        , ## context manager to preserve the current value for RooRealVar
    'var_from_name' , ## "convert" name/expression into variable/formula
    ) 
# =============================================================================
import ROOT, random
from   ostap.core.core              import Ostap, VE
from   ostap.fitting.variables      import SETVAR 
import ostap.fitting.roocollections
import ostap.fitting.roofitresult
import ostap.fitting.printable
import ostap.fitting.roocmdarg   
from   ostap.fitting.dataset        import setStorage, useStorage
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger , allright,  attention
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.rootfit' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug( 'Some useful decorations for RooFit objects')
# =============================================================================
_new_methods_ = []

# =============================================================================
## product of two PDFs 
def _pdf_mul_ ( pdf1 , pdf2 ) :
    """Easy contruct for the product of two PDFs:
    
    >>> pdf1 = ...
    >>> pdf2 = ...
    
    >>> product = pdf1 * pdf2 
    """
    return Ostap.Models.Product (
        '%s*%s'             % ( pdf1.GetName  () , pdf2.GetName  () ) ,
        'Product: %s & %s ' % ( pdf1.GetTitle () , pdf2.GetTitle () ) ,
        pdf1 , pdf2 )
# ============================================================================
ROOT.RooAbsPdf . __mul__  = _pdf_mul_ 

_new_methods_ += [
    ROOT.RooAbsPdf.__mul__  , 
    ]

# ============================================================================
## "convert" name/expression into variable/formula
def var_from_name ( w , varset ) :
    """ Convert name/expression into variable/formula
    """
    w = w.strip() 
    if    0 <= w.find('(') < what.find(')') : pass
    elif  0 <  w.find('*')                  : pass
    elif  0 <  w.find('/')                  : pass
    elif  0 <  w.find('%')                  : pass 
    elif  0 <  w.find('+')                  : pass
    elif  0 <  w.find('-')                  : pass
    else :
        v = varset[w]
        return v
    ##
    
    vlst = ROOT.RooArgList()
    for s in varset : vlst.add ( s )
    #
    f = ROOT.RooFormulaVar( w , w , vlst )
    return f 
    
# =============================================================================
## @class PDF_fun
#  Helper class to wrap PDF as 'function'
#  can be helpful for some pure math-operations
#  @code
#  pdf,var = ....
#  fun     = PDF( fun , var , xmin=0 , xmax=1 )
#  from ostap.stats.moments import mean, mode, median, CL   
#  print 'MEAN    : %s' % mean    ( fun , 0 , 1 )
#  print 'MODE    : %s' % mode    ( fun , 0 , 1 )
#  print 'MEDIAN  : %s' % median  ( fun , 0 , 1 )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-03-29
class PDF_fun(object):
    """ Helper class to wrap PDF as 'function'
    >>> pdf,var = ....
    >>> fun     = PDF( pdf , var , xmin=0 , xmax=1 )
    >>> print fun(0.1),fun(0.5) 
    >>> from ostap.stats.moments import mean, mode, median  
    >>> print 'MEAN    : %s' % mean    ( fun , 0 , 1 )
    >>> print 'MODE    : %s' % mode    ( fun , 0 , 1 )
    >>> print 'MEDIAN  : %s' % median  ( fun , 0 , 1 )
    """
    ##
    def __init__ ( self , pdf , xvar , xmin = None , xmax = None ) :
        
        self.pdf     = pdf

        ## ostap stuff: 
        if not isinstance ( pdf , ROOT.RooAbsPdf ) :
            if hasattr ( self.pdf , 'pdf' ) :
                self.pdf_ = pdf 
                self.pdf  = pdf.pdf

        self.xvar    = xvar

        self._xmin   = None 
        self._xmax   = None
        
        if not xmin is None : self._xmin = xmin 
        if not xmax is None : self._xmax = xmax

        if hasattr ( xvar , 'getMin' ) :
            if self._xmin is None : self._xmin = xvar.getMin()
            else                  : self._xmin = max ( self._xmin , xvar.getMin() )
            
        if hasattr ( xvar , 'getMax' ) :
            if self._xmax is None : self._xmax = xvar.getMax()
            else                  : self._xmax = min ( self._xmax , xvar.getMax() )
            
        if self._xmin is None :
            raise AttributeError ( "xmin can't be deduced from  input arguments" )
        if self._xmax is None :
            raise AttributeError ( "xmax can't be deduced from  input arguments" )
        
        if self._xmin > self._xmax :
            self._xmin , self._xmax = self._xmax , self._xmin
            
    def xmin     ( self ) : return self._xmin
    def xmax     ( self ) : return self._xmax
    
    ## the main method 
    def __call__ ( self , x , pars = [] ) :

        ## for ROOT.TF1
        if   hasattr ( x , '__len__' ) and hasattr( x , 'getitem' ) and not len( x )   : x = x[0]
        elif hasattr ( x , '__len__' ) and hasattr( x , 'getitem' ) and 0 != x.size()  : x = x[0]
        
        ## try to be efficient 
        if not self._xmin <= x <= self._xmax : return 0 
        
        with SETVAR( self.xvar ) :
            self.xvar.setVal ( x )
            return self.pdf.getVal()



# =============================================================================
_decorated_classes_ = ( )

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
