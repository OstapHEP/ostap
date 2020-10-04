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
    'FIXVAR'        , ## context manager to fix/unfix the variable 
    'var_from_name' , ## "convert" name/expression into variable/formula
    ) 
# =============================================================================
import ROOT, random
from   ostap.core.core              import Ostap, VE
from   ostap.fitting.variables      import SETVAR, FIXVAR  
import ostap.fitting.roocollections
import ostap.fitting.roofitresult
import ostap.fitting.printable
import ostap.fitting.roocmdarg   
from   ostap.fitting.dataset        import setStorage, useStorage
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
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
    def __init__ ( self , pdf , xvar , xmin = None , xmax = None , norm = 1 ) :
        
        self.__pdf     = pdf
        self.__norm    = float ( norm )  

        ## ostap stuff: 
        if not isinstance ( pdf , ROOT.RooAbsPdf ) :
            if hasattr ( self.pdf , 'pdf' ) :
                self.__pdf_ =   pdf 
                self.__pdf  = pdf.pdf

        self.__xvar   = xvar
        self.__xmin   = None 
        self.__xmax   = None
        
        if not xmin is None : self.__xmin = xmin 
        if not xmax is None : self.__xmax = xmax

        if hasattr ( xvar , 'getMin' ) :
            if self.__xmin is None : self.__xmin = xvar.getMin()
            else                   : self.__xmin = max ( self.__xmin , xvar.getMin() )
            
        if hasattr ( xvar , 'getMax' ) :
            if self.__xmax is None : self.__xmax = xvar.getMax()
            else                   : self.__xmax = min ( self.__xmax , xvar.getMax() )
            
        if self.__xmin is None :
            raise AttributeError ( "xmin can't be deduced from  input arguments" )

        if self.__xmax is None :
            raise AttributeError ( "xmax can't be deduced from  input arguments" )
        
        if self.__xmin > self.__xmax :
            self.__xmin , self.__xmax = self.__xmax , self.__xmin

    @property
    def pdf ( self ) :
        """``pdf'' : get the actual RooAbsPdf
        """
        return self.__pdf
    @property
    def norm ( self ) :
        """``norm'' : additional normalization factor
        """
        return self.__norm

    @property
    def xvar ( self ) :
        """``xvar'': x-variable
        """
        return self.__xvar
    
    def xmin     ( self ) : return self.__xmin
    def xmax     ( self ) : return self.__xmax
    
    ## the main method 
    def __call__ ( self , x , pars = [] ) :

        ## for ROOT.TF1
        if   hasattr ( x , '__len__' ) and hasattr( x , 'getitem' ) and not len( x )   : x = x [ 0 ]
        elif hasattr ( x , '__len__' ) and hasattr( x , 'getitem' ) and 0 != x.size()  : x = x [ 0 ]
        
        ## try to be efficient 
        if not self.__xmin <= x <= self.__xmax : return 0 
        
        with SETVAR ( self.__xvar ) :
            self.__xvar.setVal ( x )
            return self.__pdf.getVal() * self.__norm 



# =============================================================================
_decorated_classes_ = ( )

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                     The END 
# =============================================================================
