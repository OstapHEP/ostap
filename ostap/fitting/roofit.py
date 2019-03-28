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
## print RooCmdArg object 
def _rca_print_ ( self ) :
    """Print RooCmdArg object 
    """
    name = self.GetName()
    if   'NumCPU'               == name : return 'NumCPU(%d,%d)'        % ( self.getInt ( 0 ) , self.getInt ( 1 ) ) 
    elif 'Verbose'              == name : return 'Verbose(%r)'          %   self.getBool () 
    elif 'Strategy'             == name : return 'Strategy(%d)'         %   self.getInt ( 0 ) 
    elif 'PrintLevel'           == name : return 'PrintLevel(%d)'       %   self.getInt ( 0 ) 
    elif 'PrintEvalErrors'      == name : return 'PrintEvalErrors(%d)'  %   self.getInt ( 0 ) 
    elif 'Timer'                == name : return 'Timer(%r)'            %   self.getBool () 
    elif 'Warnings'             == name : return 'Warnings(%r)'         %   self.getBool () 
    elif 'SumW2Error'           == name : return 'SumW2Error(%r)'       %   self.getBool () 
    elif 'Extended'             == name : return 'Extended(%r)'         %   self.getBool () 
    elif 'Range'                == name : return 'Range(%s,%s,%r)'      % ( self.getDouble ( 0 ) ,
                                                                            self.getDouble ( 1 ) ,
                                                                            True if self.getInt ( 0 ) else False )
    elif 'RangeWithName'        == name : return "Range('%s',%r)"       % ( self.getString ( 0 ) , self.getBool() )
    elif 'Hesse'                == name : return 'Hesse(%r)'            %   self.getBool () 
    elif 'InitialHesse'         == name : return 'InitialHesse(%r)'     %   self.getBool () 
    elif 'Optimize'             == name : return 'Optimize(%r)'         %   self.getBool () 
    elif 'Minos'                == name : return 'Minos(%r)'            %   self.getBool () if self.getSet(0) else 'Minos({.})'      
    elif 'Save'                 == name : return 'Save(%r)'             %   self.getBool () 
    elif 'FitOptions'           == name : return "FitOptions('%s')"     % ( self.getString ( 0 ) )
    elif 'ExternalConstraints'  == name : return 'ExternalConstraints({.})'
    elif 'DataError'            == name : return 'DataError(%d)'        % ( self.getInt    ( 0 ) )
    elif 'Minimizer'            == name : return "Minimizer('%s','%s')" % ( self.getString ( 0 ) , self.getString( 1 ) )

    elif 'ProjectedObservables' == name : return "ProjectedObservables({.})" 
    elif 'CutRange'             == name : return "CutRange('%s')"       %   self.getString ( 0 ) 
    elif 'LineColor'            == name : return 'LineColor(%d)'        %   self.getInt    ( 0 ) 
    elif 'LineStyle'            == name : return 'LineStyle(%d)'        %   self.getInt    ( 0 ) 
    elif 'LineWidth'            == name : return 'LineWidth(%d)'        %   self.getInt    ( 0 )     
    elif 'FillColor'            == name : return 'FillColor(%d)'        %   self.getInt    ( 0 ) 
    elif 'FillStyle'            == name : return 'FillStyle(%d)'        %   self.getInt    ( 0 ) 
    elif 'MarkerColor'          == name : return 'MarkerColor(%d)'      %   self.getInt    ( 0 ) 
    elif 'MarkerStyle'          == name : return 'MarkerStyle(%d)'      %   self.getInt    ( 0 ) 
    elif 'MarkerSize'           == name : return 'MarkerSize (%s)'      %   self.getDouble ( 0 ) 
    elif 'DrawOption'           == name : return "DrawOptions('%s')"    %   self.getString ( 0 ) 
    elif 'VLines'               == name : return "VLines()" 
    elif 'VisualizeError'       == name : return "VisializeError({.})" 
    elif 'VisualizeErrorData'   == name : return "VisializeError({.})" 
    elif 'ShowProgress'         == name : return "ShowProgress()"
    
    elif 'CutSpec'              == name : return "Cut('%s')"            %   self.getString ( 0 ) 
    elif 'BinningName'          == name : return "Binning('%s')"        %   self.getString ( 0 ) 
    elif 'BinningSpec'          == name : return 'Binning(%d,%s,%s)'    % ( self.getInt    ( 0 ) , self.getDouble ( 0 ) , self.getDouble ( 1 ) )
    elif 'Normalization'        == name : return 'Normalization(%s,%d)' % ( self.getDouble ( 0 ) , self.getInt    ( 0 ) )
    elif 'SelectCompSpec'       == name : return "Component('%s')"      %   self.getString ( 0 ) 

    
    return name

def _rca_bool_ ( self ) :
    """Get boolean value"""
    return True if self.getInt ( 0 ) else False

ROOT.RooCmdArg .__str__  = _rca_print_
ROOT.RooCmdArg .__repr__ = _rca_print_
ROOT.RooCmdArg .getBool  = _rca_bool_ 

_new_methods_ += [
    ROOT.RooCmdArg.__repr__ , 
    ROOT.RooCmdArg.__str__  , 
    ROOT.RooCmdArg.getBool  , 
    ]

# =============================================================================
_decorated_classes_ = ()

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
