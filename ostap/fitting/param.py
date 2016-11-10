#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file
#  Auxillary utilities for parameterization of histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
#  
#                    $Revision$
#  Last modification $Date$
#  by                $Author$
# =============================================================================
"""Auxillary utilities for parameterization of histograms 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = () ## nothing to import 
# =============================================================================
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.param' )
else                       : logger = getLogger( __name__              )
# =============================================================================
logger.debug ( 'Auxillary utilities for Histogram parameterisation')
# =============================================================================
import ostap.fitting.funcs
import ostap.fitting.fitresult 
from   ostap.core.core import Ostap, funID

# =============================================================================
## @class H_fit
#  simple function to fit/represent the histogram with bernstein/spline
#  and other types of function expansion 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-05-09
class H_fit(object) :
    """Simple helper function to fit/represent the histogram with sum of
    Bernstein/b-spline/legendre/chebyshev, etc functions 
    """
    def __init__ ( self ,  hfit ) :
        self._hfit = hfit
        self.fun = ROOT.TF1 ( funID() , self , hfit.xmin() , hfit.xmax() , hfit.npars() )
    #
    def norm     ( self ) : return False 
    def npars    ( self ) : return self._hfit.npars () 
    def pars     ( self ) : return self._hfit.pars  ()
    #
    def draw     ( self , *args ) : return sef.fun.Draw( *args ) 
    def Draw     ( self , *args ) : return sef.fun.Draw( *args )
    
    def fit      ( self , h , opts = 'S' , *args ) : return h.Fit( self.fun , opts , *args ) 
    def Fit      ( self , h , opts = 'S' , *args ) : return h.Fit( self.fun , opts , *args ) 
    #   
    ## the major method 
    def __call__ ( self , x , pars = [] ) :
        
        x0 = x if isinstance ( x ,  ( int , long , float ) ) else x[0]
        
        if pars :
            np = self._hfit.npars() 
            for i in range ( 0 , np ) :    
                self._hfit.setPar ( i , pars[i] )
                
        return self._hfit( x0 )
    
# =============================================================================
## @class H_Nfit
#  simple function to fit/represent the histogram with normalized
#  functions, in particular:
#  - positive bernstein polynomial,
#  - positive B-spline expansion 
#  - positive monothonic B-spline expansion, etc...
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-05-09
class H_Nfit (object) :
    """ Simple helper function to fit/represent the histogram with
    the sum of bernstein positive polynominals
    """
    def __init__ ( self , hfit ) :
        self._hfit = hfit
        self.fun   = ROOT.TF1 ( funID() , self , hfit.xmin() , hfit.xmax() , hfit.npars() + 1 )
        self.fun.SetParameter ( 0 , 1 ) 
        
    def norm     ( self ) : return True  
    def npars    ( self ) : return self._hfit.npars () + 1  ## NB: normalization!  
    def pars     ( self ) : return self._hfit.pars  ()
    
    def draw     ( self , *args ) : return sef.fun.Draw( *args ) 
    def Draw     ( self , *args ) : return sef.fun.Draw( *args )
    
    def fit      ( self , h , opts = 'S' , *args ) : return h.Fit( self.fun , opts , *args ) 
    def Fit      ( self , h , opts = 'S' , *args ) : return h.Fit( self.fun , opts , *args ) 

    ## the major method 
    def __call__ ( self , x , pars = [] ) :

        norm = 1.0
        x0   = x if isinstance ( x ,  ( int , long , float ) ) else x[0]
        
        if pars :

            norm = pars[0]
            
            np   = self._hfit.npars() 
            for i in range ( 0 , np ) :    
                self._hfit.setPar ( i , pars[i+1] )
                
        return norm * self._hfit( x0 )


# =============================================================================
## helper class to wrap 1D-histogram as function
#  Optionally  normalization, bias and scale are applied
#  Seful e.g. for using a histogram as function fitting 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
class H1Func(object) :
    """Helper class to wrap 1D-histogram as function

    >>> histo =
    >>> func  = H1Func ( histo )
    >>> tf1   = ROOT.TF1('f1', func , low , high , 3 )

    Three parameters: NORM, BIAS and SCALE
    f( x )  = NORM  * histo (  ( x - BIAS ) / scale ) 
    
    >>> histo2 = ...
    >>> histo2.Fit ( tf1 , 'S' )
    
    """
    def __init__ ( self                               ,
                   histo                              ,
                   func        = lambda s : s.value() ,
                   interpolate = 2                    ,
                   edges       = True                 ) :
        
        self._histo  = histo
        self._func   = func
        self._interp = interpolate 
        self._edges  = edges  
        
    ## evaluate the function 
    def __call__ ( self , x , par = [ 1 , 0 , 1 ] ) :
        """ Evaluate the function 
        """
        #
        x0 = x if isinstance ( x , ( int , long , float ) ) else x[0]
        #
        norm  = float ( par[0] )   ## NORM 
        bias  = float ( par[1] )   ## BIAS 
        scale = float ( par[2] )   ## SCALE 
        #
        x0    = ( x0 - bias ) / scale
        # 
        return norm * self._func ( self._histo ( x0 , interpolate = self._interp , edges = self._edges ) )

    ## get corresponsing ROOT.TF1 object 
    def tf1  ( self ) :
        """Get corresponsing ROOT.TF1 object 
        """
        if not hasattr ( self , '_tf1' ) : 
            
            mn = self._histo.xmin ()
            mx = self._histo.xmax ()
            self._tf1 =  ROOT.TF1 ( funID() , self , mn , mx , 3 )
            self._tf1.SetParNames (
                'Normalization' ,
                'Bias'          ,
                'Scale'
                )
            self._tf1.FixParameter ( 0 , 1 )
            self._tf1.FixParameter ( 1 , 0 )
            self._tf1.FixParameter ( 2 , 1 )

        return self._tf1

    def Draw ( self , *args ) : return self.draw ( *args ) 
    def draw ( self , *args ) :
        t = self.tf1()
        return t.Draw( *args )


# ==============================================================================
## helper class to wrap 2D-histogram as function 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
class H2Func(object) :
    """ Helper class to Wrap 2D-histogram as function 
    """
    def __init__ ( self , histo , func = lambda s : s.value() , interpolate = True ) :
        self._histo  = histo    
        self._func   = func
        self._interp = interpolate 
        
    ## evaluate the function 
    def __call__ ( self , x ) :
        """ Evaluate the function 
        """
        x0 = x[0]
        y0 = x[1]
        return self._func ( self._histo ( x0 , y0 , interpolate = self._interp ) )




## create function object 
def  _funobj0_ ( self ) :
    """Create function object 
    """
    if hasattr ( self , '_bfit' ) : return self._bfit
    self._bfit = H_fit( self )
    return self._bfit

## create function object 
def  _funobjN_ ( self ) :
    """Create function object 
    """
    if hasattr ( self , '_bfit' ) : return self._bfit
    self._bfit = H_Nfit( self )
    return self._bfit

## draw spline object
def _sp_draw_   ( self , opts = '' ) :
    """Draw spline object 
    """
    bf = self.funobj () 
    return bf.fun.Draw( opts ) 
    

for t in (  Ostap.Math.Bernstein        ,
            Ostap.Math.BernsteinEven    ,
            Ostap.Math.ChebyshevSum     ,
            Ostap.Math.LegendreSum      ,
            Ostap.Math.FourierSum       ,
            Ostap.Math.CosineSum        ,
            Ostap.Math.Polynomial       ,
            Ostap.Math.BSpline          ) :  t.funobj = _funobj0_

for t in (  Ostap.Math.Positive         ,
            Ostap.Math.PositiveEven     , 
            Ostap.Math.Monothonic       ,
            Ostap.Math.Convex           ,
            Ostap.Math.ConvexOnly       ,            
            Ostap.Math.PositiveSpline   ,
            Ostap.Math.MonothonicSpline ,
            Ostap.Math.ConvexSpline     , 
            Ostap.Math.ConvexOnlySpline ) : t.funobj = _funobjN_
            
# =============================================================================
## construct helper class 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_as_fun_ ( self , func = lambda s : s.value () , *args , **kwargs ) :
    """Construct the function from the histogram
    """
    return H1Func ( self , func , *args , **kwargs )

# =============================================================================
## construct helper class 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-03-27
def _h1_as_spline_ ( self , func = lambda s : s.value () , *args , **kwargs ) :
    """Construct the function/spline from the histogram 
    """
    return H1Spline ( self , func , *args , **kwargs )

# =============================================================================
## construct helper class 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_as_fun_ ( self , func = lambda s : s.value () , *args , **kwargs ) :
    """Construct the helper function object from the histogram
    """
    return H2Func ( self , func , *args , **kwargs )
# =============================================================================
## construct function 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_as_tf1_ ( self                           ,
                  func   = lambda s : s.value () ,
                  spline = False , *args , **kwargs ) :
    """Construct the function from the 1D-histogram
    >>> histo = ...
    >>> fun1  = histo.asTF1 ( spline = False )
    >>> fun2  = histo.asTF1 ( spline = True  ) 
    """
    #
    if spline : fun = _h1_as_spline_ ( self , func , *args , **kwargs )
    else      : fun = _h1_as_fun_    ( self , func , *args , **kwargs )
    #
    f1 = fun .tf1  ()
    nb = self.nbins()
    
    if f1.GetNpx() <  1.2 * nb :
        f1.SetNpx ( max ( 100 , 10 * nb ) ) 
    
    return f1 


# =============================================================================
## calculate the integral of TH1
#  via conversion of ROOT::TH1 to ROOT::TF1 and using ROOT::TF1::Integral
#  @code
#  histo = ...
#  i     = histo.integral ( 0.2 , 0.16 ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-08-31
def _h1_integral_ ( histo , xlow , xhigh , *args , **kwargs ) :
    """Calculate the integral of TH1
    (convert ROOT::TH1 to ROOT::TF1 and use ROOT::TF1::Integral)    
    >>> histo = ...
    >>> i1 = histo.integral ( 0.2 , 0.16 )
    >>> i2 = histo.integral ( 0.2 , 0.16 , interpolate = 3 ) ## use cubic interpolation    
    """
    xlow  = max ( xlow  , histo.xmin () )
    xhigh = min ( xhigh , histo.xmax () )
    fun   = _h1_as_fun_ ( histo , *args , **kwargs )
    f1    = fun.tf1()
    return f1.Integral ( max ( xlow  , histo.xmin () ) ,
                         min ( xhigh , histo.xmax () ) )

# =============================================================================
## construct function 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_as_tf2_ ( self , func = lambda s : s.value () , *args , **kwargs ) :
    """Construct the function from the histogram
    >>> fun = h2.asFunc()
    >>> fun = h2.asTF  ()  ## ditto 
    >>> fun = h2.asTF1 ()  ## ditto 
    >>> fun = h2.asTF  ( interpolate = (3,3) )  ## use bi-cubic interpolation 
    """
    ax  = self.GetXaxis()
    ay  = self.GetYaxis()
    #
    fun = _h2_as_fun_ ( self , func , *args , **kwargs )
    #
    f2  = ROOT.TF2  ( funID()       ,
                      fun           ,
                      ax.GetXmin () ,
                      ax.GetXmax () ,
                      ay.GetXmin () ,
                      ay.GetXmax () ) 
    
    f2.SetNpx  ( 10 * ax.GetNbins() )
    f2.SetNpy  ( 10 * ay.GetNbins() )
    
    f2._funobj = fun  
    f2._histo  = fun._histo
    f2._func   = fun._func
    
    return f2

# =============================================================================
## calculate the integral of TH2
#  @code
#  histo = ...
#  i     = histo.integral ( 0.2 , 0.16 , 0.1 , 1.0 ) 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-08-31
def _h2_integral_ ( histo , xlow , xhigh , ylow , yhigh , *args , **kwargs ) :
    """Calculate the integral of TH2
    (convert ROOT::TH2 to ROOT::TF2 and use ROOT::TF2::Integral)
    >>> histo = ...
    >>> i     = histo.integral ( 0.2 , 0.16 , 0.1 , 1.0 ) 
    """
    f2  = _h2_as_tf2_  ( histo , *args , **kwargs )
    return f2.Integral ( xlow , xhigh , ylow , yhigh )

# =============================================================================
    
ROOT.TH1F . asTF     = _h1_as_tf1_ 
ROOT.TH1D . asTF     = _h1_as_tf1_ 
ROOT.TH2F . asTF     = _h2_as_tf2_ 
ROOT.TH2D . asTF     = _h2_as_tf2_ 
ROOT.TH1F . asTF1    = _h1_as_tf1_ 
ROOT.TH1D . asTF1    = _h1_as_tf1_ 
ROOT.TH2F . asTF2    = _h2_as_tf2_ 
ROOT.TH2D . asTF2    = _h2_as_tf2_ 
ROOT.TH1F . asFunc   = _h1_as_fun_ 
ROOT.TH1D . asFunc   = _h1_as_fun_ 
ROOT.TH2F . asFunc   = _h2_as_fun_ 
ROOT.TH2D . asFunc   = _h2_as_fun_ 
ROOT.TH1F . asSpline = _h1_as_spline_ 
ROOT.TH1D . asSpline = _h1_as_spline_ 

ROOT.TH1F . integral = _h1_integral_ 
ROOT.TH1D . integral = _h1_integral_
ROOT.TH2F . integral = _h2_integral_ 
ROOT.TH2D . integral = _h2_integral_


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.logger.line import line 
    logger.info ( __file__  + '\n' + line  ) 
    logger.info ( 80*'*'   )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' ) 

# =============================================================================
# The END 
# =============================================================================
