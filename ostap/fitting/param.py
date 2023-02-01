#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/param.py
#  Auxillary utilities for parameterization of histograms 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Auxillary utilities for parameterization of histograms 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'H_fit'  ,
    'H_Nfit' ,
    ) 
# =============================================================================
import ostap.histos.histos 
import ostap.fitting.fitresult 
from   ostap.core.core        import Ostap, funID
from   ostap.core.ostap_types import num_types, integer_types 
import ROOT, abc 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.param' )
else                       : logger = getLogger( __name__              )
# =============================================================================
logger.debug ( 'Auxillary utilities for Histogram parameterisation')
# =============================================================================
# @class C1Fun
# Helper wrapper for callable to TF1 object for fitting 
class C1Fun(object) :
    """Helper wrapper for  callable to  TF1 object for fitting
    """
    def __init__ ( self , fun , xmin , xmax ) :

        assert callable ( fun ) , "C1Fun: 'fun' must be callable!"
        
        self.__fun  = fun 
        self.__xmin = min ( xmin , xmax )
        self.__xmax = max ( xmin , xmax )

        self.__tf1  = ROOT.TF1 ( funID () , self , self.xmin , self.xmax , 3 )
        self.__tf1.SetParNames (
            'Norm'  ,
            'Bias'  ,
            'Scale' 
            )
        self.__tf1.SetParameter ( 0 , 1 )
        self.__tf1.FixParameter ( 1 , 0 )
        self.__tf1.FixParameter ( 2 , 1 )

    ## the actual call 
    def __call__ ( self , x , pars = [ 1 , 0 , 1 ] ) :
        """Call method"""
        
        x0    = x if isinstance ( x , num_types ) else x [ 0 ]
        #
        norm  = float ( pars [ 0 ] ) ## NORM 
        bias  = float ( pars [ 1 ] ) ## BIAS 
        scale = float ( pars [ 2 ] ) ## SCALE 
        #
        xx    = ( x0 - bias ) / scale
        # 
        fun = self.__fun 
        return norm * fun ( xx ) 
    
    ## make fit 
    def Fit ( self , histo , opts = 'S' , gopts = '' , *args ) :
        """Make a fit
        >>> obj   = ...
        >>> histo = ... 
        >>> obj.Fit  ( histo , 'S0Q' )
        """
        assert isinstance ( histo , ROOT.TH1 ) and 1 == histo.dim() , 'Invalid histo-type!'
        ## make a fit 
        return histo.Fit ( self.__tf1 , opts , gopts , *args )


    ## fix parameter 
    def fix ( self , index , value ) :
        """Fix parameter 
        """
        assert isinstance ( index , integer_types ) and 0 <= index <= 2  , 'Invalid index %s' % index
        value = float ( value ) 
        self.__tf1.FixParameter ( index , value )

    ## set parameter 
    def set ( self , index , value ) :
        """Set parameter 
        """
        assert isinstance ( index , integer_types ) and 0 <= index <= 2 , 'Invalid index %s' % index
        value = float ( value ) 
        self.__tf1.SetParameter ( index , value )

    ## set parameter 
    def __setitem__  ( self , index , valie ) :
        if not isinstance ( index , integer_types ) : raise IndexError ("Invalid index %s" % index )
        if not 0 <= index <= 2                      : raise IndexError ("Invalid index %s" % index )
        self.__tf1.SetParameter ( index , value )

    ## release parameter
    def release ( self , index ) :
        """Release parameter 
        """
        assert isinstance ( index , integer_types ) and 0 <= index < 2 , 'Invalid index %s' % index
        self.__tf1.ReleaseParameter ( index )

    ## 
    @property   
    def tf1  ( self ) :
        """Get corresponding ROOT.TF1 object 
        """
        return self.__tf1 

    @property
    def xmin ( self ) :
        """'xmin' : low edge"""
        return self.__xmin
    
    @property    
    def xmax ( self ) :
        """'xmax' : high edge"""
        return self.__xmax
    
    @property
    def fun  ( self ) :
        """'fun' : actual callable to   use"""
        return self.__fun

    def Draw ( self , *args , **kwargs ) : return self.draw ( *args , **kwargs ) 
    def draw ( self , *args , **kwargs ) :
        t = self.__tf1
        return t.draw( *args , **kwargs )


# =============================================================================
# @class HFit
class HFIT(object) :
    __metaclass__ = abc.ABCMeta

    ## constructor from hfit an dTF1 objects 
    def __init__ ( self , hfit , tf1 ) :
        
        assert hfit and callable   ( hfit )              , "HFIT: 'hfit' must be callable"
        assert tf1  and isinstance ( tf1   , ROOT.TF1 )  , "HFIT: 'tf1'  must be ROOT.TF1"
        
        self.__hfit = hfit
        self.__fun  = tf1
        
    @abc.abstractmethod
    def norm ( self ) :
        """Is this object represent normalised function?"""
        return True 
    @abc.abstractmethod
    def npars ( self ) :
        """number of parameters"""
        return 0
    @abc.abstractmethod
    def __call__ ( self , x , pars = [] ) :
        """The main call method"""
        return None

    @property
    def fun  ( self ) :
        """'fun' : actual ROOT.TF1  object"""
        return self.__fun 

    @property
    def hfit ( self ) :
        """'hfin' : actual `hfit' object"""
        return self.__hfit 

    ## get list of parameters 
    def pars     ( self ) : return self.hfit.pars  ()
    
    def draw ( self , *args , **kwargs )           : return self.fun.draw ( *args , **kwargs ) 
    def Draw ( self , *args , **kwargs )           : return self.fun.draw ( *args , **kwargs )    
    def fit  ( self , histo , opts = 'S' , *args ) : return h.Fit ( self.fun , opts , *args ) 
    def Fit  ( self , histo , opts = 'S' , *args ) : return h.Fit ( self.fun , opts , *args ) 
    
    
# =============================================================================
## @class H_fit
#  simple function to fit/represent the histogram with bernstein/spline
#  and other types of function expansion 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-05-09
class H_fit(HFIT) :
    """Simple helper function to fit/represent the histogram with sum of
    Bernstein/b-spline/legendre/chebyshev, etc functions 
    """
    def __init__ ( self ,  hfit , xmin = None , xmax = None ) :

        if xmin is None and hasattr ( hfit , 'xmin' ) : xmin = hfit.xmin ()
        if xmax is None and hasattr ( hfit , 'xmax' ) : xmax = hfit.xmax ()
        
        ## create the function
        tf1 = ROOT.TF1 ( funID() , self , xmin , xmax , hfit.npars() )

        ## initialize the base 
        super(H_fit,self).__init__ ( hfit , tf1 ) 
    #
    def norm     ( self ) : return False 
    def npars    ( self ) : return self.hfit.npars () 
    #
    ## the major method 
    def __call__ ( self , x , pars = [] ) :
        
        x0 = float( x  ) if isinstance ( x , num_types ) else x [ 0 ]
        
        if pars :
            np = self.hfit.npars () 
            for i in range ( 0 , np ) :    
                self.hfit.setPar ( i , pars[i] )
                
        return self.hfit( x0 )
    
# =============================================================================
## @class H_Nfit
#  simple function to fit/represent the histogram with normalized
#  functions, in particular:
#  - positive bernstein polynomial,
#  - positive B-spline expansion 
#  - positive monotonic B-spline expansion, etc...
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-05-09
class H_Nfit (HFIT) :
    """ Simple helper function to fit/represent the histogram with
    the sum of bernstein positive polynominals
    """
    def __init__ ( self , hfit , xmin = None , xmax = None ) :

        if xmin is None and hasattr ( hfit , 'xmin' ) : xmin = hfit.xmin ()
        if xmax is None and hasattr ( hfit , 'xmax' ) : xmax = hfit.xmax ()
        
        ## create function 
        tf1 = ROOT.TF1 ( funID() , self , xmin , xmax , hfit.npars() + 1 )
        
        ## initialize the base 
        super(H_Nfit,self).__init__ ( hfit , tf1 )

        ## set normalization parameter 
        self.fun.SetParameter ( 0 , 1 ) 

    def norm     ( self ) : return True  
    def npars    ( self ) : return self.hfit.npars () + 1  ## NB: normalization!  
    
    ## the major method 
    def __call__ ( self , x , pars = [] ) :

        norm = 1.0
        x0   = float ( x ) if isinstance ( x , num_types ) else x [ 0 ]
        
        if pars :

            norm = float ( pars [ 0 ] ) 
            
            np   = self.hfit.npars() 
            for i in range ( 0 , np ) :    
                self.hfit.setPar ( i , pars[i+1] )
                
        return norm * self.hfit( x0 )


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
    def __call__ ( self , x , pars = [ 1 , 0 , 1 ] ) :
        """ Evaluate the function 
        """
        #
        x0 = x if isinstance ( x , num_types ) else x [ 0 ]
        #
        norm  = float ( pars [ 0 ] ) ## NORM 
        bias  = float ( pars [ 1 ] ) ## BIAS 
        scale = float ( pars [ 2 ] ) ## SCALE 
        #
        x0    = ( x0 - bias ) / scale
        # 
        return norm * self._func ( self._histo ( x0 , interpolate = self._interp , edges = self._edges ) )

    ## get corresponding ROOT.TF1 object 
    def tf1  ( self ) :
        """Get corresponsing ROOT.TF1 object 
        """
        if not hasattr ( self , '_tf1' ) : 
            
            mn = self._histo.xmin  ()
            mx = self._histo.xmax  ()
            self._tf1 =  ROOT.TF1  ( funID() , self , mn , mx , 3 )
            self._tf1.SetParNames  ( 'Norm'  , 'Bias'  , 'Scale' )
            self._tf1.FixParameter ( 0 , 1 )
            self._tf1.FixParameter ( 1 , 0 )
            self._tf1.FixParameter ( 2 , 1 )

        return self._tf1

    def Draw ( self , *args , **kwargs ) : return self.draw ( *args , **kwargs ) 
    def draw ( self , *args , **kwargs ) :
        t = self.tf1()
        return t.draw( *args , **kwargs )


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


# ==============================================================================
## create function object 
def  _funobj0_ ( self ) :
    """Create function object 
    """
    if hasattr ( self , '_bfit' ) : return self._bfit
    self._bfit = H_fit( self )
    return self._bfit

# ==============================================================================
## create function object 
def  _funobjN_ ( self ) :
    """Create function object 
    """
    if hasattr ( self , '_bfit' ) : return self._bfit
    self._bfit = H_Nfit( self )
    return self._bfit

# ==============================================================================
## draw spline object
def _sp_draw_   ( self , opts = '' ) :
    """Draw spline object
    >>> spline = ...
    >>> spline.draw() 
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
            Ostap.Math.Monotonic       ,
            Ostap.Math.Convex           ,
            Ostap.Math.ConvexOnly       ,            
            Ostap.Math.PositiveSpline   ,
            Ostap.Math.MonotonicSpline ,
            Ostap.Math.ConvexSpline     , 
            Ostap.Math.ConvexOnlySpline ) : t.funobj = _funobjN_
            
# =============================================================================
## construct helper class 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_as_fun_ ( self , func = lambda s : s.value () , *args , **kwargs ) :
    """Construct the function from the histogram
    >>> histo = ...
    >>> fun   = histo.asFunc() 
    """
    return H1Func ( self , func , *args , **kwargs )

# =============================================================================
## construct helper class 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h2_as_fun_ ( self , func = lambda s : s.value () , *args , **kwargs ) :
    """Construct the helper function object from the histogram
    >>> histo = ...
    >>> fun   = histo.asFunc() 
    """
    return H2Func ( self , func , *args , **kwargs )

# =============================================================================
## construct function 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
def _h1_as_tf1_ ( self                           ,
                  func   = lambda s : s.value () , *args , **kwargs ) :
    """Construct the TF1-function from the 1D-histogram
    >>> histo = ...
    >>> fun1  = histo.asTF1 ()
    """
    #
    fun = _h1_as_fun_    ( self , func , *args , **kwargs )
    f1  = fun .tf1  ()
    nb  = self.nbins()
    f1._tmp_fun = fun
    
    if f1.GetNpx() <  1.2 * nb :
        f1.SetNpx ( max ( 100 , 10 * nb ) ) 
    
    f1._funobj = fun  
    f1._histo  = fun._histo
    f1._func   = fun._func
    
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

ROOT.TH1F . integral = _h1_integral_ 
ROOT.TH1D . integral = _h1_integral_
ROOT.TH2F . integral = _h2_integral_ 
ROOT.TH2D . integral = _h2_integral_

# =============================================================================
_decorated_classes_ = (
    Ostap.Math.Bernstein        ,
    Ostap.Math.BernsteinEven    ,
    Ostap.Math.ChebyshevSum     ,
    Ostap.Math.LegendreSum      ,
    Ostap.Math.FourierSum       ,
    Ostap.Math.CosineSum        ,
    Ostap.Math.Polynomial       ,
    Ostap.Math.BSpline          ,
    #
    Ostap.Math.Positive         ,
    Ostap.Math.PositiveEven     , 
    Ostap.Math.Monotonic       ,
    Ostap.Math.Convex           ,
    Ostap.Math.ConvexOnly       ,            
    Ostap.Math.PositiveSpline   ,
    Ostap.Math.MonotonicSpline ,
    Ostap.Math.ConvexSpline     , 
    Ostap.Math.ConvexOnlySpline ,
    #
    ROOT.TH1F ,
    ROOT.TH1D ,
    ROOT.TH2F ,
    ROOT.TH2D ,
    )

_new_methods_ = ( 
    ROOT.TH1F . asTF     ,
    ROOT.TH1D . asTF     ,
    ROOT.TH2F . asTF     ,
    ROOT.TH2D . asTF     ,
    ROOT.TH1F . asTF1    ,
    ROOT.TH1D . asTF1    ,
    ROOT.TH2F . asTF2    ,
    ROOT.TH2D . asTF2    ,
    ROOT.TH1F . asFunc   ,
    ROOT.TH1D . asFunc   ,
    ROOT.TH2F . asFunc   ,
    ROOT.TH2D . asFunc   ,
    #
    ROOT.TH1F . integral ,
    ROOT.TH1D . integral ,
    ROOT.TH2F . integral ,
    ROOT.TH2D . integral ,
    )
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
