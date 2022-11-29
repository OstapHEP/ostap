#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/models.py
#  Module with some useful utilities for simple functions and fit models.
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-12-01
# =============================================================================
"""Module with some useful fit-models"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = (
    'tf1'     , ## convert model/function to TF1 ,
    'tf2'     , ## convert model/function to TF2 ,
    'tf3'     , ## convert model/function to TF3 ,
    'f1_draw' , ## draw 1D-function via conversion to TF1 
    'f2_draw' , ## draw 1D-function via conversion to TF2 
    'f3_draw' , ## draw 1D-function via conversion to TF3 
    )
# =============================================================================
import  ROOT 
from    ostap.core.core        import cpp, Ostap, funID
from    ostap.core.ostap_types import num_types, integer_types
import  ostap.math.polynomials 
import  ostap.math.reduce   
import  ostap.math.derivative  as     D  
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.models' )
else                       : logger = getLogger ( __name__            )
# =============================================================================
# helper adapter for 1D-functions 
class _WO1_ (object)  :
    "Helper adapter for 1D-functions"
    def __init__ ( self , o              ) :        self.__callme = o
    def __call__ ( self , x , pars  = [] ) : return self.__callme ( x [0] )
    @property
    def callme   ( self ) : return self.__callme 
# =============================================================================
# helper adapter for 2D-functions 
class _WO2_ (object)  :
    "Helper adapter for 2D-functions"
    def __init__ ( self , o              ) :        self.__callme =  o 
    def __call__ ( self , x , pars  = [] ) : return self.__callme ( x [0] , x[1] )
    @property
    def callme   ( self ) : return self.__callme 
# =============================================================================
# helper adapter for 3D-functions 
class _WO3_ (object)  :
    "Helper adapter for 3D-functions"
    def __init__ ( self , o              ) :        self.__callme =  o 
    def __call__ ( self , x , pars  = [] ) : return self.__callme ( x [0] , x [1] , x [2] )
    @property
    def callme   ( self ) : return self.__callme 
# =============================================================================
pos_infinity = float('+inf')
neg_infinity = float('-inf')
# =============================================================================
## convert the model/function into TF1
def tf1  ( self                 ,
           xmin  = neg_infinity ,
           xmax  = pos_infinity ,
           **kwargs             ) :
    """Convert the function to TF1    
    >>> obj = ...
    >>> fun = obj.tf1 ( 3.0 , 3.2 )
    >>> fun.Draw() 
    """
    
    npars    = kwargs.pop ( 'npars'    , 0    )
    args     = kwargs.pop ( 'args'     , ()   )
    npx      = kwargs.pop ( 'npx'      , 250  )
    npoints  = kwargs.pop ( 'npoints'  , 250  )
    callme   = kwargs.pop ( 'callable' , self ) 
    title    = kwargs.pop ( 'title'    , None )
    
    if hasattr ( self , '_wo1' ) and callme is not self._wo1.callme :
        del self._wo1 
        
    if not hasattr ( self , '_wo1' ) : self._wo1 = _WO1_ ( callme )
    if not self._wo1                 : self._wo1 = _WO1_ ( callme )
    #
    xmin = float ( xmin )
    xmax = float ( xmax )
    #
    if hasattr ( self , 'xmin'  ) :
        xmn   = self.xmin
        xmin  = max ( float ( xmin ) , float ( xmn () ) if callable ( xmn ) else float ( xmn  ) )
    if hasattr ( self , 'xmax'  ) :
        xmx   = self.xmax
        xmax  = min ( float ( xmax ) , float ( xmx () ) if callable ( xmx ) else float ( xmx ) ) 
    if hasattr ( self , 'npars' ) :
        nps   = self.npars
        npars = max ( npars , nps () if callable ( nps ) else nps )
    #
    assert xmin > neg_infinity, \
          "``xmin''-parameter needs to be specified %s" % xmin
    assert xmax < pos_infinity, \
          "``xmax''-parameter needs to be specified %s" % xmax
    
    _wo = self._wo1 
    fun = ROOT.TF1 ( funID()  , _wo , xmin , xmax , npars, *args )

    if   isinstance ( npx     , integer_types ) and 1 < npx     : fun.SetNpx ( npx     ) 
    elif isinstance ( npoints , integer_types ) and 1 < npoints : fun.SetNpx ( npoints ) 

    if title is None : title = str ( self )            
    fun.SetTitle ( title ) 
    ##
    return fun 

# =============================================================================
## convert the model into TF2
def tf2 ( self ,
          xmin  = neg_infinity ,
          xmax  = pos_infinity ,
          ymin  = neg_infinity ,
          ymax  = pos_infinity ,
          npars = 0            ,
          args  = ()           ,
          npx   = 50           ,
          npy   = 50           , **kwargs ) :
    """Convert the function to TF2
    >>> obj = ...    
    >>> fun = obj.tf2 ( 3.0 , 3.2 , 3.0 , 3.2 )    
    >>> fun.Draw() 
    """
    ##
    callme   = kwargs.pop ( 'callable' , self ) 
    title    = kwargs.pop ( 'title'    , None ) 
    ##
    if not hasattr ( self , '_wo2' ) : self._wo2 = _WO2_ ( callme )
    if not self._wo2                 : self._wo2 = _WO2_ ( callme )
    ##
    xmin = float ( xmin )
    xmax = float ( xmax )
    ymin = float ( ymin )
    ymax = float ( ymax )
    #
    if hasattr ( self , 'xmin'  ) :
        xmn   = self.xmin
        xmin  = max ( float ( xmin ) , float ( xmn () ) if callable ( xmn ) else float ( xmn  ) )
    if hasattr ( self , 'xmax'  ) :
        xmx   = self.xmax
        xmax  = min ( float ( xmax ) , float ( xmx () ) if callable ( xmx ) else float ( xmx ) ) 
    if hasattr ( self , 'ymin'  ) :
        ymn   = self.ymin
        ymin  = max ( float ( ymin ) , float ( ymn () ) if callable ( ymn ) else float ( ymn  ) )
    if hasattr ( self , 'ymax'  ) :
        ymx   = self.ymax
        ymax  = min ( float ( ymax ) , float ( ymx () ) if callable ( ymx ) else float ( ymx  ) )
    if hasattr ( self , 'npars' ) :
        nps   = self.npars
        npars = max ( npars , nps () if callable ( nps ) else nps )

    ##
    assert xmin > neg_infinity, \
           "``xmin''-parameter needs to be specified %s" % xmin
    assert xmax < pos_infinity, \
           "``xmax''-parameter needs to be specified %s" % xmax
    assert ymin > neg_infinity, \
           "``ymin''-parameter needs to be specified %s" % ymin
    assert ymax < pos_infinity, \
           "``ymax''-parameter needs to be specified %s" % ymax
    ##
    _wo = self._wo2
    fun = ROOT.TF2 ( funID ()  , _wo , xmin , xmax , ymin , ymax , npars , *args )
    fun.SetNpx ( npx ) 
    fun.SetNpy ( npy ) 
    #
    if title is None : title = str ( self )
    fun.SetTitle ( title ) 
    #
    if kwargs : logger.warning ("Unused arguments %s" % [ a for a in kwargs ] )  
    return fun 

# =============================================================================
## convert the model into TF3
def tf3 ( self ,
          xmin  = neg_infinity ,
          xmax  = pos_infinity ,
          ymin  = neg_infinity ,
          ymax  = pos_infinity ,
          zmin  = neg_infinity ,
          zmax  = pos_infinity ,
          npars = 0            ,
          args  = ()           , 
          npx   = 25           ,
          npy   = 25           ,
          npz   = 25           , **kwargs ) :
    """Convert the function to TF3
    >>> obj = ...    
    >>> fun = obj.tf3 ( 3.0 , 3.2 , 3.0 , 3.2 , 1 , 2 )    
    >>> fun.Draw() 
    """
    ##
    ##
    callme   = kwargs.pop ( 'callable' , self ) 
    ##
    if not hasattr ( self , '_wo3' ) : self._wo3 = _WO3_ ( callme )
    if not self._wo3                 : self._wo3 = _WO3_ ( callme )
    ##
    xmin = float ( xmin )
    xmax = float ( xmax )
    ymin = float ( ymin )
    ymax = float ( ymax )
    zmin = float ( zmin )
    zmax = float ( zmax )
    ##
    if hasattr ( self , 'xmin'  ) :
        xmn   = self.xmin
        xmin  = max ( float ( xmin ) , float ( xmn () ) if callable ( xmn ) else float ( xmn  ) )
    if hasattr ( self , 'xmax'  ) :
        xmx   = self.xmax
        xmax  = min ( float ( xmax ) , float ( xmx () ) if callable ( xmx ) else float ( xmx ) ) 
    if hasattr ( self , 'ymin'  ) :
        ymn   = self.ymin
        ymin  = max ( float ( ymin ) , float ( ymn () ) if callable ( ymn ) else float ( ymn  ) )
    if hasattr ( self , 'ymax'  ) :
        ymx   = self.ymax
        ymax  = min ( float ( ymax ) , float ( ymx () ) if callable ( ymx ) else float ( ymx  ) )
    if hasattr ( self , 'zmin'  ) :
        zmn   = self.zmin
        zmin  = max ( float ( zmin ) , float ( zmn () ) if callable ( zmn ) else float ( zmn  ) )
    if hasattr ( self , 'zmax'  ) :
        zmx   = self.zmax
        zmax  = min ( float ( zmax ) , float ( zmx () ) if callable ( zmx ) else float ( zmx  ) )
    if hasattr ( self , 'npars' ) :
        nps   = self.npars
        npars = max ( npars , nps () if callable ( nps ) else nps )


    #
    assert xmin > neg_infinity, \
           "``xmin''-parameter needs to be specified %s" % xmin
    assert xmax < pos_infinity, \
           "``xmax''-parameter needs to be specified %s" % xmax
    assert ymin > neg_infinity, \
           "``ymin''-parameter needs to be specified %s" % ymin
    assert ymax < pos_infinity, \
           "``ymax''-parameter needs to be specified %s" % ymax
    assert zmin > neg_infinity, \
           "``zmin''-parameter needs to be specified %s" % zmin
    assert zmax < pos_infinity, \
           "``zmax''-parameter needs to be specified %s" % zmax
    #
    _wo = self._wo3
    fun = ROOT.TF3 ( funID ()  , _wo , xmin , xmax , ymin , ymax , zmin ,  zmax , npars , *args )
    fun.SetNpx ( npx ) 
    fun.SetNpy ( npy ) 
    fun.SetNpy ( npz ) 
    #
    
    if kwargs : logger.warning ("Unused arguments %s" % [ a for a in kwargs ] )
    
    return fun 

positives = ( Ostap.Math.Positive          ,
              Ostap.Math.PositiveEven      , 
              Ostap.Math.Monotonic         , 
              Ostap.Math.Convex            , 
              Ostap.Math.ConvexOnly        , 
              Ostap.Math.PositiveSpline    , 
              Ostap.Math.MonotonicSpline   , 
              Ostap.Math.ConvexSpline      ,
              Ostap.Math.ConvexOnlySpline  ,
              Ostap.Math.ExpoPositive      ,
              Ostap.Math.TwoExpoPositive   ) 

# =============================================================================
## draw the function 
def f1_draw ( self , opts ='' , **kwargs ) :
    """Drawing the function object through conversion to ROOT.TF1    
    >>> fun = ...
    >>> fun.draw()    
    """

    if hasattr ( self , '_tf1' ) and 'callable' in kwargs : del self._tf1
    
    if hasattr ( self , '_tf1' ) and 'xmin'     in kwargs :        
        xmin    = kwargs.get ( 'xmin'    , None )
        if isinstance ( xmin    , num_types     ) and float ( xmin ) != self._tf1.GetXmin () : del self._tf1 

    if hasattr ( self , '_tf1' ) and 'xmax'     in kwargs :        
        xmax    = kwargs.get ( 'xmax'    , None )
        if isinstance ( xmax    , num_types     ) and float ( xmax ) != self._tf1.GetXmax () : del self._tf1 
        
    if hasattr ( self , '_tf1' ) and 'npx'       in kwargs :        
        npx    = kwargs.get ( 'npx'    , None )
        if isinstance ( npx     , integer_types ) and 1 < npx and npx != self._tf1.GetNpx () : del self._tf1
        
    if hasattr ( self , '_tf1' ) and 'npoints'    in kwargs :        
        npx    = kwargs.get ( 'npoints'    , None )
        if isinstance ( npx     , integer_types ) and 1 < npx and npx != self._tf1.GetNpx () : del self._tf1 
        
        
    if not hasattr ( self , '_tf1'  ) :
        
        self._tf1        =  tf1 ( self , **kwargs )
        
        if type ( self ) in positives and not 'xmin' in kwargs :
            self._tf1.SetMinimum(0)
            
    kwargs.pop ( 'xmin'     , None )
    kwargs.pop ( 'xmax'     , None )
    kwargs.pop ( 'npars'    , None ) 
    kwargs.pop ( 'args'     , None )
    kwargs.pop ( 'npx'      , None )
    kwargs.pop ( 'npoints'  , None )
    kwargs.pop ( 'callable' , None ) 

    return self._tf1.draw ( opts , **kwargs )

# =============================================================================
## draw the function 
def f2_draw ( self , opts ='' , **kwargs ) :
    """Drawing the function object through conversion to ROOT.TF2    
    >>> fun = ...
    >>> fun.draw()    
    """
    
    if hasattr ( self , '_tf2' ) and 'xmin'     in kwargs :        
        xmin    = kwargs.get ( 'xmin'    , None )
        if isinstance ( xmin    , num_types     ) and float ( xmin ) != self._tf2.GetXmin () : del self._tf2 

    if hasattr ( self , '_tf2' ) and 'xmax'     in kwargs :        
        xmax    = kwargs.get ( 'xmax'    , None )
        if isinstance ( xmax    , num_types     ) and float ( xmax ) != self._tf2.GetXmax () : del self._tf2 

    if hasattr ( self , '_tf2' ) and 'ymin'     in kwargs :        
        ymin    = kwargs.get ( 'ymin'    , None )
        if isinstance ( ymin    , num_types     ) and float ( ymin ) != self._tf2.GetYmin () : del self._tf2 

    if hasattr ( self , '_tf2' ) and 'ymax'     in kwargs :        
        ymax    = kwargs.get ( 'ymax'    , None )
        if isinstance ( ymax    , num_types     ) and float ( ymax ) != self._tf2.GetYmax () : del self._tf2 

    if not hasattr ( self , '_tf2'  ) :

        self._tf2        =  tf2 ( self , **kwargs )
        
    xmin  = kwargs.pop ( 'xmin'  , neg_infinity )
    xmax  = kwargs.pop ( 'xmax'  , pos_infinity )
    ymin  = kwargs.pop ( 'ymin'  , neg_infinity )
    ymax  = kwargs.pop ( 'ymax'  , pos_infinity )
    npars = kwargs.pop ( 'npars' , 0  ) 
    args  = kwargs.pop ( 'args'  , () )
    
    return self._tf2.draw ( opts , **kwargs )

# =============================================================================
## draw the function 
def f3_draw ( self , opts ='' , **kwargs ) :
    """Drawing the function object through conversion to ROOT.TF3    
    >>> fun = ...
    >>> fun.draw()    
    """
    
    if hasattr ( self , '_tf3' ) and 'xmin'     in kwargs :        
        xmin    = kwargs.get ( 'xmin'    , None )
        if isinstance ( xmin    , num_types     ) and float ( xmin ) != self._tf3.GetXmin () : del self._tf3 

    if hasattr ( self , '_tf3' ) and 'xmax'     in kwargs :        
        xmax    = kwargs.get ( 'xmax'    , None )
        if isinstance ( xmax    , num_types     ) and float ( xmax ) != self._tf3.GetXmax () : del self._tf3 

    if hasattr ( self , '_tf3' ) and 'ymin'     in kwargs :        
        ymin    = kwargs.get ( 'ymin'    , None )
        if isinstance ( ymin    , num_types     ) and float ( ymin ) != self._tf3.GetYmin () : del self._tf3 

    if hasattr ( self , '_tf3' ) and 'ymax'     in kwargs :        
        ymax    = kwargs.get ( 'ymax'    , None )
        if isinstance ( ymax    , num_types     ) and float ( ymax ) != self._tf3.GetYmax () : del self._tf3 

    if hasattr ( self , '_tf3' ) and 'zmin'     in kwargs :        
        zmin    = kwargs.get ( 'zmin'    , None )
        if isinstance ( zmin    , num_types     ) and float ( zmin ) != self._tf3.GetZmin () : del self._tf3 

    if hasattr ( self , '_tf3' ) and 'ymax'     in kwargs :        
        zmax    = kwargs.get ( 'zmax'    , None )
        if isinstance ( zmax    , num_types     ) and float ( zmax ) != self._tf3.GetZmax () : del self._tf3 


    if not hasattr ( self , '_tf3'  ) :

        self._tf3        = tf3 ( self , **kwargs )
    
    xmin  = kwargs.pop ( 'xmin' , neg_infinity )
    xmax  = kwargs.pop ( 'xmax' , pos_infinity )
    ymin  = kwargs.pop ( 'ymin' , neg_infinity )
    ymax  = kwargs.pop ( 'ymax' , pos_infinity )
    zmin  = kwargs.pop ( 'zmin' , neg_infinity )
    zmax  = kwargs.pop ( 'zmax' , pos_infinity )
    npars = kwargs.pop ( 'npars' , 0  ) 
    args  = kwargs.pop ( 'args' , () )
    
    return self._tf3.draw ( opts , **kwargs )

# =============================================================================
## get the regular complex value for amplitude 
def _amp_ ( self , x , *a ) :
    """ Get the complex value for amplitude
    >>> fun
    >>> a = fun.amp ( x )    
    """
    v = self.amplitude ( x , *a )
    return complex( v ) 


# ==============================================================================
## get the real part of amplitude
#  @code
#  f = ...
#  r = f.amp_real ( x ) 
# @endcode
def _amp_re_ ( self, x , *a ) :
    """Get   the real part of amplitude
    >>> f = ...
    >>> r = f.amp_real ( x )     
    """
    return complex ( self.amp ( x , *a ) ).real 
    
# ==============================================================================
## get the imaginary part of amplitude
#  @code
#  f = ...
#  r = f.amp_imag ( x ) 
# @endcode
def _amp_im_ ( self, x , *a ) :
    """Get   the imaginary part of amplitude
    >>> f = ...
    >>> r = f.amp_imag ( x )     
    """    
    return complex ( self.amp ( x , *a ) ).imag

# ==============================================================================
## get the phase of amplitude
#  @code
#  f = ...
#  r = f.amp_phase ( x ) 
# @endcode
def _amp_phase_ ( self, x , *a ) :
    """Get   the phase of amplitude
    >>> f = ...
    >>> r = f.amp_phase ( x )     
    """    
    return cmath.phase (  complex ( self.amp ( x , *a ) ) )

# =============================================================================
## build the Argand diagram/graph 
#  @code
#  f = ...
#  g = f.argand ( x , xmin = ...  , xmax = ... ) 
#  @@endcode 
def _amp_argand_ ( self , xmin , xmax , npx =  500 , args = () ) : 
    """Build the Argand diagram/graph 
    >>> f = ...
    >>> g = f.argand ( x , xmin = ...  , xmax = ... ) 
    """
    import ostap.histos.graphs
    from ostap.utils.utils import vrange
    g = ROOT.TGraph  ( npx + 1 ) 
    for i , x in enumerate ( vrange ( xmin , xmax , npx ) ) :
        a = complex ( self.amp ( x , *args ) )        
        g[i] = a.real, a.imag
    return g


## Ostap.Math.LASS            . amp = _amp_
## Ostap.Math.LASS23L         . amp = _amp_
## Ostap.Math.Bugg23L         . amp = _amp_
## Ostap.Math.Flatte          . amp = _amp_
## Ostap.Math.Flatte2         . amp = _amp_
## Ostap.Math.Flatte23L       . amp = _amp_
## Ostap.Math.BW              . amp = _amp_
## Ostap.Math.Swanson         . amp = _amp_

for m in  ( Ostap.Math.Flatte    ,
            Ostap.Math.BW        ,
            Ostap.Math.BWPS      ,
            Ostap.Math.BW3L      ,
            Ostap.Math.ChannelBW ) :
    m.amp       = _amp_ 
    m.amp_real  = _amp_re_ 
    m.amp_imag  = _amp_im_ 
    m.amp_phase = _amp_phase_ 
    m.argand    = _amp_argand_ 
    
# =============================================================================
## get min/max values for bernstein polynomials
#  @code
#  p = ...
#  mn,mx = p.minmax()
#  @endcode
#  The values are guaranteed that
#  mn <= p(x) <= mx for all   xmin <= x <= xmax 
def _b_minmax_ ( bp ) :
    """Get min/max values for bernstein polynomials
    
    >>> p = ...
    >>> mn,mx = p.minmax()

    The values are such that: mn <= p(x) <= mx for all x_min<=x<x_max
    """
    b    = bp.bernstein() 
    pars = b .pars()
    mn   = min ( pars )
    mx   = max ( pars )
    return mn , mx 

# ==============================================================================
## get the maximal value for bernstein polynomial:
#  @code
#  p  = ...
#  mx = p.max()
#  @endcode
#  The values are guaranteed that
#  p(x) <= mx for all   xmin <= x <= xmax 
def _b_max_ ( bp ) :
    """Get max values for bernstein polynomials
    
    >>> p = ...
    >>> mx = p.max()

    The value is such that: p(x) <= mx  for all x_min<=x<x_max
    """
    b    = bp.bernstein() 
    pars = b .pars()
    return max ( pars )

# ==============================================================================
## get the minimal value for bernstein polynomial:
#  @code
#  p  = ...
#  mn = p.min()
#  @endcode
#  The values are guaranteed that
#  mn <= p(x) for all   xmin <= x <= xmax 
def _b_min_ ( bp ) :
    """Get min values for bernstein polynomials
    
    >>> p  = ...
    >>> mn = p.min()

    The  value is such that: mn <= p(x)  for all x_min<=x<x_max
    """
    b    = bp.bernstein() 
    pars = b .pars()
    return min ( pars )

for t in ( Ostap.Math.Bernstein     ,
           Ostap.Math.BernsteinEven ) :
    if not hasattr ( t , 'min'    ) : t.min    = _b_min_
    if not hasattr ( t , 'max'    ) : t.max    = _b_max_ 
    if not hasattr ( t , 'minmax' ) : t.minmax = _b_minmax_

# =============================================================================

# =============================================================================
## get min/max values for derived bernstein polynomials
#  @code
#  p = ...
#  mn,mx = p.minmax()
#  @endcode
#  The values are guaranteed that
#  mn <= p(x) <= mx for all   xmin <= x <= xmax 
def _p_minmax_ ( p ) :
    """Get min/max values for derived bernstein polynomials
    
    >>> p = ...
    >>> mn,mx = p.minmax()

    The values are such that: mn <= p(x) <= mx  for all x_min<=x<x_max
    """
    b    = p .bernstein() 
    pars = b .pars()
    mn   = min ( pars )
    mx   = max ( pars )
    return  max ( mn , 0 ) , mx 


# ==============================================================================
## get the minimal value for derived bernstein polynomial:
#  @code
#  p  = ...
#  mn = p.min()
#  @endcode
#  The values are guaranteed that
#  mn <= p(x) for all   xmin <= x <= xmax 
def _p_min_ ( bp ) :
    """Get min values for derived bernstein polynomials
    
    >>> p  = ...
    >>> mn = p.min()

    The  value is such that: mn <= p(x)  for all x_min<=x<x_max
    """
    b    = bp.bernstein() 
    pars = b .pars()
    return max ( min ( pars ) , 0 )

for t in ( Ostap.Math.Positive      ,
           Ostap.Math.Monotonic    ,
           Ostap.Math.Convex        ,
           Ostap.Math.ConvexOnly    ) :
    
    if not hasattr ( t , 'min'    ) : t.min    = _p_min_
    if not hasattr ( t , 'max'    ) : t.max    = _b_max_   ## ATTENTION: "b" is here
    if not hasattr ( t , 'minmax' ) : t.minmax = _p_minmax_


# =============================================================================
## try to get max-value using existing mode function
#  @code
#  f = ...
#  mx = f.max()
#  @endcode
def _f_max_mode_ ( f ) :
    """Get the max-value (using exising mode function)
    >>> f = ...
    >>> mx = f.max()    
    """
    return f ( f.mode() ) 
 
# =============================================================================
## make 1D- numerical integration
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_1D ( func , xmin , xmax , *args , **kwargs ) :
    """Make 1D numerical integration 
    
    >>> func = ...
    >>> print func.sp_integrate ( -10 , 10 )    
    """    
    from ostap.math.integral import integral as _integral 
    return _integral ( func , xmin , xmax , *args , **kwargs )

# =============================================================================
## make 2D numerical integration 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_2D ( func  ,
                      xmin  , xmax ,
                      ymin  , ymax , *args , **kwargs ) :
    """Make 2D numerical integration

    >>> func = ...  ## func ( x , y )
    ##                            xmin , xmax , ymin , ymax 
    >>> print func.sp_integrate ( -10  , 10   , -20  , 20   ) 

    """
    from ostap.math.integral import integral2 as _integral2 
    return _integral2 ( func  ,
                        xmin  , xmax ,
                        ymin  , ymax ,
                        *args , **kwargs )

# =============================================================================
## make 1D numerical integration 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_2Dx ( func  , y    ,
                       xmin  , xmax , *args , **kwargs ) :
    """Make 1D numerical integration over x-axis 
    
    >>> func = ...  ## func ( x , y )
    ##                              y   , xmin , xmax 
    >>> print func.sp_integrate_x ( 0.5 , -20  , 20   ) 
    
    """
    def _func_ ( p , *args ) :
        return func ( p , y , *args )
    
    from ostap.math.integral import integral as _integral 
    return _integral ( _func_ ,
                       xmin   , xmax ,
                       *args  , **kwargs )

# =============================================================================
## make 1D numerical integration 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_2Dy ( func  , x    ,
                       ymin  , ymax , *args , **kwargs ) :
    """Make 1D numerical integration over y-axis 
    
    >>> func = ...  ## func ( x , y )
    ##                              x   , ymin , ymax 
    >>> print func.sp_integrate_y ( 0.5 , -20  , 20   ) 
    
    """
    def _func_ ( p , *args ) :
        return func ( x , p , *args )
    
    from ostap.math.integral import integral as _integral 
    return _integral ( _func_ ,
                       ymin   , ymax ,
                       *args  , **kwargs )


# =============================================================================
## make 3D numerical integration 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_3D ( func  ,
                      xmin  , xmax ,
                      ymin  , ymax ,
                      zmin  , zmax , *args , **kwargs ) :
    """Make 3D numerical integration
    
    >>> func = ...  ## func ( x , y , z )
    ##                            xmin , xmax , ymin , ymax   zmin zmax 
    >>> print func.sp_integrate ( -10  , 10   , -20  , 20   , -1 ,   1 ) 
    """
    from ostap.math.integral import integral2 as _integral3 
    return _integral3 ( func  ,
                        xmin  , xmax ,
                        ymin  , ymax ,
                        zmin  , zmax ,
                        *args , **kwargs )

# =============================================================================
## make 1D numerical integration 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_3Dx ( func  ,
                       y     , z   ,  
                       xmin  , xmax , *args , **kwargs ) :
    """Make 1D numerical integration over x-axis 
    
    >>> func = ...  ## func ( x , y , z )
    ##                              y   ,   z , xmin , xmax 
    >>> print func.sp_integrate_x ( 0.5 , 0.1 , -20  , 20   ) 
    
    """
    def _func_ ( p , *args ) :
        return func ( p , y , z , *args )
    
    from ostap.math.integral import integral as _integral 
    return _integral ( _func_ ,
                       xmin   , xmax ,
                       *args  , **kwargs )

# =============================================================================
## make 1D numerical integration 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_3Dy ( func  ,
                       x     ,  z   , 
                       ymin  , ymax , *args , **kwargs ) :
    """Make 1D numerical integration over y-axis 
    
    >>> func = ...  ## func ( x , y , z )
    ##                              x   , z   , ymin , ymax 
    >>> print func.sp_integrate_y ( 0.5 , 0.1 , -20  , 20   ) 
    
    """
    def _func_ ( p , *args ) :
        return func ( x , p , z , *args )
    
    from ostap.math.integral import integral as _integral 
    return _integral ( _func_ ,
                       ymin   , ymax ,
                       *args  , **kwargs )

# =============================================================================
## make 1D numerical integration 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_3Dz ( func  ,
                       x     , y    , 
                       zmin  , zmax , *args , **kwargs ) :
    """Make 1D numerical integration over z-axis 
    
    >>> func = ...  ## func ( x , y , z )
    ##                              x   , y   , zmin , zmax 
    >>> print func.sp_integrate_y ( 0.5 , 0.1 , -20  , 20   ) 
    
    """
    def _func_ ( p , *args ) :
        return func ( x , y , p , *args )
    
    from ostap.math.integral import integral as _integral 
    return _integral ( _func_ ,
                       zmin   , zmax ,
                       *args  , **kwargs )


# =============================================================================
## make 2D numerical integration 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_3Dxy ( func  ,
                        z     ,
                        xmin  , xmax ,
                        ymin  , ymax , *args , **kwargs ) :
    """Make 2D numerical integration

    >>> func = ...  ## func ( x , y , z )
    ##                            z , xmin , xmax , ymin , ymax 
    >>> print func.sp_integrate_xy ( 0.5 , -10  , 10   , -20  , 20   ) 

    """
    def _func_ ( p1 , p2 , *args ) :
        return  func ( p1 , p2 , z , *args )
    
    from ostap.math.integral import integral2 as _integral2 
    return _integral2 ( func  ,
                        xmin  , xmax ,
                        ymin  , ymax ,
                        *args , **kwargs )


# =============================================================================
## make 2D numerical integration 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_3Dxz ( func  ,
                        y     ,
                        xmin  , xmax ,
                        zmin  , zmax , *args , **kwargs ) :
    """Make 2D numerical integration

    >>> func = ...  ## func ( x , y , z )
    ##                                 y , xmin , xmax , zmin , zmax 
    >>> print func.sp_integrate_xz ( 0.5 , -10  , 10   , -20  , 20   ) 

    """
    def _func_ ( p1 , p2 , *args ) :
        return  func ( p1 , y , p2 , *args )
    
    from ostap.math.integral import integral2 as _integral2 
    return _integral2 ( func  ,
                        xmin  , xmax ,
                        ymin  , ymax ,
                        *args , **kwargs )


# =============================================================================
## make 2D numerical integration 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_3Dyz ( func  ,
                        x     ,
                        ymin  , ymax ,
                        zmin  , zmax , *args , **kwargs ) :
    """Make 2D numerical integration

    >>> func = ...  ## func ( x , y , z )
    ##                                 x , ymin , ymax , zmin , zmax 
    >>> print func.sp_integrate_yz ( 0.5 , -10  , 10   , -20  , 20   ) 

    """
    def _func_ ( p1 , p2 , *args ) :
        return  func ( x , p1 , p2 , *args )
    
    from ostap.math.integral import integral2 as _integral2 
    return _integral2 ( func  ,
                        xmin  , xmax ,
                        ymin  , ymax ,
                        *args , **kwargs )

# =============================================================================
## make 1D numerical integration
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_1D_ ( pdf , xmin , xmax , *args , **kwargs ) :
    """Make 1D numerical integration over the PDF using SciPy
    """
    if hasattr ( pdf , 'setPars' ) : pdf.setPars() 
    func = pdf.function()
    return func.sp_integrate_1D ( xmin , xmax , *args , **kwargs ) 

# =============================================================================
## make 2D numerical integration
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_2D_ ( pdf   ,
                      xmin  , xmax ,
                      ymin  , ymax , *args , **kwargs ) :
    """ Make 3D numerical integration over the PDF
    """
    if hasattr ( pdf , 'setPars' ) : pdf.setPars() 
    func = pdf.function()
    return func.sp_integrate_2D ( xmin , xmax , ymin , ymax , *args , **kwargs ) 

# =============================================================================
## make 3D numerical integration 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_3D_  ( pdf   ,
                        xmin  , xmax ,
                        ymin  , ymax ,
                        zmin  , zmax ,
                        *args , **kwargs ) :
    """ Make 3D numerical integration over the PDF 
    """
    if hasattr ( pdf , 'setPars' ) : pdf.setPars() 
    func = pdf.function()
    return func.sp_integrate_3D ( xmin , xmax , ymin , ymax , zmin , zmax , *args , **kwargs ) 


from ostap.stats.moments import moment           as sp_moment
from ostap.stats.moments import central_moment   as sp_central_moment
from ostap.stats.moments import mean             as sp_mean
from ostap.stats.moments import variance         as sp_variance
from ostap.stats.moments import rms              as sp_rms 
from ostap.stats.moments import median           as sp_median
from ostap.stats.moments import quantile         as sp_quantile
from ostap.stats.moments import mode             as sp_mode 
from ostap.stats.moments import width            as sp_width
from ostap.stats.moments import fwhm             as sp_fwhm
from ostap.stats.moments import cl_symm          as sp_cl_symm
from ostap.stats.moments import cl_asymm         as sp_cl_asymm

# =============================================================================
## helper function to delegate some methods/attributes to TF1
#  @code
#  f = ...
#  f.SetLineColor(4) ## delegate to TF1
#  f.SetLineWidth(2) ## delegate to TF1
#  @endcode 
def _tf1_getattr_ ( self , attr ) :
    """Delegate some methods/attributes to TF1
    >>> f = ...
    >>> f.SetLineColor(4) ## delegate to TF1
    >>> f.SetLineWidth(2) ## delegate to TF1
    """
    if  hasattr ( ROOT.TF1 , attr ) and hasattr  ( self , '_tf1' ) and self._tf1  :
        return getattr  ( self._tf1 , attr   )
    
    raise AttributeError("Can't get attribute: %s" % attr )

# =============================================================================
## helper function to delegate some methods/attributes to TF2
#  @code
#  f = ...
#  f.SetLineColor(4) ## delegate to TF2
#  f.SetLineWidth(2) ## delegate to TF2
#  @endcode 
def _tf2_getattr_ ( self , attr ) :
    """Delegate some methods/attributes to TF2
    >>> f = ...
    >>> f.SetLineColor(4) ## delegate to TF2
    >>> f.SetLineWidth(2) ## delegate to TF2
    """
    if  hasattr ( ROOT.TF2 , attr ) and hasattr  ( self , '_tf2' ) and self._tf2 :
        return getattr  ( self._tf2 , attr   )
    
    raise AttributeError("Can't get attribute: %s" % attr )

# =============================================================================
## helper function to delegate some methods/attributes to TF3
#  @code
#  f = ...
#  f.SetLineColor(4) ## delegate to TF3
#  f.SetLineWidth(2) ## delegate to TF3
#  @endcode 
def _tf3_getattr_ ( self , attr ) :
    """Delegate some methods/attributes to TF2
    >>> f = ...
    >>> f.SetLineColor(4) ## delegate to TF2
    >>> f.SetLineWidth(2) ## delegate to TF2
    """
    if  hasattr ( ROOT.TF3 , attr ) and hasattr  ( self , '_tf3' ) and self._tf3 :
        return getattr  ( self._tf3 , attr   )
    
    raise AttributeError("Can't get attribute: %s" % attr )


from ostap.math.minimize   import sp_minimum_1D, sp_maximum_1D 
from ostap.math.rootfinder import sp_solve 

# =============================================================================
## decorate 1D-models/functions 
# =============================================================================
for model in ( Ostap.Math.Chebyshev              ,
               Ostap.Math.ChebyshevU             ,
               Ostap.Math.Legendre               ,
               Ostap.Math.Hermite                ,
               Ostap.Math.Bernstein              ,
               Ostap.Math.BernsteinEven          ,
               Ostap.Math.ChebyshevSum           ,
               Ostap.Math.LegendreSum            ,
               Ostap.Math.HermiteSum             ,
               Ostap.Math.FourierSum             ,
               Ostap.Math.CosineSum              ,
               Ostap.Math.Polynomial             ,               
               Ostap.Math.Positive               ,
               Ostap.Math.PositiveEven           ,
               Ostap.Math.Monotonic             ,
               Ostap.Math.Convex                 ,
               Ostap.Math.ConvexOnly             ,
               Ostap.Math.BifurcatedGauss        ,
               Ostap.Math.DoubleGauss            ,
               Ostap.Math.Bukin                  ,
               Ostap.Math.Novosibirsk            ,
               Ostap.Math.SinhAsinh              ,
               Ostap.Math.CrystalBall            ,
               Ostap.Math.Needham                ,
               Ostap.Math.CrystalBallDoubleSided ,
               Ostap.Math.GramCharlierA          ,
               Ostap.Math.PhaseSpace2            ,
               Ostap.Math.sPhaseSpace2           ,
               Ostap.Math.PhaseSpace3            ,
               Ostap.Math.PhaseSpace3s           ,
               Ostap.Math.PhaseSpaceLeft         ,
               Ostap.Math.PhaseSpaceRight        ,
               Ostap.Math.PSDalitz               ,
               Ostap.Math.PhaseSpaceNL           ,
               Ostap.Math.PhaseSpace23L          ,
               Ostap.Math.BreitWigner            ,
               Ostap.Math.BW                     ,
               Ostap.Math.BreitWignerMC          ,
               Ostap.Math.Rho0                   ,
               Ostap.Math.Kstar0                 ,
               Ostap.Math.Phi0                   , 
               Ostap.Math.BWPS                  ,
               Ostap.Math.BW3L                  ,
               ## Ostap.Math.Rho0FromEtaPrime       ,
               Ostap.Math.Flatte                 ,
               ## Ostap.Math.Flatte2                ,
               ## Ostap.Math.LASS                   ,
               ## Ostap.Math.LASS23L                ,
               ## Ostap.Math.Bugg23L                ,
               ## Ostap.Math.BW23L                  ,
               ## Ostap.Math.Flatte23L              ,
               ## Ostap.Math.Gounaris23L            ,
               Ostap.Math.StudentT               ,
               Ostap.Math.BifurcatedStudentT     ,
               Ostap.Math.Voigt                  ,
               Ostap.Math.PseudoVoigt            ,
               Ostap.Math.Logistic               ,
               #
               Ostap.Math.GenGaussV1             ,
               Ostap.Math.GenGaussV2             ,
               Ostap.Math.SkewGauss              ,
               Ostap.Math.ExGauss                ,
               Ostap.Math.NormalLaplace          ,
               Ostap.Math.GammaDist              ,
               Ostap.Math.GenGammaDist           ,
               Ostap.Math.Amoroso                ,
               Ostap.Math.LogGammaDist           ,
               Ostap.Math.Log10GammaDist         ,
               Ostap.Math.LogGamma               ,
               Ostap.Math.BetaPrime              ,
               Ostap.Math.Landau                 ,
               Ostap.Math.JohnsonSU              ,
               Ostap.Math.Atlas                  ,
               Ostap.Math.Sech                   ,
               Ostap.Math.Losev                  ,
               ## Ostap.Math.Swanson                ,
               Ostap.Math.Argus                  ,
               Ostap.Math.GenArgus               ,
               Ostap.Math.Slash                  ,
               Ostap.Math.AsymmetricLaplace      ,
               Ostap.Math.Tsallis                ,
               Ostap.Math.QGSM                   ,
               Ostap.Math.TwoExpos               ,
               Ostap.Math.DoubleGauss            ,
               Ostap.Math.Gumbel                 ,
               Ostap.Math.Weibull                ,
               Ostap.Math.QGaussian              ,
               Ostap.Math.KGaussian              ,
               Ostap.Math.RaisingCosine          ,
               Ostap.Math.Sigmoid                ,
               Ostap.Math.Hyperbolic             ,
               Ostap.Math.GenHyperbolic          ,
               Ostap.Math.Das                    ,
               Ostap.Math.Rice                   ,
               Ostap.Math.GenInvGauss            ,
               Ostap.Math.PearsonIV              ,
               Ostap.Math.SkewGenT               ,
               #
               Ostap.Math.BSpline                , 
               Ostap.Math.PositiveSpline         ,
               Ostap.Math.MonotonicSpline        ,
               Ostap.Math.ConvexSpline           ,
               Ostap.Math.ConvexOnlySpline       ,
               #
               Ostap.Math.BernsteinDualBasis     ,
               ## interpolators  
               Ostap.Math.Neville                ,
               Ostap.Math.Lagrange               ,
               Ostap.Math.Newton                 ,
               Ostap.Math.Berrut1st              ,
               Ostap.Math.Berrut2nd              ,
               Ostap.Math.FloaterHormann         ,
               Ostap.Math.Barycentric            ,
               Ostap.Math.Thiele                 ,
               ##
               Ostap.Math.GammaBW3               , 
               ## helper stufff
               ## Ostap.Functions.PyCallable     , 
               Ostap.Math.Piecewise              , 
               Ostap.Math.ChebyshevApproximation ,
               D.Derivative                      ,
               D.Derivative1                     ,
               D.Derivative2                     ,
               D.Derivative3                     ,
               D.Derivative4                     ,
               D.Derivative5                     ,
               D.Derivative6                     ,
               ##
               Ostap.Math.Histo1D                ,
               ## 
               Ostap.Math.Multiply               ,
               Ostap.Math.Divide                 ,
               Ostap.Math.Const                  ,
               Ostap.Math.Id                     ,
               Ostap.Math.Linear                 ,
               Ostap.Math.Compose                ,
               Ostap.Math.Multiply               ,
               Ostap.Math.Divide                 ,
               Ostap.Math.Sum                    ,
               Ostap.Math.Moebius                ,
               Ostap.Math.Step                   ,
               Ostap.Math.Min                    ,
               Ostap.Math.Max                    ,
               Ostap.Math.Apply                  ,
               Ostap.Math.KramersKronig          ,
               Ostap.Math.Abs                    ,
               Ostap.Math.Sqrt                   ,
               Ostap.Math.Cbrt                   ,
               Ostap.Math.Exp                    ,
               Ostap.Math.Log                    ,
               Ostap.Math.Log10                  ,
               Ostap.Math.Erf                    ,
               Ostap.Math.Erfc                   ,
               Ostap.Math.TGamma                 ,
               Ostap.Math.LGamma                 ,
               Ostap.Math.Sin                    ,
               Ostap.Math.Cos                    ,
               Ostap.Math.Tan                    ,
               Ostap.Math.ASin                   ,
               Ostap.Math.ACos                   ,
               Ostap.Math.ATan                   ,
               Ostap.Math.Sinh                   ,
               Ostap.Math.Cosh                   ,
               Ostap.Math.Tanh                   ,
               Ostap.Math.ASinh                  ,
               Ostap.Math.ACosh                  ,
               Ostap.Math.ATanh                  ,
               Ostap.Math.Pow                    ,
               ##
               Ostap.Math.CutOffGauss            ,
               Ostap.Math.CutOffStudent          ,
               ##
               ) :
    model.tf1          =  tf1 
    model.sp_integrate = sp_integrate_1D
    model.__getattr__  = _tf1_getattr_
    model.draw         =  f1_draw
    
    if not hasattr ( model , 'max' ) : 
        if hasattr ( model , 'mode' ) :
            model.max = _f_max_mode_ 

    if not hasattr ( model , 'mean'             ) : model.mean             = sp_mean 
    if not hasattr ( model , 'variance'         ) : model.variance         = sp_variance 
    if not hasattr ( model , 'rms'              ) : model.rms              = sp_rms  
    if not hasattr ( model , 'median'           ) : model.median           = sp_median
    if not hasattr ( model , 'mode'             ) : model.mode             = sp_mode 
    if not hasattr ( model , 'width'            ) : model.width            = sp_width
    if not hasattr ( model , 'get_width'        ) : model.get_width        = sp_width
    if not hasattr ( model , 'fwhm'             ) : model.fwhm             = sp_fwhm
    if not hasattr ( model , 'get_fwhm'         ) : model.get_fwhm         = sp_fwhm
    if not hasattr ( model , 'FWHM'             ) : model.FWHM             = sp_fwhm
    if not hasattr ( model , 'get_FWHM'         ) : model.get_FWHM         = sp_fwhm
    if not hasattr ( model , 'moment'           ) : model.moment           = sp_moment
    if not hasattr ( model , 'central_moment'   ) : model.central_moment   = sp_central_moment
    if not hasattr ( model , 'quantile'         ) : model.quantile         = sp_quantile
    if not hasattr ( model , 'cl_symm'          ) : model.cl_symm          = sp_cl_symm
    if not hasattr ( model , 'cl_asymm'         ) : model.cl_asymm         = sp_cl_asymm
    
    if sp_minimum_1D and not hasattr ( model , 'minimum' ) : model.minimum = sp_minimum_1D
    if sp_maximum_1D and not hasattr ( model , 'maximum' ) : model.maximum = sp_maximum_1D
    if sp_solve      and not hasattr ( model , 'solve'   ) : model.solve   = sp_solve


# =======================================================================================
## Special ``getattr'' for Bernstein dual basis functions: delegate the stuff to
#  the underlying bernstein polynomial
def _bdb_getattr_ ( self ,  attr ) :
    """Special ``getattr'' for Bernstein dual basis functions:
    - delegate the stuff to the underlying Bernstein polynomial
    """
    b = self.bernstein()
    return getattr ( b , attr )
Ostap.Math.BernsteinDualBasis.__getattr__ = _bdb_getattr_


## add some drawing method for some shapes 
for model in ( Ostap.Math.Bernstein         ,
               Ostap.Math.BernsteinEven     , 
               Ostap.Math.Positive          ,
               Ostap.Math.PositiveEven      ,
               Ostap.Math.Monotonic         ,
               Ostap.Math.Convex            ,
               Ostap.Math.ConvexOnly        ,
               Ostap.Math.ChebyshevSum      ,               
               Ostap.Math.LegendreSum       ,               
               Ostap.Math.HermiteSum        ,               
               Ostap.Math.FourierSum        ,               
               Ostap.Math.CosineSum         ,               
               Ostap.Math.Polynomial        ,               
               Ostap.Math.ExpoPositive      , 
               Ostap.Math.TwoExpoPositive   , 
               ##
               Ostap.Math.BSpline           ,
               Ostap.Math.MonotonicSpline   ,
               Ostap.Math.PositiveSpline    , 
               Ostap.Math.ConvexSpline      , 
               Ostap.Math.ConvexOnlySpline  ,
               ##
               Ostap.Math.Const             ,
               Ostap.Math.Id                ,
               
               ) : 
    
    model.draw = f1_draw
    model.Draw = f1_draw


# =============================================================================
def _f_print_ ( self , typ = '' ) :
    if not typ : typ = str(type(self))
    return '%s(%s,%s,%s)' % ( typ ,  self.pars() , self.xmin() , self.xmax() )

Ostap.Math.LegendreSum   .__str__  = lambda s : _f_print_ ( s , 'LegendreSum'   )
Ostap.Math.ChebyshevSum  .__str__  = lambda s : _f_print_ ( s , 'ChebyshevSum'  )
Ostap.Math.Polynomial    .__str__  = lambda s : _f_print_ ( s , 'Polynomial'    )
Ostap.Math.Bernstein     .__str__  = lambda s : _f_print_ ( s , 'Bernstein'     )
Ostap.Math.BernsteinEven .__str__  = lambda s : _f_print_ ( s , 'BernsteinEven' )
Ostap.Math.Positive      .__str__  = lambda s : _f_print_ ( s , 'Positive'      )
Ostap.Math.PositiveEven  .__str__  = lambda s : _f_print_ ( s , 'PositiveEven'  )
Ostap.Math.Convex        .__str__  = lambda s : _f_print_ ( s , 'Convex'        ) 
Ostap.Math.ConvexOnly    .__str__  = lambda s : _f_print_ ( s , 'ConvexOnly'    )
Ostap.Math.Monotonic     .__str__  = lambda s : _f_print_ ( s , 'Monotonic'    )
Ostap.Math.FourierSum    .__str__  = lambda s : _f_print_ ( s , 'FourierSum'    )
Ostap.Math.CosineSum     .__str__  = lambda s : _f_print_ ( s , 'CosineSum'     )

Ostap.Math.LegendreSum   .__repr__ = lambda s : _f_print_ ( s , 'LegendreSum'   )
Ostap.Math.ChebyshevSum  .__repr__ = lambda s : _f_print_ ( s , 'ChebyshevSum'  )
Ostap.Math.HermiteSum    .__repr__ = lambda s : _f_print_ ( s , 'HermiteSum'    )
Ostap.Math.Polynomial    .__repr__ = lambda s : _f_print_ ( s , 'Polynomial'    )
Ostap.Math.Bernstein     .__repr__ = lambda s : _f_print_ ( s , 'Bernstein'     )
Ostap.Math.BernsteinEven .__repr__ = lambda s : _f_print_ ( s , 'BernsteinEven' )
Ostap.Math.Positive      .__repr__ = lambda s : _f_print_ ( s , 'Positive'      )
Ostap.Math.PositiveEven  .__repr__ = lambda s : _f_print_ ( s , 'PositiveEven'  )
Ostap.Math.Convex        .__repr__ = lambda s : _f_print_ ( s , 'Convex'        ) 
Ostap.Math.ConvexOnly    .__repr__ = lambda s : _f_print_ ( s , 'ConvexOnly'    )
Ostap.Math.Monotonic     .__repr__ = lambda s : _f_print_ ( s , 'Monotonic'    )
Ostap.Math.FourierSum    .__repr__ = lambda s : _f_print_ ( s , 'FourierSum'    )
Ostap.Math.CosineSum     .__repr__ = lambda s : _f_print_ ( s , 'CosineSum'     )


Ostap.Math.LegendreSum2  .__repr__ = lambda s : "LegendreSum2(%d,%d)"       % ( s.nx() , s.ny() )
Ostap.Math.LegendreSum3  .__repr__ = lambda s : "LegendreSum3(%d,%d,%d)"    % ( s.nx() , s.ny() , s.nz() )
Ostap.Math.LegendreSum4  .__repr__ = lambda s : "LegendreSum4(%d,%d,%d,%d)" % ( s.nx() , s.ny() , s.nz() , s.nu() )

Ostap.Math.LegendreSum2  .__str__  =  Ostap.Math.LegendreSum2  .__repr__
Ostap.Math.LegendreSum3  .__str__  =  Ostap.Math.LegendreSum3  .__repr__
Ostap.Math.LegendreSum4  .__str__  =  Ostap.Math.LegendreSum4  .__repr__



# =============================================================================
## print function for splines 
def _sp_print_ ( self ,   typ = 'BSpline' ) :
    return '%s(%s,%s)' % ( typ, self.knots() , self.pars() )

Ostap.Math.BSpline         .__str__  = lambda s : 'BSpline(%s,%s)'             % ( s.knots ()      ,
                                                                                   s.pars  ()      )
Ostap.Math.PositiveSpline  .__str__  = lambda s : 'PositiveSpline(%s,%s)'      % ( s.knots ()      ,
                                                                                   s.pars  ()      )
Ostap.Math.ConvexOnlySpline.__str__  = lambda s : 'ConvexOnlySpline(%s,%s,%s)' % ( s.knots ()      ,
                                                                                   s.pars  ()      ,
                                                                                   s.convex()      )
Ostap.Math.MonotonicSpline.__str__  = lambda s : 'MonotonicSpline(%s,%s,%s)' % ( s.knots ()      ,
                                                                                   s.pars  ()      ,
                                                                                   s.increasing () )
Ostap.Math.ConvexSpline    .__str__  = lambda s : 'ConvexSpline(%s,%s,%s,%s)'  % ( s.knots ()      ,
                                                                                   s.pars  ()      ,
                                                                                   s.increasing () ,
                                                                                   s.convex()      )
for t in ( Ostap.Math.BSpline          ,
           Ostap.Math.PositiveSpline   ,
           Ostap.Math.ConvexOnlySpline ,
           Ostap.Math.MonotonicSpline ,
           Ostap.Math.ConvexSpline     ) : t.__repr__ = t.__str__
    
# =============================================================================
## decorate 2D-models/functions 
# =============================================================================
Ostap.Math.Spline2D    = Ostap.Math.PositiveSpline2D   
Ostap.Math.Spline2DSym = Ostap.Math.PositiveSpline2DSym

from ostap.math.minimize import sp_minimum_2D, sp_maximum_2D 

for model in ( Ostap.Math.BSpline2D           ,
               Ostap.Math.BSpline2DSym        , 
               Ostap.Math.PositiveSpline2D    ,
               Ostap.Math.PositiveSpline2DSym , 
               Ostap.Math.Bernstein2D         ,
               Ostap.Math.Positive2D          ,
               Ostap.Math.Bernstein2DSym      ,
               Ostap.Math.Positive2DSym       ,
               Ostap.Math.PS2DPol             ,
               Ostap.Math.PS2DPolSym          ,
               Ostap.Math.PS2DPol2            ,
               Ostap.Math.PS2DPol2Sym         ,
               Ostap.Math.PS2DPol3            ,
               Ostap.Math.PS2DPol3Sym         ,
               Ostap.Math.ExpoPS2DPol         ,
               Ostap.Math.Expo2DPol           ,
               Ostap.Math.Expo2DPolSym        ,
               Ostap.Math.LegendreSum2        ,
               Ostap.Math.Histo2D             ,
               Ostap.Math.Gauss2D             ) :
    
    model . tf2  =  tf2 
    model . tf   =  tf2
    model . draw =  f2_draw
    model.sp_integrate    = sp_integrate_2D
    model.sp_integrate_2D = sp_integrate_2D
    model.sp_integrate_x  = sp_integrate_2Dx
    model.sp_integrate_y  = sp_integrate_2Dy
    
    model.__getattr__     = _tf2_getattr_

    if sp_minimum_2D and not hasattr ( model , 'minimum' ) : model.minimum = sp_minimum_2D
    if sp_maximum_2D and not hasattr ( model , 'maximum' ) : model.maximum = sp_maximum_2D
    

from ostap.math.minimize import sp_minimum_3D, sp_maximum_3D 

# =============================================================================
## Decorate 3D models
# ============================================================================= 
for model in ( Ostap.Math.Bernstein3D    ,
               Ostap.Math.Bernstein3DSym ,
               Ostap.Math.Bernstein3DMix ,
               Ostap.Math.Positive3D     ,
               Ostap.Math.Positive3DSym  ,
               Ostap.Math.Positive3DMix  ,
               Ostap.Math.LegendreSum3   ,
               Ostap.Math.Histo3D        ,
               Ostap.Math.Gauss3D        ) :
    
    model . tf3  =  tf3 
    model . tf   =  tf3 
    model . draw =  f3_draw
    model.sp_integrate    = sp_integrate_3D
    model.sp_integrate_x  = sp_integrate_3Dx
    model.sp_integrate_y  = sp_integrate_3Dy
    model.sp_integrate_xy = sp_integrate_3Dxy
    model.sp_integrate_xz = sp_integrate_3Dxz
    model.sp_integrate_yz = sp_integrate_3Dyz
    model.__getattr__     = _tf3_getattr_

    if sp_minimum_2D and not hasattr ( model , 'minimum' ) : model.minimum = sp_minimum_2D
    if sp_maximum_2D and not hasattr ( model , 'maximum' ) : model.maximum = sp_maximum_2D
    
# ===============================================================================
def sp_minimum_1D_ ( pdf , xmin , xmax , x0 , *args ) :
    if hasattr ( pdf , 'setPars' ) : pdf.setPars()
    fun = pdf.function()
    return  sp_minimum_1D ( fun , xmin , xmax , x0 , *args )

def sp_maximum_1D_ ( pdf , xmin , xmax , x0 , *args ) :
    if hasattr ( pdf , 'setPars' ) : pdf.setPars()
    fun = pdf.function()
    return  sp_maximum_1D ( fun , xmin , xmax , x0 , *args ) 


# ============================================================================
## specfic ptints
# ============================================================================
Ostap.Math.NormalLaplace .__str__  = lambda s : "NormalLaplace(mu=%.4g,varsigma=%.4g,kL=%.4g,kR=%.4g)" % ( s.mu       () ,
                                                                                                       s.varsigma () ,
                                                                                                       s.kL       () ,
                                                                                                       s.kR       () )

# =============================================================================
## decorate 1D-PDFs
# =============================================================================

"""                 

for pdf in ( Ostap.Models.BreitWigner        ,
             Ostap.Models.BreitWignerMC      , 
             Ostap.Models.BWI                , 
             Ostap.Models.Flatte             ,
             Ostap.Models.Bukin              ,
             Ostap.Models.PhaseSpace2        ,
             Ostap.Models.PhaseSpaceNL       ,
             Ostap.Models.PhaseSpace23L      ,
             Ostap.Models.PhaseSpaceLeft     ,
             Ostap.Models.PhaseSpaceRight    ,
             Ostap.Models.PhaseSpacePol         ,
             Ostap.Models.PhaseSpaceLeftExpoPol ,
             Ostap.Models.Needham            ,
             Ostap.Models.CrystalBall        ,
             Ostap.Models.CrystalBallRS      ,
             Ostap.Models.CrystalBallDS      , 
             Ostap.Models.Apollonios          ,
             Ostap.Models.Apollonios2         , 
             Ostap.Models.GramCharlierA      , 
             Ostap.Models.Voigt              ,
             Ostap.Models.PseudoVoigt        ,
             Ostap.Models.Logistic           ,
             ## Ostap.Models.LASS               ,
             ## Ostap.Models.Bugg               ,
             ## Ostap.Models.LASS23L            ,
             ## Ostap.Models.Bugg23L            , 
             Ostap.Models.BWPS               , 
             Ostap.Models.BW3L               , 
             Ostap.Models.PolyPositive       ,
             Ostap.Models.ExpoPositive       ,
             Ostap.Models.TwoExpoPositive    ,
             Ostap.Models.PositiveSpline     ,
             Ostap.Models.MonotonicSpline   , 
             
             Ostap.Models.StudentT           ,
             Ostap.Models.BifurcatedStudentT , 
             Ostap.Models.GammaDist          , 
             Ostap.Models.GenGammaDist       , 
             Ostap.Models.Amoroso            ,
             Ostap.Models.LogGammaDist       ,
             Ostap.Models.Log10GammaDist     ,
             Ostap.Models.LogGamma           ,
             Ostap.Models.BetaPrime          ,
             Ostap.Models.Landau             ,
             Ostap.Models.SinhAsinh          , 
             Ostap.Models.JohnsonSU          ,
             Ostap.Models.Atlas              ,
             Ostap.Models.Sech               ,
             Ostap.Models.Losev              ,
             ## Ostap.Models.Swanson            ,
             Ostap.Models.Argus              ,
             Ostap.Models.GenArgus           ,
             Ostap.Models.Slash              ,
             Ostap.Models.AsymmetricLaplace  ,
             Ostap.Models.DoubleGauss        ,
             Ostap.Models.Gumbel             ,
             Ostap.Models.Weibull            ,
             Ostap.Models.RaisingCosine      ,
             Ostap.Models.QGaussian          ,
             Ostap.Models.KGaussian          ,
             Ostap.Models.Hyperbolic         ,
             Ostap.Models.Tsallis            ,
             Ostap.Models.QGSM               ,
             Ostap.Models.BifurcatedGauss    ,
             Ostap.Models.DoubleGauss        ,
             Ostap.Models.GenGaussV1         , 
             Ostap.Models.GenGaussV2         , 
             Ostap.Models.SkewGauss          ,
             Ostap.Models.PearsonIV          ,
             Ostap.Models.SkewGenT           
             ) :

    pdf.sp_integrate = sp_integrate_1D_
    if sp_minimum_1D and not hasattr ( pdf , 'minimum' ) : pdf . minimum = sp_minimum_1D_
    if sp_maximum_1D and not hasattr ( pdf , 'maximum' ) : pdf . maximum = sp_maximum_1D_
    
# ===============================================================================
def sp_minimum_2D_ ( pdf  ,
                     xmin , xmax ,
                     ymin , ymax ,
                     x0  = () , *args ) :
    if hasattr ( pdf , 'setPars' ) : pdf.setPars()
    fun = pdf.function()
    return  sp_minimum_2D ( fun ,
                            xmin , xmax ,
                            ymin , ymax , x0 , *args )
# ===============================================================================
def sp_maximum_2D_ ( pdf  ,
                     xmin , xmax ,
                     ymin , ymax , x0  = () , *args ) :
    if hasattr ( pdf , 'setPars' ) : pdf.setPars()
    fun = pdf.function()
    return  sp_maximum_2D ( fun  ,
                            xmin , xmax ,
                            ymin , ymax , x0 , *args ) 

# =============================================================================
## decorate 2D-PDFs
# =============================================================================

for pdf in ( Ostap.Models.Poly2DPositive     ,
             Ostap.Models.Poly2DSymPositive  , 
             Ostap.Models.PS2DPol            ,
             Ostap.Models.PS2DPolSym         , 
             Ostap.Models.PS2DPol2           ,
             Ostap.Models.PS2DPol2Sym        , 
             Ostap.Models.PS2DPol3           ,
             Ostap.Models.PS2DPol3Sym        , 
             Ostap.Models.ExpoPS2DPol        , 
             Ostap.Models.Expo2DPol          ,
             Ostap.Models.Expo2DPolSym       , 
             Ostap.Models.Spline2D           ,
             Ostap.Models.Gauss2D            ,
             Ostap.Models.Spline2DSym        ) :
    
    pdf.sp_integrate = sp_integrate_2D_
    if sp_minimum_2D and not hasattr ( pdf , 'minimum' ) : pdf.minimum = sp_minimum_2D_
    if sp_maximum_2D and not hasattr ( pdf , 'maximum' ) : pdf.maximum = sp_maximum_2D_


# ===============================================================================
def sp_minimum_3D_ ( pdf  ,
                     xmin , xmax ,
                     ymin , ymax ,
                     zmin , zmax , x0  = () , *args ) :
    if hasattr ( pdf , 'setPars' ) : pdf.setPars()
    fun = pdf.function()
    return  sp_minimum_3D ( fun ,
                            xmin , xmax ,
                            ymin , ymax ,
                            zmin , zmax , x0 , *args )
# ===============================================================================
def sp_maximum_3D_ ( pdf  ,
                     xmin , xmax ,
                     ymin , ymax ,
                     zmin , zmax , x0  = () , *args ) :
    if hasattr ( pdf , 'setPars' ) : pdf.setPars()
    fun = pdf.function()
    return  sp_maximum_3D ( fun  ,
                            xmin , xmax ,
                            ymin , ymax ,
                            zmin , zmax , x0 , *args ) 

# =============================================================================
## decorate 3D-PDFs
# =============================================================================

for pdf in ( Ostap.Models.Poly3DPositive    ,
             Ostap.Models.Poly3DSymPositive ,
             Ostap.Models.Poly3DMixPositive ) :
    
    pdf.sp_integrate = sp_integrate_3D_
    if sp_minimum_3D and not hasattr ( pdf , 'minimum' ) : pdf.minimum = sp_minimum_3D_
    if sp_maximum_3D and not hasattr ( pdf , 'maximum' ) : pdf.maximum = sp_maximum_3D_


"""

# =============================================================================
## set, get & iterator
from ostap.math.bernstein import _p_set_par_ , _p_get_par_, _p_iter_ 

for f in ( Ostap.Math.Bernstein2D    ,
           Ostap.Math.Positive2D     ,
           Ostap.Math.Bernstein2DSym ,
           Ostap.Math.Positive2DSym  ,
           ##
           ##
           Ostap.Math.BSpline2D           ,
           Ostap.Math.BSpline2DSym        ,
           Ostap.Math.PositiveSpline2D    ,
           Ostap.Math.PositiveSpline2DSym ,
           ##
           Ostap.Math.Bernstein3D    ,
           Ostap.Math.Bernstein3DSym ,
           Ostap.Math.Bernstein3DMix ,
           Ostap.Math.Positive3D     ,
           Ostap.Math.Positive3DSym  ,
           Ostap.Math.Positive3DMix  ,

           Ostap.Math.PolySum        ,
           ## 
           Ostap.Math.NSphere        ) :

    f.__setitem__  = _p_set_par_
    f.__getitem__  = _p_get_par_
    f.__len__      = lambda s     : s.npars() 
    f.__iter__     = _p_iter_
    f.__contains__ = lambda s , i : 0<=i<len(s)


# =============================================================================
## random generators 
# =============================================================================
from random import uniform as _uniform_


# =============================================================================
## generate random numbers from 2D bernstein-like distribuitions
#  @code
#  >>> func = ...
#  >>> for x,y in func.generate( 1000 ) : print x,y 
#  @endcode
def _random_generate_bernstein2D_ ( fun , num = 1 ) :
    """Generate random numbers from 2D bernstein-like distribuitions
    >>> func = ...
    >>> for x,y in func.generate( 1000 ) : print x,y 
    """
    xmn = fun.xmin ()
    xmx = fun.xmax ()
    ymn = fun.ymin ()
    ymx = fun.ymax ()
    vmx = max ( fun.bernstein().pars() )
    i   = 0 
    while i < num : 
        x = _uniform_ ( xmn , xmx ) 
        y = _uniform_ ( ymn , ymx ) 
        if fun ( x , y ) >= _uniform_ (   0 , vmx ) : 
            i+= 1 
            yield x,y

# =============================================================================
## generate random numbers from 3D bernstein-like distribuitions
#  @code
#  >>> func = ...
#  >>> for x,y,z in func.generate( 1000 ) : print x,y,z 
#  @endcode
def _random_generate_bernstein3D_ ( fun , num = 1 ) :
    """Generate random numbers from 2D bernstein-like distribuitions
    >>> func = ...
    >>> for x,y,z in func.generate( 1000 ) : print x,y,z 
    """
    xmn = fun.xmin ()
    xmx = fun.xmax ()
    ymn = fun.ymin ()
    ymx = fun.ymax ()
    zmn = fun.zmin ()
    zmx = fun.zmax ()
    vmx = max ( fun.bernstein().pars() )
    i   = 0 
    while i < num : 
        x = _uniform_ ( xmn , xmx ) 
        y = _uniform_ ( ymn , ymx ) 
        z = _uniform_ ( zmn , zmx )
        if fun ( x , y, z ) >= _uniform_ ( 0 , vmx ) :
            i+= 1 
            yield x , y , z
            
# =============================================================================
## Get random number from 2D bernstein-like distribuitions
#  @code
#  >>> func = ...
#  >>> print fun.shoot() 
#  @endcode
def _random_shoot_bernstein2D_ ( fun ) :
    """Get random number from 2D bernstein-like distribuitions
    >>> func = ...
    >>> print func.shoot()  
    """
    xmn = fun.xmin ()
    xmx = fun.xmax ()
    ymn = fun.ymin ()
    ymx = fun.ymax ()
    
    vmx = max ( fun.bernstein().pars() )
    while True : 
        x = _uniform_ ( xmn , xmx ) 
        y = _uniform_ ( ymn , ymx ) 
        if fun ( x , y ) >= _uniform_ (   0 , vmx ) : 
            return x,y

# =============================================================================
## Get random number from 3D bernstein-like distributions
#  @code
#  >>> func = ...
#  >>> print fun.shoot() 
#  @endcode
def _random_shoot_bernstein3D_ ( fun ) :
    """Get random number from 3D bernstein-like distribuitions
    >>> func = ...
    >>> print func.shoot()  
    """
    xmn = fun.xmin ()
    xmx = fun.xmax ()
    ymn = fun.ymin ()
    ymx = fun.ymax ()
    zmn = fun.zmin ()
    zmx = fun.zmax ()
    
    vmx = max ( fun.bernstein().pars() )
    while True : 
        x = _uniform_ ( xmn , xmx ) 
        y = _uniform_ ( ymn , ymx ) 
        z = _uniform_ ( zmn , zmx ) 
        if fun ( x , y , z ) >= _uniform_ ( 0 , vmx ) : 
            return x,y,z 



Ostap.Math.Positive2D    .generate = _random_generate_bernstein2D_
Ostap.Math.Positive2D    .shoot    = _random_shoot_bernstein2D_

Ostap.Math.Positive2DSym .generate = _random_generate_bernstein2D_
Ostap.Math.Positive2DSym .shoot    = _random_shoot_bernstein2D_

for p in ( Ostap.Math.Positive3D    ,
           Ostap.Math.Positive3DSym ,
           Ostap.Math.Positive3DMix ) :
    p.generate = _random_generate_bernstein3D_
    p.shoot    = _random_shoot_bernstein3D_
         
# =============================================================================
# 
# =============================================================================
## generate random numbers from single-mode distribution 
#  @code
#  >>> fun = ...
#  >>> for x in fun.generate( 1000 ) : print x 
#  @endcode
def _random_generate_mode_ ( fun , num = 1 , vmax = None , mode = None ) :
    """Generate random numbers from single-mode distribution 
    >>> fun = ...
    >>> for x in fun.generate( 1000 ) : print x 
    """
    xmn  = fun.xmin ()
    xmx  = fun.xmax ()

    if vmax is None or vmax <= 0 : 

        if mode is None or not xmn <= mode <= xmx :
            mode = fun.mode ()

        vmax  = fun ( mode )

    while i < num : 
        x = _uniform_ ( xmn , xmx ) 
        if fun ( x ) >= _uniform_ ( 0 , vmax ) :
            i += 1 
            yield x 

# =============================================================================
## Get random number from single-mode ditributions 
#  @code
#  >>> fun = ...
#  >>> print fun.shoot() 
#  @endcode
def _random_shoot_mode_ ( fun , vmax = None , mode  = None ) :
    """Get random number from single-mode distribuitions
    >>> fun = ...
    >>> print fun.shoot()  
    """
    xmn  = fun.xmin ()
    xmx  = fun.xmax ()

    if vmax is None or vmax <= 0 :

        if mode is None or not xmn <= mode <= xmx :
            mode = fun.mode ()
            
        vmax = fun ( mode )
        
    while True : 
        x = _uniform_ ( xmn , xmx ) 
        if fun ( x ) >= _uniform_ (   0 , vmax ) : 
            return x 


for m in ( Ostap.Math.PhaseSpaceNL , ) :
    m.generate = _random_generate_mode_
    m.shoot    = _random_shoot_mode_
    
# =============================================================================
## add complex amplitudes 
# =============================================================================
## Ostap.Math.LASS            . amp = _amp_
## Ostap.Math.LASS23L         . amp = _amp_
## Ostap.Math.Bugg23L         . amp = _amp_
Ostap.Math.Flatte          . amp = _amp_
## Ostap.Math.Flatte2         . amp = _amp_
## Ostap.Math.Flatte23L       . amp = _amp_
Ostap.Math.BW                 . amp = _amp_
    
# =============================================================================

import ostap.math.derivative as _D1
import ostap.math.integral   as _D2
for i in ( _D1.Derivative , _D2.Integral , _D2.IntegralCache ) :
    if not hasattr ( i , 'tf1' ) : i.tf1 = tf1

# =============================================================================
def _ff_str_ ( ff ) :
    """Self-printout for FormFactor"""
    return ff.describe()
Ostap.Math.FormFactor.__str__  = _ff_str_
Ostap.Math.FormFactor.__repr__ = _ff_str_
Ostap.Math.Channel   .__str__  = _ff_str_
Ostap.Math.Channel   .__repr__ = _ff_str_

def _bw_str_   ( bw ) :
    """Self-printout for Breit-Wigner function"""
    return "BreitWigner  (%s,%s)" % ( bw.m0() , bw.channel  () )
def _bwmc_str_ ( bw ) :
    """Self-printout for multi-channel Breit-Wigner function"""
    channels = [ bw.channel( i ) for i in range ( bw.nChannels() ) ]
    return "BreitWignerMC(%s,%s)" % ( bw.m0() , channels )

Ostap.Math.BreitWigner  .__str__  = _bw_str_
Ostap.Math.BreitWigner  .__repr__ = _bw_str_
Ostap.Math.BreitWignerMC.__str__  = _bwmc_str_
Ostap.Math.BreitWignerMC.__repr__ = _bwmc_str_


## trick ...
Ostap.Math.GammaBW3.xmin =  lambda s : s.dalitz().s_min() 


# =============================================================================
## more opreations for primitives 
# =============================================================================
## absolute value for primitive functions
#  \f$  g(x)  = \left|  f(x) \right| \f$
# @see Ostap::Magth::Abs 
def _prims_abs_ ( f ) :
    """Absolute value for primitive functions: g(x)  = |f(x)||
    - see `Ostap.Magth.Abs`
    """
    return Ostap.Math.Abs.create ( f ) 

# =============================================================================
## add two primitive functions
#  @see Ostap::Math::Sum 
def  _prims_add_ ( f1  , f2  ) :
    """Add two primitive functions
    - see `Ostap.Math.Sum`
    """
    return Ostap.Math.Sum.create ( f1 , f2 )

# =============================================================================
## add two primitive functions
#  @see Ostap::Math::Sum 
def  _prims_radd_ ( f1  , f2  ) :
    """Add two primitive functions
    - see `Ostap.Math.Sum``
    """
    return Ostap.Math.Sum.create ( f2 , f1 )

# =============================================================================
## multiply two primitive functions
#  @see Ostap::Math::Multiply 
def  _prims_mul_ ( f1  , f2  ) :
    """Multiply two primitive functions
    - see `Ostap.Math.Multiply`
    """
    return Ostap.Math.Multiply.create ( f1 , f2 )

# =============================================================================
## multiply two primitive functions
#  @see Ostap::Math::Multiply 
def  _prims_rmul_ ( f1  , f2  ) :
    """Multiply two primitive functions
    - see `Ostap.Math.Multiply`
    """
    return Ostap.Math.Multiply.create ( f2 , f1 )


# =============================================================================
## subtract two primitive functions
#  @see Ostap::Math::Subtract 
def  _prims_sub_ ( f1  , f2  ) :
    """Subtract two primitive functions
    - see `Ostap.Math.Subtract`
    """
    return Ostap.Math.Subtract.create ( f1 , f2 )

# =============================================================================
## subtract two primitive functions
#  @see Ostap::Math::Subtract 
def  _prims_rsub_ ( f1  , f2  ) :
    """Subtract two primitive functions
    - see `Ostap.Math.Subtract``
    """
    return Ostap.Math.Subtract.create ( f2 , f1 )

# =============================================================================
## division two primitive functions
#  @see Ostap::Math::Divide 
def  _prims_div_ ( f1  , f2  ) :
    """Division for two primitive functions
    - see `Ostap.Math.Divide`
    """
    return Ostap.Math.Divide.create ( f1 , f2 )

# =============================================================================
## division two primitive functions
#  @see Ostap::Math::Divide 
def  _prims_rdiv_ ( f1  , f2  ) :
    """Division for two primitive functions
    - see `Ostap.Math.Divide`
    """
    return Ostap.Math.Divide.create ( f2 , f1 )

# =============================================================================
## pow for two primitive functions
#  @see Ostap::Math::Pow 
def  _prims_pow_ ( f1  , f2  ) :
    """Pow for two primitive functions
    - see `Ostap.Math.Pow`
    """
    return Ostap.Math.Pow.create ( f1 , f2 )

# =============================================================================
## pow for two primitive functions
#  @see Ostap::Math::Pow 
def  _prims_rpow_ ( f1  , f2  ) :
    """Pow for two primitive functions
    - see `Ostap.Math.Pow`
    """
    return Ostap.Math.Pow.create ( f2 , f1 )


# =============================================================================
## exp for primitive functions
#  @see Ostap::Math::Exp
def _prims_exp_ ( f ):
    """Exp for primitive functions
    - see `Ostap.Math.Exp`
    """
    return Ostap.Math.Exp.create ( f )  

# =============================================================================
## log for primitive functions
#  @see Ostap::Math::Log
def _prims_log_ ( f ):
    """Exp for primitive functions
    - see `Ostap.Math.Log`
    """
    return Ostap.Math.Log.create ( f )  

# =============================================================================
## log10 for primitive functions
#  @see Ostap::Math::Log10
def _prims_log10_ ( f ):
    """Exp for primitive functions
    - see `Ostap.Math.Log10`
    """
    return Ostap.Math.Log10.create ( f )  

# =============================================================================
## sqrt for primitive functions
#  @see Ostap::Math::Sqrt
def _prims_sqrt_ ( f ):
    """Sqrt for primitive functions
    - see `Ostap.Math.Sqrt`
    """
    return Ostap.Math.Sqrt.create ( f )  

# =============================================================================
## Cbrt for primitive functions
#  @see Ostap::Math::Cbrt
def _prims_cbrt_ ( f ):
    """Cbrt for primitive functions
    - see `Ostap.Math.Cbrt`
    """
    return Ostap.Math.Cbrt.create ( f )  


# =============================================================================
## Erf for primitive functions
#  @see Ostap::Math::Erf
def _prims_erf_ ( f ):
    """Erf for primitive functions
    - see `Ostap.Math.Erf`
    """
    return Ostap.Math.Erf.create ( f )  

# =============================================================================
## Erfc for primitive functions
#  @see Ostap::Math::Erfc
def _prims_erfc_ ( f ):
    """Erfc for primitive functions
    - see `Ostap.Math.Erfc`
    """
    return Ostap.Math.Erfc.create ( f )  

# =============================================================================
## Gamma function for primitive functions
#  @see Ostap::Math::TGamma
def _prims_tgamma_ ( f ):
    """Gamma function for primitive functions
    - see `Ostap.Math.TGamma`
    """
    return Ostap.Math.TGamma.create ( f )  

# =============================================================================
## log Gamma function for primitive functions
#  @see Ostap::Math::LGamma
def _prims_lgamma_ ( f ):
    """log-Gamma function for primitive functions
    - see `Ostap.Math.LGamma`
    """
    return Ostap.Math.LGamma.create ( f )  

# =============================================================================
## Sin function for primitive functions
#  @see Ostap::Math::Sin
def _prims_sin_ ( f ):
    """sin function for primitive functions
    - see `Ostap.Math.Sin`
    """
    return Ostap.Math.Sin.create ( f )  

# =============================================================================
## Cos function for primitive functions
#  @see Ostap::Math::Cos
def _prims_cos_ ( f ):
    """cos function for primitive functions
    - see `Ostap.Math.Cos`
    """
    return Ostap.Math.Cos.create ( f )  

# =============================================================================
## Tan function for primitive functions
#  @see Ostap::Math::Tan
def _prims_tan_ ( f ):
    """tan function for primitive functions
    - see `Ostap.Math.Tan`
    """
    return Ostap.Math.Tan.create ( f )  

# =============================================================================
## ASin function for primitive functions
#  @see Ostap::Math::ASin
def _prims_asin_ ( f ):
    """asin function for primitive functions
    - see `Ostap.Math.ASin`
    """
    return Ostap.Math.ASin.create ( f )  

# =============================================================================
## ACos function for primitive functions
#  @see Ostap::Math::ACos
def _prims_acos_ ( f ):
    """acos function for primitive functions
    - see `Ostap.Math.ACos`
    """
    return Ostap.Math.ACos.create ( f )  

# =============================================================================
## ATan function for primitive functions
#  @see Ostap::Math::ATan
def _prims_atan_ ( f ):
    """atan function for primitive functions
    - see `Ostap.Math.ATan`
    """
    return Ostap.Math.ATan.create ( f )  

# =============================================================================
## Sinh function for primitive functions
#  @see Ostap::Math::Sinh
def _prims_sinh_ ( f ):
    """sinh function for primitive functions
    - see `Ostap.Math.Sinh`
    """
    return Ostap.Math.Sinh.create ( f )  

# =============================================================================
## Cosh function for primitive functions
#  @see Ostap::Math::Cosh
def _prims_cosh_ ( f ):
    """cosh function for primitive functions
    - see `Ostap.Math.Cosh`
    """
    return Ostap.Math.Cosh.create ( f )  

# =============================================================================
## Tanh function for primitive functions
#  @see Ostap::Math::Tanh
def _prims_tanh_ ( f ):
    """tanh function for primitive functions
    - see `Ostap.Math.Tanh`
    """
    return Ostap.Math.Tanh.create ( f )  

# =============================================================================
## ASinh function for primitive functions
#  @see Ostap::Math::ASinh
def _prims_asinh_ ( f ):
    """asinh function for primitive functions
    - see `Ostap.Math.ASinh`
    """
    return Ostap.Math.ASinh.create ( f )  

# =============================================================================
## ACosh function for primitive functions
#  @see Ostap::Math::ACosh
def _prims_acosh_ ( f ):
    """acosh function for primitive functions
    - see `Ostap.Math.ACosh`
    """
    return Ostap.Math.ACosh.create ( f )  

# =============================================================================
## ATanh function for primitive functions
#  @see Ostap::Math::ATanh
def _prims_atanh_ ( f ):
    """atan function for primitive functions
    - see `Ostap.Math.ATan`
    """
    return Ostap.Math.ATanh.create ( f )  


for p in ( Ostap.Math.Const    ,
           Ostap.Math.Id       ,
           Ostap.Math.Linear   ,
           Ostap.Math.Compose  ,
           Ostap.Math.Multiply ,
           Ostap.Math.Divide   ,
           Ostap.Math.Sum      ,
           Ostap.Math.Subtract ,
           Ostap.Math.Moebius  ,
           Ostap.Math.Step     ,
           Ostap.Math.Max      ,
           Ostap.Math.Min      ,
           Ostap.Math.Apply    ,
           Ostap.Math.Abs      ,
           Ostap.Math.Sqrt     ,
           Ostap.Math.Cbrt     ,
           Ostap.Math.Exp      ,
           Ostap.Math.Log      ,
           Ostap.Math.Log10    ,
           Ostap.Math.Erf      ,
           Ostap.Math.Erfc     ,
           Ostap.Math.TGamma   ,
           Ostap.Math.LGamma   ,
           Ostap.Math.Sin      ,
           Ostap.Math.Cos      ,
           Ostap.Math.Tan      ,
           Ostap.Math.ASin     ,
           Ostap.Math.ACos     ,
           Ostap.Math.ATan     ,
           Ostap.Math.Sinh     ,
           Ostap.Math.Cosh     ,
           Ostap.Math.Tanh     ,
           Ostap.Math.ASinh    ,
           Ostap.Math.ACosh    ,
           Ostap.Math.ATanh    ,
           Ostap.Math.Pow      ) :

           p. __abs__      = _prims_abs_
           
           p. __add__      = _prims_add_
           p. __radd__     = _prims_radd_

           p. __sub__      = _prims_sub_
           p. __rsub__     = _prims_rsub_

           p. __mul__      = _prims_mul_
           p. __rmul__     = _prims_rmul_

           p. __div__      = _prims_div_
           p. __rdiv__     = _prims_rdiv_

           p. __truediv__  = _prims_div_
           p. __rtruediv__ = _prims_rdiv_

           p. __pow__      = _prims_pow_
           p. __rpow__     = _prims_rpow_

           p. __exp__       = _prims_exp_
           p. __log__       = _prims_log_
           p. __log10__     = _prims_log10_
           p. __sqrt__      = _prims_sqrt_
           p. __cbrt__      = _prims_cbrt_
           p. __erf__       = _prims_erf_
           p. __erfc__      = _prims_erfc_
           p. __tgamma__    = _prims_tgamma_
           p. __gamma__     = _prims_tgamma_
           p. __lgamma__    = _prims_lgamma_           
           p. __sin__       = _prims_sin_
           p. __cos__       = _prims_cos_
           p. __tan__       = _prims_tan_
           p. __asin__      = _prims_asin_
           p. __acos__      = _prims_acos_
           p. __atan__      = _prims_atan_
           p. __sinh__      = _prims_sinh_
           p. __cosh__      = _prims_cosh_
           p. __tanh__      = _prims_tanh_
           p. __asinh__     = _prims_asinh_
           p. __acosh__     = _prims_acosh_
           p. __atanh__     = _prims_atanh_

    
# =============================================================================
_decorated_classes_ = set( [
    ##
    Ostap.Math.Positive               ,
    Ostap.Math.PositiveEven           , 
    Ostap.Math.Monotonic              , 
    Ostap.Math.Convex                 , 
    Ostap.Math.ConvexOnly             , 
    Ostap.Math.PositiveSpline         , 
    Ostap.Math.MonotonicSpline        , 
    Ostap.Math.ConvexSpline           ,
    Ostap.Math.ConvexOnlySpline       ,
    Ostap.Math.ExpoPositive           ,
    Ostap.Math.TwoExpoPositive        ,
    ## Ostap.Math.LASS                   , 
    ## Ostap.Math.LASS23L                , 
    ## Ostap.Math.Bugg23L                , 
    Ostap.Math.Flatte                 , 
    ## Ostap.Math.Flatte2                , 
    ## Ostap.Math.Flatte23L              ,
    Ostap.Math.BreitWigner            ,
    Ostap.Math.BW                     ,
    Ostap.Math.BreitWignerMC          ,
    Ostap.Math.BWPS                   ,
    Ostap.Math.BW3L                   ,
    ## Ostap.Math.Swanson                ,
    ##
    Ostap.Math.Chebyshev              ,
    Ostap.Math.ChebyshevU             ,
    Ostap.Math.Legendre               ,
    Ostap.Math.Hermite                ,
    Ostap.Math.Bernstein              ,
    Ostap.Math.BernsteinEven          ,
    Ostap.Math.ChebyshevSum           ,
    Ostap.Math.LegendreSum            ,
    Ostap.Math.HermiteSum             ,
    Ostap.Math.FourierSum             ,
    Ostap.Math.CosineSum              ,
    Ostap.Math.Polynomial             ,               
    Ostap.Math.Positive               ,
    Ostap.Math.PositiveEven           ,
    Ostap.Math.Monotonic              ,
    Ostap.Math.Convex                 ,
    Ostap.Math.ConvexOnly             ,
    Ostap.Math.BifurcatedGauss        ,
    Ostap.Math.Bukin                  ,
    Ostap.Math.Novosibirsk            ,
    Ostap.Math.SinhAsinh              ,
    Ostap.Math.CrystalBall            ,
    Ostap.Math.Needham                ,
    Ostap.Math.CrystalBallDoubleSided ,
    Ostap.Math.GramCharlierA          ,
    Ostap.Math.PhaseSpace2            ,
    Ostap.Math.sPhaseSpace2           ,
    Ostap.Math.PhaseSpace3            ,
    Ostap.Math.PhaseSpace3s           ,
    Ostap.Math.PhaseSpaceLeft         ,
    Ostap.Math.PhaseSpaceRight        ,
    Ostap.Math.PSDalitz               ,
    Ostap.Math.PhaseSpaceNL           ,
    Ostap.Math.PhaseSpace23L          ,
    Ostap.Math.BreitWigner            ,
    Ostap.Math.BW                     ,
    Ostap.Math.BreitWignerMC          ,
    Ostap.Math.Rho0                   ,
    Ostap.Math.Kstar0                 ,
    Ostap.Math.Phi0                   ,
    ## Ostap.Math.Rho0FromEtaPrime       ,
    Ostap.Math.Flatte                 ,
    ## Ostap.Math.Flatte2                ,
    ## Ostap.Math.LASS                   ,
    ## Ostap.Math.LASS23L                ,
    ## Ostap.Math.Bugg23L                ,
    ## Ostap.Math.BW23L                  ,
    ## Ostap.Math.Flatte23L              ,
    ## Ostap.Math.Gounaris23L            ,
    Ostap.Math.StudentT               ,
    Ostap.Math.BifurcatedStudentT     ,
    Ostap.Math.Voigt                  ,
    Ostap.Math.PseudoVoigt            ,
    Ostap.Math.Logistic               ,
    #
    Ostap.Math.GenGaussV1             ,
    Ostap.Math.GenGaussV2             ,
    Ostap.Math.SkewGauss              , 
    Ostap.Math.ExGauss                , 
    Ostap.Math.NormalLaplace          , 
    Ostap.Math.GammaDist              ,
    Ostap.Math.GenGammaDist           ,
    Ostap.Math.Amoroso                ,
    Ostap.Math.LogGammaDist           ,
    Ostap.Math.Log10GammaDist         ,
    Ostap.Math.LogGamma               ,
    Ostap.Math.BetaPrime              ,
    Ostap.Math.Landau                 ,
    Ostap.Math.JohnsonSU              ,
    Ostap.Math.Atlas                  ,
    Ostap.Math.Sech                   ,
    Ostap.Math.Losev                  ,
    ## Ostap.Math.Swanson                ,
    Ostap.Math.Argus                  ,
    Ostap.Math.GenArgus               ,
    Ostap.Math.Slash                  ,
    Ostap.Math.AsymmetricLaplace      ,
    Ostap.Math.DoubleGauss            ,
    Ostap.Math.Gumbel                 ,
    Ostap.Math.Weibull                ,
    Ostap.Math.QGaussian              ,
    Ostap.Math.KGaussian              ,
    Ostap.Math.RaisingCosine          ,
    Ostap.Math.Tsallis                ,
    Ostap.Math.QGSM                   ,
    Ostap.Math.TwoExpos               ,
    Ostap.Math.Sigmoid                ,
    #
    Ostap.Math.BSpline                , 
    Ostap.Math.PositiveSpline         ,
    Ostap.Math.MonotonicSpline       ,
    Ostap.Math.ConvexSpline           ,
    Ostap.Math.ConvexOnlySpline       ,
    #
    Ostap.Math.BernsteinDualBasis     ,
    ##
    Ostap.Math.Bernstein         ,
    Ostap.Math.BernsteinEven     , 
    Ostap.Math.Positive          ,
    Ostap.Math.PositiveEven      ,
    Ostap.Math.Monotonic        ,
    Ostap.Math.Convex            ,
    Ostap.Math.ConvexOnly        ,
    Ostap.Math.ChebyshevSum      ,               
    Ostap.Math.LegendreSum       ,               
    Ostap.Math.HermiteSum        ,               
    Ostap.Math.FourierSum        ,               
    Ostap.Math.CosineSum         ,               
    Ostap.Math.Polynomial        ,               
    Ostap.Math.ExpoPositive      , 
    Ostap.Math.TwoExpoPositive   , 
    #
    Ostap.Math.BSpline           ,
    Ostap.Math.MonotonicSpline  ,
    Ostap.Math.PositiveSpline    , 
    Ostap.Math.ConvexSpline      , 
    Ostap.Math.ConvexOnlySpline  ,
    #
    Ostap.Math.LegendreSum    , 
    Ostap.Math.ChebyshevSum  , 
    Ostap.Math.Polynomial    , 
    Ostap.Math.Bernstein     ,
    Ostap.Math.BernsteinEven ,
    Ostap.Math.Positive      ,
    Ostap.Math.PositiveEven  ,
    Ostap.Math.FourierSum    ,
    Ostap.Math.CosineSum     ,
    Ostap.Math.LegendreSum   ,
    Ostap.Math.ChebyshevSum  ,
    Ostap.Math.HermiteSum    ,
    Ostap.Math.Polynomial    ,
    Ostap.Math.Bernstein     ,
    Ostap.Math.BernsteinEven ,
    Ostap.Math.Positive      ,
    Ostap.Math.PositiveEven  ,
    Ostap.Math.Convex        ,
    Ostap.Math.ConvexOnly    ,
    Ostap.Math.Monotonic    ,
    Ostap.Math.FourierSum    ,
    Ostap.Math.CosineSum     ,
    ##
    Ostap.Math.BSpline2D       ,
    Ostap.Math.BSpline2DSym    , 
    Ostap.Math.PositiveSpline2D    ,
    Ostap.Math.PositiveSpline2DSym , 
    Ostap.Math.Bernstein2D    ,
    Ostap.Math.Positive2D     ,
    Ostap.Math.Bernstein2DSym ,
    Ostap.Math.Positive2DSym  ,
    Ostap.Math.PS2DPol        ,
    Ostap.Math.PS2DPolSym     ,
    Ostap.Math.PS2DPol2       ,
    Ostap.Math.PS2DPol2Sym    ,
    Ostap.Math.PS2DPol3       ,
    Ostap.Math.PS2DPol3Sym    ,
    Ostap.Math.ExpoPS2DPol    ,
    Ostap.Math.Expo2DPol      ,
    Ostap.Math.Expo2DPolSym   ,
    ##
    ## Ostap.Models.BreitWigner        , 
    ## Ostap.Models.BreitWignerMC      , 
    ## Ostap.Models.BWI                , 
    ## Ostap.Models.Flatte             ,
    ## Ostap.Models.Bukin              ,
    ## Ostap.Models.PhaseSpace2        ,
    ## Ostap.Models.PhaseSpaceNL       ,
    ## Ostap.Models.PhaseSpace23L      ,
    ## Ostap.Models.PhaseSpaceLeft     ,
    ## Ostap.Models.PhaseSpaceRight    ,
    ## Ostap.Models.PhaseSpacePol      ,
    ## Ostap.Models.PhaseSpaceLeftExpoPol ,
    ## Ostap.Models.Needham            ,
    ## Ostap.Models.CrystalBall        ,
    ## Ostap.Models.CrystalBallRS      ,
    ## Ostap.Models.CrystalBallDS      , 
    ## Ostap.Models.Apollonios         ,
    ## Ostap.Models.Apollonios2        , 
    ## Ostap.Models.GramCharlierA      , 
    ## Ostap.Models.Voigt              ,
    ## Ostap.Models.PseudoVoigt        ,
    ## Ostap.Models.Logistic           ,
    ## ## Ostap.Models.LASS               ,
    ## ## Ostap.Models.Bugg               ,
    ## ## Ostap.Models.LASS23L            ,
    ## ## Ostap.Models.Bugg23L            , 
    ## ## Ostap.Models.BW23L              , 
    ## Ostap.Models.PolyPositive       ,
    ## Ostap.Models.ExpoPositive       ,
    ## Ostap.Models.TwoExpoPositive    ,
    ## Ostap.Models.PositiveSpline     ,
    ## Ostap.Models.MonotonicSpline   , 
    ## ##
    ## Ostap.Models.StudentT           ,
    ## Ostap.Models.BifurcatedStudentT , 
    ## Ostap.Models.GammaDist          , 
    ## Ostap.Models.GenGammaDist       , 
    ## Ostap.Models.Amoroso            ,
    ## Ostap.Models.LogGammaDist       ,
    ## Ostap.Models.Log10GammaDist     ,
    ## Ostap.Models.LogGamma           ,
    ## Ostap.Models.BetaPrime          ,
    ## Ostap.Models.Landau             ,
    ## Ostap.Models.SinhAsinh          , 
    ## Ostap.Models.JohnsonSU          ,
    ## Ostap.Models.Atlas              ,
    ## Ostap.Models.Sech               ,
    ## Ostap.Models.Losev              ,
    ## ## Ostap.Models.Swanson            ,
    ## Ostap.Models.Argus              ,
    ## Ostap.Models.GenArgus           ,
    ## Ostap.Models.Slash              ,
    ## Ostap.Models.AsymmetricLaplace  ,
    ## Ostap.Models.Tsallis            ,
    ## Ostap.Models.QGSM               ,
    ## Ostap.Models.BifurcatedGauss    ,
    ## Ostap.Models.GenGaussV1         , 
    ## Ostap.Models.GenGaussV2         , 
    ## ##
    ## Ostap.Models.Poly2DPositive     ,
    ## Ostap.Models.Poly2DSymPositive  , 
    ## Ostap.Models.PS2DPol            ,
    ## Ostap.Models.PS2DPolSym         , 
    ## Ostap.Models.PS2DPol2           ,
    ## Ostap.Models.PS2DPol2Sym        , 
    ## Ostap.Models.PS2DPol3           ,
    ## Ostap.Models.PS2DPol3Sym        , 
    ## Ostap.Models.ExpoPS2DPol        , 
    ## Ostap.Models.Expo2DPol          ,
    ## Ostap.Models.Expo2DPolSym       , 
    ## Ostap.Models.Spline2D           ,
    ## Ostap.Models.Spline2DSym        ,
    ##
    Ostap.Math.Positive       ,
    Ostap.Math.PositiveEven   ,  
    Ostap.Math.Bernstein      , 
    Ostap.Math.BernsteinEven  , 
    Ostap.Math.Bernstein2D    ,
    Ostap.Math.Positive2D     ,
    Ostap.Math.Bernstein2DSym ,
    Ostap.Math.Positive2DSym  ,
    ##
    Ostap.Math.BSpline        ,
    Ostap.Math.PositiveSpline ,
    Ostap.Math.Spline2D       ,
    Ostap.Math.Spline2DSym    ,
    Ostap.Math.Gauss2D        ,
    ## 
    Ostap.Math.PolySum        ,
    ## 
    Ostap.Math.NSphere        ,
    ##
    ## Ostap.Math.LASS          , 
    ## Ostap.Math.LASS23L       ,
    ## Ostap.Math.Bugg23L       , 
    Ostap.Math.Flatte        ,
    ## Ostap.Math.Flatte2       , 
    ## Ostap.Math.Flatte23L     , 
    Ostap.Math.BreitWigner     ,
    Ostap.Math.BW              ,
    Ostap.Math.BreitWignerMC   ,
    ##
    Ostap.Math.Bernstein3D    ,
    Ostap.Math.Bernstein3DSym ,
    Ostap.Math.Bernstein3DMix ,
    Ostap.Math.Positive3D     ,
    Ostap.Math.Positive3DSym  ,
    Ostap.Math.Positive3DMix  ,
    ##
    Ostap.Functions.PyCallable        , 
    Ostap.Math.Piecewise              , 
    Ostap.Math.ChebyshevApproximation ,
    D.Derivative                      ,
    
    Ostap.Math.Multiply               ,
    Ostap.Math.Divide                 ,
    Ostap.Math.Const                  ,
    Ostap.Math.Linear                 ,
    Ostap.Math.Sum                    ,
    Ostap.Math.Moebius                ,
    Ostap.Math.Compose                ,
    Ostap.Math.Step                   ,
    Ostap.Math.Min                    ,
    Ostap.Math.Max                    ,
    Ostap.Math.Apply                  ,
    Ostap.Math.KramersKronig          ,
    ])

# ============================================================================
import ostap.math.bernstein 
import ostap.math.bspline

     
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================
