#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file models.py
#
#  Module with some useful utilities for simple functions and fit models.
#
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2011-12-01
#
# =============================================================================
"""Module with some useful fit-models"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = ()
# =============================================================================
import  ROOT 
from    ostap.core.core import cpp, Ostap, funID
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
    def __init__ ( self , o              ) :        self._o   =  o 
    def __call__ ( self , x , pars  = [] ) : return self._o ( x [0]        )
# helper adapter for 1D-functions 
class _WO2_ (object)  :
    "Helper adapter for 2D-functions"
    def __init__ ( self , o              ) :        self._o   =  o 
    def __call__ ( self , x , pars  = [] ) : return self._o ( x [0] , x[1] )
# =============================================================================
pos_infinity = float('+inf')
neg_infinity = float('-inf')
# =============================================================================
## convert the model into TF1
def _tf1_ ( self                 ,
            xmin  = neg_infinity ,
            xmax  = pos_infinity ,
            npars = 0            , *args ) :
    """Convert the function to TF1    
    >>> obj = ...
    >>> fun = obj.tf1 ( 3.0 , 3.2 )
    >>> fun.Draw() 
    """
    #
    if not hasattr ( self , '_wo1' ) : self._wo1 = _WO1_ ( self )
    if not self._wo1                 : self._wo1 = _WO1_ ( self )
    #
    if hasattr ( self , 'xmin'  ) : xmin  = max ( xmin  , self.xmin () )
    if hasattr ( self , 'xmax'  ) : xmax  = min ( xmax  , self.xmax () )
    if hasattr ( self , 'npars' ) : npars = max ( npars , self.npars() )
    #
    _wo = self._wo1 
    fun = ROOT.TF1 ( funID()  , _wo , xmin , xmax , npars, *args )
    fun.SetNpx ( 500 ) 
    #
    return fun 

# =============================================================================
## convert the model into TF2
def _tf2_ ( self ,
            xmin  = neg_infinity ,
            xmax  = pos_infinity ,
            ymin  = neg_infinity ,
            ymax  = pos_infinity ,
            npars = 0            , *args ) :
    """Convert the function to TF2
    >>> obj = ...    
    >>> fun = obj.tf2 ( 3.0 , 3.2 , 3.0 , 3.2 )    
    >>> fun.Draw() 
    """
    ##
    if not hasattr ( self , '_wo2' ) : self._wo2 = _WO2_ ( self )
    if not self._wo2                 : self._wo2 = _WO2_ ( self )
    ## 
    if hasattr ( self , 'xmin'  ) : xmin  = max ( xmin  , self.xmin () )
    if hasattr ( self , 'xmax'  ) : xmax  = min ( xmax  , self.xmax () )
    if hasattr ( self , 'ymin'  ) : ymin  = max ( ymin  , self.ymin () )
    if hasattr ( self , 'ymax'  ) : ymax  = min ( ymax  , self.ymax () )
    if hasattr ( self , 'npars' ) : npars = max ( npars , self.npars() )
    #
    _wo = self._wo2
    fun = ROOT.TF2 ( funID ()  , _wo , xmin , xmax , ymin , ymax , npars , *args )
    fun.SetNpx ( 100 ) 
    fun.SetNpy ( 100 ) 
    #
    return fun 

# =============================================================================
## draw the function 
def _f1_draw_ ( self , *opts ) :
    """Drawing the function object through conversion to ROOT.TF1    
    >>> fun = ...
    >>> fun.draw()    
    """
    if not hasattr ( self , '_tf1'  ) :

        self._tf1 = _tf1_ ( self )
        if type(self) in ( Ostap.Math.Positive          ,
                           Ostap.Math.PositiveEven      , 
                           Ostap.Math.Monothonic        , 
                           Ostap.Math.Convex            , 
                           Ostap.Math.ConvexOnly        , 
                           Ostap.Math.PositiveSpline    , 
                           Ostap.Math.MonothonicSpline  , 
                           Ostap.Math.ConvexSpline      ,
                           Ostap.Math.ConvexOnlySpline  ,
                           Ostap.Math.ExpoPositive      ,
                           Ostap.Math.TwoExpoPositive   ) :                                
            self._tf1.SetMinimum(0)
            
    return self._tf1.Draw ( *opts )
    
    
# =============================================================================
## get the regular complex value for amplitude 
def _amp_ ( self , x ) :
    """ Get the complex value for amplitude
    >>> fun
    >>> a = fun.amp ( x )    
    """
    v = self.amplitude ( x )
    #
    return complex( v.real () , v.imag () ) 

Ostap.Math.LASS        . amp = _amp_
Ostap.Math.LASS23L     . amp = _amp_
Ostap.Math.Bugg23L     . amp = _amp_
Ostap.Math.Flatte      . amp = _amp_
Ostap.Math.Flatte2     . amp = _amp_
Ostap.Math.Flatte23L   . amp = _amp_
Ostap.Math.BreitWigner . amp = _amp_
Ostap.Math.Swanson     . amp = _amp_

# =============================================================================
## make 1D-integration using SciPy
#  @see http://www.scipy.org/
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_1D ( func , xmin , xmax , *args , **kwargs ) :
    """Make 1D-integration using SciPy
    
    >>> func = ...
    >>> print func.sp_integrate ( -10 , 10 )    
    """    
    from scipy import integrate
    result = integrate.quad ( func , xmin , xmax , *args , **kwargs )
    return result[0]

# =============================================================================
## make 2D-integration using SciPy
#  @see http://www.scipy.org/
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_2D ( func  ,
                      xmin  , xmax ,
                      ymin  , ymax , *args , **kwargs ) :
    """Make 2D-integration using SciPy

    >>> func = ...  ## func ( x , y )
    ##                            xmin , xmax , ymin , ymax 
    >>> print func.sp_integrate ( -10  , 10   , -20  , 20   ) 
    
    Note different naming with respect to SciPy:
    - SciPy first integrates over 2nd variable
    """
    from scipy import integrate
    result = integrate.dblquad ( func ,
                                 ymin ,
                                 ymax ,
                                 lambda x : xmin ,
                                 lambda x : xmax , 
                                 *args , **kwargs )
    return result[0]

# =============================================================================
## make 1D-integration using SciPy
#  @see http://www.scipy.org/
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_1D_ ( pdf , xmin , xmax , *args , **kwargs ) :
    """Make 1D integration over the PDF using SciPy
    """
    if hasattr ( pdf , 'setPars' ) : pdf.setPars() 
    func = pdf.function()
    return func.sp_integrate_1D ( xmin , xmax , *args , **kwargs ) 

# =============================================================================
## make 2D-integration using SciPy
#  @see http://www.scipy.org/
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-12-01
def sp_integrate_2D_ ( pdf   ,
                      xmin  , xmax ,
                      ymin  , ymax , *args , **kwargs ) :
    """ Make 3D integration over the PDF using SciPy
    """
    if hasattr ( pdf , 'setPars' ) : pdf.setPars() 
    func = pdf.function()
    return func.sp_integrate_2D ( xmin , xmax , ymin , ymax , *args , **kwargs ) 


from ostap.stats.moments import moment   as sp_moment
from ostap.stats.moments import mean     as sp_mean
from ostap.stats.moments import variance as sp_variance
from ostap.stats.moments import rms      as sp_rms 
from ostap.stats.moments import median   as sp_median
from ostap.stats.moments import quantile as sp_quantile
from ostap.stats.moments import mode     as sp_mode 


# =============================================================================
## helper function to delegate some methods/attributes to TF1
#  @code
#  f = ...
#  f.SetLineColor(4) ## delegate to TF1
#  f.SetLineWidth(2) ## delegate to TF1
#  @endcode 
def _f1_getattr_ ( self , attr ) :
    """Delegate some methods/attributes to TF1
    >>> f = ...
    >>> f.SetLineColor(4) ## delegate to TF1
    >>> f.SetLineWidth(2) ## delegate to TF1
    """
    if  hasattr ( ROOT.TF1 , attr ) :
        if not hasattr  ( self      , '_tf1' ) : self._tf1 = self.tf1 ( ) 
        return getattr  ( self._tf1 , attr   )
    
    raise AttributeError

    
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
               Ostap.Math.Monothonic             ,
               Ostap.Math.Convex                 ,
               Ostap.Math.ConvexOnly             ,
               Ostap.Math.BifurcatedGauss        ,
               Ostap.Math.Bukin                  ,
               Ostap.Math.Novosibirsk            ,
               Ostap.Math.CrystalBall            ,
               Ostap.Math.Needham                ,
               Ostap.Math.CrystalBallDoubleSided ,
               Ostap.Math.GramCharlierA          ,
               Ostap.Math.PhaseSpace2            ,
               Ostap.Math.PhaseSpaceLeft         ,
               Ostap.Math.PhaseSpaceRight        ,
               Ostap.Math.PhaseSpaceNL           ,
               Ostap.Math.PhaseSpace23L          ,
               Ostap.Math.BreitWigner            ,
               Ostap.Math.Rho0                   ,
               Ostap.Math.Rho0FromEtaPrime       ,
               Ostap.Math.Flatte                 ,
               Ostap.Math.Flatte2                ,
               Ostap.Math.LASS                   ,
               Ostap.Math.LASS23L                ,
               Ostap.Math.Bugg23L                ,
               Ostap.Math.BW23L                  ,
               Ostap.Math.Flatte23L              ,
               Ostap.Math.Gounaris23L            ,
               Ostap.Math.StudentT               ,
               Ostap.Math.BifurcatedStudentT     ,
               Ostap.Math.Voigt                  ,
               Ostap.Math.PseudoVoigt            ,
               Ostap.Math.Logistic               ,
               #
               Ostap.Math.GenGaussV1             ,
               Ostap.Math.GenGaussV2             ,
               ## Ostap.Math.SkewGauss              , ## (temporarily removed)
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
               Ostap.Math.Swanson                ,
               Ostap.Math.Argus                  ,
               Ostap.Math.Tsallis                ,
               Ostap.Math.QGSM                   ,
               Ostap.Math.TwoExpos               ,
               Ostap.Math.Sigmoid                ,
               #
               Ostap.Math.BSpline                , 
               Ostap.Math.PositiveSpline         ,
               Ostap.Math.MonothonicSpline       ,
               Ostap.Math.ConvexSpline           ,
               Ostap.Math.ConvexOnlySpline       ,
               #
               Ostap.Math.BernsteinDualBasis     ,

               ) :
    model.tf1          = _tf1_ 
    model.sp_integrate = sp_integrate_1D
    model.__getattr__  = _f1_getattr_
    
    if not hasattr ( model , 'mean'     ) : model.mean     = sp_mean 
    if not hasattr ( model , 'variance' ) : model.variance = sp_variance 
    if not hasattr ( model , 'rms'      ) : model.rms      = sp_rms  
    if not hasattr ( model , 'median'   ) : model.median   = sp_median
    if not hasattr ( model , 'mode'     ) : model.mode     = sp_mode 
    if not hasattr ( model , 'moment'   ) : model.moment   = sp_moment
    if not hasattr ( model , 'quantile' ) : model.quantile = sp_quantile 
        
## add some drawing method for some shapes 
for model in ( Ostap.Math.Bernstein         ,
               Ostap.Math.BernsteinEven     , 
               Ostap.Math.Positive          ,
               Ostap.Math.PositiveEven      ,
               Ostap.Math.Monothonic        ,
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
               Ostap.Math.MonothonicSpline  ,
               Ostap.Math.PositiveSpline    , 
               Ostap.Math.ConvexSpline      , 
               Ostap.Math.ConvexOnlySpline  ) : 
    
    model.draw = _f1_draw_
    model.Draw = _f1_draw_

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
Ostap.Math.Monothonic    .__repr__ = lambda s : _f_print_ ( s , 'Monothonic'    )
Ostap.Math.FourierSum    .__repr__ = lambda s : _f_print_ ( s , 'FourierSum'    )
Ostap.Math.CosineSum     .__repr__ = lambda s : _f_print_ ( s , 'CosineSum'     )


# =============================================================================
## decorate 2D-models/functions 
# =============================================================================
for model in ( Ostap.Math.Spline2D       ,
               Ostap.Math.Spline2DSym    , 
               Ostap.Math.Bernstein2D    ,
               Ostap.Math.Positive2D     ,
               Ostap.Math.Bernstein2DSym ,
               Ostap.Math.Positive2DSym  ,
               Ostap.Math.PS2DPol        ,
               Ostap.Math.PS2DPolSym     ,
               Ostap.Math.ExpoPS2DPol    ,
               Ostap.Math.Expo2DPol      ,
               Ostap.Math.Expo2DPolSym   ) :
    
    model . tf2 = _tf2_ 
    model.sp_integrate = sp_integrate_2D

for pdf in ( Ostap.Models.BreitWigner          , 
             Ostap.Models.Rho0               ,
             Ostap.Models.Kstar              ,
             Ostap.Models.Phi                , 
             Ostap.Models.Flatte             ,
             Ostap.Models.Flatte2            ,  
             Ostap.Models.Bukin              ,
             Ostap.Models.PhaseSpace2        ,
             Ostap.Models.PhaseSpaceNL       ,
             Ostap.Models.PhaseSpace23L      ,
             Ostap.Models.PhaseSpaceLeft     ,
             Ostap.Models.PhaseSpaceRight    ,
             Ostap.Models.PhaseSpacePol      ,
             Ostap.Models.Needham            ,
             Ostap.Models.CrystalBall        ,
             Ostap.Models.CrystalBallRS      ,
             Ostap.Models.CrystalBallDS      , 
             Ostap.Models.Apolonios          ,
             Ostap.Models.Apolonios2         , 
             Ostap.Models.GramCharlierA      , 
             Ostap.Models.Voigt              ,
             Ostap.Models.PseudoVoigt        ,
             Ostap.Models.Logistic           ,
             Ostap.Models.LASS               ,
             Ostap.Models.Bugg               ,
             Ostap.Models.LASS23L            ,
             Ostap.Models.Bugg23L            , 
             Ostap.Models.BW23L              , 
             Ostap.Models.PolyPositive       ,
             Ostap.Models.ExpoPositive       ,
             Ostap.Models.TwoExpoPositive    ,
             Ostap.Models.PositiveSpline     ,
             Ostap.Models.MonothonicSpline   , 
             
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
             Ostap.Models.Swanson            ,
             Ostap.Models.Argus              ,
             Ostap.Models.Tsallis            ,
             Ostap.Models.QGSM               ,
             Ostap.Models.BifurcatedGauss    ,
             Ostap.Models.GenGaussV1         , 
             Ostap.Models.GenGaussV2         , 
             ## Ostap.Models.SkewGauss          ## (temporarily removed) 
             ) :

    pdf.sp_integrate = sp_integrate_1D_

for pdf in ( Ostap.Models.Poly2DPositive     ,
             Ostap.Models.Poly2DSymPositive  , 
             Ostap.Models.PS2DPol            ,
             Ostap.Models.PS2DPolSym         , 
             Ostap.Models.ExpoPS2DPol        , 
             Ostap.Models.Expo2DPol          ,
             Ostap.Models.Expo2DPolSym       , 
             Ostap.Models.Spline2D           ,
             Ostap.Models.Spline2DSym        ) :
    
    pdf.sp_integrate = sp_integrate_2D_

    
# =============================================================================
## set parameter for polynomial/spline functions
#  @code
#  fun = ...
#  fun[1] = 10.0
#  @endcode 
def _p_set_par_ ( o , index , value ) :
    """Set parameter for polynomial/spline function
    >>> fun = ...
    >>> fun[1] = 10.0
    """
    return o.setPar ( index , value )

# =============================================================================
## get parameter from polynomial/spline functions
#  @code
#  fun = ...
#  print fun[1]
#  @endcode 
def _p_get_par_ ( o , index ) :
    """Get parameter from polynomial/spline function
    >>> fun = ...
    >>> print fun[1]
    """
    return o.par ( index )

# =============================================================================
## iterator over parameters of polynomial function
#  @code
#  fun = ...
#  for i in fun  : print fun[i] 
#  @endcode
def _p_iter_ ( o ) :
    """Iterator over parameters of polynomial function
    >>> fun = ...
    >>> for i in fun  : print fun[i] 
    """
    np = o.npars()
    for i in  range(np) :
        yield i
        
    
for f in ( Ostap.Math.Positive       ,
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
           ## 
           Ostap.Math.PolySum        ,
           ## 
           Ostap.Math.NSphere        ) :

    f.__setitem__  = _p_set_par_
    f.__getitem__  = _p_get_par_
    f.__len__      = lambda s     : s.npars() 
    f.__iter__     = _p_iter_
    f.__contains__ = lambda s , i : 0<=i<len(s)

# =============================================================================
##  Long polynomial division
#   f(x) = q(x)*g(x)+r(x), where  deg(f)=m >= def(g)=n, and
#   deg(q)<=m-n, deg(r)<n
#   @code
#   f = ... # the first polynom
#   g = ... # the second polynom
#   q,r = divmod ( f , g ) 
#   @endcode
def _b_divmod_ ( b , p ) :
    """Long polynomial division
    f(x) = q(x)*g(x)+r(x), where  deg(f)=m >= def(g)=n, and
    deg(q)<=m-n, deg(r)<n
    
    >>> f = ... # the first polynom
    >>> g = ... # the second polynom
    >>> q,r = divmod ( f , g ) 
    """
    rr = b.divmod(p)
    return rr.first, rr.second

# =============================================================================
##  Long polynomial division
#   f(x) = q(x)*g(x)+r(x), where  deg(f)=m >= def(g)=n, and
#   deg(q)<=m-n, deg(r)<n
#   @code
#   f = ... # the first polynom
#   g = ... # thr second  polynom
#   q = f // g # get quotient 
#   @endcode
def _b_floordiv_ ( b , p ) :
    """Long polynomial division
    f(x) = q(x)*g(x)+r(x), where  deg(f)=m >= def(g)=n, and
    deg(q)<=m-n, deg(r)<n
    
    >>> f = ... # the first polynom
    >>> g = ... # the second polynom
    >>> q = f // g 
    """
    rr = b.divmod(p)
    return rr.first

# =============================================================================
##  Long polynomial division
#   f(x) = q(x)*g(x)+r(x), where  deg(f)=m >= def(g)=n, and
#   deg(q)<=m-n, deg(r)<n
#   @code
#   f = ... # the first polynom
#   g = ... # the second  polynom
#   r = f % g # get reminder 
#   @endcode
def _b_mod_ ( b , p ) :
    """Long polynomial division
    f(x) = q(x)*g(x)+r(x), where  deg(f)=m >= def(g)=n, and
    deg(q)<=m-n, deg(r)<n
    
    >>> f = ... # the first polynom
    >>> g = ... # the second polynom
    >>> r = f % g 
    """
    rr = b.divmod(p)
    return rr.second

# =============================================================================
## power function for polynomials
#  @code
#  fun1 = ..
#  fun2 = fun1**3  
#  @endcode
def _b_pow_ ( b , n ) :
    """Power function for polynomials
    >>> fun1 = ..
    >>> fun2 = fun1**3  
    """
    if n != int(n) : raise ValueError('Illegal non-integer  exponent:%s' % n )
    n = int(n) 
    if 0  > n      : raise ValueError('Illegal negative     exponent:%s' % n )
    return b.pow(n)

Ostap.Math.Bernstein. __divmod__   = _b_divmod_
Ostap.Math.Bernstein. __floordiv__ = _b_floordiv_
Ostap.Math.Bernstein. __mod__      = _b_mod_
Ostap.Math.Bernstein. __pow__      = _b_pow_
# =============================================================================


# =============================================================================
## random generators 
# =============================================================================
from random import uniform as _uniform_

# =============================================================================
## generate random numbers from bernstein-like distribuitions
#  @code
#  >>> func = ...
#  >>> for x in func.generate( 1000 ) : print x 
#  @endcode
def _random_generate_bernstein_ ( fun , num ) :
    """Generate random numbers from bernstein-like distribuitions
    >>> func = ...
    >>> for x in func.generate( 1000 ) : print x 
    """
    xmn = fun.xmin ()
    xmx = fun.xmax ()
    ymx = max ( fun.bernstein().pars() )
    i   = 0 
    while i < num : 
        x = _uniform_ ( xmn , xmx ) 
        y = _uniform_ (   0 , ymx )
        v = fun ( x )
        if v >= y :
            i+= 1 
            yield x
            
# =============================================================================
## Get random number from bernstein-like distribuitions
#  @code
#  >>> func = ...
#  >>> print fun.shoot() 
#  @endcode
def _random_shoot_bernstein_ ( fun ) :
    """Get random number from bernstein-like distribuitions
    >>> func = ...
    >>> print func.shoot()  
    """
    xmn = fun.xmin ()
    xmx = fun.xmax ()
    ymx = max ( fun.bernstein().pars() )
    i   = 0 
    while True : 
        x = _uniform_ ( xmn , xmx ) 
        y = _uniform_ (   0 , ymx )
        v = fun ( x )
        if v >= y : return x 


# =============================================================================
## generate random numbers from b-spline-distribuitions
#  @code
#  >>> func = ...
#  >>> for x in func.generate( 1000 ) : print x 
#  @endcode
def _random_generate_bspline_ ( fun , num ) :
    """Generate random numbers from bspline-like distribuitions
    >>> func = ...
    >>> for x in func.generate( 1000 ) : print x 
    """
    xmn = fun.xmin ()
    xmx = fun.xmax ()
    ymx = max ( fun.bspline().pars() )
    i   = 0 
    while i < num : 
        x = _uniform_ ( xmn , xmx ) 
        y = _uniform_ (   0 , ymx )
        v = fun ( x )
        if v >= y :
            i+= 1 
            yield x
            
# =============================================================================
## Get random number from bspline-like distribuitions
#  @code
#  >>> func = ...
#  >>> print fun.shoot() 
#  @endcode
def _random_shoot_bspline_ ( fun ) :
    """Get random number from bspline-like distribuitions
    >>> func = ...
    >>> print func.shoot()  
    """
    xmn = fun.xmin ()
    xmx = fun.xmax ()
    ymx = max ( fun.bspline().pars() )
    i   = 0 
    while True : 
        x = _uniform_ ( xmn , xmx ) 
        y = _uniform_ (   0 , ymx )
        v = fun ( x )
        if v >= y : return x 

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
    zmx = max ( fun.bernstein().pars() )
    i   = 0 
    while i < num : 
        x = _uniform_ ( xmn , xmx ) 
        y = _uniform_ ( ymn , ymx ) 
        x = _uniform_ (   0 , zmx )
        v = fun ( x , y )
        if v >= z :
            i+= 1 
            yield x,y

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
    zmx = max ( fun.bernstein().pars() )
    i   = 0 
    while True : 
        x = _uniform_ ( xmn , xmx ) 
        y = _uniform_ ( ymn , ymx ) 
        z = _uniform_ (   0 , zmx )
        v = fun ( x , y )
        if v >= z : return x,y

Ostap.Math.Bernstein     .generate = _random_generate_bernstein_
Ostap.Math.Bernstein     .shoot    = _random_shoot_bernstein_

Ostap.Math.BernsteinEven .generate = _random_generate_bernstein_
Ostap.Math.BernsteinEven .shoot    = _random_shoot_bernstein_

Ostap.Math.Positive      .generate = _random_generate_bernstein_
Ostap.Math.Positive      .shoot    = _random_shoot_bernstein_

Ostap.Math.PositiveEven  .generate = _random_generate_bernstein_
Ostap.Math.PositiveEven  .shoot    = _random_shoot_bernstein_

Ostap.Math.PositiveSpline.generate = _random_generate_bspline_
Ostap.Math.PositiveSpline.shoot    = _random_shoot_bspline_

Ostap.Math.Positive2D    .generate = _random_generate_bernstein2D_
Ostap.Math.Positive2D    .shoot    = _random_shoot_bernstein2D_

Ostap.Math.Positive2DSym .generate = _random_generate_bernstein2D_
Ostap.Math.Positive2DSym .shoot    = _random_shoot_bernstein2D_

# =============================================================================
## add complex amplitudes 
# =============================================================================
Ostap.Math.LASS        . amp = _amp_
Ostap.Math.LASS23L     . amp = _amp_
Ostap.Math.Bugg23L     . amp = _amp_
Ostap.Math.Flatte      . amp = _amp_
Ostap.Math.Flatte2     . amp = _amp_
Ostap.Math.Flatte23L   . amp = _amp_
Ostap.Math.BreitWigner . amp = _amp_
    
# =============================================================================


import ostap.math.derivative as _D1
import ostap.math.integral   as _D2
for i in ( _D1.Derivative , _D2.Integral , _D2.IntegralCache ) :
    if not hasattr ( i , 'tf1' ) : i.tf1 = _tf1_
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
