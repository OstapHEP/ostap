#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/roofuncs.py
#  Set of useful functions for manipulations with RooAbsReal objects 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2020-03-08
# =============================================================================
""" Set of useful functions for manipulations with RooAbsReal objects 
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2020-03-08"
# =============================================================================
__all__     = (
    'BernsteinPoly'  , ## generic polynomial in Bernstein form     (RooAbsReal)
    'MonotonicPoly'  , ## monotonic polynomial                     (RooAbsReal)
    'ConvexPoly'     , ## monotonic convex/concave polynomial      (RooAbsReal)
    'ConvexOnlyPoly' , ## convex/concave polynomial                (RooAbsReal)
    'ScaleAndShift'  , ## scale and shift                          (RooAbsReal)
    ##
    'var_sum'        , ## sum                          for RooAbsReal objects           
    'var_mul'        , ## product                      for RooAbsReal objects           
    'var_sub'        , ## subtraction                  for RooAbsReal objects           
    'var_div'        , ## division                     for RooAbsReal objects           
    'var_fraction'   , ## fraction                     for RooAbsReal objects           
    'var_asymmetry'  , ## asymmetry                    for RooAbsReal objects           
    'var_pow'        , ## pow                 function for RooAbsReal objects           
    'var_abs'        , ## absolutevalue       function for RooAbsReal objects           
    'var_exp'        , ## exponent            function for RooAbsReal objects           
    'var_log'        , ## logarithm           function for RooAbsReal objects           
    'var_log10'      , ## logarithm           function for RooAbsReal objects           
    'var_erf'        , ## error               function for RooAbsReal objects           
    'var_sin'        , ## sine                function for RooAbsReal objects           
    'var_cos'        , ## cosine              function for RooAbsReal objects           
    'var_tan'        , ## tangent             function for RooAbsReal objects           
    'var_tanh'       , ## hyperbolic tangent  function for RooAbsReal objects           
    'var_atan2'      , ## inverse tangent     function for RooAbsReal objects           
    'var_gamma'      , ## gamma               function for RooAbsReal objects           
    'var_lgamma'     , ## logarithm of gamma  function for RooAbsReal objects           
    'var_igamma'     , ## 1/gamma             function for RooAbsReal objects
    'scale_var'      , ## var_mul
    'add_var'        , ## var_sum
    'sum_var'        , ## var_sum
    'ratio_var'      , ## var_div
    'fraction_var'   , ## var_fraction
    'asymmetry_var'  , ## var_asymmetry
   )
# =============================================================================
import ROOT, math
# =============================================================================
from   ostap.logger.logger          import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.roofuncs' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
from ostap.core.core                import Ostap 
from ostap.core.ostap_types         import num_types
from ostap.fitting.utils            import ParamsPoly , ShiftScalePoly, MakeVar  
from ostap.fitting.funbasic         import FUNC, Fun1D, Fun2D, Fun3D, make_fun 
# =============================================================================
## @class BernsteinPoly
#  Polynomial in Bernstein form
#  \f[ p(x) = \sum_{k=0}^{n} a_i B_n^k(x)  \f],
#  where \f$ B_n^k(x)\f$ is basic Bernstein polynomial   
#  @code
#  p1 =  BernsteinPoly ( 'P1' , xvar = (0,1) , power = 3 ) 
#  p2 =  BernsteinPoly ( 'P2' , xvar = (0,1) , pars  = ... ) 
#  @endcode 
#  @see Ostap::MoreRooFit::Bernstein
#  @see Ostap::Math::Bernstein
class BernsteinPoly ( FUNC , ParamsPoly ) :
    """Polynomial in Bernstein form
    >>> p1 =  BernsteinPoly ( 'P1' , xvar = (0,1) , power = 3   ) 
    >>> p2 =  BernsteinPoly ( 'P2' , xvar = (0,1) , pars  = ... ) 
    - see Ostap.MoreRooFit.Bernstein
    - see Ostap.Math.Bernstein
    """
    def __init__ ( self , name , xvar , power = 1 , pars = None ) :
        
        ## initialize the base class 
        FUNC      .__init__ ( self  , name  , xvar = xvar )
        ParamsPoly.__init__ ( self          ,
                              power = power ,
                              pars  = pars  )
        
        xmin , xmax = self.xminmax ()
        
        ## create the function
        self.fun    = Ostap.MoreRooFit.Bernstein (
            'bpol_%s'       % self.name ,
            'Bernstein(%s)' % self.name ,
            self.xvar                   ,
            xmin                        ,
            xmax                        ,
            self.pars_lst               ) 

        self.tricks = True 
        self.config = {
            'name' : self.name ,
            'xvar' : self.xvar ,
            'pars' : self.pars }

# =============================================================================
## @class MonotonicPoly
#  Simple monotonical polynomial 
#  \f[ p(x) = a + b P(x) \f], 
#  where \f$ P(x)\f$ is normalized positive monotonic polynomial:
#  - \f$ \int_{x_{\mathrm{min}}}^{x_{\mathrm{max}}} P(x) dx = 1 \f$ 
#  - \f$ P (x) \ge              0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$ 
#  - \f$ P^{\prime}(x) \ge(\le) 0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$
#  @code
#  p1 =  MonotonicPoly ( 'P1' , xvar = (0,1) , increasing = True  , power = 3 ) 
#  p2 =  MonotinicPoly ( 'P2' , xvar = (0,1) , increasing = False , pars  = ... ) 
#  @endcode 
#  @see Ostap::MoreRooFit::Monotonic
#  @see Ostap::Math::Monotonic
class MonotonicPoly(FUNC,ShiftScalePoly) :
    """Monotonic polynomial in Bernstein form
    >>> p1 =  MonotonicPoly ( 'P1' , xvar = (0,1) , increasing = True  , power = 3 ) 
    >>> p2 =  MonotonicPoly ( 'P2' , xvar = (0,1) , increasing = False , pars  = ... ) 
    - see Ostap.MoreRooFit.Monotonic
    """
    def __init__ ( self              ,
                   name              ,
                   xvar              ,
                   increasing = True ,
                   a          = 0.0  ,
                   b          = 1.0  ,  
                   power      = 1    ,
                   pars       = None ) :
        
        ## initialize the base class
        FUNC          .__init__ ( self  , name  , xvar = xvar ) 
        ShiftScalePoly.__init__ ( self          ,
                                  a     = a     ,
                                  b     = b     ,
                                  power = power ,
                                  pars  = pars  )
        
        self.__increasing = True if increasing else False

        xmin , xmax = self.xminmax ()

        ## create the function
        self.fun    = Ostap.MoreRooFit.Monotonic (
            'mpol_%s'       % self.name ,
            'Bernstein(%s)' % self.name ,
            self.xvar                   ,
            self.increasing             ,
            xmin                        ,
            xmax                        ,
            self.a                      ,
            self.b                      ,            
            self.pars_lst               ) 

        self.tricks = True 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'pars'       : self.pars       ,
            'a'          : self.a          ,
            'b'          : self.b          ,
            'increasing' : self.increasing }

    @property
    def increasing ( self ) :
        """``increasing'' : increasing polynomial?"""
        return self.__increasing
    
    @property
    def decreasing ( self ) :
        """``decreasing'' : decreasing polynomial?"""
        return not self.__increasing 

# =============================================================================
## @class ConvexPoly
#  Simple convex/concave polynomial 
#  \f[ p(x) = a + b P(x) \f], 
#  where \f$ P(x)\f$ is normalized positive polynomial:
#  - \f$ \int_{x_{\mathrm{min}}}^{x_{\mathrm{max}}} P(x) dx = 1 \f$ 
#  - \f$ P (x) \ge              0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$ 
#  - \f$ P^{\prime}(x) \ge(\le) 0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$
#  - \f$ P^{\prime\prime}(x) \ge(\le) 0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$ 
#  @code
#  p1 = ConvexPoly ( 'P1' , xvar = (0,1) , increasing = True  , convex = True  , power = 3   ) 
#  p2 = ConvexPoly ( 'P2' , xvar = (0,1) , increasing = False , convex = False , pars  = ... ) 
#  @endcode 
#  @see Ostap::MoreRooFit::Convex
#  @see Ostap::Math::Convex
class ConvexPoly(FUNC,ShiftScalePoly) :
    """Monotonic polynomial in Bernstein form
    >>> p1 = ConvexPoly ( 'P1' , xvar = (0,1) , increasing = True  , convex = True  , power = 3 ) 
    >>> p2 = ConvexPoly ( 'P2' , xvar = (0,1) , increasing = False , convex = False , pars  = ... ) 
    - see Ostap.MoreRooFit.Monotonic
    """
    def __init__ ( self              ,
                   name              ,
                   xvar              ,
                   increasing = True ,
                   convex     = True ,
                   a          = 0.0  ,
                   b          = 1.0  ,  
                   power      = 1    ,
                   pars       = None ) :
        
        ## initialize the base class
        FUNC          .__init__ ( self , name   , xvar = xvar ) 
        ShiftScalePoly.__init__ ( self          ,
                                  a     = a     ,
                                  b     = b     ,
                                  power = power ,
                                  pars  = pars  )
        
        self.__increasing = True if increasing else False
        self.__convex     = True if convex     else False

        xmin , xmax = self.xminmax ()
        
        ## create the function
        self.fun    = Ostap.MoreRooFit.Convex (
            'cpol_%s'       % self.name ,
            'Bernstein(%s)' % self.name ,
            self.xvar                   ,
            self.increasing             ,
            self.iconvex                ,
            xmin                        ,
            xmax                        ,
            self.a                      ,
            self.b                      ,            
            self.pars_lst               ) 

        self.tricks = True 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'pars'       : self.pars       ,
            'a'          : self.a          ,
            'b'          : self.b          ,
            'convex'     : self.convex     ,
            'increasing' : self.increasing }

    @property
    def increasing ( self ) :
        """``increasing'' : increasing polynomial?"""
        return self.__increasing
    
    @property
    def decreasing ( self ) :
        """``decreasing'' : decreasing polynomial?"""
        return not self.__increasing 
    
    @property
    def convex ( self ) :
        """``convex'' : convex polynomial?"""
        return self.__convex 

    @property
    def concave ( self ) :
        """``concave'' : concave polynomial?"""
        return not self.__convex 

# =============================================================================
## @class ConvexOnlyPoly
#  Simple convex/concave polynomial 
#  \f[ p(x) = a + b P(x) \f], 
#  where \f$ P(x)\f$ is normalized positive polynomial:
#  - \f$ \int_{x_{\mathrm{min}}}^{x_{\mathrm{max}}} P(x) dx = 1 \f$ 
#  - \f$ P (x) \ge              0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$ 
#  - \f$ P^{\prime\prime}(x) \ge(\le) 0 \f$ for \f$ {x_{\mathrm{min}}}\lex\le {x_{\mathrm{max}}}\f$ 
#  @code
#  p1 =  ConvexOnlyPoly ( 'P1' , xvar = (0,1) , convex = True  , power = 3 ) 
#  p2 =  ConvexOnlyPoly ( 'P2' , xvar = (0,1) , convex = False , pars  = ... ) 
#  @endcode 
#  @see Ostap::MoreRooFit::Monotonic
#  @see Ostap::Math::Monotonic
class ConvexOnlyPoly(FUNC,ShiftScalePoly) :
    """Monotonic polynomial in Bernstein form
    >>> p1 =  ConvexOnlyPoly ( 'P1' , xvar = (0,1) , convex = True  , power = 3 ) 
    >>> p2 =  ConvexOnlyPoly ( 'P2' , xvar = (0,1) , convex = False , pars  = ... ) 
    - see Ostap.MoreRooFit.Monotonic
    """
    def __init__ ( self              ,
                   name              ,
                   xvar              ,
                   convex     = True ,
                   a          = 0.0  ,
                   b          = 1.0  ,  
                   power      = 1    ,
                   pars       = None ) :
        
        ## initialize the base class
        FUNC          .__init__ ( self  , name  , xvar = xvar ) 
        ShiftScalePoly.__init__ ( self          ,
                                  a     = a     ,
                                  b     = b     ,
                                  power = power ,
                                  pars  = pars  )
        
        self.__convex     = True if convex     else False

        xmin , xmax = self.xminmax ()
        
        ## create the function
        self.fun    = Ostap.MoreRooFit.ConvexOnly (
            'cpol_%s'       % self.name ,
            'Bernstein(%s)' % self.name ,
            self.xvar                   ,
            self.iconvex                ,
            xmin                        ,
            xmax                        ,
            self.a                      ,
            self.b                      ,            
            self.pars_lst               ) 

        self.tricks = True 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'pars'       : self.pars       ,
            'a'          : self.a          ,
            'b'          : self.b          ,
            'convex'     : self.convex     }

    @property
    def convex ( self ) :
        """``convex'' : convex polynomial?"""
        return self.__convex 

    @property
    def concave ( self ) :
        """``concave'' : concave polynomial?"""
        return not self.__convex 

# =============================================================================
# @class ScaleAndShift
# Scale and shift another function/PDF 
# \f[ f = a + b c \f\
# @see Ostap::MoreRooFit::ScaleAndShift 
class ScaleAndShift ( FUNC ) :
    """ Scale and shift another function/PDF: 
    f = a+b*func
    - see Ostap.MoreRooFit.ScaleAndShift 
    """
    def __init__ ( self     ,
                   name     ,
                   xvar     ,
                   func     ,   ## c 
                   a  = 0.0 ,   ## a 
                   b  = 1.0 ) : ## b 

        FUNC.__init__ ( self , name , xvar = xvar )

        assert isinstance ( func , ( ROOT.RooAbsReal, FUNC ) ) ,\
               "Invalid type of ``func'' %s/%s" % ( func , type ( func ) )

        self.__func0 = func

        self.__c  = func.fun if isnstance ( func , FUNC ) else func 

        self.__a  = self.make_var ( a ,
                                    "a_%s"              % self.name ,
                                    "shift/bias for %s" % self.name , False , a ) 
        self.__b  = self.make_var ( b ,
                                    "b_%s"              % self.name ,
                                    "scale for %s"      % self.name , False , b ) 


        if not self.xvar in self.a.getParameters ( 0 ) and \
           not self.xvar in self.b.getParameters ( 0 ) and \
           not self.xvar in self.c.getParameters ( 0 ) :
            self.warning ( "Function does not depend on xvar=%s" % self.xvar.name )
           
        self.__bc = Ostap.MoreRooFit.Product  ( self.b , self.c )
        self.fun  = Ostap.MoreRooFit.Addition (
            "scaleshift_%s"   % self.name ,
            "Scale&Shift(%s)" % self.name ,
            self.a      ,
            self.__bc  
            )
        
        self.config = {
            'name' : self.name ,
            'xvar' : self.xvar ,
            'a'    : self.a    ,
            'b'    : self.b    ,
            'func' : self.c    }

    @property
    def a ( self ) :
        """``a'' : bias parameter for polynomial:  f(x) = a + b*c"""
        return self.__a
    @a.setter
    def a ( self , value ) :
        vv = float ( value )
        if self.__a.minmax () and not vv in self.__a  :
            self.error ("Value %s is outside the allowed region %s"  % ( vv , self.__a.minmax() ) )
        self.__a.setVal ( vv )

    @property
    def b ( self ) :
        """``scale'' : bias parameter for polynomial:  f(x) = a + b*c"""
        return self.__b
    @b.setter
    def b ( self , value ) :
        vv = float ( value )
        if self.__b.minmax () and not vv in self.__b  :
            self.error ("Value %s is outside the allowed region %s"  % ( vv , self.__b.minmax() ) )
        self.__b.setVal ( vv )
    
    @property
    def c ( self ) :
        """``a'' : bias parameter for polynomial:  f(x) = a + b*c"""
        return self.__c

# =============================================================================
# @class FAB
# helper Base class 
class FAB(FUNC):
    """ Helper base class 
    """
    def __init__ ( self                  ,
                   a                     ,  ## a                   
                   xvar    = None        ,
                   b       = 1.0         ,  ## b 
                   name    = ''          ,
                   pattern = 'FAB_%s_%s' ) :


        self.__a0 = a
        self.__b0 = b
        
        if instance ( a , FUNC ) :
            a    = a.fun
            if not xvar : xvar = a.xvar
            
        if not name :
            if hasattr ( b , 'name' ) : name = pattern % ( a.name , b.name  )
            else                      : name = pattern % ( a.name , b       )

        FUNC.__init__ ( self , name , xvar = xvar )

        if isinstance ( b , FUNC ) : b = b.fun
                
        self.__a  = a
        self.__b  = self.make_var ( b ,
                                    "b_%s"      % self.name ,
                                    "B for %s"  % self.name , False , b ) 
        
        if not self.xvar in self.a.getParameters ( 0 ) and \
           not self.xvar in self.b.getParameters ( 0 ) and \
           not self.xvar is self.a                     and \
           not self.xvar is self.b :                    
            self.warning ( "Function does not depend on xvar=%s" % self.xvar.name )
        
        self.config = {
            'name' : self.name ,
            'xvar' : self.xvar ,
            'a'    : self.a    ,
            'b'    : self.b    }

    @property
    def a ( self ) :
        """``a'' : f(x) = F(ab) """
        return self.__a

    @property
    def b ( self ) :
        """``b''  : f(x) = F(ab)"""
        return self.__b


_b_types = num_types + ( FUNC , ROOT.RooAbsReal ) 
# ===============================================================================
def _fn_make_fun_ ( afun         ,
                    b            ,
                    cpptyp       ,
                    namepat      ,
                    swap = False ) :
    
    
    assert isinstance  ( b , _b_types ) , 'Invalid argument type %s' % type ( b )
    
    vars = [ v for v in afun.variables ] 
    
    a = afun.fun
    
    if   isinstance ( b , num_types ) :
        b = ROOT.RooRealConstant.value ( b )
    elif insistance ( b , FUNC ) :
        vars += [ v for v in b.variables if not v in vars ]
        b = b.fun

        
    afun.aux_keep.append ( b )  

    if swap :
        a , b = b , a
        
    ## construct the function
    result = cpptyp ( a , b )

    if namepat : name = namepat % ( a.name , b.name )
    else       : name = result.name 

    if 3 < len ( vars ) :
        afun.warning ( "Skip extra variables: %s" % vars[3:] )
        vars = vars[:3]

    return make_fun ( result , vars , name ) 

# ==============================================================================
## Absolute value for the function  \f$ f = \left| ab \right| \f$
#  @code
#  f =
#  a = f.abs (   )
#  a = f.abs ( b )
#  a =   abs ( f )  
#  @endcode 
def  _fn_abs_ ( self , b = 1 ) :
    """Absolute value for the function:: f = abs(ab)
    >>> f =
    >>> a = f.abs (   )
    >>> a = f.abs ( b )
    >>> a =   abs ( f )  
    """
    return _fn_make_fun_ ( self                  ,
                           b                     ,
                           Ostap.MoreRooFit.Abs  , 
                           'abs_%s_%s'           )

# ==============================================================================
## Exponent \f$ f = {\mathrm{e}}^{ab} \f$
#  @code
#  f =
#  a = f.exp (   )
#  a = f.exp ( b )
#  a =   exp ( f )  
#  @endcode 
def  _fn_exp_ ( self , b = 1 ) :
    """ Exponent f = exp(ab)
    >>> f =
    >>> a = f.exp (   )
    >>> a = f.exp ( b )
    >>> a =   exp ( f )  
    """
    return _fn_make_fun_ ( self                  ,
                           b                     ,
                           Ostap.MoreRooFit.Exp  , 
                           'exp_%s_%s'           )

# ==============================================================================
## Natural logarithm  \f$ f = \log ab  \f$
#  @code
#  f =
#  a = f.log (   )
#  a = f.log ( b )
#  a =   log ( f )  
#  @endcode 
def  _fn_log_ ( self , b = 1 ) :
    """ Natural logarithm  f = log(ab)
    >>> f =
    >>> a = f.log (   )
    >>> a = f.log ( b )
    >>> a =   log ( f )  
    """
    return _fn_make_fun_ ( self                  ,
                           b                     ,
                           Ostap.MoreRooFit.Log  , 
                           'log_%s_%s'           )

# ==============================================================================
## Decimal logarithm  \f$ f = \log_{10} ab  \f$
#  @code
#  f =
#  a = f.log10 (   )
#  a = f.log10 ( b )
#  a =   log10 ( f )  
#  @endcode 
def  _fn_log10_ ( self , b = 1 ) :
    """ Decimal logarithm  f = log10(ab)
    >>> f =
    >>> a = f.log10 (   )
    >>> a = f.log10 ( b )
    >>> a =   log10 ( f )  
    """
    return _fn_make_fun_ ( self                    ,
                           b                       ,
                           Ostap.MoreRooFit.Log10  , 
                           'log10_%s_%s'           )

# ==============================================================================
## Error function \f$ f = erf (ab)  \f$
#  @code
#  f =
#  a = f.erf (   )
#  a = f.erf ( b )
#  a =   erc ( f )  
#  @endcode 
def  _fn_erf_ ( self , b = 1 ) :
    """ Error function f = erf(ab)
    >>> f =
    >>> a = f.erf (   )
    >>> a = f.erf ( b )
    >>> a =   erf ( f )  
    """
    return _fn_make_fun_ ( self                  ,
                           b                     ,
                           Ostap.MoreRooFit.Erf  , 
                           'erf_%s_%s'           )

# ==============================================================================
## Sine function: \f$ f = sin (ab)  \f$
#  @code
#  f =
#  a = f.sin (   )
#  a = f.sin ( b )
#  a =   sin ( f )  
#  @endcode 
def  _fn_sin_ ( self , b = 1 ) :
    """ Sine function: f = sin(ab)
    >>> f =
    >>> a = f.sin (   )
    >>> a = f.sin ( b )
    >>> a =   sin ( f )  
    """
    return _fn_make_fun_ ( self                  ,
                           b                     ,
                           Ostap.MoreRooFit.Sin  , 
                           'sin_%s_%s'           )

# ==============================================================================
## Cosine function: \f$ f = cos (ab)  \f$
#  @code
#  f =
#  a = f.cos (   )
#  a = f.cos ( b )
#  a =   cos ( f )  
#  @endcode 
def  _fn_cos_ ( self , b = 1 ) :
    """ Cosine function: f = cos(ab)
    >>> f =
    >>> a = f.cos (   )
    >>> a = f.cos ( b )
    >>> a =   cos ( f )  
    """
    return _fn_make_fun_ ( self                  ,
                           b                     ,
                           Ostap.MoreRooFit.Cos  , 
                           'cos_%s_%s'           )


# ==============================================================================
## Tagent function: \f$ f = tan (ab)  \f$
#  @code
#  f =
#  a = f.tan (   )
#  a = f.tan ( b )
#  a =   tan ( f )  
#  @endcode 
def  _fn_tan_ ( self , b = 1 ) :
    """ Tangent function: f = tan(ab)
    >>> f =
    >>> a = f.tan (   )
    >>> a = f.tan ( b )
    >>> a =   tan ( f )  
    """
    return _fn_make_fun_ ( self                  ,
                           b                     ,
                           Ostap.MoreRooFit.Tan  , 
                           'tan_%s_%s'           )


# ==============================================================================
## Hyperbolic Tagent function: \f$ f = tanh (ab)  \f$
#  @code
#  f =
#  a = f.tanh (   )
#  a = f.tanh ( b )
#  a =   tanh ( f )  
#  @endcode 
def  _fn_tanh_ ( self , b = 1 ) :
    """ Hyperbolic tangent function: f = tan(ab)
    >>> f =
    >>> a = f.tanh (   )
    >>> a = f.tanh ( b )
    >>> a =   tanh ( f )  
    """
    return _fn_make_fun_ ( self                   ,
                           b                      ,
                           Ostap.MoreRooFit.Tanh  , 
                           'tanh_%s_%s'           )


# ==============================================================================
## Inverse agent function: \f$ f = atan2 ( a , b)  \f$
#  @code
#  f =
#  a = f.atan2 (   )
#  a = f.atan2 ( b )
#  a =   atan2 ( f )  
#  @endcode 
def  _fn_atan2_ ( self , b = 1 ) :
    """ Inverse tangent function: f = atan2(a,b)
    >>> f =
    >>> a = f.atan2 (   )
    >>> a = f.atan2 ( b )
    >>> a =   atan2 ( f )  
    """
    return _fn_make_fun_ ( self                    ,
                           b                       ,
                           Ostap.MoreRooFit.Atan2  , 
                           'atan2_%s_%s'           )

# ==============================================================================
## Gamma function: \f$ f = \Gamma  ( ab )  \f$
#  @code
#  f =
#  a = f.tgamma (   )
#  a = f.tgamma ( b )
#  a =   tgamma ( f )  
#  @endcode 
def  _fn_tgamma_ ( self , b = 1 ) :
    """ Gamma function: f = Gamma(a,b)
    >>> f =
    >>> a = f.tgamma (   )
    >>> a = f.tgamma ( b )
    >>> a =    gamma  ( f )  
    """
    return _fn_make_fun_ ( self                    ,
                           b                       ,
                           Ostap.MoreRooFit.Gamma  , 
                           'gamma_%s_%s'           )


# ==============================================================================
## log-Gamma function: \f$ f = \log \Gamma  ( ab )  \f$
#  @code
#  f =
#  a = f.lgamma (   )
#  a = f.lgamma ( b )
#  a =   lgamma ( f )  
#  @endcode 
def  _fn_lgamma_ ( self , b = 1 ) :
    """ Gamma function: f = log(Gamma(ab))
    >>> f =
    >>> a = f.lgamma (   )
    >>> a = f.lgamma ( b )
    >>> a =   lgamma  ( f )  
    """
    return _fn_make_fun_ ( self                     ,
                           b                        ,
                           Ostap.MoreRooFit.LGamma  , 
                           'lgamma_%s_%s'           )

# ==============================================================================
## 1/Gamma function: \f$ f = 1/\Gamma  ( ab )  \f$
#  @code
#  f =
#  a = f.igamma (   )
#  a = f.igamma ( b )
#  a =   igamma ( f )  
#  @endcode 
def  _fn_igamma_ ( self , b = 1 ) :
    """ 1/Gamma function: f = 1/Gamma(ab))
    >>> f =
    >>> a = f.igamma (   )
    >>> a = f.igamma ( b )
    >>> a =   igamma  ( f )  
    """
    return _fn_make_fun_ ( self                     ,
                           b                        ,
                           Ostap.MoreRooFit.IGamma  , 
                           'igamma_%s_%s'           )


# ==============================================================================
## Power function: \f$ f = pow( a , b )  \f$
#  @code
#  f =
#  a = f.pow ( b )
#  a =  a ** b   
#  @endcode 
def  _fn_pow_ ( self , b ) :
    """ Power function: f = pow( a, b )  )
    >>> f =
    >>> a = f.pow ( b )
    >>> a =  f ** b 
    """
    return _fn_make_fun_ ( self                    ,
                           b                       ,
                           Ostap.MoreRooFit.Power  , 
                           'pow_%s_%s'              )

# ==============================================================================
## Right power function: \f$ f = pow( b , a )  \f$
#  @code
#  f =
#  a = f.rpow ( b )
#  a =   b**f    
#  @endcode 
def  _fn_rpow_ ( self , b ) :
    """ Righ power function: f = pow( b, a )  )
    >>> f =
    >>> a = f.rpow  ( b )
    >>> a =   b ** f   
    """
    return _fn_make_fun_ ( self                    ,
                           b                       ,
                           Ostap.MoreRooFit.Power  , 
                           'rpow_%s_%s'             , swap = True  )

# ==============================================================================
## Right power function: \f$ f = pow( b , a )  \f$
#  @code
#  f =
#  a = f.rpow ( b )
#  a =   b**f    
#  @endcode 
def  _fn_rpow2_ ( self , b  ) :
    """ Righ power function: f = pow( b, a )  )
    >>> f =
    >>> a = f.rpow  ( b )
    >>> a =   b ** f   
    """

    if not isinstance  ( b , _b_types ) : return NotImplemented 

    return _fn_rpow_ ( self , b  )


# ==============================================================================
## Fraction  \f$ f =  a / ( a + b ) \f$
#  @code
#  f =
#  a = f.fraction ( b )
#  @endcode 
def  _fn_fraction_ ( self , b , swap = False ) :
    """ Fraction : f = a / ( a+ b )   )
    >>> f =
    >>> a = f.fraction( b )
    """
    return _fn_make_fun_ ( self                      ,
                           b                         ,
                           Ostap.MoreRooFit.Fraction , 
                           'frac_%s_%s'              , swap = swap )


# ==============================================================================
## Asymmetry  \f$ f =  ( a - b ) / ( a + b ) \f$
#  @code
#  f =
#  a = f.asymmetry ( b )
#  @endcode 
def  _fn_asymmetry_ ( self , b , swap = False ) :
    """ Asymmetry  : f = ( a - b ) / ( a+ b )   )
    >>> f =
    >>> a = f.fraction( b )
    """
    return _fn_make_fun_ ( self                       ,
                           b                          ,
                           Ostap.MoreRooFit.Asymmetry , 
                           'asymm_%s_%s'              , swap = swap )



# =============================================================================
FUNC.__abs__     = _fn_abs_
FUNC.__exp__     = _fn_exp_
FUNC.__log__     = _fn_log_
FUNC.__log10__   = _fn_log10_
FUNC.__erf__     = _fn_erf_
FUNC.__sin__     = _fn_sin_
FUNC.__cos__     = _fn_cos_
FUNC.__tan__     = _fn_tan_
FUNC.__tanh__    = _fn_tanh_
FUNC.__atan2__   = _fn_atan2_
FUNC.__tgamma__  = _fn_tgamma_
FUNC.__lgamma__  = _fn_lgamma_
FUNC.__igamma__  = _fn_igamma_
FUNC.__pow__     = _fn_pow_
FUNC.__rpow__    = _fn_rpow2_

FUNC.  abs       = _fn_abs_
FUNC.  exp       = _fn_exp_
FUNC.  log       = _fn_log_
FUNC.  log10     = _fn_log10_
FUNC.  erf       = _fn_erf_
FUNC.  sin       = _fn_sin_
FUNC.  cos       = _fn_cos_
FUNC.  tan       = _fn_tan_
FUNC.  tanh      = _fn_tanh_
FUNC.  atan2     = _fn_atan2_
FUNC.  tgamma    = _fn_tgamma_
FUNC.  lgamma    = _fn_lgamma_
FUNC.  igamma    = _fn_igamma_
FUNC.  pow       = _fn_pow_
FUNC.  rpow      = _fn_rpow_

FUNC.  fraction  = _fn_fraction_
FUNC.  asymmetry = _fn_asymmetry_


# ==============================================================================
## absolute value   \f$ f = abs{ab}\f$
#  @code
#  var = ...
#  e   = var_abs ( var ) 
#  @endcode 
def var_abs ( a , b = 1 , name = '' , title = '' ) :
    """Absolute value: f(x) = abs(ab)
    >>> var = ...
    >>> e   = var_abs ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.abs ( float ( a ) * float ( b ) )
        return ROOT.RooRealConstant.value ( ab )     
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.Abs( a, b , name , title )
    #

# ==============================================================================
## exponent  \f$ f = \mathrm{e}^{ab}\f$
#  @code
#  var = ...
#  e   = var_exp ( var ) 
#  @endcode 
def var_exp ( a , b = 1 , name = '' , title = '' ) :
    """Exponent: f(x) = exp(ab)
    >>> var = ...
    >>> e   = var_exp ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.exp ( float ( a ) * float ( b ) )
        return ROOT.RooRealConstant.value ( ab )          ## RETURN
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.Exp ( a, b , name , title )

# ==============================================================================
## logarithm  \f$ f = \log ab \f$
#  @code
#  var = ...
#  e   = var_log ( var ) 
#  @endcode 
def var_log ( a , b = 1 , name = '' , title = '' ) :
    """logarithm f(x) = log(ab)
    >>> var = ...
    >>> e   = var_log ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.log ( float ( a ) * float ( b ) )       ## RETURN
        return ROOT.RooRealConstant.value ( ab )
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.Log ( a, b , name , title ) 

# ==============================================================================
## logarithm  \f$ f = \log10 ab \f$
#  @code
#  var = ...
#  e   = var_log10 ( var ) 
#  @endcode 
def var_log10 ( a , b = 1 , name = '' , title = '' ) :
    """logarithm f(x) = log10(ab)
    >>> var = ...
    >>> e   = var_log10 ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.log10 ( float ( a ) * float ( b ) )       ## RETURN
        return ROOT.RooRealConstant.value ( ab )
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.Log10 ( a, b , name , title ) 

# ==============================================================================
## error function \f$ f = erf ( ab) \f$
#  @code
#  var = ...
#  e   = var_erf ( var ) 
#  @endcode 
def var_erf ( a , b = 1 , name = '' , title = '' ) :
    """Error function f(x) = erf(ab)
    >>> var = ...
    >>> e   = var_erf ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.erf ( float ( a ) * float ( b ) )
        return ROOT.RooRealConstant.value ( ab )          ## RETURN
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.Erf ( a, b , name , title ) 

# ==============================================================================
## Sine \f$ f = \sin ab\f$
#  @code
#  var = ...
#  e   = var_sin ( var ) 
#  @endcode 
def var_sin ( a , b = 1 , name = '' , title = '' ) :
    """Sine  f(x) = sin(ab)
    >>> var = ...
    >>> e   = var_sin ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.sin ( float ( a ) * float ( b ) )
        return ROOT.RooRealConstant.value ( ab )          ## RETURN 
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.Sin ( a, b , name , title ) 
    

# ==============================================================================
## Cosine \f$ f = \cos ab\f$
#  @code
#  var = ...
#  e   = var_cos ( var ) 
#  @endcode 
def var_cos ( a , b = 1 , name = '' , title = '' ) :
    """Cosine  f(x) = cos(ab)
    >>> var = ...
    >>> e   = var_cos ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.cos ( float ( a ) * float ( b ) )
        return ROOT.RooRealConstant.value ( ab )       ## RETURN
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.Cos ( a, b , name , title )

# ==============================================================================
## Tangent\f$ f = \tan ab\f$
#  @code
#  var = ...
#  e   = var_tan ( var ) 
#  @endcode 
def var_tan ( a , b = 1 , name = '' , title = '' ) :
    """Tangent  f(x) = tan(ab)
    >>> var = ...
    >>> e   = var_tan ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.tan ( float ( a ) * float ( b ) )   ##  RETURN
        return ROOT.RooRealConstant.value ( ab )
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.Tan ( a, b , name , title )

# ==============================================================================
## Hyperboilic tangent\f$ f = \tanh ab\f$
#  @code
#  var = ...
#  e   = var_tanh ( var ) 
#  @endcode 
def var_tanh ( a , b = 1 , name = '' , title = '' ) :
    """Hyperbolic tangent  f(x) = tan(ab)
    >>> var = ...
    >>> e   = var_tanh ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.tanh ( float ( a ) * float ( b ) ) 
        return ROOT.RooRealConstant.value ( ab )      ## RETURN
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.Tanh ( a, b , name , title )

# ==============================================================================
## arctangent\f$ f = atan2 (a,b)\f$
#  @code
#  var = ...
#  e   = var_atan2 ( var ) 
#  @endcode 
def var_atan2 ( a , b = 1 , name = '' , title = '' ) :
    """Inverse tangent  f(x) = atan2(a,b)
    >>> var = ...
    >>> e   = var_atan2 ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.tanh ( float ( a ) * float ( b ) ) 
        return ROOT.RooRealConstant.value ( ab )      ## RETURN
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.Atan2 ( a, b , name , title )

# ==============================================================================
## Gamma function \f$ f =    \Gamma(ab) \f$
#  @code
#  a = ...
#  e = var_gamma ( a ) 
#  @endcode 
def var_gamma ( a , b = 1 , name = '' , title = '' ) :
    """Gamma function  f = Gamma(ab)
    >>> a = ...
    >>> e = var_gamma ( a ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.gamma ( float ( a ) * float ( b ) )      
        return ROOT.RooRealConstant.value ( ab )           ## RETURN 
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.Gamma ( a, b , name , title ) 

# ==============================================================================
## logarithm of Gamma function \f$ f = \log \Gamma(ab) \f$
#  @code
#  a = ...
#  e = var_lgamma ( a ) 
#  @endcode 
def var_lgamma ( a , b = 1 , name = '' , title = '' ) :
    """logarithm of Gamma function  f = log Gamma(ab)
    >>> a = ...
    >>> e = var_lgamma ( a ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.lgamma ( float ( a ) * float ( b ) )      
        return ROOT.RooRealConstant.value ( ab )           ## RETURN 
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.LGamma ( a, b , name , title ) 

# ==============================================================================
## 1/Gamma function \f$ f = \frac{1}{\Gamma(ab)} \f$
#  @code
#  a = ...
#  e = var_igamma ( a ) 
#  @endcode 
def var_igamma ( a , b = 1 , name = '' , title = '' ) :
    """1/Gamma  f = 1/Gamma(ab)
    >>> a = ...
    >>> e = var_igamma ( a ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = Ostap.Math.igamma ( float ( a ) * float ( b ) )
        return ROOT.RooRealConstant.value ( ab )              ## RETURN 
    elif fa : a = ROOT.RooRealConstant.value ( a )
    elif fb : b = ROOT.RooRealConstant.value ( b )
    #
    return Ostap.MoreRooFit.IGamma ( a, b , name , title ) 

# =============================================================================
## Sum of two RooAbsReal objects
#  @code
#  v1 = ...
#  v2 = ...
#  v  = var_sum ( v1 ,  v2 )
#  @endcode
def var_sum ( v1 , v2 , name = '' , title = '' ) :
    """Sum of two RooAbsReal objects
    >>> v1 = ...
    >>> v2 = ...
    >>> v  = var_sum ( v1 , v2 )
    """
    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )
    
    if f1 and f2 :
        r = float ( v1 ) + float ( v2 ) 
        return ROOT.RooRealConstant.value ( r )                 ## RETURN
    
    elif f1 : v1 = ROOT.RooRealConstant.value ( float ( v1 ) )        
    elif f2 : v2 = ROOT.RooRealConstant.value ( float ( v2 ) ) 
    #
    return Ostap.MoreRooFit.Addition ( v1 , v2 , name , title )

# =============================================================================
## Product of two RooAbsReal objects
#  @code
#  v1 = ...
#  v2 = ...
#  v  = var_mul ( v1 ,  v2 )
#  @endcode
def var_mul ( v1 , v2 , name = '' , title = '' ) :
    """Product of two RooAbsReal objects
    >>> v1 = ...
    >>> v2 = ...
    >>> v  = var_mul ( v1 ,  v2 )
    """

    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )
    
    if f1 and f2 :
        r = float ( v1 ) * float ( v2 ) 
        return ROOT.RooRealConstant.value ( r )                 ## RETURN
    
    elif f1 : v1 = ROOT.RooRealConstant.value ( float ( v1 ) )        
    elif f2 : v2 = ROOT.RooRealConstant.value ( float ( v2 ) ) 
    #
    return Ostap.MoreRooFit.Product ( v1 , v2 , name , title ) 

# =============================================================================
## Subtraction of two RooAbsReal objects
#  @code
#  v1 = ...
#  v2 = ...
#  v  = var_sub ( v1 , v2 )
#  @endcode
def var_sub ( v1 , v2 , name = '' , title = '' ) :
    """Subraction of two RooAbsReal objects
    >>> v1 = ...
    >>> v2 = ...
    >>> v  = var_sub ( v1 , v2 )  
    """

    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )
    
    if f1 and f2 :
        r = float ( v1 ) - float ( v2 ) 
        return ROOT.RooRealConstant.value ( r )                 ## RETURN
    
    elif f1 : v1 = ROOT.RooRealConstant.value ( float ( v1 ) )        
    elif f2 : v2 = ROOT.RooRealConstant.value ( float ( v2 ) ) 
    # 
    return Ostap.MoreRooFit.Subtraction ( v1 , v2 , name , title ) 

# =============================================================================
## Division of two RooAbsReal objects
#  @code
#  v1 = ...
#  v2 = ...
#  v  = var_div ( v1 - v2 )
#  @endcode
def var_div ( v1 , v2 , name = '' , title = '' ) :
    """Division of two RooAbsReal objects
    >>> v1 = ...
    >>> v2 = ...
    >>> v  = var_div ( v1 , v2 )  
    """
    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )
    
    if f1 and f2 :
        r = float ( v1 ) / float ( v2 ) 
        return ROOT.RooRealConstant.value ( r )                 ## RETURN
    
    elif f1 : v1 = ROOT.RooRealConstant.value ( float ( v1 ) )        
    elif f2 : v2 = ROOT.RooRealConstant.value ( float ( v2 ) ) 
    #
    return Ostap.MoreRooFit.Division ( v1 , v2 , name , title ) 

# ==============================================================================
## "Fraction" of two RooAbsReal objects: f = a/(a+b)
#  @code
#  a = ...
#  b = ...
#  e   = var_fraction ( a , b ) 
#  @endcode 
def var_fraction ( a , b , name = '' , title = '' ) :
    """``Fraction''  f(x) = a/(a+b)
    >>> a = ...
    >>> b = ...
    >>> e = var_fraction ( a , b  ) 
    """
    
    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )
    
    if f1 and f2 :
        r = float ( v1) / ( float ( v1 ) + float ( v1 ) )  
        return ROOT.RooRealConstant.value ( r )                 ## RETURN
    
    elif f1 : v1 = ROOT.RooRealConstant.value ( float ( v1 ) )        
    elif f2 : v2 = ROOT.RooRealConstant.value ( float ( v2 ) ) 
    #
    return Ostap.MoreRooFit.Fraction ( v1 , v2 , name , title ) 


# ==============================================================================
## "Asymmetry" of two RooAbsReal objects: f = (a-b)/(a+b)
#  @code
#  a = ...
#  b = ...
#  e = var_asymmetry ( a , b ) 
#  @endcode 
def var_asymmetry ( a , b , name = '' , title = '' ) :
    """``Asymmetry''  f(x) = (a-b)/(a+b)
    >>> a = ...
    >>> b = ...
    >>> e = var_asymmetry ( a , b  ) 
    """

    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )
    
    if f1 and f2 :
        r = ( float ( v1 ) - float ( v2 ) ) / ( float ( v1 ) + float ( v2 ) )  ## 
        return ROOT.RooRealConstant.value ( r )                 ## RETURN
    
    elif f1 : v1 = ROOT.RooRealConstant.value ( float ( v1 ) )        
    elif f2 : v2 = ROOT.RooRealConstant.value ( float ( v2 ) ) 
    #
    return Ostap.MoreRooFit.Asymmetry ( v1 , v2 , name , title ) 


# =============================================================================
## pow for two RooAbsReal objects
#  @code
#  v1 = ...
#  v2 = ...
#  v  = var_pow ( v1 , v2 )
#  @endcode
def var_pow ( v1 , v2 , name = '' , title = '' ) :
    """pow for two RooAbsReal objects
    >>> v1 = ...
    >>> v2 = ...
    >>> v  = var_pow ( v1 ,  v2 ) 
    """

    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )
    
    if f1 and f2 :
        r = float ( v1 ) ** float ( v2 ) 
        return ROOT.RooRealConstant.value ( r )                 ## RETURN
    
    elif f1 : v1 = ROOT.RooRealConstant.value ( float ( v1 ) )        
    elif f2 : v2 = ROOT.RooRealConstant.value ( float ( v2 ) ) 
    #
    return Ostap.MoreRooFit.Power ( v1 , v2 , name , title ) 


scale_var     = var_mul
add_var       = var_sum
sum_var       = var_sum
ratio_var     = var_div
fraction_var  = var_fraction
asymmetry_var = var_asymmetry


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
 
# =============================================================================
##                                                                      The END 
# =============================================================================
