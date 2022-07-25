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
    'var_sech'       , ## hyperbolic secant   function for RooAbsReal objects           
    'var_atan2'      , ## inverse tangent     function for RooAbsReal objects           
    'var_min'        , ## minimal             function for RooAbsReal objects           
    'var_max'        , ## minimal             function for RooAbsReal objects           
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
from   ostap.core.core                import Ostap 
from   ostap.core.ostap_types         import num_types
from   ostap.math.base                import iszero, isequal 
from   ostap.fitting.fithelpers       import ParamsPoly , ShiftScalePoly
import ostap.fitting.variables 
from   ostap.fitting.funbasic         import FUN1, Fun1D, Fun2D, Fun3D
# =============================================================================
from   ostap.logger.logger          import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.roofuncs' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
## is valuer equal to 1?
isone = lambda x : isequal ( float ( x ) , 1 )
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
class BernsteinPoly(FUN1,ParamsPoly) :
    """Polynomial in Bernstein form
    >>> p1 =  BernsteinPoly ( 'P1' , xvar = (0,1) , power = 3   ) 
    >>> p2 =  BernsteinPoly ( 'P2' , xvar = (0,1) , pars  = ... ) 
    - see Ostap.MoreRooFit.Bernstein
    - see Ostap.Math.Bernstein
    """
    def __init__ ( self , name , xvar , power = 1 , pars = None ) :
        
        ## initialize the base class 
        FUN1      .__init__ ( self  , name  , xvar = xvar )
        ParamsPoly.__init__ ( self          ,
                              power = power ,
                              pars  = pars  )
        
        xmin , xmax = self.xminmax ()
        
        ## create the function
        self.fun    = Ostap.MoreRooFit.Bernstein (
            self.roo_name ( 'bpol_' )   ,
            'Bernstein %s' % self.name  ,
            self.xvar                   ,
            self.pars_lst               ,            
            xmin                        ,
            xmax                        )

        self.tricks = True 
        self.config = {
            'name' : self.name ,
            'xvar' : self.xvar ,
            'pars' : self.pars
            }
    
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
class MonotonicPoly(FUN1,ShiftScalePoly) :
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
        FUN1          .__init__ ( self  , name  , xvar = xvar ) 
        ShiftScalePoly.__init__ ( self          ,
                                  a     = a     ,
                                  b     = b     ,
                                  power = power ,
                                  pars  = pars  )
        
        self.__increasing = True if increasing else False

        xmin , xmax = self.xminmax ()

        ## create the function
        self.fun    = Ostap.MoreRooFit.Monotonic (
            self.roo_name ( 'mpol_' )   ,
            '%s Bernstein %s' % ( "increasing" if increasing else "decreasing" , self.name ) , 
            self.xvar                   ,
            self.pars_lst               ,
            self.increasing             ,
            xmin                        ,
            xmax                        ,
            self.a                      ,
            self.b                      )            

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
        """'increasing' : increasing polynomial?"""
        return self.__increasing
    
    @property
    def decreasing ( self ) :
        """'decreasing' : decreasing polynomial?"""
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
class ConvexPoly(FUN1,ShiftScalePoly) :
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
        FUN1          .__init__ ( self , name   , xvar = xvar ) 
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
            self.roo_name ( 'cpol_' )   ,
            '%s %s Bernstein %s' % ( "convex"     if convex     else "concave"    ,
                                     "increasing" if increasing else "decreasing" ,
                                     self.name ) , 
            self.xvar                   ,
            self.pars_lst               ,
            self.increasing             ,
            self.convex                 ,
            xmin                        ,
            xmax                        ,
            self.a                      ,
            self.b                      )            

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
        """'increasing' : increasing polynomial?"""
        return self.__increasing
    
    @property
    def decreasing ( self ) :
        """'decreasing' : decreasing polynomial?"""
        return not self.__increasing 
    
    @property
    def convex ( self ) :
        """'convex' : convex polynomial?"""
        return self.__convex 

    @property
    def concave ( self ) :
        """'concave' : concave polynomial?"""
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
class ConvexOnlyPoly(FUN1,ShiftScalePoly) :
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
        FUN1          .__init__ ( self  , name  , xvar = xvar ) 
        ShiftScalePoly.__init__ ( self          ,
                                  a     = a     ,
                                  b     = b     ,
                                  power = power ,
                                  pars  = pars  )
        
        self.__convex     = True if convex     else False

        xmin , xmax = self.xminmax ()
        
        ## create the function
        self.fun    = Ostap.MoreRooFit.ConvexOnly (
            self.roo_name ( 'copol_' )   ,
            '%s Bernstein %s' % ( "convex" if convex else "concave" , self.name ) , 
            self.xvar                   ,
            self.pars_lst               ,
            self.convex                 ,
            xmin                        ,
            xmax                        ,
            self.a                      ,
            self.b                      )            

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
        """'convex' : convex polynomial?"""
        return self.__convex 

    @property
    def concave ( self ) :
        """'concave' : concave polynomial?"""
        return not self.__convex 

# =============================================================================
# @class ScaleAndShift
# Scale and shift another FUN/PDF 
# \f[ f = a + b c \f\
# @see Ostap::MoreRooFit::Product 
# @see Ostap::MoreRooFit::Addition
class ScaleAndShift (FUN1) :
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

        FUN1.__init__ ( self , name , xvar = xvar )

        assert isinstance ( func , ( ROOT.RooAbsReal , AFUN1 ) ) ,\
               "Invalid type of 'func' %s/%s" % ( func , type ( func ) )

        self.__func0 = func

        self.__c  = func.fun if isnstance ( func , AFUN1 ) else func 
        
        self.__a  = self.make_var ( a ,
                                    "a_%s"              % self.name ,
                                    "shift/bias for %s" % self.name , False ) 
        self.__b  = self.make_var ( b ,
                                    "b_%s"              % self.name ,
                                    "scale for %s"      % self.name , False ) 

        if not self.xvar in self.a.getParameters ( 0 ) and \
           not self.xvar in self.b.getParameters ( 0 ) and \
           not self.xvar in self.c.getParameters ( 0 ) :
            self.warning ( "Function does not depend on xvar=%s" % self.xvar.name )
           
        self.__bc = Ostap.MoreRooFit.Product  ( self.b , self.c )
        self.fun  = Ostap.MoreRooFit.Addition (
            self.roo_name ( 'ss_' )      ,
            "Scale&Shift %s" % self.name ,
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
        """'a' : bias parameter for function :  f(x) = a + b*c"""
        return self.__a
    @a.setter
    def a ( self , value ) :
        self.set_value ( self.__a , value )

    @property
    def b ( self ) :
        """'b' : scale parameter for function:  f(x) = a + b*c"""
        return self.__b
    @b.setter
    def b ( self , value ) :
        self.set_value ( self.__b , value )
    
    @property
    def c ( self ) :
        """'c' : the function:  f(x) = a + b * c"""
        return self.__c
    @c.setter
    def c ( self , value ) :
        self.set_value ( self.__c , value )

# ==============================================================================
# primitive functions for RooAbsReal objects 
# ==============================================================================

# ==============================================================================
## absolute value   \f$ f = abs{ab} \f$
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
        return ROOT.RooFit.RooConst ( ab )     
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 0 ) 
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
        return ROOT.RooFit.RooConst ( ab )          ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 1 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 1 ) 
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
        return ROOT.RooFit.RooConst ( ab )
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
        return ROOT.RooFit.RooConst ( ab )
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
        return ROOT.RooFit.RooConst ( ab )          ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 0 ) 
    #
    return Ostap.MoreRooFit.Erf ( a, b , name , title ) 

# ==============================================================================
## complementary error function \f$ f = erfc ( ab) \f$
#  @code
#  var = ...
#  e   = var_erfc ( var ) 
#  @endcode 
def var_erfc ( a , b = 1 , name = '' , title = '' ) :
    """Error function f(x) = erfc(ab)
    >>> var = ...
    >>> e   = var_erfc ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.erfc ( float ( a ) * float ( b ) )
        return ROOT.RooFit.RooConst ( ab )          ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 1 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 1 ) 
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
        return ROOT.RooFit.RooConst ( ab )          ## RETURN 
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 0 ) 
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
        return ROOT.RooFit.RooConst ( ab )       ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 1 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 1 ) 
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
        return ROOT.RooFit.RooConst ( ab )
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 0 ) 
    #
    return Ostap.MoreRooFit.Tan ( a, b , name , title )


# ==============================================================================
## Hyprbolic sine \f$ f = \sinh ab\f$
#  @code
#  var = ...
#  e   = var_sinh ( var ) 
#  @endcode 
def var_sinh ( a , b = 1 , name = '' , title = '' ) :
    """Hyperbolic sine  f(x) = sinh(ab)
    >>> var = ...
    >>> e   = var_sinh ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.sinh ( float ( a ) * float ( b ) )
        return ROOT.RooFit.RooConst ( ab )          ## RETURN 
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 0 ) 
    #
    return Ostap.MoreRooFit.Sinh ( a, b , name , title ) 
    

# ==============================================================================
## Hyperbolic cosine \f$ f = \cosh ab\f$
#  @code
#  var = ...
#  e   = var_cosh ( var ) 
#  @endcode 
def var_cosh ( a , b = 1 , name = '' , title = '' ) :
    """Hyperbolic cosine  f(x) = cos(ab)
    >>> var = ...
    >>> e   = var_cosh ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = math.cosh ( float ( a ) * float ( b ) )
        return ROOT.RooFit.RooConst ( ab )       ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 1 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 1 ) 
    #
    return Ostap.MoreRooFit.Cosh ( a, b , name , title )

# ==============================================================================
## Hyperboilic tangent\f$ f = \tanh ab\f$
#  @code
#  var = ...
#  e   = var_tanh ( var ) 
#  @endcode 
def var_tanh ( a , b = 1 , name = '' , title = '' ) :
    """Hyperbolic tangent  f(x) = tanh(ab)
    >>> var = ...
    >>> e   = var_tanh ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :        
        ab = math.tanh ( float ( a ) * float ( b ) ) 
        return ROOT.RooFit.RooConst ( ab )      ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 0 ) 
    #
    return Ostap.MoreRooFit.Tanh ( a, b , name , title )

# ==============================================================================
## Hyperboilic secant \f$ f = \sech ab\f$
#  @code
#  var = ...
#  e   = var_sech ( var ) 
#  @endcode 
def var_sech ( a , b = 1 , name = '' , title = '' ) :
    """Hyperbolic tangent  f(x) = tanh(ab)
    >>> var = ...
    >>> e   = var_tanh ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :        
        ab = Ostap.Math.sech ( float ( a ) * float ( b ) ) 
        return ROOT.RooFit.RooConst ( ab )      ## RETURN
    elif fa and iszero ( fa ) : return ROOT.RooFit.RooConst ( 1 ) 
    elif fb and iszero ( fb ) : return ROOT.RooFit.RooConst ( 1 ) 
    #
    return Ostap.MoreRooFit.Sech ( a, b , name , title )

# ==============================================================================
## arctangent \f$ f = atan2 (a,b)\f$
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
        ab = math.atan2 ( float ( a ) , float ( b ) ) 
        return ROOT.RooFit.RooConst ( ab )      ## RETURN
    return Ostap.MoreRooFit.Atan2 ( a, b , name , title )


# ==============================================================================
## maximal \f$ f = max (a,b)\f$
#  @code
#  var1 = ...
#  var2 = ...
#  var  = var_max ( var1 , var2 ) 
#  @endcode 
def var_max ( a , b = 1 , name = '' , title = '' ) :
    """Maximal from two fnuctions f(x) = max(a,b)
    >>> var = ...
    >>> e   = var_max ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = max ( float ( a ) , float ( b ) ) 
        return ROOT.RooFit.RooConst ( ab )      ## RETURN
    return Ostap.MoreRooFit.MaxV ( a, b , name , title )


# ==============================================================================
## minimal \f$ f = min (a,b)\f$
#  @code
#  var1 = ...
#  var2 = ...
#  var  = var_min ( var1 , var2 ) 
#  @endcode 
def var_min ( a , b = 1 , name = '' , title = '' ) :
    """Minimal from two fnuctions f(x) = min(a,b)
    >>> var = ...
    >>> e   = var_min ( var ) 
    """
    fa = isinstance ( a , num_types )
    fb = isinstance ( b , num_types )
    if fa and fb :
        ab = min ( float ( a ) , float ( b ) ) 
        return ROOT.RooFit.RooConst ( ab )      ## RETURN
    return Ostap.MoreRooFit.MinV ( a, b , name , title )


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
        return ROOT.RooFit.RooConst ( ab )           ## RETURN 
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
        return ROOT.RooFit.RooConst ( ab )           ## RETURN 
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
        return ROOT.RooFit.RooConst ( ab )              ## RETURN 
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
        return ROOT.RooFit.RooConst ( r )                 ## RETURN    
    elif f1 and iszero ( v1 ) : return v2
    elif f2 and iszero ( v2 ) : return v1
    #
    return Ostap.MoreRooFit.Addition ( v1 , v2 , name , title )

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
        return ROOT.RooFit.RooConst ( r )                 ## RETURN    
    elif f1 and iszero ( v1 ) : return var_mul ( v2 , -1 , name , title ) 
    elif f2 and iszero ( v2 ) : return v1 
    # 
    return Ostap.MoreRooFit.Subtraction ( v1 , v2 , name , title ) 


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
        return ROOT.RooFit.RooConst ( r )                 ## RETURN
    elif f1 and iszero ( v1 ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif f2 and iszero ( v2 ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif f1 and isone  ( v1 ) : return v2 
    elif f2 and iseone ( v2 ) : return v1 
    #
    return Ostap.MoreRooFit.Product ( v1 , v2 , name , title ) 

# =============================================================================
## Division of two RooAbsReal objects
#  @code
#  v1 = ...
#  v2 = ...
#  v  = var_div ( v1 , v2 )
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
        return ROOT.RooFit.RooConst ( r )                 ## RETURN
    elif f1 and iszero ( v1 ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif f2 and isone  ( v2 ) : return v1
    #
    return Ostap.MoreRooFit.Division ( v1 , v2 , name , title ) 


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
        return ROOT.RooFit.RooConst ( r )                 ## RETURN    
    elif f2 and iszero  ( v2 ) : return ROOT.RooFit.RooConst ( 1 )
    elif f2 and isone   ( v2 ) : return v1 
    elif f1 and iszero  ( v1 ) : return ROOT.RooFit.RooConst ( 0 )
    elif f1 and isone   ( v1 ) : return ROOT.RooFit.RooConst ( 1 )
    #
    return Ostap.MoreRooFit.Power ( v1 , v2 , name , title ) 

# ==============================================================================
## "Fraction" of two RooAbsReal objects: f = a/(a+b)
#  @code
#  a = ...
#  b = ...
#  e   = var_fraction ( a , b ) 
#  @endcode 
def var_fraction ( a , b , name = '' , title = '' ) :
    """'Fraction'  f(x) = a/(a+b)
    >>> a = ...
    >>> b = ...
    >>> e = var_fraction ( a , b  ) 
    """
    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )    
    if f1 and f2 :
        r = float ( v1) / ( float ( v1 ) + float ( v1 ) )  
        return ROOT.RooFit.RooConst ( r )                 ## RETURN    
    elif f1 and iszero  ( v1 ) : return ROOT.RooFit.RooConst ( 0 ) 
    elif f2 and iszero  ( v2 ) : return ROOT.RooFit.RooConst ( 1 ) 
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
    """'Asymmetry'  f(x) = (a-b)/(a+b)
    >>> a = ...
    >>> b = ...
    >>> e = var_asymmetry ( a , b  ) 
    """
    f1 = isinstance ( v1 , num_types )
    f2 = isinstance ( v2 , num_types )    
    if f1 and f2 :
        r = ( float ( v1 ) - float ( v2 ) ) / ( float ( v1 ) + float ( v2 ) )  ## 
        return ROOT.RooFit.RooConst ( r )                                ## RETURN    
    elif f1 and iszero  ( v1 ) : return ROOT.RooFit.RooConst ( -1 ) 
    elif f2 and iszero  ( v2 ) : return ROOT.RooFit.RooConst (  1 ) 
    #
    return Ostap.MoreRooFit.Asymmetry ( v1 , v2 , name , title ) 


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
