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
    'BernsteinPoly'        , ## generic polynomial in Bernstein form     (RooAbsReal)
    'MonotonicPoly'        , ## monotonic polynomial                     (RooAbsReal)
    'ConvexPoly'           , ## monotonic convex/concave polynomial      (RooAbsReal)
    'ConvexOnlyPoly'       , ## convex/concave polynomial                (RooAbsReal)
    'RooPoly'              , ## simple wrapper for RooPolyVar            (RooAbsReal)
    'ScaleAndShift'        , ## scale and shift                          (RooAbsReal)
    'BSplineFun'           , ## BSpline                                  (RooAbsReal)
    'RationalFun'          , ## Ratioal function                         (RooAbsReal)
    'RationalBernsteinFun' , ## Ratioal function                         (RooAbsReal)
    'Shape1D_fun'          , ## arbitrary fixed shape                    (RooAbsReal)
    'Histo1D_fun'          , ## fixed shap form historgam                (RooAbsReal)
    'Histo1DErr_fun'       , ## fixed shap from historgam errors         (RooAbsReal)
    ##
    ## 'var_sum'        , ## sum                          for RooAbsReal objects           
    ## 'var_mul'        , ## product                      for RooAbsReal objects           
    ## 'var_sub'        , ## subtraction                  for RooAbsReal objects           
    ## 'var_div'        , ## division                     for RooAbsReal objects           
    ## 'var_fraction'   , ## fraction                     for RooAbsReal objects           
    ## 'var_asymmetry'  , ## asymmetry                    for RooAbsReal objects
    ## 'var_pow'        , ## pow                 function for RooAbsReal objects           
    ## 'var_abs'        , ## absolutevalue       function for RooAbsReal objects           
    ## 'var_exp'        , ## exponent            function for RooAbsReal objects           
    ## 'var_log'        , ## logarithm           function for RooAbsReal objects           
    ## 'var_log10'      , ## logarithm           function for RooAbsReal objects           
    ## 'var_erf'        , ## error               function for RooAbsReal objects           
    ## 'var_sin'        , ## sine                function for RooAbsReal objects           
    ## 'var_cos'        , ## cosine              function for RooAbsReal objects           
    ## 'var_tan'        , ## tangent             function for RooAbsReal objects           
    ## 'var_tanh'       , ## hyperbolic tangent  function for RooAbsReal objects           
    ## 'var_sech'       , ## hyperbolic secant   function for RooAbsReal objects           
    ## 'var_atan2'      , ## inverse tangent     function for RooAbsReal objects           
    ## 'var_min'        , ## minimal             function for RooAbsReal objects           
    ## 'var_max'        , ## minimal             function for RooAbsReal objects           
    ## 'var_gamma'      , ## gamma               function for RooAbsReal objects           
    ## 'var_lgamma'     , ## logarithm of gamma  function for RooAbsReal objects           
    ## 'var_igamma'     , ## 1/gamma             function for RooAbsReal objects
    ## ##
    ## 'scale_var'      , ## var_mul
    ## 'add_var'        , ## var_sum
    ## 'sum_var'        , ## var_sum
    ## 'ratio_var'      , ## var_div
    ## 'fraction_var'   , ## var_fraction
    ## 'asymmetry_var'  , ## var_asymmetry
   )
# =============================================================================
from   ostap.core.core                import Ostap, VE  
from   ostap.core.ostap_types         import num_types, integer_types
from   ostap.fitting.fithelpers       import ParamsPoly , ShiftScalePoly
import ostap.fitting.variables 
from   ostap.fitting.funbasic         import FUN1, Fun1D, Fun2D, Fun3D
import ostap.histos.histos 
import ROOT, math, array
# =============================================================================
from   ostap.logger.logger          import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.roofuncs' )
else                       : logger = getLogger ( __name__                 )
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
        ParamsPoly.__init__ ( self              ,
                              npars = power + 1 ,
                              pars  = pars      )

        assert 1 <= self.npars , 'Invalid number of parameters! '
        
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
            'name'  : self.name      ,
            'xvar'  : self.xvar      ,
            'pars'  : self.pars      ,
            'power' : self.npars - 1 
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


# ============================================================================
## @class RooPoly
#  Simple wrapper for class RooPolyVar
#  @see RoPolyVar
#  @code
#  p1 =  RooPoly ( 'P1' , xvar = (0,1) , power = 3 ) 
#  p2 =  RooPoly ( 'P2' , xvar = (0,1) , pars  = ... ) 
#  @endcode 
#  @see RooPolyVar 
class RooPoly(FUN1,ParamsPoly) :
    """ Simple wrapper for class RooPolyVar
    >>> p1 =  RooPoly ( 'P1' , xvar = (0,1) , power = 3 ) 
    >>> p2 =  RooPoly ( 'P2' , xvar = (0,1) , pars  = ... ) 
    - see `ROOT.RooPolyVar`
    """
    def __init__ ( self , name , xvar , power = 1 , pars = None ) :
        
        ## initialize the base class 
        FUN1      .__init__ ( self  , name  , xvar = xvar )
        ParamsPoly.__init__ ( self              ,
                              npars = power + 1 ,
                              pars  = pars      )

        assert 1 <= self.npars , 'Invalid number of parameters! '
        
        xmin , xmax = self.xminmax ()
        
        ## create the function
        self.fun    = ROOT.RooPolyVar (
            self.roo_name ( 'rpol_' )   ,
            'RooPoly %s' % self.name    ,
            self.xvar                   ,
            self.pars_lst               ) 
        
        self.tricks = True 
        self.config = {
            'name'  : self.name      ,
            'xvar'  : self.xvar      ,
            'pars'  : self.pars      ,
            'power' : self.npars - 1  
            }
        
# =============================================================================
## @class BSplineFun 
#  BSpline as RooAbsReal
#  @code
#  knots = ...
#  p1 =  Bspline ( 'S1' , xvar = (0,1) , knots = knots , power = 3  )
#  p2 =  Bspline ( 'S1' , xvar = (0,1) , knots = knots , pars = ... )
#  @endcode 
#  @see Ostap::MoreRooFit::BSpline
#  @see Ostap::Math::BSpline
class BSplineFun(FUN1,ParamsPoly) :
    """ BSpline as RooAbsReal
    >>> knots = ...
    >>> p1 =  Bspline ( 'S1' , xvar = (0,1) , knots = knots , power = 3  )
    - see Ostap.MoreRooFit.BSpline
    - see Ostap.Math.BSpline
    """
    def __init__ ( self , name , xvar , knots , power  , pars = None ) :
        
        from ostap.math.base import isint as _isint 
        assert pars or ( ( isinstance ( power, integer_types ) or _isint ( power ) ) and 0 <= power ) ,\
               "BSplineFun: Inconsistent 'power' setting!"
 
        ## initialize the first base class 
        FUN1      .__init__ ( self  , name  , xvar = xvar )

        from ostap.math.base import doubles 
        knots = doubles ( knots ) 
        assert 2 <= len ( knots ) , "BSplineFun: invalid size of 'knots'!"

        ## number of inner knots 
        inner = len ( knots ) - 2 

        ## number of parameters 
        npars = power + inner + 1 
        
        ## initialize the second base class        
        ParamsPoly.__init__ ( self          ,
                              npars = npars ,
                              pars  = pars  )
        
        ## create the function
        self.fun    = Ostap.MoreRooFit.BSpline (
            self.roo_name ( 'bspline_' )   ,
            'BSpline %s' % self.name  ,
            self.xvar                 ,
            knots                     ,            
            self.pars_lst             )


        self.__knots = array.array ( 'd' , self.fun.knots () )
        
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'knots' : self.knots ,
            'power' : self.power , 
            'pars'  : self.pars  ,
            }

    @property 
    def knots ( self ) :
        """'knots' : array of knots for BSpline object""" 
        return self.__knots
    
    @property 
    def power ( self ) :
        """'power' : degree for BSpline object""" 
        return self.fun.degree() 


# =============================================================================
## @class RationalFun 
#  A simple pole-free rational function at interval \f$ x_{min} \le x \le x_{max}\f$
#  \f[ F(x) = \frac{p(x)}{q(x)} \f]
#  Actually internally it uses 
#  the Floater-Hormann rational barycentric interpolant 
#  and parameters are the function values at Chebyshev's nodes  
#   
#  @see Ostap::MoreRooFit::Rational
#  @see Ostap::Math::Rational
#  @see Ostap::Math::FloaterHormann
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2023-09-21
class RationalFun(FUN1,ParamsPoly) :
    r"""A simple pole-free rational function at interval \f$ x_{min} \le x \le x_{max}\f$
    >>> p1 =  BernsteinPoly ( 'P1' , xvar = (0,1) , power = 3   ) 
    >>> p2 =  BernsteinPoly ( 'P2' , xvar = (0,1) , pars  = ... ) 
    - see Ostap.MoreRooFit.Bernstein
    - see Ostap.Math.Bernstein
    """
    def __init__ ( self        ,
                   name        ,
                   xvar        ,
                   n    = 3    ,
                   d    = 1    ,
                   pars = None ) :
        
        ## initialize the base class 
        FUN1      .__init__ ( self  , name  , xvar = xvar )
        ParamsPoly.__init__ ( self          ,
                              npars = n + 1 ,
                              pars  = pars  )
        
        assert isinstance ( n , integer_types ) and 0 <= n      , 'Invalid parameter n'
        assert isinstance ( d , integer_types ) and 0 <= d <= n , 'Invalid parameter d'
        
        xmin , xmax = self.xminmax ()
        
        ## create the function
        self.fun    = Ostap.MoreRooFit.Rational (
            self.roo_name ( 'rfun_' )   ,
            'Rational %s' % self.name   ,
            self.xvar                   ,
            self.pars_lst               ,
            d                           , 
            xmin                        ,
            xmax                        )
        
        ## self.tricks = True 
        self.config = {
            'name'  : self.name      ,
            'xvar'  : self.xvar      ,
            'n'     : self.fun.n ()  ,
            'd'     : self.fun.d ()  ,
            'pars'  : self.pars      ,
            }    

    @property 
    def n ( self ) :
        """'n' : n-parameter of rational function""" 
        return self.fun.n ()
    @property 
    def d ( self ) :
        """'d' : d-parameter of rational function""" 
        return self.fun.d ()



# =============================================================================
## @class RationalBernsteinFun 
#  A simple rational function at interval \f$ x_{min} \le x \le x_{max}\f$ as ratio
#  of Bernstein and positive Bernstein polynomials 
#  \f[ F(x) = \frac{p(x)}{q(x)} \f]
#  @see Ostap::MoreRooFit::RationalBernstein
#  @see Ostap::Math::RationalBernstein
#  @see Ostap::Math::Bernstein
#  @see Ostap::Math::Positive
#  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
#  @date   2023-09-21
class RationalBernsteinFun(FUN1,ParamsPoly) :
    r""" A simple rational function at interval \f$ x_{min} \le x \le x_{max}\f$ as ratio
    of Bernstein and positive Bernstein polynomials 
    - see Ostap.MoreRooFit.RationalBernstein
    - see Ostap.Math.RationalBernstein
    - see Ostap.Math.Bernstein
    - see Ostap.Math.Positive
    """
    def __init__ ( self        ,
                   name        ,
                   xvar        ,
                   p    = 3    ,
                   q    = 3    ,
                   pars = None ) :
        
        assert isinstance ( p , integer_types ) and 0 <= p , 'Invalid parameter p'
        assert isinstance ( q , integer_types ) and 0 <= q , 'Invalid parameter q'
        
        ## initialize the base class 
        FUN1      .__init__ ( self  , name  , xvar = xvar )
        ParamsPoly.__init__ ( self          ,
                              npars = p + q + 1 ,
                              pars  = pars  )

        assert self.pars , 'Invalid number of parameters!'
        p = min ( p , self.npars  - 1 )
        
        if not pars :
            for i, v in enumerate ( self.pars ) :
                if i < p + 1 : v.setVal ( 1 ) 
                else : 
                    v.setMin ( -5 * math.pi ) 
                    v.setMax ( +5 * math.pi )
                    
        xmin , xmax = self.xminmax ()

        ## create the function
        self.fun    = Ostap.MoreRooFit.RationalBernstein (
            self.roo_name ( 'rbf_' )           ,
            'RationalBernstein %s' % self.name ,
            self.xvar                   ,
            self.pars_lst               ,
            p                           , 
            xmin                        ,
            xmax                        )
        
        ## self.tricks = True 
        self.config = {
            'name'  : self.name      ,
            'xvar'  : self.xvar      ,
            'p'     : self.fun.p ()  ,
            'q'     : self.fun.q ()  ,
            'pars'  : self.pars      ,
            }    

    @property 
    def p ( self ) :
        """'p' : p-parameter of rational function""" 
        return self.fun.p ()
    @property 
    def q ( self ) :
        """'q' : q-parameter of rational function""" 
        return self.fun.q ()

    
# =============================================================================
## Generic 1D-shape from C++ callable
#  @see Ostap::MoreRooFit::Shape1D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Shape1D_fun(FUN1) :
    """ Generic 1D-shape from C++ callable
    - see Ostap::MoreRooFit:Shape1D
    """
    
    def __init__ ( self , name , shape , xvar ) :

        if isinstance ( shape , ROOT.TH1 ) and not isinstance ( shape , ROOT.TH2 ) and not xvar :
            xvar = shape.xminmax() 

        if isinstance ( shape , ROOT.TH1 ) and not isinstance ( shape , ROOT.TH2 ) :
            self.histo = shape
            shape      = Ostap.Math.Histo1D ( shape )

        ##  initialize the base 
        FUN1.__init__ ( self , name , xvar = xvar ) 

        
        self.__shape = shape

        if isinstance ( self.shape , Ostap.Math.Histo1D ) :
        
            ## create the actual function
            self.fun = Ostap.MoreRooFit.Histo1D (
                self.roo_name ( 'histo1fun_' ) , 
                "Histo-1D-fun %s" % self.name  ,
                self.xvar                      ,
                self.shape                     )
            
        else :
            
            ## create the actual pdf
            self.fun = Ostap.MoreRooFit.Shape1D.create  (
                self.roo_name ( 'shape1fun_' ) , 
                "Shape-1D-fun %s" % self.name  ,
                self.xvar                      ,
                self.shape                     ) 
            
        ## save the configuration
        self.config = {
            'name'    : self.name    , 
            'shape'   : self.shape   , 
            'xvar'    : self.xvar    , 
            }
        
    @property
    def shape  ( self ) :
        """'shape': the actual C++ callable shape"""
        return self.__shape 

# =============================================================================
## Generic 1D-shape from histogram 
#  @see Ostap::MoreRooFit::Histo1D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Histo1D_fun(Shape1D_fun) :
    """ Generic 1D-shape from C++ callable
    - see Ostap::MoreRooFit:Shape1D
    """
    
    def __init__ ( self , name , histo , xvar ) :

        assert ( isinstance ( histo , ROOT.TH1 ) and 1 == histo.dim() ) or \
               isinstance ( histo , Ostap.Math.Histo1D ) , \
               "Invalid 'histo' argument!"

        ##  initialize the base class 
        Shape1D_fun.__init__ ( self , name , xvar = xvar , shape = histo )

        ## save the configuration
        self.config = {
            'name'    : self.name  , 
            'histo'   : self.shape , 
            'xvar'    : self.xvar  , 
            }
        
    @property
    def histo  ( self ) :
        """'histo': the actual histogram (same as 'shape')"""
        return self.shape 

# =============================================================================
## Generic 1D-shape from histogram errors 
#  @see Ostap::MoreRooFit::Histo1D
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2020-07-20
class Histo1DErr_fun(Shape1D_fun) :
    """ Generic 1D-shape from hisgram errors 
    - see Ostap::MoreRooFit:Shape1D
    """
    
    def __init__ ( self , name , histo , xvar ) :

        assert ( isinstance ( histo , ROOT.TH1 ) and 1 == histo.dim() ) or \
               isinstance ( histo , Ostap.Math.Histo1D ) , \
               "Invalid 'histo' argument!"

        ## get the histogram
        h    = histo.h    () if isinstance ( histo , Ostap.Math.Histo1D ) else histo

        self.__init_histo = histo 

        ## clone it and store uncertainties
        herr = h.clone ()
        for i, x, y in herr.items() : herr [ i ] = VE ( y.error() , 0 ) 

        if isinstance ( histo , Ostap.Math.Histo1D ) :
            he = Ostap<Math.Histo1D ( he , histo ) 
            
        ##  initialize the base class 
        Histo1D_fun.__init__ ( self , name , xvar = xvar , histo = herr )

        ## save the configuration
        self.config = {
            'name'    : self.name       , 
            'histo'   : self.init_histo , 
            'xvar'    : self.xvar       , 
            }
        
    @property
    def init_histo  ( self ) :
        """'init_histo': initial histogram"""
        return self.__init_histo

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
 
# =============================================================================
##                                                                      The END 
# =============================================================================
