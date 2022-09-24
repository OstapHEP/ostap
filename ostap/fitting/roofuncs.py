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
    'RooPoly'        , ## simple wrapper for RooPolyVar            (RooAbsReal)
    'ScaleAndShift'  , ## scale and shift                          (RooAbsReal)
    'BSplineFun'     , ## BSpline                                  (RooAbsReal)
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
from   ostap.core.core                import Ostap 
from   ostap.core.ostap_types         import num_types
from   ostap.fitting.fithelpers       import ParamsPoly , ShiftScalePoly
import ostap.fitting.variables 
from   ostap.fitting.funbasic         import FUN1, Fun1D, Fun2D, Fun3D
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
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
 
# =============================================================================
##                                                                      The END 
# =============================================================================
