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
    'MonotonicPoly'  , ## monitonic polynomial                     (RooAbsReal)
    'ConvexPoly'     , ## monotonic convex/concave polynomial      (RooAbsReal)
    'ConvexOnlyPoly' , ## convex/concave polynomial                (RooAbsReal)
    'ScaleAndShift'  , ## scale and shift                          (RooAbsReal)
    'Exp'            , ## exponent                                 (RooAbsReal)
    'Log'            , ## logarithm                                (RooAbsReal)
    'Sin'            , ## Sine                                     (RooAbsReal)
    'Cos'            , ## Cosine                                   (RooAbsReal)
    'Tan'            , ## Tangent                                  (RooAbsReal)
    'Tanh'           , ## hyperbolic tangent                       (RooAbsReal)
    #
    'var_sum'        , ## sum                         for RooAbsReal objects           
    'var_mul'        , ## product                     for RooAbsReal objects           
    'var_sub'        , ## subtraction                 for RooAbsReal objects           
    'var_div'        , ## division                    for RooAbsReal objects           
    'var_fraction'   , ## fraction                    for RooAbsReal objects           
    'var_asymmetry'  , ## asymmetry                   for RooAbsReal objects           
    'var_pow'        , ## pow                 function for RooAbsReal objects           
    'var_exp'        , ## exponent            function for RooAbsReal objects           
    'var_log'        , ## logarithm           function for RooAbsReal objects           
    'var_erf'        , ## error               function for RooAbsReal objects           
    'var_sin'        , ## sine                function for RooAbsReal objects           
    'var_cos'        , ## cosine              function for RooAbsReal objects           
    'var_tan'        , ## tangent             function for RooAbsReal objects           
    'var_tanh'       , ## hyperbolic tangent  function for RooAbsReal objects           
    'var_atan2'      , ## inverse tangent     function for RooAbsReal objects           
    'var_gamma'      , ## gamma               function for RooAbsReal objects           
    'var_lgamma'     , ## logaruithm of gamma function for RooAbsReal objects           
    'var_igamma'     , ## 1/gamma             function for RooAbsReal objects           
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
from ostap.fitting.funbasic         import FUNC
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
                   xvar                  ,
                   b       = 1.0         ,  ## b 
                   name    = ''          ,
                   pattern = 'FAB_%s_%s' ) :
                   
        assert isinstance ( a , ( ROOT.RooAbsReal, FUNC ) ) ,\
               "Invalid type of ``func'' %s/%s" % ( a , type ( a ) )

        if not name :
            if hasattr ( b , 'name' ) : name = pattern % ( a.name , b.name  )
            else                      : name = pattern % ( a.name , b       )
                    
        FUNC.__init__ ( self , name , xvar = xvar )

        self.__func0 = a 

        self.__a  = a.fun if isinstance ( a , FUNC ) else a 
        self.__b  = self.make_var ( b ,
                                    "b_%s"              % self.name ,
                                    "scale for %s"      % self.name , False , b ) 
        
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
    @a.setter
    def a ( self , value ) :
        vv = float ( value )
        if self.__a.minmax () and not vv in self.__a  :
            self.error ("Value %s is outside the allowed region %s"  % ( vv , self.__a.minmax() ) )
        self.__a.setVal ( vv )

    @property
    def b ( self ) :
        """``b''  : f(x) = F(ab)"""
        return self.__b
    @b.setter
    def b ( self , value ) :
        vv = float ( value )
        if self.__b.minmax () and not vv in self.__b  :
            self.error ("Value %s is outside the allowed region %s"  % ( vv , self.__b.minmax() ) )
        self.__b.setVal ( vv )


# =============================================================================
# @class Exp
# \f[ f = \mathrm{e}^{ab} \f] 
# @see Ostap::MoreRooFit::Exp
class Exp (FAB) :
    """ Exponent   
    f = exp( ab ) 
    - see Ostap.MoreRooFit.Exp
    """
    def __init__ ( self       ,
                   a          ,   ## a 
                   xvar       ,
                   b    = 1.0 ,  ## b 
                   name = ''  ) :

        FAB.__init__ ( self , a = a , xvar = xvar , b = b , name = name , pattern = "Exp_%s_%s" )
            
        self.fun  = Ostap.MoreRooFit.Exp (
            "exp_%s"   % self.name ,
            "Exp(%s|%s,%s)" % ( self.name , self.a.name , self.b.name ) , 
            self.a      ,
            self.b  
            )
        
# =============================================================================
# @class Log
# \f[ f = \log ab \f] 
# @see Ostap::MoreRooFit::Log
class Log (FAB) :
    """Logarithm
    f = log(ab) 
    - see Ostap.MoreRooFit.Log
    """ 
    def __init__ ( self       ,
                   a          ,   ## a 
                   xvar       ,
                   b    = 1.0 ,  ## b 
                   name = ''  ) :

        FAB.__init__ ( self , a = a , xvar = xvar , b = b , name = name , pattern = "Log_%s_%s" )
            
        self.fun  = Ostap.MoreRooFit.Log (
            "log_%s"        % self.name ,
            "Log(%s|%s,%s)" % ( self.name , self.a.name , self.b.name ) , 
            self.a      ,
            self.b  
            )

# =============================================================================
# @class Sin
# \f[ f = \sin ab \f] 
# @see Ostap::MoreRooFit::Sin
class Sin (FAB) :
    """Sine
    f = sin(ab) 
    - see Ostap.MoreRooFit.Sin
    """
    def __init__ ( self       ,
                   a          ,   ## a 
                   xvar       ,
                   b    = 1.0 ,  ## b 
                   name = ''  ) :

        FAB.__init__ ( self , a = a , xvar = xvar , b = b , name = name , pattern = "Sin_%s_%s" )
            
        self.fun  = Ostap.MoreRooFit.Sin (
            "sin_%s"        % self.name ,
            "Sin(%s|%s,%s)" % ( self.name , self.a.name , self.b.name ) , 
            self.a      ,
            self.b  
            )

# =============================================================================
# @class Cos
# \f[ f = \cos ab \f] 
# @see Ostap::MoreRooFit::Cos
class Cos (FAB) :
    """Cosine
    f = cos(ab) 
    - see Ostap.MoreRooFit.Cos
    """
    def __init__ ( self       ,
                   a          ,   ## a 
                   xvar       ,
                   b    = 1.0 ,  ## b 
                   name = ''  ) :

        FAB.__init__ ( self ,  a = a , xvar = xvar , b = b , name = name , pattern = "Cos_%s_%s" )
            
        self.fun  = Ostap.MoreRooFit.Cos (
            "cos_%s"        % self.name ,
            "Cos(%s|%s,%s)" % ( self.name , self.a.name , self.b.name ) , 
            self.a      ,
            self.b  
            )

# =============================================================================
# @class Tan
# \f[ f = \tan ab \f] 
# @see Ostap::MoreRooFit::Tan
class Tan (FAB) :
    """Tangent
    f = tan(ab) 
    - see Ostap.MoreRooFit.Tan
    """
    def __init__ ( self       ,
                   a          ,   ## a 
                   xvar       ,
                   b    = 1.0 ,  ## b 
                   name = ''  ) :

        FAB.__init__ ( self , a = a , xvar = xvar , b = b , name = name , pattern = "Tan_%s_%s" )
            
        self.fun  = Ostap.MoreRooFit.Tan (
            "tan_%s"        % self.name ,
            "Tan(%s|%s,%s)" % ( self.name , self.a.name , self.b.name ) , 
            self.a      ,
            self.b  
            )

# =============================================================================
# @class Tanh
# \f[ f = \tanh ab \f] 
# @see Ostap::MoreRooFit::Tanh
class Tanh ( FUNC ) :
    """hyperbolic Tangent
    f = tanh(ab) 
    - see Ostap.MoreRooFit.Tanh
    """
    def __init__ ( self       ,
                   a          ,   ## a 
                   xvar       ,
                   b    = 1.0 ,  ## b 
                   name = ''  ) :

        FAB.__init__ ( self , a = a , xvar = xvar , b = b , name = name , pattern = "Tanh_%s_%s" )
            
        self.fun  = Ostap.MoreRooFit.Tanh (
            "tanh_%s"        % self.name ,
            "Tanh(%s|%s,%s)" % ( self.name , self.a.name , self.b.name ) , 
            self.a      ,
            self.b  
            )

# =============================================================================
## local storage of temporary variables 
KEEPER = MakeVar()
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
    result    = Ostap.MoreRooFit.Exp ( a, b , name , title )
    #
    KEEPER.aux_keep.add ( a )
    KEEPER.aux_keep.add ( b )
    KEEPEP.aux_keep.add ( result )
    #
    return result

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
    result    = Ostap.MoreRooFit.Log ( a, b , name , title ) 
    #
    KEEPER.aux_keep.add ( a )
    KEEPER.aux_keep.add ( b )
    KEEPEP.aux_keep.add ( result )
    #
    return result

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
    result    = Ostap.MoreRooFit.Erf ( a, b , name , title ) 
    #
    KEEPER.aux_keep.add ( a )
    KEEPER.aux_keep.add ( b )
    KEEPEP.aux_keep.add ( result )
    #
    return result


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
    result    = Ostap.MoreRooFit.Sin ( a, b , name , title ) 
    #
    KEEPER.aux_keep.add ( a )
    KEEPER.aux_keep.add ( b )
    KEEPEP.aux_keep.add ( result )
    #
    return result


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
    result    = Ostap.MoreRooFit.Cos ( a, b , name , title )
    #
    KEEPER.aux_keep.add ( a )
    KEEPER.aux_keep.add ( b )
    KEEPEP.aux_keep.add ( result )
    #
    return result

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
    result    = Ostap.MoreRooFit.Tan ( a, b , name , title )
    #
    KEEPER.aux_keep.add ( a )
    KEEPER.aux_keep.add ( b )
    KEEPEP.aux_keep.add ( result )
    #
    return result


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
    result    = Ostap.MoreRooFit.Tanh ( a, b , name , title )
    #
    KEEPER.aux_keep.add ( a )
    KEEPER.aux_keep.add ( b )
    KEEPEP.aux_keep.add ( result )
    #
    return result

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
    result    = Ostap.MoreRooFit.Atan2 ( a, b , name , title )
    #
    KEEPER.aux_keep.add ( a )
    KEEPER.aux_keep.add ( b )
    KEEPEP.aux_keep.add ( result )
    #
    return result

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
    result    = Ostap.MoreRooFit.Gamma ( a, b , name , title ) 
    #
    KEEPER.aux_keep.add ( a )
    KEEPER.aux_keep.add ( b )
    KEEPEP.aux_keep.add ( result )
    #
    return result

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
    result    = Ostap.MoreRooFit.LGamma ( a, b , name , title ) 
    #
    KEEPER.aux_keep.add ( a )
    KEEPER.aux_keep.add ( b )
    KEEPEP.aux_keep.add ( result )
    #
    return result

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
    result    = Ostap.MoreRooFit.IGamma ( a, b , name , title ) 
    #
    KEEPER.aux_keep.add ( a )
    KEEPER.aux_keep.add ( b )
    KEEPEP.aux_keep.add ( result )
    #
    return result

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
    result     = Ostap.MoreRooFit.Addition ( v1 , v2 , name , title )
    #
    KEEPER.aux_keep.add ( v1 )
    KEEPER.aux_keep.add ( v2 )
    KEEPEP.aux_keep.add ( result )
    #
    return result

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
    result     = Ostap.MoreRooFit.Product ( v1 , v2 , name , title ) 
    #
    KEEPER.aux_keep.add ( v1 )
    KEEPER.aux_keep.add ( v2 )
    KEEPEP.aux_keep.add ( result )
    #
    return result


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
    result     = Ostap.MoreRooFit.Subtraction ( v1 , v2 , name , title ) 
    #
    KEEPER.aux_keep.add ( v1 )
    KEEPER.aux_keep.add ( v2 )
    KEEPEP.aux_keep.add ( result )
    #
    return result


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
    result     = Ostap.MoreRooFit.Division ( v1 , v2 , name , title ) 
    #
    KEEPER.aux_keep.add ( v1 )
    KEEPER.aux_keep.add ( v2 )
    KEEPEP.aux_keep.add ( result )
    #
    return result


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
        r = float ( a ) / ( float ( a ) + float ( b ) )  
        return ROOT.RooRealConstant.value ( r )                 ## RETURN
    
    elif f1 : v1 = ROOT.RooRealConstant.value ( float ( v1 ) )        
    elif f2 : v2 = ROOT.RooRealConstant.value ( float ( v2 ) ) 
    #
    result     = Ostap.MoreRooFit.Fraction ( v1 , v2 , name , title ) 
    #
    KEEPER.aux_keep.add ( v1 )
    KEEPER.aux_keep.add ( v2 )
    KEEPEP.aux_keep.add ( result )
    #
    return result


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
        r = ( float ( a ) - float ( b ) ) / ( float ( a ) + float ( b ) )  ## 
        return ROOT.RooRealConstant.value ( r )                 ## RETURN
    
    elif f1 : v1 = ROOT.RooRealConstant.value ( float ( v1 ) )        
    elif f2 : v2 = ROOT.RooRealConstant.value ( float ( v2 ) ) 
    #
    result     = Ostap.MoreRooFit.Asymmetry ( v1 , v2 , name , title ) 
    #
    KEEPER.aux_keep.add ( v1 )
    KEEPER.aux_keep.add ( v2 )
    KEEPEP.aux_keep.add ( result )
    #
    return result

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
        r = float ( a ) ** float ( b ) 
        return ROOT.RooRealConstant.value ( r )                 ## RETURN
    
    elif f1 : v1 = ROOT.RooRealConstant.value ( float ( v1 ) )        
    elif f2 : v2 = ROOT.RooRealConstant.value ( float ( v2 ) ) 
    #
    result     = Ostap.MoreRooFit.Power ( v1 , v2 , name , title ) 
    #
    KEEPER.aux_keep.add ( v1 )
    KEEPER.aux_keep.add ( v2 )
    KEEPEP.aux_keep.add ( result )
    #
    return result


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
 
# =============================================================================
##                                                                      The END 
# =============================================================================
