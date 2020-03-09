#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/roofuncs.py
#
#  Set of useful functions for various 1D and 2D fits
#  It includes
#  - some empricial PDFs to describe narrow peaks: Gauss, CrystalBall, ....
#  - some PDF to describe "wide" peaks: BreitWigner,LASS, Bugg, Flatter, ...
#  - some useful PDFs to describe smooth background: phase space ;
#    expo times polynomial; phase space times polynomial, ...
#  - set of smooth non-facrorizeable model for 2D fits 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
# 
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
    )
# =============================================================================
import ROOT, math
# =============================================================================
from   ostap.logger.logger          import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.roofuncs' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
from ostap.core.core                import Ostap 
from ostap.fitting.utils            import ParamsPoly , ShiftScalePoly  
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
                   func     ,
                   a  = 0.0 ,
                   b  = 1.0 ) :

        FUNC.__init__ ( self , name , xvar = xvar )


        assert isinstance ( func , ( ROOT.RooAbsReal, FUNC ) ) ,\
               "Invalid type of ``fnuc'' %s/%s" % ( func , type ( func ) )

        self.__func0 = func

        self.__c  = func.fun if isnstance ( func , FUNC ) else func 

        self.__a  = self.make_var ( a ,
                                    "a_%s"              % self.name ,
                                    "shift/bias for %s" % self.name , False , a ) 
        self.__b  = self.make_var ( b ,
                                    "b_%s"              % self.name ,
                                    "scale for %s"      % self.name , False , b ) 
        
        self.fun = Ostap.MoreRooFit.ScaleAndShift (
            "scaleshift_%s"   % self.name ,
            "Scale&Shift(%s)" % self.name ,
            self.a ,
            self.b ,
            self.c ) 
        
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
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
 
# =============================================================================
##                                                                      The END 
# =============================================================================
