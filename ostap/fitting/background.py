#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/background.py
#  Set of useful smooth 1D-models to describe smooth ``background'' distribtions
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# =============================================================================
"""Set of useful smooth 1D-models to describe ``background'' distributions"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    'Bkg_pdf'           , ## An exponential function, modulated by positive polynomial
    'PSPol_pdf'         , ## A phase space  function, modulated by positive polynomial
    'PSLeftExpoPol_pdf' , ## A phase space function, modulated by positive polynomial and exponent
    'PolyPos_pdf'       , ## Bernstein positive polynomial
    'PolyEven_pdf'      , ## Bernstine positive even polynomial
    'Monotonic_pdf'     , ## Bernstein positive monotonic polynomial
    'Convex_pdf'        , ## Bernstein positive polynomial with fixed sign first and second derivatives 
    'ConvexOnly_pdf'    , ## Bernstein  positive polynomial with fixed sign second derivatives 
    'Sigmoid_pdf'       , ## Background: sigmoid modulated by positive polynom 
    'TwoExpoPoly_pdf'   , ## difference of two exponents, modulated by positive polynomial
    ##
    'Linear_pdf'        , ## positive linear polynom 
    'Parabolic_pdf'     , ## positive parabolic polynom 
    ## 
    'PSpline_pdf'       , ## positive            spline 
    'MSpline_pdf'       , ## positive monotonic spline 
    'CSpline_pdf'       , ## positive monotonic convex or concave spline 
    'CPSpline_pdf'      , ## positive convex or concave spline 
    ##
    'PS2_pdf'           , ## 2-body phase space (no parameters)
    'PSLeft_pdf'        , ## Low  edge of N-body phase space 
    'PSRight_pdf'       , ## High edge of L-body phase space from N-body decays  
    'PSNL_pdf'          , ## L-body phase space from N-body decays  
    'PS23L_pdf'         , ## 2-body phase space from 3-body decays with orbital momenta
    ##
    'PSSmear_pdf'       , ## smeared (left/right, Gaussian) PhaseSpace-based PDF 
    'PSSmear2_pdf'      , ## smeared (left, generic)        PhaseSpace-based PDF
    ##
    'KarlinShapley_pdf' , ## Kaglin-Shapley positive polynomial 
    'KarlinStudden_pdf' , ## Kaglin-Studden positive polynomial
    ##
    'Rational_pdf'      , ## Rational fuuction
    ## 
    ## get the native RooFit background shapes
    ##
    'RooPoly_pdf'       , ## wrapper for RooPolynomial 
    'RooCheb_pdf'       , ## wrapper for RooChebyshev
    'RooKeys1D_pdf'     , ## wrapper for RooNDKeysPdf 
    ##
    'make_bkg'          , ## helper function to create backgrounds 
    )
# =============================================================================
from   ostap.core.core          import Ostap
from   ostap.core.ostap_types   import integer_types , num_types 
from   ostap.math.base          import iszero
from   ostap.fitting.pdfbasic   import PDF1, Generic1D_pdf, Flat1D,  Sum1D
from   ostap.fitting.fithelpers import Phases, ParamsPoly
from   ostap.utils.ranges       import vrange 
import ROOT, math
# =============================================================================
from   ostap.logger.logger      import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.background' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
## list of "left" phase space functions 
PSL = ( Ostap.Math.PhaseSpaceLeft , Ostap.Math.PhaseSpaceNL ,
        Ostap.Math.PhaseSpace2    , Ostap.Math.PhaseSpace3  , Ostap.Math.PhaseSpace3s )
# =============================================================================
models  = []
# =============================================================================
##  @class PolyBase
#   helper base class to implement various polynomial-like shapes
class PolyBase(PDF1,Phases) :
    """Helper base class to implement various polynomial-like shapes
    """
    def __init__ ( self , name , power , xvar , the_phis = None ) :
        ## check  the arguments 
        PDF1  .__init__ ( self , name  , xvar = xvar )
        Phases.__init__ ( self , power , the_phis    )


# =============================================================================        
## @class  Bkg_pdf
#  The exponential modified with the positive polynomial 
#  @see Ostap::Models::ExpoPositive
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Bkg_pdf(PolyBase) :
    """Exponential function, modulated by the positive polynomial:
    
    f(x) ~ exp(-tau*X) * Pol_n(x)
    where Pol_n(x) is POSITIVE polynomial (Pol_n(x)>=0 over the whole range) 
    
    >>>  mass = ROOT.RooRealVar( ... ) 
    >>>  bkg  = Bkg_pdf ( 'B' , mass , power = 3 )    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   power    = 0     ,   ## degree of polynomial
                   tau      = None  ,   ## exponential slope 
                   xmin     = None  ,  ## optional x-min
                   xmax     = None  ,  ## optional x-max 
                   the_phis = None  ) : ## the phis...
        
        ##            
        PolyBase.__init__  ( self , name , power , xvar , the_phis = the_phis )
        #                
        self.__power = len ( self.phis ) 
        #
        ## xmin/xmax
        xmin , xmax = self.xmnmx ( xmin , xmax )
        #
        ## limits for |tau|
        mmax        = max ( abs ( xmin ) , abs ( xmax ) )
        limits_tau  = -500. / mmax ,  500. / mmax             
        # 
        ## the exponential slope
        #
        self.__tau  = self.make_var ( tau              ,
                                      "tau_%s"  % name ,
                                      "tau(%s)" % name ,
                                      None , 0 , *limits_tau  )
        
        self.pdf  = Ostap.Models.ExpoPositive (
            self.roo_name ( "exppol_"  ) ,
            "Exponential modulated by positive polynom %s" % self.name , 
            self.xvar            ,
            self.tau             ,
            self.phi_list        ,
            xmin                 ,
            xmax                 )

        ## save configuration 
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'power'    : self.power ,            
            'tau'      : self.tau   ,            
            'the_phis' : self.phis  ,            
            'xmin'     : xmin       ,
            'xmax'     : xmax       ,            
            }
        
    @property
    def power ( self ) :
        """``power''-parameter (polynomial order) for expo*pol function"""
        return self.__power
    @property
    def tau ( self ) :
        """``tau''-parameter (exponential slope) for expo*pol function"""
        return self.__tau
    @tau.setter
    def tau ( self , value ) :
        value = float ( value  )
        if self.xminmax() : 
            mn , mx     = self.xminmax()
            mmax        = max ( abs ( mn ) , abs ( mx ) )
            assert -500./mmax <= value <=  500/mmax, "``tau''-parameter is too large"                              
        self.__tau.setVal ( value  ) 
        return self.__tau.getVal() 
           
models.append ( Bkg_pdf ) 

# =============================================================================
## @class  PolyPos_pdf
#  A positive polynomial 
#  @see Ostap::Models::PolyPositive 
#  @see Ostap::Math::Positive
#  @see Linear_pdf 
#  @see Parabolic_pdf 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class PolyPos_pdf(PolyBase) :
    """Positive (Bernstein) polynomial: 
    
    f(x) = Pol_n(x)
    with Pol_n(x)>= 0 over the whole range 
    
    >>>  mass = ROOT.RooRealVar( ... )
    >>>  bkg  = PolyPos_pdf ( 'B' , mass , power = 2 )
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,  ## the name 
                   xvar             ,  ## the variable 
                   power = 1        ,  ## degree of the polynomial
                   xmin     = None  ,  ## optional x-min
                   xmax     = None  ,  ## optional x-max 
                   the_phis = None  ) : 
        #
        PolyBase.__init__ ( self , name , power , xvar , the_phis = the_phis )
        #
        self.__power = len ( self.phis ) 
        #
        ## xmin/xmax
        xmin , xmax = self.xmnmx ( xmin , xmax )
        #
        ## build the model
        self.pdf   = Ostap.Models.PolyPositive (
            self.roo_name ( "pol_"  ) ,
            "Positive polynom %s" % self.name , 
            self.xvar            ,
            self.phi_list        ,
            xmin , xmax )
        
        ## save configuration 
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'power'    : self.power ,            
            'the_phis' : self.phis  ,
            'xmin'     : xmin       ,
            'xmax'     : xmax       ,            
            }
                
    @property
    def power ( self ) :
        """``power''-parameter (polynomial order) for PolyPos function"""
        return self.__power



models.append ( PolyPos_pdf ) 


# =============================================================================
## @class KarlinShapley_pdf
#  A positive polynomial 
#  @see Ostap::Models::KarlinShapley
#  @see Ostap::Math::KarlinShapley
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class KarlinShapley_pdf(PolyBase) :
    """Positive Karlin-Shapley  polynomial: 
    
    f(x) = Pol_n(x)
    with Pol_n(x)>= 0 over the whole range 
    
    >>>  mass = ROOT.RooRealVar( ... )
    >>>  bkg  = KarlinShapley_pdf ( 'B' , mass , power = 2 )
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,  ## the name 
                   xvar             ,  ## the variable 
                   power    = 1     ,  ## degree of the polynomial
                   xmin     = None  ,  ## optional x-min
                   xmax     = None  ,  ## optional x-max 
                   the_phis = None  ) : 
        #
        PolyBase.__init__ ( self , name , power , xvar , the_phis = the_phis )
        #
        self.__power = len ( self.phis ) 
        #
        ## xmin/xmax
        self.__x_min , self.__x_max = self.xmnmx ( xmin , xmax )
        #

        ## build the model
        self.pdf   = Ostap.Models.KarlinShapley (
            self.roo_name ( "ksh_"  ) ,
            "Karlin-Shapley polynomial %s" % self.name , 
            self.xvar            ,
            self.phi_list        ,
            self.x_min , self.x_max )
        
        ## save configuration 
        self.config = {
            'name'     : self.name   ,
            'xvar'     : self.xvar   ,
            'power'    : self.power  ,            
            'the_phis' : self.phis   ,
            'xmin'     : self.x_min  ,
            'xmax'     : self.x_max  ,            
            }
                
    @property
    def power ( self ) :
        """``power''-parameter (polynomial order) for Karlin-Shapley function"""
        return self.__power

    @property
    def x_min ( self ) :
        """'xmin' - minimal x for Karlin-Studden polynomial"""
        return self.__x_min

    @property
    def x_max ( self ) :
        """'xmax' - maximal x for Karlin-Studden polynomial"""
        return self.__x_max

    @property
    def karlin_shapley ( self ) :
        """'karlin-shapley' : polinomial onject"""
        self.pdf.setPars()
        return self.pdf.karlin_shapley()
    
    @property
    def pars ( self ) :
        """'pars' : phases (same as 'phis')"""
        return self.phis

    @property
    def troots ( self ) :
        """'troots' : Karlin-Shapley t-roots"""
        return tuple ( p for p in self.karlin_shapley.troots() )
    
models.append ( KarlinShapley_pdf ) 

# =============================================================================
## @class KarlinStudden_pdf
#  A positive polynomial 
#  @see Ostap::Models::KarlinStudden
#  @see Ostap::Math::KarlinStudden
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class KarlinStudden_pdf(PolyBase) :
    """Positive Karlin-Studden polynomial: 
    
    f(x) = Pol_n(x)
    with Pol_n(x)>= 0 for x > xmin 
    
    >>>  mass = ROOT.RooRealVar( ... )
    >>>  bkg  = KarlinStudden_pdf ( 'B' , mass , power = 2 )
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,  ## the name 
                   xvar             ,  ## the variable 
                   power    = 1     ,  ## degree of the polynomial
                   xmin     = None  ,  ## optional x-min
                   scale    = None  ,  ## scale 
                   the_phis = None  ) : 
        #
        PolyBase.__init__ ( self , name , power , xvar , the_phis = the_phis )
        #
        self.__power = len ( self.phis ) 
        #
        assert isinstance ( xmin , num_types ) or self.xvar.hasMin() , \
               'Invalid setting for xmin!'
        
        if isinstance ( xmin , num_types ) : self.__x_min = float ( xmin )
        else                               : self.__x_min = self.xvar.getMin() 

        assert ( isinstance ( scale , num_types ) and scale ) or self.xvar.hasMax() and self.x_min < self.xvar.getMax() , \
               'Invalid setting for scale!'
        
        if isinstance ( scale , num_types ) and scale : self.__scale =       abs ( float ( scale ) )
        else                                          : self.__scale = 0.5 * abs ( self.xvar.getMax() - self.x_min )
        
        #
        ## build the model
        self.pdf   = Ostap.Models.KarlinStudden (
            self.roo_name ( "kst_"  ) ,
            "Karlin-Studden polynomial %s" % self.name , 
            self.xvar     ,
            self.phi_list ,
            self.x_min    ,
            self.scale    )
        
        ## save configuration 
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'power'    : self.power ,            
            'the_phis' : self.phis  ,
            'xmin'     : self.x_min ,
            'scale'    : self.scale ,            
            }
                
    @property
    def power ( self ) :
        """'power'-parameter (polynomial order) for Karlin-Studden function"""
        return self.__power

    @property
    def x_min ( self ) :
        """'xmin' - minimal x for Karlin-Studden polynomial"""
        return self.__x_min

    @property
    def scale ( self ) :
        """'scale' - scale for Karlin-Studden polynomial"""
        return self.__scale 

    @property
    def pars ( self ) :
        """'pars' : phases (same as 'phis')"""
        return self.phis

    @property
    def karlin_studden ( self ) :
        """'karlin-studden' : polinomial onject"""
        self.pdf.setPars()
        return self.pdf.karlin_studden()

    @property
    def troots ( self ) :
        """'troots' : Karlin-Studden t-roots"""
        self.pdf.setPars()
        return tuple ( p for p in self.pdf.karlin_studden.troots() )
    
    @property
    def zroots ( self ) :
        """'zroots' : Karlin-Studden z-roots"""
        self.pdf.setPars()
        return tuple ( p for p in self.pdf.karlin_studden.zroots() )
    

models.append ( KarlinStudden_pdf ) 


# =============================================================================
## @class Rational_pdf
#  Ratinal function: ratio of two positve bernstein polynomials 
#  @see Ostap::Models::Rational
#  @see Ostap::Math::RationalPositive
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2023-09-15
class Rational_pdf(PolyBase) :
    """Rational fnuction: ratio of two positive Bernstein polynomials  
    - see Ostap.Models.Rational
    - see Ostap.Math.RationalPositive
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,  ## the name 
                   xvar             ,  ## the variable 
                   np       = 1     ,  ## degree if numerator 
                   nq       = 1     ,  ## degree of denumerator 
                   xmin     = None  ,  ## optional x-min
                   xmax     = None  ,  ## optional x-max 
                   the_phis = None  ) :

        assert isinstance ( np , int ) and 0 <= np , "Invalid `np'-parameter"
        assert isinstance ( nq , int ) and 0 <= nq , "Invalid `nq'-parameter"
        
        ## initialize the base 
        PolyBase.__init__ ( self , name , np + nq , xvar , the_phis = the_phis )
        #
        n  = len ( self.phis )
        np = min ( n , np )
        nq =       n - np 
        
        #
        ## xmin/xmax
        self.__x_min , self.__x_max = self.xmnmx ( xmin , xmax )
        #
        
        ## build the model
        self.pdf   = Ostap.Models.Rational (
            self.roo_name ( "rat_"  ) ,
            "Rational function %s" % self.name , 
            self.xvar            ,
            np                   , 
            self.phi_list        ,
            self.x_min           ,
            self.x_max           )
        
        ## save configuration 
        self.config = {
            'name'     : self.name     ,
            'xvar'     : self.xvar     ,
            'np'       : self.pdf.p () ,            
            'nq'       : self.pdf.q () ,            
            'the_phis' : self.phis     ,
            'xmin'     : self.x_min    ,
            'xmax'     : self.x_max    ,            
            }
                
    @property
    def np ( self ) :
        """`np'-parameter - degree of numerator"""
        return self.pdf.p ()
    @property
    def nq ( self ) :
        """`nq'-parameter - degree of denomerator"""
        return self.pdf.q()
    
    @property
    def x_min ( self ) :
        """'x_min' - minimal x for Rational function"""
        return self.__x_min

    @property
    def x_max ( self ) :
        """'x_max' - maximal x for Rational function"""
        return self.__x_max
    
    @property
    def pars ( self ) :
        """'pars' : phases (same as 'phis')"""
        return self.phis

    @property
    def ppars ( self ) :
        """'ppars' : parameters of numerator"""
        return self.pars[ : self.np ]

    @property
    def qpars ( self ) :
        """'qpars' : parameters of denumorator"""
        return self.pars[ self.np : ]

    
models.append ( Rational_pdf ) 

# =============================================================================
## @class Linear_pdf
#  A positive linear polynomial 
#  @see Ostap::Models::PolyPositive 
#  @see Ostap::Math::Positive
#  @see PolyPos_pdf 
#  @see Parabolic_pdf 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-01-31
class Linear_pdf (PolyPos_pdf) :
    """Positive Linear (Bernstein) polynomial:     
    >>>  mass = ROOT.RooRealVar( ... )
    >>>  bkg  = Linear_pdf ( 'B' , mass )
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,  ## the name 
                   xvar             ,  ## the variable 
                   xmin     = None  ,  ## optional x-min
                   xmax     = None  ,  ## optional x-max 
                   the_phis = None  ) : 
        
        ## initialize the base class 
        PolyPos_pdf.__init__ ( self ,
                               name     = name     ,
                               xvar     = xvar     ,
                               power    = 1        ,
                               xmin     = xmin     ,
                               xmax     = xmax     ,
                               the_phis = the_phis ) 
        
        ## save configuration 
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'the_phis' : self.phis  ,
            'xmin'     : xmin       ,
            'xmax'     : xmax       ,            
            }
        
models.append ( Linear_pdf ) 


# =============================================================================
## @class Parabolic_pdf
#  A positive parabolic polynomial 
#  @see Ostap::Models::PolyPositive 
#  @see Ostap::Math::Positive
#  @see PolyPos_pdf 
#  @see Linear_pdf 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-01-31
class Parabolic_pdf (PolyPos_pdf) :
    """Positive parabolic (Bernstein) polynomial:     
    >>>  mass = ROOT.RooRealVar( ... )
    >>>  bkg  = Linear_pdf ( 'B' , mass )
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,  ## the name 
                   xvar             ,  ## the variable 
                   xmin     = None  ,  ## optional x-min
                   xmax     = None  ,  ## optional x-max 
                   the_phis = None  ) : 
        
        ## initialize the base class 
        PolyPos_pdf.__init__ ( self ,
                               name     = name     ,
                               xvar     = xvar     ,
                               power    = 2        ,
                               xmin     = xmin     ,
                               xmax     = xmax     ,
                               the_phis = the_phis ) 
        
        ## save configuration 
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'the_phis' : self.phis  ,
            'xmin'     : xmin       ,
            'xmax'     : xmax       ,            
            }

        
models.append ( Parabolic_pdf ) 

# =============================================================================
## @class  PolyEven_pdf
#  A positive even polynomial
#  \f$ f( \frac{x_{max}+x_{min}}{2} -x ) = f( \frac{x_{max}+x_{min}}{2} +x )\f$ 
#  @see Ostap::Models::PolyPositiveEven 
#  @see Ostap::Math::PositiveEven
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2016-10-03
class PolyEven_pdf(PolyBase) :
    """Positive (Bernstein) even polynomial: 
    
    f(x) = Pol_n(x)
    with Pol_n(x)>= 0 over the whole range
    and  Pol_n(0.5*(xmax+xmin)-x) = Pol_n(0.5*(xmax+xmin)+x) 
    
    >>>  mass = ROOT.RooRealVar( ... )
    >>>  bkg  = PolyPos_pdf ( 'B' , mass , power = 2 )
    """
    ## constructor
    def __init__ ( self             ,
                   name             , ## the name 
                   xvar             , ## the varibale 
                   power = 1        , ## (half)degree of the polynomial
                   xmin     = None  ,  ## optional x-min
                   xmax     = None  ,  ## optional x-max 
                   the_phis = None  ) :
        #
        PolyBase.__init__ ( self , name , power , xvar , the_phis = the_phis )
        self.__power = power
        #        
        ## xmin/xmax
        xmin , xmax = self.xmnmx ( xmin , xmax )
        #
        ## build PDF 
        self.pdf  = Ostap.Models.PolyPositiveEven (
            self.roo_name ( "epol_"  ) ,
            "Positive even polynom %s" % self.name , 
            self.xvar                     ,
            self.phi_list                 ,
            xmin ,
            xmax )
        
        ## save configuration 
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'power'    : self.power ,            
            'the_phis' : self.phis  ,            
            'xmin'     : xmin       ,
            'xmax'     : xmax       ,            
            }
                
    @property
    def power ( self ) :
        """``power''-parameter (polynomial order) for PolyEven function"""
        return self.__power
        
models.append ( PolyEven_pdf ) 

# =============================================================================
## @class  Monotonic_pdf
#  A positive monotonic polynomial 
#  @see Ostap::Models::PolyMonotonic 
#  @see Ostap::Math::Monotonic
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Monotonic_pdf(PolyBase) :
    """Positive monotonic (Bernstein) polynomial:
    
    f(x) = Pol_n(x)
    with f(x)>= 0 over the whole range and
    the derivative f' do not change the sign
    
    >>>  mass = ROOT.RooRealVar( ... )
    
    # increasing background 
    >>>  bkg_inc  = Monotonic_pdf ( 'B1' , mass , power = 2 , increasing = True  )
    
    # decreasing background 
    >>>  bkg_dec  = Monotonic_pdf ( 'B2' , mass , power = 2 , increasing = False  )
    """
    ## constructor
    def __init__ ( self              ,
                   name              ,  ## the name 
                   xvar              ,  ## the variable
                   power      = 2    ,  ## degree of the polynomial
                   increasing = True ,  ## increasing or decreasing ?
                   xmin       = None ,  ## optional x-min
                   xmax       = None ,  ## optional x-max 
                   the_phis   = None ) : 
        #
        PolyBase.__init__ ( self , name , power , xvar , the_phis = the_phis )
        #
        self.__power      = power
        self.__increasing = True if increasing else False 
        #
        ## xmin/xmax
        xmin , xmax = self.xmnmx ( xmin , xmax )
        #
        ## build PDF
        self.pdf  = Ostap.Models.PolyMonotonic (
            self.roo_name ( "mpol_"  ) ,
            "%s polynomial %s" % ( "increasing" if self.increasing else "decreasing" , self.name ) , 
            self.xvar            ,
            self.phi_list        ,
            xmin                 ,
            xmax                 ,           
            self.increasing      )
        
        ## save configuration 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'power'      : self.power      ,            
            'increasing' : self.increasing , 
            'the_phis'   : self.phis       ,             
            'xmin'       : xmin            ,
            'xmax'       : xmax            ,            
           }

    @property
    def power ( self ) :
        """``power''-parameter (polynomial order) for Monotonic function"""
        return self.__power

    @property
    def increasing ( self ) :
        """``increasing''-parameter for Monotonic function"""
        return self.__increasing

    @property
    def decreasing ( self ) :
        """``decreasing''-parameter for Monotonic function"""
        return not self.increasing
    
        
models.append ( Monotonic_pdf ) 
# =============================================================================
## @class  Convex_pdf
#  A positive polynomial with fixed signs of the first and second derivative 
#  @see Ostap::Models::PolyConvex 
#  @see Ostap::Math::Convex
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Convex_pdf(PolyBase) :
    """ Positive monotonic (Bernstein) polynomial with fixed-sign second derivative:
    
    f(x) = Pol_n(x)
    with f(x)>= 0 over the whole range and
    the both frist derivative f' and second derivative f'' do not change the sign
    
    >>>  mass = ROOT.RooRealVar( ... )
    
    # increasing convex background 
    >>>  bkg  = Convex_pdf ( 'B1' , mass , power = 2 , increasing = True  , convex = True  )

    # dereasing convex background 
    >>>  bkg  = Convex_pdf ( 'B2' , mass , power = 2 , increasing = False , convex = True  )

    # increasing concave background 
    >>>  bkg  = Convex_pdf ( 'B3' , mass , power = 2 , increasing = True  , convex = False )

    # dereasing concave background 
    >>>  bkg  = Convex_pdf ( 'B3' , mass , power = 2 , increasing = False , convex = False )   
    """
    ## constructor
    def __init__ ( self              ,
                   name              ,   ## the name 
                   xvar              ,   ## the variable
                   power = 2         ,   ## degree of the polynomial
                   increasing = True ,   ## increasing or decreasing ?
                   convex     = True ,   ## convex or concave ?
                   xmin       = None ,   ## optional x-min
                   xmax       = None ,   ## optional x-max 
                   the_phis   = None ) : ## 
        #
        PolyBase.__init__ ( self , name , power , xvar , the_phis = the_phis )
        #
        self.__power      = power
        self.__increasing = True if increasing else False 
        self.__convex     = True if convex     else False 
        #
        ## xmin/xmax
        xmin , xmax = self.xmnmx ( xmin , xmax )
        #
        ## build PDF 
        self.pdf  = Ostap.Models.PolyConvex (
            self.roo_name ( "cpol_"  ) ,
            "%s %s polynomial %s" % ( "convex"     if self.convex     else "concave"    ,
                                      "increasing" if self.increasing else "decreasing" ,
                                      self.name ) , 
            self.xvar            ,
            self.phi_list        ,
            xmin ,
            xmax , 
            self.increasing      ,
            self.convex          )
        
        ## save configuration 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'power'      : self.power      ,            
            'increasing' : self.increasing , 
            'convex'     : self.convex     , 
            'the_phis'   : self.phis       ,            
            'xmin'       : xmin            ,
            'xmax'       : xmax            ,            
            }

    @property
    def power ( self ) :
        """``power''-parameter (polynomial order) for Convex function"""
        return self.__power

    @property
    def increasing ( self ) :
        """``increasing''-parameter for Convex function"""
        return self.__increasing

    @property
    def decreasing ( self ) :
        """``decreasing''-parameter for Convex function"""
        return not self.increasing
    
    @property
    def convex ( self ) :
        """``convex''-parameter for Convex function"""
        return self.__convex
    
    @property
    def concave ( self ) :
        """``concave''-parameter for Convex function"""
        return not self.convex


models.append ( Convex_pdf ) 
# =============================================================================
## @class  ConvexOnly_pdf
#  A positive polynomial with fixed signs of the second derivative 
#  @see Ostap::Models::PolyConvexOnly 
#  @see Ostap::Math::ConvexOnly
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class ConvexOnly_pdf(PolyBase) :
    """ Positive (Bernstein) polynomial with fixed-sign second derivative:
    
    f(x) = Pol_n(x)
    with f(x)>= 0 over the whole range and
    the second derivative f'' do not change the sign
    
    >>>  mass = ROOT.RooRealVar( ... )
    
    >>>  bkg  = ConvexOnly_pdf ( 'B1' , mass , power = 2 , convex = True  )
    >>>  bkg  = ConvexOnly_pdf ( 'B1' , mass , power = 2 , convex = False  )
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   power    = 2     ,   ## degree of the polynomial
                   convex   = True  ,   ## convex or concave ?
                   xmin     = None  ,   ## optional x-min
                   xmax     = None  ,   ## optional x-max 
                   the_phis = None  ) :
        #
        PolyBase.__init__ ( self , name , power , xvar , the_phis = the_phis )
        #
        self.__power      = power
        self.__convex     = True if convex else False 
        #
        ## xmin/xmax
        xmin , xmax = self.xmnmx ( xmin , xmax )
        #
        ## build PDF 
        self.pdf  = Ostap.Models.PolyConvexOnly (
            self.roo_name ( "copol_"  ) ,
            "%s polynomial %s" % ( "convex"     if self.convex     else "concave"    ,
                                   self.name ) , 
            self.xvar            ,
            self.phi_list        ,
            xmin                 ,
            xmax                 ,
            self.convex          )
        
        ## save configuration 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'power'      : self.power      ,            
            'convex'     : self.convex     , 
            'the_phis'   : self.phis       ,            
            'xmin'       : xmin            ,
            'xmax'       : xmax            ,            
            }
    
    @property
    def power ( self ) :
        """``power''-parameter (polynomial order) for ConvexOnly function"""
        return self.__power

    @property
    def convex ( self ) :
        """``convex''-parameter for ConvexOnly function"""
        return self.__convex
    
    @property
    def concave ( self ) :
        """``concave''-parameter for ConvexOnly function"""
        return not self.convex


models.append ( ConvexOnly_pdf ) 

# =============================================================================
## @class  PSPol_pdf
#  Phase space function modulated with the positive polynomial 
#  @see Ostap::Models::PhaseSpacePol
#  @see Ostap::Math::PhaseSpacePol
#  @see Ostap::Math::PhaseSpaceNL 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class PSPol_pdf(PolyBase) :
    """ The phase space function modified with positive polynomial 
    
    f(x) ~ PS( x ) * Pol_n(x)
    where
    - PS       is a phase-space function
    - Pol_n(x) is POSITIVE polynomial (Pol_n(x)>=0 over the whole range) 

    >>>  mass = ROOT.RooRealVar( ... )
    >>>  ps   = Ostap.Math.PhaseSpaceNL ( low , high , 3 , 4 )  
    >>>  bkg  = PSPol_pdf ( 'B' , mass , ps , power = 2 )    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable 
                   phasespace       ,   ## Ostap::Math::PhaseSpaceNL 
                   power    = 1     ,   ## degree of the polynomial
                   the_phis = None  ) : 
        
        #
        #
        PolyBase.__init__  ( self , name , power , xvar , the_phis = the_phis )
        #
        if   isinstance ( phasespace , Ostap.Math.PhaseSpaceNL ) : pass
        elif isinstance ( phasespace , tuple ) and phasespace :
            phasespace = Ostap.Math.PhaseSpaceNL ( *phasespace )
            logger.debug ('Phase space created: PhaseSpaceNL(%s,%s,%d,%d)' % ( phasespace.lowEdge  () ,
                                                                               phasespace.highEdge () ,
                                                                               phasespace.L        () ,
                                                                               phasespace.N        () ) ) 
        else :
            raise AttributeError ('Illegal type of "phasespace" parameter')
        
        self.__ps    = phasespace  ## Ostap::Math::PhaseSpaceNL
        self.__power = power
        # 
        self.pdf  = Ostap.Models.PhaseSpacePol (
            self.roo_name ( "pspol_"  ) ,
            "Phase space modulated by polynomial %s" % self.name  , 
            self.xvar            ,
            self.phasespace      ,  ## Ostap::Math::PhaseSpaceNL 
            self.phi_list        )
        
        ## save configuration 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'phasespace' : self.phasespace ,            
            'power'      : self.power      ,            
            'the_phis'   : self.phis       ,            
            }

    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (alias for ``x'' or ``xvar'')"""
        return self.xvar
    
    @property
    def phasespace ( self ) :
        """``phasespace''-function for PS*pol function"""
        return self.__ps 

    @property
    def power ( self ) :
        """``power''-parameter (polynomial order) for polynomial function"""
        return self.__power

models.append ( PSPol_pdf ) 

# =============================================================================
## @class  PSLeftExpoPol_pdf
#  Left Phase space factor, exponent and the positive polynomial 
#   \f[ f(x) \propto 
#      \Phi_{l}(x;x_{low}) \mathrm{e}^{-\left|\tau\right| x } P_{N}(x) \f]
#  where :
#   -  \f$  \Phi_{l}(x;x_{low}) \f$  is a phase space of 
#     l-particles near the threshold 
#   -  \f$ P_{N}(x) \f$ is a positive polynomial of degree N
#  @see Ostap::Models::PhaseSpaceLeftExpoPol
#  @see Ostap::Math::PhaseSpaceLeftExpoPol
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class PSLeftExpoPol_pdf(PolyBase) :
    r"""The phase space function modified with positive polynomial 
    
    # f(x) ~ Phi_{l}(x;x_{low}) \mathrm{e}^{-\left|\tau\right| x } P_{N}(x) 
    where :
    - \Phi_{l}(x;x_{low}) is a phase space of l-particles near the threshold 
    -  P_{N}(x) is a positive polynomial of degree N

    >>>  mass = ROOT.RooRealVar( ... )
    >>>  ps   = Ostap.Math.PhaseSpaceLeft ( low , high , 3 , 4 )  
    >>>  bkg  = PSPol_pdf ( 'B' , mass , ps , power = 2 )    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,  ## the name 
                   xvar             ,  ## the variable
                   phasespace       ,  ## phase space object 
                   power    = 2     ,  ## degree of the polynomial
                   tau      = None  ,  ## the exponent 
                   scale    = 1     ,  ## the exponent 
                   the_phis = None  ) :

        #
        ## check the type of phasespace 
        assert isinstance ( phasespace , PSL ) or ( isinstance ( phasespce , integer_types ) and 2 <= ps ) , \
               "Invalid type of ``phasespace'' %s/%s" % ( phasespace , type ( phasespace ) )

        ## initialize the base 
        PolyBase.__init__  ( self , name , power , xvar , the_phis = the_phis )
        #

        self.__phasespace = phasespace 
        self.__power      = power
        
        limits_tau = () 
        if self.xminmax() : 
            mn , mx     = self.xminmax()
            delta       = 0.5 * ( mx - mn ) 
            limits_tau  = -500.0/delta , 500.0/delta

        #
        ## the exponential slope
        #
        self.__tau  = self.make_var ( tau              ,
                                      "tau_%s"  % name ,
                                      "tau(%s)" % name ,
                                      None , 0 , *limits_tau  )

        #
        ## the scale factor 
        #
        self.__scale = self.make_var ( scale             ,
                                      "scale_%s"  % name ,
                                      "scale(%s)" % name ,
                                       None , 1 , 1.e-3 , 1.e+6 )

        self.pdf  = Ostap.Models.PhaseSpaceLeftExpoPol (
            self.roo_name ( "psepol_"  ) ,
            "Phase space and exponent modulated by polynomial %s" % self.name  , 
            self.xvar            ,
            self.phasespace      ,  
            self.tau             , 
            self.scale           , 
            self.phi_list        )
        
        ## save configuration 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'phasespace' : self.phasespace ,            
            'power'      : self.power      ,            
            'tau'        : self.tau        ,            
            'scale'      : self.scale      ,            
            'the_phis'   : self.phis       ,            
            }

    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (alias for ``x'' or ``xvar'')"""
        return self.xvar
    
    @property
    def power ( self ) :
        """``power''-parameter (polynomial order) for PSLeftExpoPol-function"""
        return self.__power
    
    @property
    def phasespace ( self ) :
        """``phasespace''-function for PS*exp*pol function"""
        return self.__phasespace
    
    @property
    def tau ( self ) :
        """``tau''-parameter - exponential slope"""
        return self.__tau
    @tau.setter
    def tau ( self , value ) :
        self.set_value ( self.__tau , float ( value ) ) 
        
    @property
    def scale( self ) :
        """``scale''-parameter - scale factor for phase space function"""
        return self.__scale
    @scale.setter
    def scale ( self , value ) :
        self.set_value ( self.__scale , float ( value ) ) 
    
models.append ( PSLeftExpoPol_pdf ) 

# =============================================================================
## @class  TwoExpoPoly_pdf
#  Difference of two exponents, modulated by positive polynomial 
#  @see Ostap::Models::TwoExpoPositive
#  @see Ostap::Math::TwoExpoPositive
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-03-26
class TwoExpoPoly_pdf(PolyBase) :
    """Difference of two exponential function, modulated by the positive polynomial:
    
    f(x) ~ ( exp(-alpha*x) - exp(-(alpha_delta)*x) *  Pol_n(x)
    where Pol_n(x) is POSITIVE polynomial (Pol_n(x)>=0 over the whole range) 
    
    >>>  mass = ROOT.RooRealVar( ... ) 
    >>>  bkg  = TwoExpoPoly_pdf ( 'B' , mass , alpah = 1 , delta = 1 , x0 = 0 , power = 3 )
    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   alpha    = None  ,   ## the slope of the first exponent 
                   delta    = None  ,   ## (alpha+delta) is the slope of the first exponent
                   x0       = 0     ,   ## f(x)=0 for x<x0 
                   power    = 0     ,   ## degree of polynomial
                   xmin     = None  ,   ## optional x-min
                   xmax     = None  ,   ## optional x-max                    
                   the_phis = None  ) : 
        #
        PolyBase.__init__  ( self , name , power , xvar , the_phis = the_phis  )
        #                
        self.__power = power
        #
        ## xmin/xmax
        xmin , xmax = self.xmnmx ( xmin , xmax )
        #
        ## the exponential slope
        #
        mn , mx      = xmin , xmax 
        dm           = mx - mn 
        mmax         = max ( abs ( mn ) , abs ( mx ) )
        limits_alpha = 1.e-6 , 1.e-16 , 300. / mmax             
        limits_delta = 1.e-6 , 1.e-16 , 300. / mmax             
        limits_x0    = mn , mn - 10 * dm , mx + 10 * dm 


        self.__alpha  = self.make_var ( alpha               ,
                                        "alpha_%s"   % name ,
                                        "#alpha(%s)" % name ,
                                        None , *limits_alpha ) 
        
        self.__delta  = self.make_var ( delta               ,
                                        "delta_%s"   % name ,
                                        "#delta(%s)" % name ,
                                        None , *limits_delta )
        
        self.__x0     = self.make_var ( x0                  ,
                                        "x0_%s"     % name  ,
                                        "x_{0}(%s)" % name  ,
                                        None ,  *limits_x0 )
        
        #
        xmin , xmax = self.xminmax() 
        self.pdf  = Ostap.Models.TwoExpoPositive (
            self.roo_name ( "pol2e_"  ) ,
            "Two exponents  modulated by polynomial %s" % self.name  , 
            self.xvar             ,
            self.alpha            ,
            self.delta            ,
            self.x0               ,
            self.phi_list         ,
            xmin , xmax )
        
        ## save configuration 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'alpha'      : self.alpha      , 
            'delta'      : self.delta      , 
            'x0'         : self.x0         , 
            'power'      : self.power      ,            
            'the_phis'   : self.phis       ,
            'xmin'       : xmin            , 
            'xmax'       : xmax            , 
            }
        
    @property
    def power ( self ) :
        """``power''-parameter (polynomial order) for TwoExpo function"""
        return self.__power

    @property
    def alpha ( self ) :
        """``alpha''-parameter (slope of the leading exponenet) of the TwoExpo function"""
        return self.__alpha
    @alpha.setter
    def alpha ( self , value ) :
        value = float ( value )
        assert 0 < value, "``alpha''-parameter must be positive "
        self.alpha.setVal ( value ) 
        return self.__alpha.getVal()
    
    @property
    def delta ( self ) :
        """``delta''-parameter (second exponent slope is ``alpha+delta'') of the TwoExpo function"""
        return self.__delta
    @delta.setter
    def delta ( self , value ) :
        value = float ( value )
        assert 0 < value, "``delta''-parameter must be positive"
        self.delta.setVal ( value ) 
        return self.__delta.getVal()
    
    @property
    def x0 ( self ) :
        """``x0''-parameter (shift) of the TwoExpo function  (f(x)=0 for x<x0)"""
        return self.__x0
    @x0.setter
    def x0 ( self , value ) :
        value = float ( value )
        self.x0.setVal ( value ) 
        return self.__x0.getVal()
 
models.append ( TwoExpoPoly_pdf ) 


# =============================================================================
## @class  Sigmoid_pdf
#  Sigmoid function modulated wit hpositive polynomial
#  @see Ostap::Models::PolySigmoid
#  @see Ostap::Math::Sigmoid
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Sigmoid_pdf(PolyBase) :
    """A ``sigmoid'' function modulated by the positive (Bernstein) polynomial 
    f(x) = 0.5*(1+tahn(alpha*(x-x0))*Pol_n(x)
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,  ## the name 
                   xvar             ,  ## the variable
                   power    = 2     ,  ## degree of the polynomial
                   alpha    = None  ,  ## 
                   x0       = None  ,  ## x0
                   xmin     = None  ,  ## optional x-min
                   xmax     = None  ,  ## optional x-max 
                   the_phis = None  ) :
        #
        #
        PolyBase.__init__ ( self , name , power , xvar , the_phis = the_phis )
        #
        self.__power      = power
        #
        ## xmin/xmax
        xmin , xmax = self.xmnmx ( xmin , xmax )
        ## 
        mn, mx       = xmin , xmax 
        dx           = mx - mn 
        alpmx        = 20
        limits_alpha = -1000./dx, +1000./dx
        limits_x0    = 0.5 * ( mn + mx ) , mn - 0.2 *  dx , mx + 0.2 * dx

        ## alpha 
        self.__alpha  = self.make_var ( alpha               ,
                                        'alpha_%s'  % name  ,
                                        'alpha(%s)' % name  ,
                                        None , *limits_alpha )
        
        ## x0 
        self.__x0    = self.make_var  ( x0                  ,
                                        'x0_%s'     % name  ,
                                        'x0(%s)'    % name  ,
                                        None , *limits_x0    )

        xmin,xmax = self.xminmax() 
        self.pdf  = Ostap.Models.PolySigmoid (
            self.roo_name ( "sigmoid_"  ) ,
            "Sigmoid modulated by polynomial %s" % self.name  , 
            self.xvar            ,
            self.phi_list        ,
            xmin                 ,
            xmax                 ,
            self.alpha           ,
            self.x0              )
        
        ## save configuration 
        self.config = {
            'name'       : self.name  ,
            'xvar'       : self.xvar  ,
            'power'      : self.power ,            
            'alpha'      : self.alpha , 
            'x0'         : self.x0    , 
            'the_phis'   : self.phis  ,
            'xmin'       : xmin       ,
            'xmax'       : xmax       ,                        
            }

    @property
    def power ( self ) :
        """``power''-parameter (polynomial order) for Sigmoid function"""
        return self.__power

    @property
    def alpha ( self ) :
        """``alpha''-parameter for Sigmoid function"""
        return self.__alpha
    @alpha.setter
    def alpha ( self , value ) :
        value = float ( value )
        self.alpha.setVal ( value ) 
        return self.__alpha.getVal()
    
    @property
    def x0 ( self ) :
        """``x0''-parameter for Sigmoid fuction"""
        return self.__x0
    @x0.setter
    def x0 ( self , value ) :
        value = float ( value )
        self.x0.setVal ( value ) 
        return self.__x0.getVal()
      
        
models.append ( Sigmoid_pdf ) 
# =============================================================================
## @class  PSpline_pdf
#  The special spline for non-negative function
#  Actually it is a sum of M-splines with non-negative coefficients 
#  \f$ f(x) = \sum_i \alpha_i * M_i^k(x) \f$,
#  with constraints  \f$  \sum_i \alpha_i=1\f$ and \f$ 0 \le \alpha_i\f$.
#  @see Ostap::Models::PositiveSpline 
#  @see Ostap::Math::PositiveSpline 
#  @see http://en.wikipedia.org/wiki/M-spline
#  @see http://en.wikipedia.org/wiki/B-spline
#  The constrains are implemented as N-sphere
#  @see Ostap::Math::NSphere
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class PSpline_pdf(PolyBase) :
    """A positive spline, a compositon of M-splines with non-negative coefficients

    >>> mass   = ... ## the variable
    >>> order  = 3   ## spline order
    
    ## create uniform spline 
    >>> inner  = 3   ## number of inner knots between min and max 
    >>> spline = Ostap.Math.PositiveSpline( mass.xmin() , mass.xmax() , inner , order )

    ## create non-uniform spline with
    >>> knots = std.vector('double)()
    >>> knots.push_back ( mass.xmin() )
    >>> knots.push_back ( mass.xmax() )
    >>> knots.push_back ( ... )
    >>> spline = Ostap.Math.PositiveSpline( knots , order )

    >>> bkg = PSpline_pdf ( 'Spline' , mass , spline ) 
    
    """
    ## constructor
    def __init__ ( self             ,
                   name             , ## the name 
                   xvar             , ## the variable
                   spline           , ## the spline object Ostap::Math::PositiveSpline
                   the_phis = None  ) : 
        #
        if   isinstance ( spline , Ostap.Math.PositiveSpline ) :  pass
        elif isinstance ( spline , tuple ) :
            spline = Ostap.Math.PositiveSpline ( *spline )
            logger.debug ( 'Spline created %s' % spline ) 
            
        PolyBase.__init__ ( self , name , spline.npars() , xvar , the_phis = the_phis )
        #
        self.__spline = spline
        #
        
        self.pdf  = Ostap.Models.PositiveSpline (
            self.roo_name ( "pspline_"  ) ,
            "Positive spline %s" % self.name  , 
            self.xvar            ,
            self.spline          , 
            self.phi_list        )

        ## save configuration 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'spline'     : self.spline     ,            
            'the_phis'   : self.phis       ,            
            }

    @property
    def spline( self ) :
        """``spline''-function for PSpline PDF"""
        return self.__spline
        
models.append ( PSpline_pdf ) 
# =============================================================================
## @class  MSpline_pdf
#  The special spline for non-negative monotonic function
#  @see Ostap::Models::MonotonicSpline 
#  @see Ostap::Math::MonotonicSpline 
#  @see http://en.wikipedia.org/wiki/I-spline
#  @see http://en.wikipedia.org/wiki/M-spline
#  @see http://en.wikipedia.org/wiki/B-spline
#  @see http://link.springer.com/chapter/10.1007%2F978-3-0348-7692-6_6
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class MSpline_pdf(PolyBase) :
    """A positive monotonic spline

    >>> mass   = ... ## the variable
    >>> order  = 3   ## spline order
    
    ## create uniform spline 
    >>> inner  = 3   ## number of inner knots between min and max 
    >>> spline = Ostap.Math.MonotonicSpline( mass.xmin() , mass.xmax() , inner , order , True )

    ## create non-uniform spline with
    >>> knots = std.vector('double)()
    >>> knots.push_back ( mass.xmin() )
    >>> knots.push_back ( mass.xmax() )
    >>> knots.push_back ( ... )
    >>> spline = Ostap.Math.MonotonicSpline( knots , order , True )

    >>> bkg = MSpline_pdf ( 'Spline' , mass , spline ) 
    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   spline           ,   ## the spline object Ostap::Math::MonotonicSpline
                   the_phis = None  ) : 
        #
        PolyBase.__init__ ( self , name , spline.npars() , xvar , the_phis = the_phis  )
        #
        self.__spline = spline
        # 
        self.pdf  = Ostap.Models.MonotonicSpline (
            self.roo_name ( "mspline_"  )      ,
            "Monotonic spline %s" % self.name  , 
            self.xvar                          ,
            self.spline                        , 
            self.phi_list                      )

        ## save configuration 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'spline'     : self.spline     ,            
            'the_phis'   : self.phis       ,            
            }

    @property
    def spline( self ) :
        """``spline''-function for MSpline PDF"""
        return self.__spline

models.append ( MSpline_pdf )
# =============================================================================
## @class  CSpline_pdf
#  The special spline for non-negative monotonic convex/concave function
#  @see Ostap::Models::ConvexSpline 
#  @see Ostap::Math::ConvexSpline 
#  @see http://en.wikipedia.org/wiki/I-spline
#  @see http://en.wikipedia.org/wiki/M-spline
#  @see http://en.wikipedia.org/wiki/B-spline
#  @see http://link.springer.com/chapter/10.1007%2F978-3-0348-7692-6_6
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class CSpline_pdf(PolyBase) :
    """A positive monotonic convex/concave spline

    >>> mass   = ... ## the variable
    >>> order  = 3   ## spline order
    
    ## create uniform spline 
    >>> inner  = 3   ## number of inner knots between min and max 
    >>> spline = Ostap.Math.ConvexSpline( mass.xmin() , mass.xmax() , inner , order , True , True )

    ## create non-uniform spline with
    >>> knots = std.vector('double)()
    >>> knots.push_back ( mass.xmin() )
    >>> knots.push_back ( mass.xmax() )
    >>> knots.push_back ( ... )
    >>> spline = Ostap.Math.ConvexSpline( knots , order, True  , True )

    >>> bkg = CSpline_pdf ( 'Spline' , mass , spline ) 
    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   spline           ,   ## the spline object Ostap::Math::ConvexSpline
                   the_phis = None  ) : ##
        #
        PolyBase.__init__ ( self , name , spline.npars() , xvar , the_phis = the_phis )
        self.__spline = spline
        #
        self.pdf  = Ostap.Models.ConvexSpline (
            self.roo_name ( "cspline_"  )   ,
            "Convex spline %s" % self.name  , 
            self.xvar                ,
            self.spline              , 
            self.phi_list            )
        
        ## save configuration 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'spline'     : self.spline     ,            
            'the_phis'   : self.phis       ,            
            }

    @property
    def spline( self ) :
        """``spline''-function for CSpline PDF"""
        return self.__spline

models.append ( CSpline_pdf ) 


# =============================================================================
## @class  CPSpline_pdf
#  The special spline for non-negative convex/concave function
#  @see Ostap::Models::ConvexOnlySpline 
#  @see Ostap::Math::ConvexOnlySpline 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class CPSpline_pdf(PolyBase) :
    """A positive convex/concave spline

    >>> mass   = ... ## the variable
    >>> order  = 3   ## spline order
    
    ## create uniform spline 
    >>> inner  = 3   ## number of inner knots between min and max 
    >>> spline = Ostap.Math.ConvexOnlySpline( mass.xmin() , mass.xmax() , inner , order , True )

    ## create non-uniform spline with
    >>> knots = std.vector('double)()
    >>> knots.push_back ( mass.xmin() )
    >>> knots.push_back ( mass.xmax() )
    >>> knots.push_back ( ... )
    >>> spline = Ostap.Math.ConvexonlySpline( knots , order, True )

    >>> bkg = CPSpline_pdf ( 'Spline' , mass , spline ) 
    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   spline           ,   ## the spline object Ostap::Math::ConvexOnlySpline
                   the_phis = None  ) : 
        #
        PolyBase.__init__ ( self , name , spline.npars() , xvar , the_phis = the_phis )
        self.__spline = spline
        #
        self.pdf  = Ostap.Models.ConvexOnlySpline (
            self.roo_name ( "cospline_"  )   ,
            "Convex spline %s" % self.name  , 
            self.xvar                ,
            self.spline              , 
            self.phi_list            )
        
        ## save configuration 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'spline'     : self.spline     ,            
            'the_phis'   : self.phis       ,            
            }

    @property
    def spline( self ) :
        """``spline''-function for CPSpline PDF"""
        return self.__spline

models.append ( CPSpline_pdf ) 

# =============================================================================
## @class  PS2_pdf
#  Primitive 2-body phase space function
#  @code 
#  mass  = ... ## mass variable
#  m_pi  = 139
#  m_K   = 493 
#  model = PS2_pdf ( 'PDF' , mass , m_K , m_pi )
#  @endcode 
#  @see Ostap::Models::PhaseSpace2
#  @see Ostap::Math::PhaseSpace2
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class PS2_pdf(PDF1) :
    """ Primitive 2-body phase space function
    >>> mass  = ... ## mass variable
    >>> m_pi  = 139
    >>> m_K   = 493 
    >>> model = PS2_pdf ( 'PDF' , mass , m_K , m_pi )   
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   m1               ,   ## the first  mass (constant)
                   m2               ) : ## the second mass (constant)
        
        ## initialize the base 
        PDF1.__init__ ( self , name , xvar = xvar )
        #
        self.__am1 = abs ( float ( m1 ) )
        self.__am2 = abs ( float ( m2 ) )

        am1 = self.m1
        am2 = self.m2
        
        if self.xminmax() :
            mn, mx =  self.xminmax()
            assert am1  + am2 < mx , 'PS2: senseless setting of edges/threshold %s,%s vs %s' % ( mn , mx , am1+am2 )   
            
        self.pdf  = Ostap.Models.PhaseSpace2 (
            self.roo_name ( "ps2_"  )   ,
            "Two body phase space  %s" % self.name  , 
            self.xvar                , self.m1  , self.m2 )
        
        ## save configuration 
        self.config = {
            'name'       : self.name ,
            'xvar'       : self.xvar ,
            'm1'         : self.m1   ,            
            'm2'         : self.m2   ,            
            }
    
    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar

    @property
    def m1 ( self ) :
        """Tthe mass of the first particle for PS2 function"""
        return self.__am1
    @property
    def m2 ( self ) :
        """Tthe mass of the first particle for PS2 function"""
        return self.__am2
        
models.append ( PS2_pdf ) 
# =============================================================================
## @class  PSLeft_pdf
#  Left edge of N-body phase space (with possible scaling)
#  @code 
#  mass  = ... ## mass variable
#  low   = 139 + 139 + 139  ## 3 pion mass
#  model = PSLeft_pdf ( 'PDF' , mass , 3 , low )
#  @endcode 
#  @see Ostap::Models::PhaseSpaceLeft
#  @see Ostap::Math::PhaseSpaceLeft
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class PSLeft_pdf(PDF1) :
    """Left edge of N-body phase space function (with possible scaling)
    >>> mass  = ... ## mass variable
    >>> low   = 139 + 139 + 139  ## 3 pion mass
    >>> model = PSLeft_pdf ( 'PDF' , mass , 3 , low )   
    """
    ## constructor
    def __init__ ( self           ,
                   name           , ## the name 
                   xvar           , ## the variable
                   phasespace     , ## phase space object 
                   left   = None  ,
                   scale  = None  ) : 
        
        ## check the type of phasespace 
        assert isinstance ( phasespace , PSL ) or ( isinstance ( phasespce , integer_types ) and 2 <= ps ) , \
               "Invalid type of ``phasespace'' %s/%s" % ( phasespace , type ( phasespace ) )

        # 
        ## initialize the base 
        PDF1.__init__ ( self , name , xvar = xvar )
        #
        
        self.__phasespace = phasespace

        if self.xminmax() :
            mn , mx = self.xminmax()
            mn = max ( 0 , mn - 0.25 * ( mx - mn ) )
            lminmax = mn , mx 
        else :
            lminmax = ()
            
        ## left threshold variable 
        self.__left = self.make_var ( left                     ,
                                      'left_%s'    % self.name ,
                                      'm_left(%s)' % self.name ,
                                      None , *lminmax ) 
        
        ## scale variable 
        if scale is None :
            self.__scale = ROOT.RooFit.RooConst ( 1.0 )
        else             :
            self.__scale = self.make_var ( scale                   ,
                                           'scale_%s'  % sefl.name ,
                                           'scale(%s)' % self.name ,
                                           None , 0.1 , 1.e+6 )
            
        if  self.xminmax() and self.left.minmax() : 
            mn  , mx  = self.xminmax()
            lmn , lmx = self.left.minmax()            
            assert lmn < mx, "PSLeft_pdf: senseless setting of edges/thresholds: %s,%s vs %s,%s"  % (  mn, mx , lmn, lmx ) 


        self.pdf  = Ostap.Models.PhaseSpaceLeft (
            self.roo_name ( "psl_"  )   ,
            "Left edge of the phase space %s" % self.name  , 
            self.xvar       ,
            self.left       ,
            self.scale      ,
            self.phasespace )
        
        ## save configuration 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'phasespace' : self.phasespace ,            
            'left'       : self.left       ,            
            'scale'      : self.scale      ,            
            }
        
    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar
    
    @property
    def scale ( self ) : 
        """``scale'' : scale-factor for left-threshold"""
        return self.__scale
    @scale.setter 
    def scale ( self , value ) :
        self.set_value ( self.__scale , float ( value ) , ok = lambda a , b : 0 < b ) 
        
    @property
    def left( self ) :
        """(left) threshold for N-body phase space"""
        return self.__left
    @left.setter
    def left ( self , value ) :
        self.set_value ( self.__left , float ( value ) )
    
    @property
    def threshold ( self ) :
        """``threshold'' : threshold position, same as ``left''"""
        return self.left
    @threshold.setter
    def threshold ( self , value ) :
        self.left = value
        
    @property
    def phasespace ( self ) :
        """``phasespace'' : phase space object (or number of particles)"""
        return self.__phasespace

models.append ( PSLeft_pdf ) 
# =============================================================================
## @class  PSRight_pdf
#  Right edge of L-body phase space from N-body decay 
#  @see Ostap::Models::PhaseSpaceRight
#  @see Ostap::Math::PhaseSpaceRight
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class PSRight_pdf(PDF1) :
    """ Right edge of L-body phase space for N-body decay 
    >>> mass  = ... ## mass variable
    >>> high  = 5.278-3.096 ## m(B)-m(J/psi)
    >>> model = PSRight_pdf ( 'PDF' , mass , high , 3  , 4 )   ## e.g. 3 pions from B->J/psi 3pi
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   L                ,   ## L  
                   N                ,   ## N
                   right   = None   ) : 
        #
        PDF1.__init__ ( self , name , xvar = xvar )
        #
        self.__L = L 
        self.__N = N

        limits = ()  
        if self.xminmax() :
            mn , mx  = self.xminmax()
            mx = mx + ( mx  - mn )
            limits = mn , mx
            
        self.__right = self.make_var ( right ,
                                       'right_%s'      % name  ,
                                       'm_{right}(%s)' % name  ,
                                       None , *limits )

        if self.xminmax() and self.right.minmax() :
            mn  , mx  = self      .xminmax()
            rmn , rmx = self.right. minmax()            
            assert rmx > mn, "PSRight_pdf: senseless setting of edges/thresholds: %s,%s vs %s,%s"  % (  mn, mx , rmn, rmx ) 
            
        self.pdf  = Ostap.Models.PhaseSpaceRight (
            self.roo_name ( "psr_"  )   ,
            "Right edge of phase space %s" % self.name  , 
            self.xvar  ,
            self.right ,
            self.L     , 
            self.N     ) 

        ## save configuration 
        self.config = {
            'name'       : self.name  ,
            'xvar'       : self.xvar  ,
            'L'          : self.L     ,            
            'N'          : self.N     ,            
            'right'      : self.right ,            
            }

    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar
    
    @property
    def N ( self ) :
        """define ``L-from-N''phase space"""
        return self.__N
    @property
    def L ( self ) :
        """define ``L-from-N''phase space"""
        return self.__L

    @property
    def right( self ) :
        """(Right) threshold for ``L-from-N''-body phase space"""
        return self.__right
    @right.setter
    def right ( self , value ) :
        value = float ( value )
        self.__right.setVal ( value  )
        return self.__right.getVal()
    
models.append ( PSRight_pdf ) 
# =============================================================================
## @class  PSNL_pdf
#  L-body phase space from N-body decay
#  @code 
#  mass  = ... ## mass variable
#  low   = 3*139       ## 3*m(pi)
#  high  = 5.278-3.096 ## m(B)-m(J/psi)
#  ## e.g. 3 pions from B->J/psi 3pi    
#  model = PSNL_pdf ( 'PDF' , mass , high , 3  , 4 , low , high )
#  @endcode 
#  @see Ostap::Models::PhaseSpaceNL
#  @see Ostap::Math::PhaseSpaceNL
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class PSNL_pdf(PDF1) :
    """L-body phase space from N-body decay
    >>> mass  = ... ## mass variable
    >>> low   = 3*139       ## 3*m(pi)
    >>> high  = 5.278-3.096 ## m(B)-m(J/psi)
    ## e.g. 3 pions from B->J/psi 3pi    
    >>> model = PSNL_pdf ( 'PDF' , mass , high , 3  , 4 , low , high )
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable 
                   L                ,   ## L  
                   N                ,   ## N
                   left  = None     , 
                   right = None     ) : 

        assert isinstance ( L , int ) and 2<=L , 'PSNL: invalid L=%s (must be L>=2)' % L 
        assert isinstance ( N , int ) and L< N , 'PSNL: invalid N=%s (must be N>L)'  % N

        ## initialize the base 
        PDF1.__init__ ( self , name , xvar = xvar )
        #
        self.__L = L
        self.__N = N
        #
        limits_left  = ()
        limits_right = ()
        if self.xminmax() :
            mn, mx = self.xminmax()
            dx     = mx - mn 
            limits_left  = 0.95 * mn + 0.05 * mx , mn - dx , mx + dx 
            limits_right = 0.05 * mn + 0.95 * mx , mn - dx , mx + dx 
            
        self.__left  = self.make_var ( left ,
                                       'left_%s'        % name ,
                                       'm_{left}(%s)'   % name ,
                                       None , *limits_left )
        
        self.__right = self.make_var ( right ,
                                       'right_%s'       % name ,
                                       'm_{right}(%s)'  % name ,
                                       None , *limits_right )

        ## pdf 
        self.pdf  = Ostap.Models.PhaseSpaceNL (
            self.roo_name ( "ps_"  )   ,
            "Phase space %s" % self.name  , 
            self.xvar  ,
            self.left  ,
            self.right ,
            self.L     , 
            self.N     )
       
        ## save configuration 
        self.config = {
            'name'       : self.name  ,
            'xvar'       : self.xvar  ,
            'L'          : self.L     ,            
            'N'          : self.N     ,            
            'left'       : self.left  ,            
            'right'      : self.right ,            
            }

    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar
    
    @property
    def N ( self ) :
        """``N'':define ``L-from-N''phase space"""
        return self.__N
    @property
    def L ( self ) :
        """``L'':define ``L-from-N''phase space"""
        return self.__L

    @property
    def left( self ) :
        """(Left) threshold for ``L-from-N''-body phase space"""
        return self.__left
    @left.setter
    def left ( self , value ) :
        value = float ( value )
        self.__left.setVal ( value  )
        return self.__left.getVal()
    @property
    def right( self ) :
        """(Right) threshold for ``L-from-N''-body phase space"""
        return self.__right
    @right.setter
    def right ( self , value ) :
        value = float ( value )
        self.__right.setVal ( value  )
        return self.__right.getVal()      

models.append ( PSNL_pdf )

# =============================================================================
## @class  PS23L_pdf
#  2-body phase space from 3-body decay with orbital momenta 
#  e.g. m(2pi) from Bs -> J/psi pi pi decay
#  @code 
#  mass  = ... # mass variable for 2 pions
#  m_pi  = 139
#  m_psi = 3096
#  m_Bs  = 5366
#  dalitz = Ostap.Kinematics.Dalitz ( m_Bs , m_psi , m_pi , m_pi ) 
#  ## S-wave 
#  model = PS23L_pdf ( 'S' , mass , dalitz , 1 , 0 )
#  @endcode
#  @see Ostap::Models::PhaseSpace23L
#  @see Ostap::Math::PhaseSpace23L
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class PS23L_pdf(PDF1) :
    """2-body phase space from 3-body decay with orbital momenta
    e.g. m(2pi) from Bs -> J/psi pi pi decay
    >>> mass  = ... # mass variable for 2 pions 
    >>> m_pi  = 139
    >>> m_psi = 3096
    >>> m_Bs  = 5366
    >>> dalitz = Ostap.Kinematics.Dalitz ( m_Bs , m_psi , m_pi , m_pi )     
    ## S-wave 
    >>> model = PS23L_pdf ( 'S' , mass , dalitz , 1 , 0 ) 
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   dalitz           ,   ## Dalizt configurtaion 
                   L                ,   ## orbital momentum between (1,2) and 3
                   l  = 0           ) : ## orbital momentum between 1 and 2
        #
        ## initialize the base 
        PDF1.__init__ ( self , name , xvar = xvar )
        #
        assert isinstance ( dalitz , Ostap.Kinematics.Dalitz ), \
               "Invalid type of ``dalitz''!"
        self.__dalitz = dalitz 

        assert  isinstance ( L , integer_types ) and 0 <= L < 10 , "Invalid ``L'' for the phase space function"
        assert  isinstance ( l , integer_types ) and 0 <= l < 10 , "Invalid ``l'' for the phase space function"

        self.__L  = L
        self.__l  = l

        self.pdf  = Ostap.Models.PhaseSpace23L (
            self.roo_name ( "ps23_"  )   ,
            "Phase space 2 from 3 %s" % self.name  , 
            self.xvar   ,
            self.dalitz , 
            self.L     ,
            self.l     )
        
        ## save configuration 
        self.config = {
            'name'       : self.name   ,
            'xvar'       : self.xvar   ,
            'dalitz'     : self.dalitz ,
            'L'          : self.L      ,            
            'l'          : self.l      ,            
            }

    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar

    @property
    def dalitz ( self ) :
        """``dalitz'' : Dalitz's plot configurartion """
        return self.__dalitz

    @property
    def l  ( self ) :
        """``l''-parameter, the orbital momentum between particles ``(1)'' and ``(2)''"""
        return self.__l

    @property
    def L  ( self ) :
        """``L''-parameter, the orbital momentum between particles ``(1,2)'' and ``3''"""
        return self.__L
        
models.append ( PS23L_pdf ) 

# ==============================================================================
## @class PSSmear_pdf
#  Usefull class to represent "smear" phase space function
#  @code
#  pspdf = PSPol_pdf ( ... )
#  smeared_pdf = PSSmear_pdf ( pspdf , sigma = 2.5 * MeV , step = 0.5 , nstep = 8 ) 
#  @endcode
#  It is very useful to smear left or right thresholds for phase-space-based PDFs.
#
#  This PDF is a weighted sum of <code>2*nstep+1</code> components
#  with fixed coefficients. Each component corresponds to small shoft in
#  left or right edge of the phase space, and the components 
#  are weighted accoring to gaussian function
#
#  It is *not* a real convolution with Gaussian function, but some approximation
#  taking <code>step</endcode> small enough and <code>nstep</code> large enough,
#  such as product of step* and step exceed 3, one can get very good approximation
#  to real convolution.
# 
#  The gaussian function for smear/convolution is evaluated at following points
#  \f$ x_0, x_0\pm s\sigma,x_0\pm 2s\sigma,...x_0\pm ns\sigma,...\f$, where
#   - \f$\sigma\f$ corresponds to <code>sigma</code> 
#   - \f$s\f$ corresponds to <code>step</code> 
#   - \f$n\f$ corresponds to <code>nsteps</code> 
# 
#  The phase-space-based PDF is required to have a method <code>phasespace</code>
#  that returns object of the type :
#  - class Ostap::Math::PhaseSpaceNL    
#  - class Ostap::Math::PhaseSpaceLeft
#  - class Ostap::Math::PhaseSpaceRight 
#
#  The phase-space-based PDF is required to have constructor with keyword
#  <code>phasespace</code>
# 
class PSSmear_pdf ( PDF1 ) :
    """Usefull class to represent ``smear'' phase space function
    >>> pspdf = PSPol_pdf ( ... )
    >>> smeared_pdf = PSSmear_pdf ( pspdf , sigma = 2.5 * MeV , step = 0.5 , nstep = 8 ) 

    This PDF is a weighted sum of ``2*nstep+1`` components
    with fixed coefficients. Each component corresponds to small shift in
    left or right edge of the phase space, and the components
    are weighted accoring to gaussian function
    
    It is *not* a real convolution with Gaussian function, but some approximation
    taking ``step`` small enough and ``nstep`` large enough,
    such as product of ``step`` and ``nstep`` exceeds 3,
    one can get very good approximation to real convolution
    
    The gaussian function for smear/convolution is evaluated at folloiwing points
    

    The phase-space-based PDF is required to have a method <code>phasespace</code>
    that returns object of the type :
    - class Ostap.Math.PhaseSpaceNL    
    - class Ostap.Math.PhaseSpaceLeft
    - class Ostap.Math.PhaseSpaceRight 
    
    The phase-space-based PDF is required to have constructor with keyword ``phasespace``    
    """
    def __init__ ( self          ,
                   pdf0          ,
                   sigma         ,
                   nsteps = 10   , 
                   step   = 0.35 , 
                   name   = ''   ) :
        
        assert sigma      and isinstance ( sigma      , num_types     ) ,\
               "Invalid value for ``sigma'' : %s" % sigma
        assert nsteps and isinstance ( nsteps, integer_types ) and 0 < nsteps ,\
               "Invalid value for ``nsteps'': %s" % nsteps
        assert step       and isinstance ( step       , num_types     ) and 0 < step < 1.5 ,\
               "Invalid value for ``step''  : %s" % step 

        self.__sigma  = sigma 
        self.__step   = step 
        self.__nsteps = nsteps
        
        ## initialize the base
        name =  name if name else pdf0.name + '_smear'
        PDF1.__init__ ( self , name , xvar = pdf0.xvar )

        ## keep the template 
        self.__pdf0 = pdf0

        ##
        ps0    = self.pdf0.phasespace
        asigma = abs ( sigma ) 

        left   = False
        rright = False
        
        if   isinstance  ( ps0 , Ostap.Math.PhaseSpaceNL    ) :
            
            xmn , xmx , L , N = ps0.lowEdge () , ps0.highEdge () , ps0.L() , ps0.N() 
            
            ## ## 1.75 is a safety factor here 
            ## assert xmx - xmn > 1.75 * 0.5 * asigma * nhalfsigma, \
            ##        'Interval is too small (%s,%s) for such large smearing %s*%s' % (
            ##    xmn , xmx , sigma , 0.5 * nhalfsigma )
            
            if   0 < sigma : PS = lambda x : Ostap.Math.PhaseSpaceNL ( xmn + x , xmx , L , N )
            elif 0 > sigma : PS = lambda x : Ostap.Math.PhaseSpaceNL ( xmn , xmx + x , L , N )

        elif isinstance  ( ps0 , Ostap.Math.PhaseSpaceLeft  ) :
            
            N = ps0.N ()
            if 2 <= N :
                t , s = ps0.threshold() , ps0.scale() 
                PS = lambda x : Ostap.Math.PhaseSpaceLeft ( t + x , N , s )
            else      :
                ps2 = ps0.ps2() 
                m1 , m2 , s = ps2.m1() , ps2.m2() , ps0.scale() 
                PS2 = lambda x : Ostap.Math.PhaseSpace2    ( m1 + x    , m2 )
                PS  = lambda x : Ostap.Math.PhaseSpaceLeft ( PS2 ( x ) , s  )
            
        elif isinstance  ( ps0 , Ostap.Math.PhaseSpaceRight ) :

            t , L , N = ps0.thresholsd() , ps0.L() , ps0.N() 
            PS = lambda x : Ostap.Math.PhaseSpaceRight ( t + x , L , N ) 

        else :
            raise AttributeError ('Illegal type of "phasespace"!')
        
        from ostap.math.math_ve import gauss_pdf as gpdf 
        from ostap.math.math_ve import gauss_cdf as gcdf 

        self.__pdfs = []
        self.__ws   = [] 
        for i in range ( 1 , self.nsteps + 1 ) :

            delta = 1.0 * self.step * i * asigma 

            psp = PS (  delta )
            psm = PS ( -delta )

            pdfp = self.pdf0.clone ( name = self.generate_name ( name + '%dp' % i ) , phasespace = psp )
            pdfm = self.pdf0.clone ( name = self.generate_name ( name + '%dm' % i ) , phasespace = psm )
            
            fraction = ROOT.RooConstVar ( 'fF%s%dfix' % ( name , i ) , '' , 0.5 )
            pdfi     = Sum1D ( pdfs = ( pdfp , pdfm ) , name = self.generate_name ( name + '%d' % i  ) , fraction = fraction )
            
            ## wp = gpdf ( step * i ) 
            wp = gcdf ( step * ( i + 0.5 ) ) - gcdf ( step * ( i - 0.5 ) ) 
            wm = wp
            wi = wp + wm
            
            self.__pdfs.append ( pdfi )
            self.__ws  .append ( wi   )

        self.__pdfs.reverse ()
        self.__ws  .reverse ()
        
        self.__pdfs.append  ( self.pdf0  )
        
        ## self.__ws  .append  ( gpdf ( 0 ) )
        self.__ws  .append  ( gcdf ( 0.5 * step ) - gcdf ( -0.5 * step ) )
        
        self.__pdfs = tuple ( self.__pdfs )

        ## normalize fractions 
        wsum = sum ( self.__ws )
        self.__ws = [ v/wsum for v in self.__ws ]
        self.__ws = tuple ( self.__ws [ :-1 ] )   ## and skip the last!
        
        ## create fixed fractions 
        self.__ff = []
        for i , w in enumerate ( self.__ws  ) :
            f = ROOT.RooConstVar ( 'fF%spm%dfix' % ( name , i + 1 ) , '' , w )
            self.__ff.append ( f )
        self.__ff = tuple ( self.__ff )
                
        ## ## create the final PDF 
        ## self.pdf , self.__fractions , _ = self.add_pdf (
        ##     [ i.pdf for i in self.__pdfs ] ,
        ##     self.name                 ,          
        ##     'Smeared phasel space %s' % self.name , 
        ##     'some_pattern%d'          , 
        ##     'some_pattern%d'          ,
        ##     recursive = False         ,
        ##     fractions = self.__ff     )
        
        ## combine all components into single PDF 
        self.__combined_pdf = Sum1D ( pdfs      = self.pdfs      ,
                                      xvar      = self.xvar      ,
                                      name      = 'Smeared phase space %s' % self.name , 
                                      recursive = False          ,
                                      fractions = self.__ff      )
        
        self.pdf = self.combined_pdf.pdf
        
        self.config = {
            'name'   : self.name   ,            
            'pdf0'   : self.pdf0   ,
            'sigma'  : self.sigma  ,
            'step'   : self.step   ,
            'nsteps' : self.nsteps ,            
            }
        
    @property
    def sigma  ( self ):
        """``sigma'' : resolution parameter for PhaseSpace smearing"""
        return self.__sigma
    @property
    def step   ( self ) :
        """``step''  : step (fraction of sigma) for PhaseSpace smearing"""
        return self.__step
    @property
    def nsteps ( self ) :
        """``nstep's'  : number of steps for PhaseSpace smearing"""
        return self.__nsteps
    @property
    def pdf0   ( self ) :
        """``pdf0'' : non-smeared original PDF"""
        return self.__pdf0
    @property
    def pdfs   ( self ) :
        """``pdfs'' : the tuple of all PDFs"""
        return self.__pdfs
    @property
    def ws     ( self ) :
        """``ws'' : calculates fractions"""
        return self.__ws
    @property
    def ffs    ( self ) :
        """``ffs'' : calculated fractions"""
        return self.__ff 
    @property
    def fractions ( self ) :
        """``fractions'' : calcualted fractions"""
        return self.combined_pdf.fractions 
    @property
    def combined_pdf ( self ) :
        """``combinedpdf'' : the final combined PDF (e.g. for drawing)"""
        return self.__combinedpdf
    
# ==============================================================================
## @class PSSmear2_pdf
#  Usefull class to represent "smear" phase space function
#  @code
#  pspdf       = PSPol_pdf ( ... )
#  gamma       = 10 * MeV 
#  shape       = lambda x : 1.0/(x*x+0.25*gamma*gamma)
#  values      = vrange ( -5.0*gamma , 5.0 * gamma , 100 ) 
#  smeared_pdf = PSSmear2_pdf ( pspdf , profile = shape , values = values ) 
#  @endcode
#  It is very useful to smear left thresholds for phase-space-based PDFs.
#
#  This PDF is a weighted sum of <code>len(values)</code> components
#  with fixed coefficients. Each component corresponds to small shift in
#  left edge of the phase space, and the components 
#  are weighted accoring to profile function
#
#  It is *not* a real convolution with profile function, but some approximation
# 
#  The phase-space-based PDF is required to have a method <code>phasespace</code>
#  that returns object of the type :
#  - Ostap::Math::PhaseSpaceNL    
#  - Ostap::Math::PhaseSpace2    
#  - Ostap::Math::PhaseSpace3    
#  - Ostap::Math::PhaseSpace3s    
#  - Ostap::Math::PhaseSpaceLeft
#
#  The phase-space-based PDF is required to have constructor with keyword
#  <code>phasespace</code>
#  @see Ostap::Math::PhaseSpaceNL    
#  @see Ostap::Math::PhaseSpace2    
#  @see Ostap::Math::PhaseSpace3    
#  @see Ostap::Math::PhaseSpace3s    
#  @see Ostap::Math::PhaseSpaceLeft
class PSSmear2_pdf ( PDF1 ) :
    """ Usefull class to represent ``smear'' phase space function
    >>> pspdf  = PSPol_pdf ( ... )
    >>> gamma  = 10 * MeV 
    >>> shape  = lambda x : 1.0/(x*x+0.25*gamma*gamma)
    >>> values = vrange ( -5.0*gamma , 5.0 * gamma , 100 ) 
    >>> smeared_pdf = PSSmear2_pdf ( pspdf , profile = shape , values = values )
    
    - It is very useful to smear left thresholds for phase-space-based PDFs.
    
    - This PDF is a weighted sum of `len(values)` components
    with fixed coefficients. Each component corresponds to small shift in
    left edge of the phase space, and the components 
    are weighted accoring to profile function
    
    - It is *not* a real convolution with profile function, but some approximation
    
    - The phase-space-based PDF is required to have a method `phasespace`
    that returns object of the type :
    - `Ostap.Math.PhaseSpaceNL`    
    - `Ostap.Math.PhaseSpace2`    
    - `Ostap.Math.PhaseSpace3`    
    - `Ostap.Math.PhaseSpace3s`    
    - `Ostap.Math.PhaseSpaceLeft`
    """
    def __init__ ( self          ,
                   pdf0          ,
                   profile       ,
                   values        , 
                   name   = ''   ) :

        assert callable ( profile ) , "PSSmear2_pdf: ``profile'' must be callable!"

        self.__profile  = profile
        
        ## initialize the base
        name =  name if name else pdf0.name + '_smear2'
        PDF1.__init__ ( self , name , xvar = pdf0.xvar )

        ## keep the template 
        self.__pdf0 = pdf0

        vals = list ( set ( [ v for v in values ] ) )
        
        vals.sort()
        self.__values = tuple ( vals )
        
        nn = len ( self.values )        
        assert 3 < nn , 'PSSmear2_pdf: At least three different points must be specified!'
        
        ##
        ps0    = self.pdf0.phasespace

        if   isinstance  ( ps0 , Ostap.Math.PhaseSpaceNL    ) :
            
            xmn , xmx , L , N = ps0.lowEdge () , ps0.highEdge () , ps0.L() , ps0.N() 
            
            PS = lambda x : Ostap.Math.PhaseSpaceNL ( xmn + x , xmx , L , N )

        elif isinstance  ( ps0 , Ostap.Math.PhaseSpace2     ) :

            m1 , m2 = ps0.m1() , ps0.m2() 
            PS = lambda x : Ostap.Math.PhaseSpace2 ( m1 + x , m2 )
            
        elif isinstance  ( ps0 , Ostap.Math.PhaseSpace3     ) :

            m1 , m2 , m3 = ps0.m1() , ps0.m2() , ps0.m3() 
            PS = lambda x : Ostap.Math.PhaseSpace3 ( m1 + x , m2 , m3 )
            
        elif isinstance  ( ps0 , Ostap.Math.PhaseSpace3s    ) :

            m1 , m2 , m3 = ps0.m1() , ps0.m2() , ps0.m3() 
            PS = lambda x : Ostap.Math.PhaseSpace3s ( m1 + x , m2 , m3 )
            
        elif isinstance  ( ps0 , Ostap.Math.PhaseSpaceLeft  ) :
            
            N = ps0.N ()
            if 2 <= N :
                t , s = ps0.threshold() , ps0.scale() 
                PS = lambda x : Ostap.Math.PhaseSpaceLeft ( t + x , N , s )
            else      :
                ps2 = ps0.ps2() 
                m1 , m2 , s = ps2.m1() , ps2.m2() , ps0.scale() 
                PS2 = lambda x : Ostap.Math.PhaseSpace2    ( m1 + x    , m2 )
                PS  = lambda x : Ostap.Math.PhaseSpaceLeft ( PS2 ( x ) , s  )

        else :
            
            raise AttributeError ('Illegal type of "phasespace"!')

        ## get numerical integration 
        from ostap.math.integral import integral as _integral 

        self.__pdfs = []
        self.__ws   = []

        xminmax     = self.xminmax() 

        
        nn = len ( self.values )        
        for i , v in enumerate ( self.values ) : 

            first = ( i     == 0  ) 
            last  = ( i + 1 == nn )             
            
            if first  : xmin =         self.values [ i ]
            else      : xmin = 0.5 * ( self.values [ i ] + self.values [ i - 1 ] )
            
            if last   : xmax =         self.values [ i ] 
            else      : xmax = 0.5 * ( self.values [ i ] + self.values [ i + 1 ] )

            ## calculate the contribution of this component
            iint = _integral ( profile , xmin = xmin , xmax = xmax )
            ## if first or last : iint *= 2 ## Needed? 

            wv   = iint
            
            psv  = PS ( v )
            
            ## skip this component
            if xminmax and xminmax[1] <= psv.lowEdge() : continue 
            
            pdfv = self.pdf0.clone ( name = self.generate_name ( self.name + '%dv' % i ) , phasespace = psv )

            self.__pdfs.append ( pdfv )
            self.__ws  .append ( wv   )


        assert 2 <= len ( self.__pdfs ) , 'Not enough valid components!'

        self.__pdfs = tuple ( self.__pdfs )

        ## normalize fractions
        wsum = sum ( self.__ws )
        assert 0 < wsum , 'PSSmear2_pdf: illegal sum of weights: %s/%d' % ( wsum , len ( self.__ws ) )
        
        self.__ws   = [ v/wsum for v in self.__ws ]
        self.__ws   = tuple ( self.__ws [ :-1 ] )   ## and skip the last!

        ## create fixed fractions 
        self.__ff = []
        for i , w in enumerate ( self.__ws  ) :
            f = ROOT.RooConstVar ( self.generate_name ( 'f' + self.name + '_Ffix%d' % i ) , 'fixed fraction' , w )
            self.__ff.append ( f )
        self.__ff = tuple ( self.__ff )

        ## ## create the final PDF 
        ## self.pdf , self.__fractions , _ = self.add_pdf (
        ##     [ i.pdf for i in self.__pdfs ] ,
        ##     self.name                 ,          
        ##     'Smeared phase-space %s' % self.name , 
        ##     'some_pattern%d'          , 
        ##     'some_pattern%d'          ,
        ##     recursive = False         ,
        ##     fractions = self.__ff     )
        
        ## combine all components into single PDF 
        self.__combined_pdf = Sum1D ( pdfs      = self.pdfs      ,
                                      xvar      = self.xvar      ,
                                      name      = 'Smeared phase space %s' % self.name , 
                                      recursive = False          ,
                                      fractions = self.__ff      )
        
        self.pdf = self.combined_pdf.pdf
        
        self.config = {
            'name'    : self.name    ,            
            'pdf0'    : self.pdf0    ,
            'values'  : self.values  , 
            'profile' : self.profile ,  
            }
        
    @property
    def profile ( self ) :
        """``profile'' : profiel for smearing"""
        return self.__profile 
    @property
    def pdf0   ( self ) :
        """``pdf0'' : non-smeared original PDF"""
        return self.__pdf0
    @property
    def pdfs   ( self ) :
        """``pdfs'' : the tuple of all PDFs"""
        return self.__pdfs
    @property
    def ws     ( self ) :
        """``ws'' : calculated fractions"""
        return self.__ws
    @property
    def ffs    ( self ) :
        """``ffs'' : calculated fractions"""
        return self.__ff 
    @property
    def values ( self ) :
        """``values'' : list of abscissas where profile is evaluated"""
        return self.__values 
    @property
    def fractions ( self ) :
        """``fractions'' : calcualted fractions"""
        return self.combined_pdf.fractions 
    @property
    def combined_pdf ( self ) :
        """``combined_pdf'' : the final combined PDF (e.g. for drawing)"""
        return self.__combined_pdf

    
# ==============================================================================
##  @class RooPolyBase
#   helper base class to implement various polynomial-like shapes
class RooPolyBase(PDF1,ParamsPoly) :
    """Helper base class to implement various polynomial-like shapes
    """
    def __init__ ( self , name , xvar , power = 1 , pars = None ) :
        ## check  the arguments 
        PDF1       .__init__  ( self , name          , xvar = xvar )
        ParamsPoly .__init__  ( self , npars = power , pars = pars )
        
# =============================================================================        
## @class RooPoly_pdf
#  Trivial Ostap wrapper for the native <code>RooPolynomial</code> PDF from RooFit
#  @code
#  xvar = ...
#  poly = RooPoly_pdf ( 'P4' , xvar , 4 ) ;
#  poly = RooPoly_pdf ( 'P3' , xvar , coefficients = [ 1, 2 , 3 ] )
#  @endcode
#  @see RooPolynomial 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-04-27
class RooPoly_pdf(RooPolyBase) :
    """Trivial Ostap wrapper for the native RooPolynomial PDF from RooFit
    - see ROOT.RooPolynomial
    >>> xvar = ...
    >>> poly = RooPoly_pdf ( 'P4' , xvar , 4 ) ;
    >>> poly = RooPoly_pdf ( 'P3' , xvar , coefficients = [ 1, 2 , 3 ] )
    """
    ## constructor
    def __init__ ( self                ,
                   name                , ## the name 
                   xvar                , ## the variable
                   power        = 1    , ## degree of polynomial
                   coefficients = [] ) : ## the list of coefficients 

        coeffs = [ c for c in coefficients ] 
        if not coeffs :
            if isinstance ( xvar , ROOT.RooAbsReal ) and xvar.hasMin() and xvar.hasMax() : 
                from ostap.core.core import iszero 
                xmin , xmax = xvar.getMin() , xvar.getMax() 
                limits = []
                for p in range ( 1 , power + 1 ) :
                    vv   = tuple ( v ** p for v in vrange ( xmin , xmax ,  ) )
                    vmin = min ( vv )
                    vmax = max ( vv )
                    a1   = +1.e+7 if iszero ( vmin ) else -1.0/vmin
                    a2   = -1.e+7 if iszero ( vmax ) else -1.0/vmax
                    amin = min ( 0 , a1 , a2 )
                    amax = max ( 0 , a1 , a2 )
                    lims = 0 , amin , amax
                    limits.append ( lims )
                coeffs = limits
                    
        ## initialize the base class 
        RooPolyBase.__init__ (  self  , name , xvar = xvar , power = power , pars = coeffs  )

        ## create PDF
        self.pdf = ROOT.RooPolynomial (
            self.roo_name ( 'rpoly_' ) ,
            'RooPolynomial %s|%d' % ( self.name , self.power ) ,
            self.xvar , self.pars_lst )

        ## save configuration 
        self.config = {
            'name'        : self.name         ,
            'xvar'        : self.xvar         ,
            'power'       : self.power        ,            
            'coefficients': self.coefficients ,            
            }

    @property
    def power ( self ) :
        """'power' : degree/order/power of polynomial"""
        return self.npars
    
    @property
    def coefficients ( self ) :
        """'coefficients' : polynomial coeffcients"""
        return self.pars 

# =============================================================================        
## @class RooCheb_pdf
#  Trivial Ostap wrapper for the native <code>RooChebyshev</code> PDF from RooFit
#  @code
#  xvar = ...
#  poly = RooCheb_pdf ( 'P4' , xvar , 4 ) ;
#  poly = RooCheb_pdf ( 'P3' , xvar , coefficients = [ 1, 2 , 3 ] )
#  @endcode
#  @see RooChebyshev
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-04-27
class RooCheb_pdf(RooPolyBase) :
    """Trivial Ostap wrapper for the native RooChebyshev PDF from RooFit
    - see ROOT.RooChebushev 
    >>> xvar = ...
    >>> poly = RooCheb_pdf ( 'P4' , xvar , 4 ) ;
    >>> poly = RooCheb_pdf ( 'P3' , xvar , coefficients = [ 1, 2 , 3 ] )
    """
    ## constructor
    def __init__ ( self                ,
                   name                , ## the name 
                   xvar                , ## the variable
                   power        = 0    , ## degree of polynomial
                   coefficients = [] ) : ## the list of coefficients 

        coeffs = [ c for c in coefficients ] 
        if not coeffs :
            coeffs = [] 
            for n in range ( power ) : 
                limits = 0 , -1 , 1
                coeffs.append ( limits ) 

        ## initialize the base class 
        RooPolyBase.__init__ ( self , name , xvar , power ,  pars = coeffs )

        ## create PDF 
        self.pdf = ROOT.RooChebychev (
            self.roo_name ( 'rcheb_' ) ,
            'RooChebychev %s|%d' % ( self.name , self.power ) ,
            self.xvar , self.pars_lst )
        
        ## save configuration 
        self.config = {
            'name'        : self.name         ,
            'xvar'        : self.xvar         ,
            'power'       : self.power        ,            
            'coefficients': self.coefficients ,            
            }        

    @property
    def power ( self ) :
        """'power' : degree/order/power of polynomial"""
        return self.npars
    
    @property
    def coefficients ( self ) :
        """'coefficients' : polynomial coeffcients"""
        return self.pars
    
# =============================================================================        
## @class RooKeys1D_pdf
#  Trivial Ostap wrapper for the native <code>RooNDKeysPdf</code> from RooFit
#  @code
#  xvar = ...
#  pdf  = RooKeys1D_pdf ( 'P4' , xvar , data  ) ;
#  @endcode
#  @see RooNDKeysPdf
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-04-27
class RooKeys1D_pdf(PDF1) :
    """Trivial Ostap wrapper for the native RooNDKeysPdf from RooFit
    - see ROOT.RooKeysPdf
    >>> xvar = ...
    >>> pdf  = RooKeys_pdf ( 'Keys' , xvar , data  )
    """
    ## constructor
    def __init__ ( self           ,
                   name           ,   ## the name 
                   xvar           ,   ## the variable
                   data           ,   ## data set 
                   options = 'ma' ,   ## options 
                   rho     = 1    ,   ## global scale 
                   nsigma  = 3    ,   ##
                   rotate  = True ,   
                   sort    = True ) : 
        
        ## initialize the base class 
        PDF1.__init__ (  self , name , xvar = xvar )

        self.__data    = data
        self.__rho     = rho 
        self.__options = options 
        self.__nsigma  = nsigma 
        self.__rotate  = True if rotate else False 
        self.__sort    = True if sort   else False 

        self.__keys_vlst = ROOT.RooArgList()
        self.__keys_vlst.add ( self.xvar )
        
        ## create PDF
        self.pdf = ROOT.RooNDKeysPdf (
            self.roo_name ( 'keys1_' ) , 
            "Keys 1D %s"  % self.name ,
            self.__keys_vlst  ,
            self.data    ,
            self.options , 
            self.rho     , 
            self.nsigma  , 
            self.rotate  , 
            self.sort    ) 

        ## save configuration
        self.config = {
            'name'    : self.name    ,
            'xvar'    : self.xvar    ,
            'data'    : self.data    ,            
            'options' : self.options ,            
            'rho'     : self.rho     ,            
            'nsigma'  : self.nsigma  ,            
            'rotate'  : self.rotate  ,            
            'sort'    : self.sort    ,            
            }

    @property
    def data   ( self ) :
        """`data' : the actual data set for RooNDKeysPdf"""
        return self.__data
    @property
    def options ( self ) :
        """`options' : `options' string for RooNDKeysPdf"""
        return self.__options        
    @property
    def rho    ( self )  :
        """`rho' : `rho' parameter for RooNDKeysPdf"""
        return self.__rho 
    @property
    def rotate    ( self )  :
        """`rotate' : `rotate' flag for RooNDKeysPdf"""
        return self.__rotate
    @property
    def sort      ( self )  :
        """`sort' : `sort' flag for RooNDKeysPdf"""
        return self.__sort 
    @property
    def nsigma    ( self )  :
        """`nsigma' : `nsigma' parameter for RooNDKeysPdf"""
        return self.__nsigma
    
# =============================================================================
## create popular 1D `background'  function
#  @param bkg  the type of background function/PDF
#  @param name the name of background function/PDF
#  @param xvar the observable
#  Possible values for <code>bkg</code>:
#  - None or 0          : <code>Flat1D</code>
#  - positive integer N : <code>Bkg_pdf(power=N)</code> 
#  - negative integer K : <code>PolyPos_pdf(power=abs(K))</code> 
#  - any Ostap/PDF      : PDF will be copied or cloned  
#  - any RooAbsPdf    P : <code>Generic1D_pdf(pdf=P)</code> 
#  - any RooAbsReal   V : <code>Bkg_pdf(power=0,tau=V)</code> 
#  - math.exp           : <code>Bkg_pdf(power=0)</code>
#  - ''  , 'const', 'constant' , 'flat' , 'uniform' : <code>Flat1D</code>
#  - 'p0', 'pol0' , 'poly0' : <code>Flat1D</code>
#  - 'e' , 'exp'  , 'expo'  : <code>Bkg_pdf(power=0)</code>
#  - 'e+', 'exp+' , 'expo+' : <code>Bkg_pdf(power=0)</code> with tau>0
#  - 'e-', 'exp-' , 'expo-' : <code>Bkg_pdf(power=0)</code> with tau<0
#  - 'e0', 'exp0' , 'expo0' : <code>Bkg_pdf(power=0)</code>
#  - 'eN', 'expN' , 'expoN' : <code>Bkg_pdf(power=N)</code>
#  - 'pN', 'polN' , 'polyN' : <code>PolyPos_pdf(power=N)</code>
#  - 'iN', 'incN' , 'incrN','increasingN' : <code>Monotonic_pdf(power=N,increasing=True)</code>
#  - 'dN', 'decN' , 'decrN','decreasingN' : <code>Monotonic_pdf(power=N,increasing=False)</code>
#  - 'cxN' , 'convexN'                    : <code>ConvexOnly_pdf(power=N,convex=True)</code>
#  - 'cvN' , 'concaveN'                   : <code> ConvexOnly_pdf(power=N,convex=False)</code>
#  - 'roopolyN','rpN'                     : <code>RooPoly_pdf (power=N)</code>>
#  - 'roochebyshevN',roochebN'            : <code>RooCheb_pdf (power=N)</code> 
#  - 'chebyshevN', 'chebN','rcN'          : <code>RooCheb_pdf (power=N)</code>
# 
#  @see ostap.fitting.basic.PDF.make_bkg 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2015-04-03
def make_bkg ( bkg , name , xvar , logger = None , **kwargs ) :
    """Helper function to create various popular 1D-background models
    
    Possible values for ``bkg'':
    
    - None or 0                               : Flat1D
    - positive integer ``N''                  : Bkg_pdf(power=N)
    - negative integer ``K''                  : PolyPos_pdf(power=abs(K))
    - any Ostap-PDF                           : PDF will be copied or cloned  
    - RooAbsPdf      ``pdf''                  : Generic1D_pdf(pdf=pdf)
    - RooAbsReal     ``var''                  : Bkg_pdf(power=0,tau=var)
    - math.exp                                : Bkg_pdf(power=0)
    - 'const' or 'constant'                   : Flat1D
    - '' , 'flat' or 'uniform'                : Flat1D
    - 'e' , 'exp'  or 'expo'                  : Bkg_pdf(power=0)
    - 'e+', 'exp+' or 'expo+'                 : Bkg_pdf(power=0) with tau>0 (increasing)
    - 'e-', 'exp-' or 'expo-'                 : Bkg_pdf(power=0) with tau<0 (decreasing)
    - 'e0', 'exp0' or 'expo0'                 : Bkg_pdf(power=0) 
    - 'eN', 'expN' or 'expoN'                 : Bkg_pdf(power=N)
    - 'p0', 'pol0' or 'poly0'                 : Flat1D
    - 'pN', 'polN' or 'polyN'                 : PolyPos_pdf(power=N)
    - 'iN', 'incN' , 'incrN' or 'increasingN' : Monotonic_pdf(power=N,increasing=True)
    - 'dN', 'decN' , 'decrN' or 'decreasingN' : Monotonic_pdf(power=N,increasing=False)
    - 'cxN' , 'convexN'                       : ConvexOnly_pdf(power=N,convex=True)
    - 'cvN' , 'concaveN'                      : ConvexOnly_pdf(power=N,convex=False)
    - 'roopolyN','rpN'                        : RooPoly_pdf (power=N) 
    - 'roochebyshevN' or roochebyshevN'       : RooCheb_pdf (power=N) 
    - 'chebyshevN', 'chebN' or 'rcN'          : RooCheb_pdf (power=N) 
    
    >>> x =   .. ## the variable
    
    ## None or non-negative integer:  construct PDF using Bkg_pdf 
    >>> bkg1  = make_bkg ( 3      , 'B1' , x ) ## use Bkg_pdf ( 'B1' , x , power = 3 )
    
    ## generic RooAbsPdf: use this PDF
    >>> pdf   = RooPolynomial( ... )
    >>> bkg2  = make_bkg ( pdf    , 'B2' , x ) ## use Generic1D_pdf ( pdf , x , 'B2' )
    
    ## some Ostap-based model, use it as it is  
    >>> model = Convex_pdf ( ... )
    >>> bkg3  = make_bkg ( models , 'B3' , x )  
    
    ## some RooAbsReal, use is as exponenial slope for Bkg_pdf
    >>> tau  = RooRealVar( ....  )
    >>> bkg4 = make_bkg ( tau     , 'B4' , x ) 
    """
    if not logger : logger = globals()['logger'] 
    
    model = None

    prefix = kwargs.pop ( 'prefix' , '' )
    suffix = kwargs.pop ( 'suffix' , '' )
    if prefix : name = prefix + name  
    if suffix : name = name   + suffix 
    
    if   bkg is None :         
        model = Flat1D ( name = name , xvar =  xvar )
        if kwargs : logger.warning ('make_bkg: kwargs %s are ignored' % kwargs )
        
    ## regular case: use Bkg_pdf or PolyPos_pdf as baseline background shapes 
    elif isinstance ( bkg , integer_types ) :

        if   0 < bkg : model =     Bkg_pdf ( name , power =       bkg   , xvar = xvar , **kwargs )
        elif 0 > bkg : model = PolyPos_pdf ( name , power = abs ( bkg ) , xvar = xvar , **kwargs )
        else         :
            model = Flat1D      ( name = name             , xvar = xvar )
            if kwargs : logger.warning ( 'make_bkg: kwargs %s are ignored' % kwargs )

    ## native RooFit pdf ? 
    elif isinstance ( bkg , ROOT.RooAbsPdf ) :

        ## use Generic1D_pdf 
        model = Generic1D_pdf ( bkg , xvar = xvar , name = name ) 
        if kwargs : logger.warning ('make_bkg: kwargs %s are ignored' % kwargs )

    ## some Ostap-based background model ?
    elif isinstance ( bkg , PDF1 ) : 

        ## return the same model/PDF 
        if   xvar is bkg.xvar and not  kwargs :
            model = bkg  ##  use the same stuff 
            logger.debug ( 'make_bkg: %s model is copied to %s' % ( bkg , model ) )

        else :           ## make a clone : 
            model = bkg.clone ( name = name , xvar = xvar , **kwargs )
            logger.debug ( 'make_bkg: %s model is cloned to %s' % ( bkg , model ) )
        
    ## interprete it as exponential slope for Bkg-pdf 
    elif isinstance ( bkg , ROOT.RooAbsReal ) \
             or   isinstance ( bkg , float )  \
             or ( isinstance ( bkg , tuple ) and 1 < len ( bkg ) <=3 ) :
        
        model = Bkg_pdf ( name , xvar = xvar , tau = bkg , power = 0 , **kwargs )
        if kwargs : logger.warning ( 'make_bkg: kwargs %s are ignored' % kwargs )

    ## exponent 
    elif bkg is math.exp : 
        
        model = Bkg_pdf ( name , xvar = xvar , power = 0 , **kwargs )
        if kwargs : logger.warning ( 'make_bkg: kwargs %s are ignored' % kwargs )

    ## strings ....
    elif isinstance ( bkg , str ) :

        bkg = bkg.strip().lower()

        if   bkg in ( '' , 'const' , 'constant' , 'flat' , 'uniform' , 'p0' , 'pol0' , 'poly0' ) :
            return make_bkg ( 0 , name , xvar   , logger = logger , **kwargs ) 
        elif bkg in ( 'e' , 'exp' , 'expo' , 'e0' , 'exp0' , 'expo0' ) :
            return Bkg_pdf ( name , xvar = xvar , power = 0 , **kwargs )        
        elif bkg in ( 'e+' , 'exp+' , 'expo+' ) : 
            model = Bkg_pdf ( name , xvar = xvar , power = 0 , **kwargs )
            model.tau.setMin ( 0 )            
            return model 
        elif bkg in ( 'e-' , 'exp-' , 'expo-' ) :
            model = Bkg_pdf ( name , xvar = xvar , power = 0 , **kwargs )
            model.tau.setMax ( 0 ) 
            return model
        
        import re        
        poly = re.search ( r'(poly|pol|p)(( *)|(_*))(?P<degree>\d{1,2})' , bkg , re.IGNORECASE )
        if poly :
            degree = -1 * abs ( int ( poly.group ( 'degree' ) ) ) 
            return make_bkg ( degree  , name ,  xvar , logger = logger , **kwargs  )
        
        expo = re.search ( r'(expo|exp|e)(( *)|(_*))(?P<degree>\d{1,2})' , bkg , re.IGNORECASE )
        if expo :
            degree =            int ( expo.group ( 'degree' ) ) 
            return make_bkg ( degree  , name ,  xvar , logger = logger , **kwargs  )
        
        incr = re.search ( r'(increasing|increase|incr|inc|i)(( *)|(_*))(?P<degree>\d{1,2})' , bkg , re.IGNORECASE )
        if incr : 
            degree = int ( incr.group ( 'degree' ) )
            bkg    = Monotonic_pdf ( name , xvar , power = degree , increasing = True )
            return make_bkg ( bkg , name ,  xvar , logger = logger , **kwargs  )
        
        decr = re.search ( r'(decreasing|decrease|decr|dec|d)(( *)|(_*))(?P<degree>\d{1,2})' , bkg , re.IGNORECASE )
        if decr : 
            degree = int ( decr.group ( 'degree' ) )
            bkg    = Monotonic_pdf ( name , xvar , power = degree , increasing = False )
            return make_bkg ( bkg , name ,  xvar , logger = logger , **kwargs  )
        
        decr = re.search ( r'(convex|cx)(( *)|(_*))(?P<degree>\d{1,2})' , bkg , re.IGNORECASE )
        if decr :
            degree = int ( decr.group ( 'degree' ) )
            bkg    = ConvexOnly_pdf ( name , xvar , power = degree , convex = True )
            return make_bkg ( bkg , name ,  xvar , logger = logger , **kwargs  )

        decr = re.search ( r'(concave|cv)(( *)|(_*))(?P<degree>\d{1,2})' , bkg , re.IGNORECASE )
        if decr :
            degree = int ( decr.group ( 'degree' ) )
            bkg    = ConvexOnly_pdf ( name , xvar , power = degree , convex = False )
            return make_bkg ( bkg , name ,  xvar , logger = logger , **kwargs  )

        decr = re.search ( r'(convex|cx)(( *)|(_*))(decreasing|dec)(( *)|(_*))(?P<degree>\d{1,3})' , bkg , re.IGNORECASE )
        if decr :
            degree = int ( decr.group ( 'degree' ) )
            bkg    = Convex_pdf ( name , xvar , power = degree , increasing = False, convex = True )
            return make_bkg ( bkg , name ,  xvar , logger = logger , **kwargs  )
            
        decr = re.search ( r'(convex|cx)(( *)|(_*))(increasing|inc)(( *)|(_*))(?P<degree>\d{1,3})' , bkg , re.IGNORECASE )
        if decr :
            degree = int ( decr.group ( 'degree' ) )
            bkg    = Convex_pdf ( name , xvar , power = degree , increasing = True , convex = True )
            return make_bkg ( bkg , name ,  xvar , logger = logger , **kwargs  )
            
        decr = re.search ( r'(concave|cv)(( *)|(_*))(decreasing|dec)(( *)|(_*))(?P<degree>\d{1,3})' , bkg , re.IGNORECASE )
        if decr :
            degree = int ( decr.group ( 'degree' ) )
            bkg    = Convex_pdf ( name , xvar , power = degree , increasing = False, convex = False )
            return make_bkg ( bkg , name ,  xvar , logger = logger , **kwargs  )
            
        decr = re.search ( r'(concave|cv)(( *)|(_*))(increasing|inc)(( *)|(_*))(?P<degree>\d{1,3})' , bkg , re.IGNORECASE )
        if decr :
            degree = int ( decr.group ( 'degree' ) )
            bkg    = Convex_pdf ( name , xvar , power = degree , increasing = True , convex = False )
            return make_bkg ( bkg , name ,  xvar , logger = logger , **kwargs  )

        decr = re.search ( r'(roochebyshev|roocheb|chebyshev|cheb|rc)(( *)|(_*))(?P<degree>\d{1,2})' , bkg , re.IGNORECASE )
        if decr :
            degree = int ( decr.group ( 'degree' ) )
            bkg    = RooCheb_pdf ( name , xvar , power = degree )
            return make_bkg ( bkg , name ,  xvar , logger = logger , **kwargs  )
        
        decr = re.search ( r'(roopoly|rp|r)(( *)|(_*))(?P<degree>\d{1,2})' , bkg , re.IGNORECASE )
        if decr :
            degree = int ( decr.group ( 'degree' ) )
            bkg    = RooPoly_pdf ( name , xvar , power = degree )
            return make_bkg ( bkg , name ,  xvar , logger = logger , **kwargs  )
        
    if model :
        logger.debug ( 'make_bkg: created model is %s' % model ) 
        return model

    raise  TypeError("make_bkg: Wrong type of 'bkg' object: %s/%s " % ( bkg , type ( bkg ) ) ) 


# =============================================================================
if '__main__' == __name__ : 

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )

# =============================================================================
##                                                                      The END 
# =============================================================================
