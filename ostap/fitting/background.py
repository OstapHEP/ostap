#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file backgtround.py
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
    'Bkg_pdf'         , ## An exponential function, modulated by positive polynomial
    'PSPol_pdf'       , ## A phase space  function, modulated by positive polynomial
    'PolyPos_pdf'     , ## A positive polynomial
    'PolyEven_pdf'    , ## A positive even polynomial
    'Monotonic_pdf'  , ## A positive monotonic polynomial
    'Convex_pdf'      , ## A positive polynomial with fixed sign first and second derivatives 
    'ConvexOnly_pdf'  , ## A positive polynomial with fixed sign second derivatives 
    'Sigmoid_pdf'     , ## Background: sigmoid modulated by positive polynom 
    'TwoExpoPoly_pdf' , ## difference of two exponents, modulated by positive polynomial
    ##
    'PSpline_pdf'     , ## positive            spline 
    'MSpline_pdf'     , ## positive monotonic spline 
    'CSpline_pdf'     , ## positive monotonic convex or concave spline 
    'CPSpline_pdf'    , ## positive convex or concave spline 
    ##
    'PS2_pdf'         , ## 2-body phase space (no parameters)
    'PSLeft_pdf'      , ## Low  edge of N-body phase space 
    'PSRight_pdf'     , ## High edge of L-body phase space from N-body decays  
    'PSNL_pdf'        , ## L-body phase space from N-body decays  
    'PS23L_pdf'       , ## 2-body phase space from 3-body decays with orbital momenta
    ##
    'make_bkg'        , ## helper function to create backgrounds 
    )
# =============================================================================
import ROOT, math
from   ostap.core.core     import cpp, Ostap
from   ostap.math.base     import iszero
from   ostap.fitting.basic import PDF
from   ostap.fitting.utils import Phases 
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.models_bkg' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
models = []
# =============================================================================
##  @class PolyBase
#   helper base class to implement various polynomial-like shapes
class PolyBase(PDF,Phases) :
    """Helper base class to implement various polynomial-like shapes
    """
    def __init__ ( self , name , power , xvar = None , the_phis = None ) :
        ## check  the arguments 
        xvar = self.make_var  ( xvar , 'xvar' , 'x-variable' )
        PDF   .__init__ ( self , name  , xvar      )
        Phases.__init__ ( self , power , the_phis  )
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
                   the_phis = None  ) : ## the phis...
        
        ##            
        PolyBase.__init__  ( self , name , power , xvar , the_phis )
        #                
        self.__power = power
        #
        limits_tau = () 
        if self.xminmax() : 
            mn , mx     = self.xminmax()
            mmax        = max ( abs ( mn ) , abs ( mx ) )
            limits_tau  = -500. / mmax ,  500. / mmax             
        # 
        ## the exponential slope
        #
        self.__tau  = self.make_var ( tau              ,
                                "tau_%s"  % name ,
                                "tau(%s)" % name , tau , 0 , *limits_tau  )
        
        self.pdf  = Ostap.Models.ExpoPositive (
            'expopos_%s'  % name ,
            'expopos(%s)' % name ,
            self.xvar            ,
            self.tau             ,
            self.phi_list        ,
            *self.xminmax ()     )

        ## save configuration 
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'power'    : self.power ,            
            'tau'      : self.tau   ,            
            'the_phis' : self.phis  ,            
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
                   xvar             ,  ## the varibale 
                   power = 1        ,  ## degree of the polynomial
                   the_phis = None  ) : 
        #
        PolyBase.__init__ ( self , name , power , xvar , the_phis )
        #
        self.__power = power
        #
        xmin, xmax = self.xminmax () 
        self.pdf   = Ostap.Models.PolyPositive (
            'pp_%s'            % name ,
            'PolyPositive(%s)' % name ,
            self.xvar            ,
            self.phi_list        ,
            xmin , xmax )
        
        ## save configuration 
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'power'    : self.power ,            
            'the_phis' : self.phis  ,            
            }
                
    @property
    def power ( self ) :
        """``power''-parameter (polynomial order) for PolyPos function"""
        return self.__power



models.append ( PolyPos_pdf ) 

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
                   the_phis = None  ) :
        #
        PolyBase.__init__ ( self , name , power , xvar , the_phis )
        self.__power = power
        ##
        xmin , xmax = self.xminmax() 
        self.pdf  = Ostap.Models.PolyPositiveEven (
            'ppe_%s'               % name ,
            'PolyPositiveEven(%s)' % name ,
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
                   name              ,   ## the name 
                   xvar              ,   ## the variable
                   power      = 2    ,   ## degree of the polynomial
                   increasing = True ,   ## increasing or decreasing ?
                   the_phis   = None ) : 
        #
        PolyBase.__init__ ( self , name , power , xvar )
        #
        self.__power      = power
        self.__increasing = True if increasing else False 
        #
        xmin,  xmax = self.xminmax() 
        self.pdf  = Ostap.Models.PolyMonotonic (
            'pp_%s'              % name ,
            'PolyMonotonic(%s)' % name ,
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
    """Positive monotonic (Bernstein) polynomial with fixed-sign second derivative:
    
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
                   the_phis   = None ) : ## 
        #
        PolyBase.__init__ ( self , name , power , xvar )
        #
        self.__power      = power
        self.__increasing = True if increasing else False 
        self.__convex     = True if convex     else False 
        #
        xmin,xmax = self.xminmax()
        
        self.pdf  = Ostap.Models.PolyConvex (
            'pp_%s'          % name ,
            'PolyConvex(%s)' % name ,
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
    """Positive (Bernstein) polynomial with fixed-sign second derivative:
    
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
                   the_phis = None  ) :
        #
        PolyBase.__init__ ( self , name , power , xvar , the_phis )
        #
        self.__power      = power
        self.__convex     = True if convex else False 
        #
        xmin,xmax = self.xminmax() 
        self.pdf  = Ostap.Models.PolyConvexOnly (
            'pp_%s'          % name ,
            'PolyConvex(%s)' % name ,
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
    """The phase space function modified with positive polynomial 
    
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
                   xvar             ,   ## the varibale 
                   phasespace       ,   ## Ostap::Math::PhaseSpaceNL 
                   power    = 1     ,   ## degree of the polynomial
                   the_phis = None  ) : 
        
        #
        #
        PolyBase.__init__  ( self , name , power , xvar , the_phis )
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
            'pspol_%s'          % name ,
            'PhaseSpacePol(%s)' % name ,
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
    def power ( self ) :
        """``power''-parameter (polynomial order) for PS*pol function"""
        return self.__power
    @property
    def phasespace ( self ) :
        """``phasespace''-function for PS*pol function"""
        return self.__ps 

models.append ( PSPol_pdf ) 
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
    >>>  bkg  = TwoExpoPolu_pdf ( 'B' , mass , alpah = 1 , delta = 1 , x0 = 0 , power = 3 )
    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   alpha = None     ,   ## the slope of the first exponent 
                   delta = None     ,   ## (alpha+delta) is the slope of the first exponent
                   x0    = 0        ,   ## f(x)=0 for x<x0 
                   power = 0        ,   ## degree of polynomial
                   the_phis = None  ) : 
        #
        PolyBase.__init__  ( self , name , power , xvar , the_phis )
        #                
        self.__power = power
        #
        ## the exponential slope
        #
        limits_alpha = () 
        limits_delta = ()
        limits_x0    = ()
        if self.xminmax() : 
            mn , mx      = self.xminmax()
            dm           = mx - mn 
            mmax         = max ( abs ( mn ) , abs ( mx ) )
            limits_alpha = 1.e-6 , 1.e-16 , 300. / mmax             
            limits_delta = 1.e-6 , 1.e-16 , 300. / mmax             
            limits_x0    = mn , mn - 10 * dm , mx + 10 * dm 

        self.__alpha  = self.make_var ( alpha               ,
                                  "alpha_%s"   % name ,
                                  "#alpha(%s)" % name , alpha , *limits_alpha ) 
        
        self.__delta  = self.make_var ( delta               ,
                                  "delta_%s"   % name ,
                                  "#delta(%s)" % name , delta , *limits_delta )
        
        self.__x0     = self.make_var ( x0                  ,
                                  "x0_%s"     % name  ,
                                  "x_{0}(%s)" % name  , x0    ,  *limits_x0 )

        #
        xmin , xmax = self.xminmax() 
        self.pdf  = Ostap.Models.TwoExpoPositive (
            '2expopos_%s'  % name ,
            '2expopos(%s)' % name ,
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
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   power    = 2     ,   ## degree of the polynomial
                   alpha    = None  ,   ## 
                   x0       = None  ,
                   the_phis = None  ) :
        #
        #
        PolyBase.__init__ ( self , name , power , xvar , the_phis )
        #
        self.__power      = power

        limits_alpha = ()
        limits_x0    = ()
        if self.xminmax() :
            mn, mx       = self.xminmax() 
            dx           = mx - mn 
            alpmx        = 20
            limits_alpha = -1000./dx, +1000./dx
            limits_x0    = 0.5 * ( mn + mx ) , mn - 0.2 *  dx , mx + 0.2 * dx
            
        self.__alpha  = self.make_var ( alpha               ,
                                  'alpha_%s'  % name  ,
                                  'alpha(%s)' % name  , alpha , *limits_alpha )
        
        self.__x0    = self.make_var  ( x0                  ,
                                  'x0_%s'     % name  ,
                                  'x0(%s)'    % name  , x0    , *limits_x0    )

        xmin,xmax = self.xminmax() 
        self.pdf  = Ostap.Models.PolySigmoid (
            'ps_%s'           % name ,
            'PolySigmoid(%s)' % name ,
            self.xvar            ,
            self.phi_list        ,
            xmin                 ,
            xmax                 ,
            self.alpha           ,
            self.x0              )
        
        ## save configuration 
        self.config = {
            'name'       : self.name       ,
            'xvar'       : self.xvar       ,
            'power'      : self.power      ,            
            'alpha'      : self.alpha      , 
            'x0'         : self.x0         , 
            'the_phis'   : self.phis       ,            
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
            
        PolyBase.__init__ ( self , name , spline.npars() , xvar , the_phis )
        #
        self.__spline = spline
        #
        
        self.pdf  = Ostap.Models.PositiveSpline (
            'ps_%s'              % name ,
            'PositiveSpline(%s)' % name ,
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
        PolyBase.__init__ ( self , name , spline.npars() , xvar , the_phis )
        #
        self.__spline = spline
        # 
        self.pdf  = Ostap.Models.MonotonicSpline (
            'is_%s'                % name ,
            'MonotonicSpline(%s)' % name ,
            self.xvar                     ,
            self.spline                   , 
            self.phi_list                 )

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
        PolyBase.__init__ ( self , name , spline.npars() , xvar , the_phis )
        self.__spline = spline
        #
        self.pdf  = Ostap.Models.ConvexSpline (
            'is_%s'            % name ,
            'ConvexSpline(%s)' % name ,
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
        PolyBase.__init__ ( self , name , spline.npars() , xvar , the_phis )
        self.__spline = spline
        #
        self.pdf  = Ostap.Models.ConvexOnlySpline (
            'is_%s'                % name ,
            'ConvexOnlySpline(%s)' % name ,
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
class PS2_pdf(PDF) :
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
        PDF.__init__ ( self , name , xvar )
        #
        self.__am1 = abs ( float ( m1 ) )
        self.__am2 = abs ( float ( m2 ) )

        am1 = self.m1
        am2 = self.m2
        
        if self.xminmax() :
            mn, mx =  self.xminmax()
            assert am1  + am2 < mx , 'PS2: senseless setting of edges/threshold %s,%s vs %s' % ( mn , mx , am1+am2 )   
            
        self.pdf  = Ostap.Models.PhaseSpace2 (
            'ps2_%s'          % name ,
            'PhaseSpace2(%s)' % name ,
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
#  Left edge of N-body phase space
#  @code 
#  mass  = ... ## mass variable
#  low   = 139 + 139 + 139  ## 3 pion mass
#  model = PSLeft_pdf ( 'PDF' , mass , low , 3  )
#  @endcode 
#  @see Ostap::Models::PhaseSpaceLeft
#  @see Ostap::Math::PhaseSpaceLeft
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class PSLeft_pdf(PDF) :
    """Left edge of N-body phase space function
    >>> mass  = ... ## mass variable
    >>> low   = 139 + 139 + 139  ## 3 pion mass
    >>> model = PSLeft_pdf ( 'PDF' , mass , low , 3  )   
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   N                ,   ## N 
                   left   = None    ) : 
        #
        ## initialize the base 
        PDF.__init__ ( self , name , xvar )
        #
        self.__N    = N
        self.__left = self.make_var ( left                ,
                                'left_%s'    % name ,
                                'm_left(%s)' % name ,
                                None , *self.xminmax() ) 

        if  self.xminmax() and self.left.minmax() : 
            mn  , mx  = self.xminmax()
            lmn , lmx = self.left.minmax()            
            assert lmn <= mx, "PSLeft_pdf: senseless setting of edges/thresholds: %s,%s vs %s,%s"  % (  mn, mx , lmn, lmx ) 

        self.pdf  = Ostap.Models.PhaseSpaceLeft (
            'psl_%s'             % name ,
            'PhaseSpaceLeft(%s)' % name ,
            self.xvar ,
            self.left ,
            self.N    ) 
        
        ## save configuration 
        self.config = {
            'name'       : self.name ,
            'xvar'       : self.xvar ,
            'N'          : self.N    ,            
            'left'       : self.left ,            
            }

    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar
    
    @property
    def N( self ) :
        """define N-body phase space"""
        return self.__N
    @property
    def left( self ) :
        """(left) threshold for N-body phase space"""
        return self.__left
    @left.setter
    def left ( self , value ) :
        value = float ( value )
        self.__left.setVal ( value  )
        return self.__left.getVal()

models.append ( PSLeft_pdf ) 
# =============================================================================
## @class  PSRight_pdf
#  Right edge of L-body phase space fro N-body decay 
#  @see Ostap::Models::PhaseSpaceRight
#  @see Ostap::Math::PhaseSpaceRight
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class PSRight_pdf(PDF) :
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
        PDF.__init__ ( self , name , xvar )
        #
        self.__L = L 
        self.__N = N 
        self.__right = self.make_var ( right ,
                                       'right_%s'      % name ,
                                       'm_{right}(%s)' % name ,
                                       None , *self.xminmax() )
        
        if self.xminmax() and self.right.minmax() :
            mn  , mx  = self       .xminmax()
            rmn , rmx = self.rright. minmax()            
            assert rmn <= mx, "PSRight_pdf: senseless setting of edges/thresholds: %s,%s vs %s,%s"  % (  mn, mx , rmn, rmx ) 
            
        self.pdf  = Ostap.Models.PhaseSpaceRight (
            'psr_%s'              % name ,
            'PhaseSpaceRight(%s)' % name ,
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
class PSNL_pdf(PDF) :
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
        PDF.__init__ ( self , name , xvar )
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
                                 'm_{left}(%s)'   % name , left  , *limits_left )
        
        self.__right = self.make_var ( right ,
                                 'right_%s'       % name ,
                                 'm_{right}(%s)'  % name , right , *limits_right )

        ## pdf 
        self.pdf  = Ostap.Models.PhaseSpaceNL (
            'psnl_%s'          % name ,
            'PhaseSpaceNL(%s)' % name ,
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
#  ## S-wave 
#  model = PS23L_pdf ( 'S' , mass , m_pi , m_pi , m_psi , m_Bs , 1 , 0 )
#  @endcode
#  @see Ostap::Models::PhaseSpace23L
#  @see Ostap::Math::PhaseSpace23L
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class PS23L_pdf(PDF) :
    """2-body phase space from 3-body decay with orbital momenta
    e.g. m(2pi) from Bs -> J/psi pi pi decay
    >>> mass  = ... # mass variable for 2 pions 
    >>> m_pi  = 139
    >>> m_psi = 3096
    >>> m_Bs  = 5366
    ## S-wave 
    >>> model = PS23L_pdf ( 'S' , mass , m_pi , m_pi , m_psi , m_Bs , 1 , 0 ) 
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   m1               ,   ## mass the first particle  (const)
                   m2               ,   ## mass the second particle (const)
                   m3               ,   ## mass the third particle  (const)
                   mm               ,   ## mass of the whole system (const)
                   L                ,   ## orbital momenutm between (1,2) and 3
                   l  = 0           ) : ## orbital momentum between 1 and 2
        #
        ## initialize the base 
        PDF.__init__ ( self , name , xvar )
        #
        am1 = float ( m1 )
        am2 = float ( m2 )
        am3 = float ( m3 )
        amm = float ( mm )
        
        assert 0 <= am1          , "The first mass ``m1'' must be non-negative"
        assert 0 <= am2          , "The second mass ``m2'' must be non-negative"
        assert 0 <= am3          , "The third mass ``m3'' must be non-negative"
        assert am1+am2+am3 < amm , "The total mass ``mm'' is too low"

        assert  isinstance ( L , (int,long) ) and 0 <= L < 10 , "Invalid ``L'' for the phase space function"
        assert  isinstance ( l , (int,long) ) and 0 <= l < 10 , "Invalid ``l'' for the phase space function"
        self.__m1 = am1
        self.__m2 = am2
        self.__m3 = am3
        self.__mm = amm
        self.__L  = L
        self.__l  = l
        
        self.pdf  = Ostap.Models.PhaseSpace23L (
            'ps23l_%s'          % name ,
            'PhaseSpace23L(%s)' % name ,
            self.xvar  ,
            self.m1    ,
            self.m2    ,
            self.m3    ,
            self.mm    ,
            self.L     ,
            self.l     )
        
        ## save configuration 
        self.config = {
            'name'       : self.name  ,
            'xvar'       : self.xvar  ,
            'm1'         : self.m1    ,
            'm2'         : self.m2    ,
            'm3'         : self.m3    ,
            'mm'         : self.mm    ,            
            'L'          : self.L     ,            
            'l'          : self.l     ,            
            }

    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar

    @property
    def m1 ( self ) :
        """``m1''-parameter, the mass of the first particle ``(1)'' """
        return self.__m1

    @property
    def m2 ( self ) :
        """``m2''-parameter, the mass of the second particle ``(2)''"""
        return self.__m2

    @property
    def m3 ( self ) :
        """``m3''-parameter, the mass of the third particle ``(3)''"""
        return self.__m3

    @property
    def mm ( self ) :
        """``mm''-parameter, the mass of the mother particle"""
        return self.__mm

    @property
    def l  ( self ) :
        """``l''-parameter, the orbital momentum between particles ``(1)'' and ``(2)''"""
        return self.__l

    @property
    def L  ( self ) :
        """``L''-parameter, the orbital momentum between particles ``(1,2)'' and ``3''"""
        return self.__L
        
models.append ( PS23L_pdf ) 

# =============================================================================
## create popular 1D ``background''  function
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
    
    from ostap.fitting.basic import Flat1D, Generic1D_pdf
    model = None

    if   bkg is None :         
        model = Flat1D ( name = name , xvar =  xvar )
        if kwargs : logger.warning ('make_bkg: kwargs %s are ignored' % kwargs )
        
    ## regular case: use Bkg_pdf or PolyPos_pdf as baseline background shapes 
    elif isinstance ( bkg , ( int , long ) ) :

        if   0 < bkg : model =     Bkg_pdf ( name , power =       bkg  , xvar = xvar , **kwargs )
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
    elif isinstance ( bkg , PDF ) : 
        
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
             or ( isinstance ( bkg , tuple ) and 1 < len ( tuple ) <=3 ) :
        
        model = Bkg_pdf ( name , mass = xvar , tau = bkg , power = 0 , **kwargs )
        if kwargs : logger.warning ( 'make_bkg: kwargs %s are ignored' % kwargs )

    ## exponent 
    elif bkg is math.exp : 
        
        model = Bkg_pdf ( name , mass = xvar , power = 0 , **kwargs )
        if kwargs : logger.warning ( 'make_bkg: kwargs %s are ignored' % kwargs )

    ## strings ....
    elif isinstance ( bkg , str ) :

        bkg = bkg.strip().lower()

        if   bkg in ( '' , 'const' , 'constant' , 'flat' , 'uniform' , 'p0' , 'pol0' , 'poly0' ) :            
            return make_bkg ( 0 , name , xvar   , logger = logger , **kwargs ) 
        elif bkg in ( 'e' , 'exp' , 'expo' , 'e0' , 'exp0' , 'expo0' ) : 
            model = Bkg_pdf ( name , mass = xvar , power = 0 , **kwargs )
        elif bkg in ( 'e+' , 'exp+' , 'expo+' ) : 
            model = Bkg_pdf ( name , mass = xvar , power = 0 , **kwargs )
            model.tau.setMin ( 0 ) 
        elif bkg in ( 'e-' , 'exp-' , 'expo-' ) :             
            model = Bkg_pdf ( name , mass = xvar , power = 0 , **kwargs )
            model.tau.setMax ( 0 ) 

        import re
        
        poly = re.search ( r'(poly|pol|p)(( *)|(_*))(?P<degree>\d)' , bkg , re.IGNORECASE )
        if poly :
            degree = -1 * abs ( int ( poly.group ( 'degree' ) ) ) 
            return make_bkg ( degree  , name ,  xvar , logger = logger , **kwargs  )
        
        expo = re.search ( r'(expo|exp|e)(( *)|(_*))(?P<degree>\d)' , bkg , re.IGNORECASE )
        if expo :
            degree =            int ( expo.group ( 'degree' ) ) 
            return make_bkg ( degree  , name ,  xvar , logger = logger , **kwargs  )

        incr = re.search ( r'(increasing|increase|incr|inc|i)(( *)|(_*))(?P<degree>\d)' , bkg , re.IGNORECASE )
        if incr : 
            degree = int ( incr.group ( 'degree' ) )
            bkg    = Monotonic_pdf ( name , xvar , degree , True )
            return make_bkg ( bkg , name ,  xvar , logger = logger , **kwargs  )
        
        decr = re.search ( r'(decreasing|decrease|decr|dec|d)(( *)|(_*))(?P<degree>\d)' , bkg , re.IGNORECASE )
        if decr : 
            degree = int ( decr.group ( 'degree' ) )
            bkg    = Monotonic_pdf ( name , xvar , degree , False )
            return make_bkg ( bkg , name ,  xvar , logger = logger , **kwargs  )

    if model :
        logger.debug ( 'make_bkg: created model is %s' % model ) 
        return model
    
    raise  TypeError("Wrong type of bkg object: %s/%s " % ( bkg , type ( bkg ) ) ) 


# =============================================================================
if '__main__' == __name__ : 

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )

# =============================================================================
# The END 
# =============================================================================
