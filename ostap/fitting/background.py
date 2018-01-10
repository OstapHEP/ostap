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
    'Monothonic_pdf'  , ## A positive monothonic polynomial
    'Convex_pdf'      , ## A positive polynomial with fixed sign first and second derivatives 
    'ConvexOnly_pdf'  , ## A positive polynomial with fixed sign second derivatives 
    'Sigmoid_pdf'     , ## Background: sigmoid modulated by positive polynom 
    'TwoExpoPoly_pdf' , ## difference of two exponents, modulated by positive polynomial
    ##
    'PSpline_pdf'     , ## positive            spline 
    'MSpline_pdf'     , ## positive monothonic spline 
    'CSpline_pdf'     , ## positive monothonic convex or concave spline 
    'CPSpline_pdf'    , ## positive convex or concave spline 
    ##
    'PS2_pdf'         , ## 2-body phase space (no parameters)
    'PSLeft_pdf'      , ## Low  edge of N-body phase space 
    'PSRight_pdf'     , ## High edge of L-body phase space from N-body decays  
    'PSNL_pdf'        , ## L-body phase space from N-body decays  
    'PS23L_pdf'       , ## 2-body phase space from 3-body decays with orbital momenta
    )
# =============================================================================
import ROOT, math
from   ostap.core.core     import cpp, Ostap
from   ostap.math.base     import iszero
from   ostap.fitting.basic import makeVar, PDF, Phases 
# =============================================================================
from   ostap.logger.logger     import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.models_bkg' )
else                       : logger = getLogger ( __name__             )
# =============================================================================
logger.debug ( __doc__ )
models = []
# =============================================================================
##  @class PolyBase
#   helper base class to implement various polynomial-like shapes
class PolyBase(PDF,Phases) :
    """Helper base class to implement various polynomial-like shapes
    """
    def __init__ ( self , name , power , xvar = None , the_phis = None ) :
        xvar = makeVar ( xvar , 'xvar' , 'x-variable' )
        PDF   .__init__ ( self , name  , xvar      )
        Phases.__init__ ( self , power , the_phis  )
    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar
        
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
                   mass             ,   ## the variable
                   power    = 0     ,   ## degree of polynomial
                   tau      = None  ,   ## exponential slope 
                   the_phis = None  ) : ## the phis... 
        #
        PolyBase.__init__  ( self , name , power , mass )
        #                
        self.power = power
        #
        mn,mx   = self.xminmax
        mc      = 0.5 * ( mn + mx )
        taumax  = 100
        #
        if not iszero ( mn ) : taumax =                100.0 / abs ( mn ) 
        if not iszero ( mc ) : taumax = min ( taumax , 100.0 / abs ( mc ) )
        if not iszero ( mx ) : taumax = min ( taumax , 100.0 / abs ( mx ) )
        # 
        ## the exponential slope
        #
        self.__tau  = makeVar ( tau              ,
                                "tau_%s"  % name ,
                                "tau(%s)" % name , tau , 0 , -taumax, taumax )
        #
            
        if 0 >= self.power :
            
            while self.phis :
                del self.phis[-1] 
            self.phi_list.removeAll()
            
            self.pdf      = ROOT.RooExponential (
                'exp_%s' % name  , 'exp(%s)' % name , mass , self.tau )
            
        else :
            
            self.pdf  = Ostap.Models.ExpoPositive (
                'expopos_%s'  % name ,
                'expopos(%s)' % name ,
                mass                 ,
                self.tau             ,
                self.phi_list        ,
                mass.getMin()        ,
                mass.getMax()        )

    @property
    def tau ( self ) :
        """Tau-parameter (exponential slope) for expo*pol function"""
        return self.__tau
    @tau.setter
    def tau ( self , value ) :
        value = float ( value  )
        self.__tau.setVal ( value  ) 
        return self.__tau.getVal() 
    
models.append ( Bkg_pdf ) 
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
                   mass             ,   ## the varibale 
                   phasespace       ,   ## Ostap::Math::PhaseSpaceNL 
                   power = 1        ) : ## degree of the polynomial
        
        #
        PolyBase.__init__  ( self , name , power , mass )
        #
        self.ps    = phasespace  ## Ostap::Math::PhaseSpaceNL
        self.power = power

        ## make PDF 
        self.pdf  = Ostap.Models.PhaseSpacePol (
            'pspol_%s'          % name ,
            'PhaseSpacePol(%s)' % name ,
            self.mass            ,
            self.ps              ,  ## Ostap::Math::PhaseSpaceNL 
            self.phi_list        )

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
    
    f(x) ~ ( exp(-alpha*x) - exp(-(alpha+delta)*x) *  Pol_n(x)
    where Pol_n(x) is POSITIVE polynomial (Pol_n(x)>=0 over the whole range) 
    
    >>>  mass = ROOT.RooRealVar( ... ) 
    >>>  bkg  = TwoExpoPolu_pdf ( 'B' , mass , alpah = 1 , delta = 1 , x0 = 0 , power = 3 )
    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   mass             ,   ## the variable
                   alpha = None     ,   ## the slope of the first exponent 
                   delta = None     ,   ## (alpha+delta) is the slope of the first exponent
                   x0    = 0        ,   ## f(x)=0 for x<x0 
                   power = 0        ,   ## degree of polynomial
                   tau   = None     ) : ##  
        #
        PolyBase.__init__  ( self , name , power , mass )
        #                
        self.power = power
        #
        mn,mx   = self.xminmax
        mc      = 0.5 * ( mn + mx )
        taumax  = 100
        #
        if not iszero ( mn ) : taumax =                100.0 / abs ( mn ) 
        if not iszero ( mc ) : taumax = min ( taumax , 100.0 / abs ( mc ) )
        if not iszero ( mx ) : taumax = min ( taumax , 100.0 / abs ( mx ) )
        # 
        ## the exponential slope
        #
        self.__alpha  = makeVar ( alpha               ,
                                  "alpha_%s"   % name ,
                                  "#alpha(%s)" % name , alpha , 1 , 0 , taumax )
        
        self.__delta  = makeVar ( delta               ,
                                  "delta_%s"   % name ,
                                  "#delta(%s)" % name , delta , 1 , 0 , taumax )
        
        self.__x0     = makeVar ( x0                 ,
                                  "x0_%s"     % name ,
                                  "x_{0}(%s)" % name , x0  ,
                                  mn , mn-0.5*(mx-mn) , mx+0.5*(mx-mn) ) 
        #
        self.pdf  = Ostap.Models.TwoExpoPositive (
            '2expopos_%s'  % name ,
            '2expopos(%s)' % name ,
            mass                  ,
            self.alpha            ,
            self.delta            ,
            self.x0               ,
            self.phi_list         ,
            mass.getMin()         ,
            mass.getMax()         )
        
    @property
    def alpha ( self ) :
        """Alpha-parameter (slope of leading exponnet) of the TwoExpo function"""
        return self.__alpha
    @alpha.setter
    def alpha ( self , value ) :
        value = float ( value )
        assert 0 <= value, 'Alpha-parameter must be non-negative'
        self.alpha.setVal ( value ) 
        return self.__alpha.getVal()
    
    @property
    def delta ( self ) :
        """Delta-parameter (second exponent slope is ``alpha+delta'') of the TwoExpo function"""
        return self.__delta
    @delta.setter
    def delta ( self , value ) :
        value = float ( value )
        assert 0 <= value, 'Delta-parameter must be non-negative'
        self.delta.setVal ( value ) 
        return self.__delta.getVal()
    
    @property
    def x0 ( self ) :
        """x0-parameter of the TwoExpo function  (f(x)=0 for x<x0)"""
        return self.__x0
    @x0.setter
    def x0 ( self , value ) :
        value = float ( value )
        self.x0.setVal ( value ) 
        return self.__x0.getVal()
    

models.append ( TwoExpoPoly_pdf ) 
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
                   name             ,   ## the name 
                   mass             ,   ## the varibale 
                   power = 1        ) : ## degree of the polynomial
        #
        PolyBase.__init__ ( self , name , power , mass )
        #
        self.power = power

        ## make PDF
        self.pdf  = Ostap.Models.PolyPositive (
            'pp_%s'            % name ,
            'PolyPositive(%s)' % name ,
            self.mass            ,
            self.phi_list        ,
            self.mass.getMin()   ,
            self.mass.getMax()   ) 
        
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
                   name             ,   ## the name 
                   mass             ,   ## the varibale 
                   power = 1        ) : ## (half)degree of the polynomial
        #
        PolyBase.__init__ ( self , name , power , mass )
        #
        self.power = power
        
        ## make PDF
        self.pdf  = Ostap.Models.PolyPositiveEven (
            'ppe_%s'               % name ,
            'PolyPositiveEven(%s)' % name ,
            self.mass            ,
            self.phi_list        ,
            self.mass.getMin()   ,
            self.mass.getMax()   ) 
        
models.append ( PolyPos_pdf ) 

# =============================================================================
## @class  Monothonic_pdf
#  A positive monothonic polynomial 
#  @see Ostap::Models::PolyMonothonic 
#  @see Ostap::Math::Monothonic
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Monothonic_pdf(PolyBase) :
    """Positive monothonic (Bernstein) polynomial:
    
    f(x) = Pol_n(x)
    with f(x)>= 0 over the whole range and
    the derivative f' do not change the sign
    
    >>>  mass = ROOT.RooRealVar( ... )
    
    # increasing background 
    >>>  bkg_inc  = Monothonic_pdf ( 'B1' , mass , power = 2 , increasing = True  )
    
    # decreasing background 
    >>>  bkg_dec  = Monothonic_pdf ( 'B2' , mass , power = 2 , increasing = False  )
    """
    ## constructor
    def __init__ ( self              ,
                   name              ,   ## the name 
                   mass              ,   ## the variable
                   power = 2         ,   ## degree of the polynomial
                   increasing = True ) : ## increasing or decreasing ?
        #
        PolyBase.__init__ ( self , name , power , mass )
        #
        self.power      = power
        self.increasing = increasing
        # 
        self.pdf  = Ostap.Models.PolyMonothonic (
            'pp_%s'              % name ,
            'PolyMonothonic(%s)' % name ,
            self.mass            ,
            self.phi_list        ,
            self.mass.getMin()   ,
            self.mass.getMax()   ,
            self.increasing      )
        
        
models.append ( Monothonic_pdf ) 
# =============================================================================
## @class  Convex_pdf
#  A positive polynomial with fixed signs of the first and second derivative 
#  @see Ostap::Models::PolyConvex 
#  @see Ostap::Math::Convex
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Convex_pdf(PolyBase) :
    """Positive monothonic (Bernstein) polynomial with fixed-sign second derivative:
    
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
                   mass              ,   ## the variable
                   power = 2         ,   ## degree of the polynomial
                   increasing = True ,   ## increasing or decreasing ?
                   convex     = True ) : ## convex or concave ?
        #
        PolyBase.__init__ ( self , name , power , mass )
        #
        self.power      = power
        self.increasing = increasing
        self.convex     = convex  
        ## make PDF 
        self.pdf  = Ostap.Models.PolyConvex (
            'pp_%s'          % name ,
            'PolyConvex(%s)' % name ,
            self.mass            ,
            self.phi_list        ,
            self.mass.getMin()   ,
            self.mass.getMax()   ,
            self.increasing      ,
            self.convex          )

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
    def __init__ ( self              ,
                   name              ,   ## the name 
                   mass              ,   ## the variable
                   power = 2         ,   ## degree of the polynomial
                   convex     = True ) : ## convex or concave ?
        #
        PolyBase.__init__ ( self , name , power , mass )
        #
        self.power      = power
        self.convex     = convex
        
        ## make PDF 
        self.pdf  = Ostap.Models.PolyConvexOnly (
            'pp_%s'          % name ,
            'PolyConvex(%s)' % name ,
            self.mass            ,
            self.phi_list        ,
            self.mass.getMin()   ,
            self.mass.getMax()   ,
            self.convex          )

models.append ( ConvexOnly_pdf ) 
# =============================================================================
## @class  Sigmoid_pdf
#  Sigmoid function modulated wit hpositive polynomial
#  @see Ostap::Models::PolySigmoid
#  @see Ostap::Math::Sigmoid
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Sigmoid_pdf(PolyBase) :
    """A sigmoid function modulated by positive (Bernstein) polynomial 
    f(x) = 0.5*(1+tahn(alpha*(x-x0))*Pol_n(x)
    """
    ## constructor
    def __init__ ( self              ,
                   name              ,   ## the name 
                   mass              ,   ## the variable
                   power = 2         ,   ## degree of the polynomial
                   alpha = None      ,   ##
                   x0    = None      ) :
        #
        PolyBase.__init__ ( self , name , power , mass )
        #
        self.power      = power

        xmin  = mass.getMin()
        xmax  = mass.getMax()
        dx    = xmax - xmin 
        alpmx = 2000.0/dx 
        
        self.__alpha  = makeVar ( alpha               ,
                                  'alpha_%s'  % name  ,
                                  'alpha(%s)' % name  ,
                                  alpha               , 0 , -alpmx , alpmx ) 
        
        self.__x0    = makeVar  ( x0               ,
                                'x0_%s'  % name  ,
                                'x0(%s)' % name  ,
                                x0               ,
                                0.5*(xmax+xmin)  ,
                                xmin - 0.1 * dx  ,
                                xmax + 0.1 * dx  ) 
            
        self.pdf  = Ostap.Models.PolySigmoid (
            'ps_%s'           % name ,
            'PolySigmoid(%s)' % name ,
            self.mass            ,
            self.phi_list        ,
            self.mass.getMin()   ,
            self.mass.getMax()   ,
            self.alpha           ,
            self.x0              )

    @property
    def alpha ( self ) :
        """Alpha-parameter for Sigmoid function"""
        return self.__alpha
    @alpha.setter
    def alpha ( self , value ) :
        value = float ( value )
        self.alpha.setVal ( value ) 
        return self.__alpha.getVal()
    
    @property
    def x0 ( self ) :
        """x0-parameter for Sigmoid fuction"""
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
    """A positive spline, a composion of M-splines with non-negative coefficients

    >>> mass   = ... ## the variable
    >>> order  = 3   ## spline order
    
    ## create uniform spline 
    >>> inner  = 3   ## number of inner knots between min and max 
    >>> spline = Gaudi.Math.PositiveSpline( mass.xmin() , mass.xmax() , inner , order )

    ## create non-uniform spline with
    >>> knots = std.vector('double)()
    >>> knots.push_back ( mass.xmin() )
    >>> knots.push_back ( mass.xmax() )
    >>> knots.push_back ( ... )
    >>> spline = Gaudi.Math.PositiveSpline( knots , order )

    >>> bkg = PSpline_pdf ( 'Spline' , mass , spline ) 
    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   mass             ,   ## the variable
                   spline           ) : ## the spline object Ostap::Math::PositiveSpline
        #
        PolyBase.__init__ ( self , name , spline.npars() , mass )
        #
        self.spline = spline 

        ## make PDF 
        self.pdf  = Ostap.Models.PositiveSpline (
            'ps_%s'              % name ,
            'PositiveSpline(%s)' % name ,
            self.mass            ,
            self.spline          , 
            self.phi_list        )

models.append ( PSpline_pdf ) 
# =============================================================================
## @class  MSpline_pdf
#  The special spline for non-negative monothonic function
#  @see Ostap::Models::MonothonicSpline 
#  @see Ostap::Math::MonothonicSpline 
#  @see http://en.wikipedia.org/wiki/I-spline
#  @see http://en.wikipedia.org/wiki/M-spline
#  @see http://en.wikipedia.org/wiki/B-spline
#  @see http://link.springer.com/chapter/10.1007%2F978-3-0348-7692-6_6
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class MSpline_pdf(PolyBase) :
    """A positive monothonic spline

    >>> mass   = ... ## the variable
    >>> order  = 3   ## spline order
    
    ## create uniform spline 
    >>> inner  = 3   ## number of inner knots between min and max 
    >>> spline = Gaudi.Math.MonothonicSpline( mass.xmin() , mass.xmax() , inner , order , True )

    ## create non-uniform spline with
    >>> knots = std.vector('double)()
    >>> knots.push_back ( mass.xmin() )
    >>> knots.push_back ( mass.xmax() )
    >>> knots.push_back ( ... )
    >>> spline = Gaudi.Math.MonothonicSpline( knots , order , True )

    >>> bkg = MSpline_pdf ( 'Spline' , mass , spline ) 
    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   mass             ,   ## the variable
                   spline           ) : ## the spline object Ostap::Math::MonothonicSpline
        #
        PolyBase .__init__ ( self , name  , spline.npars() , mass )
        #
        self.spline = spline         
        # make PDF
        self.pdf  = Ostap.Models.MonothonicSpline (
            'is_%s'                % name ,
            'MonothonicSpline(%s)' % name ,
            self.mass                     ,
            self.spline                   , 
            self.phi_list                 )

models.append ( MSpline_pdf )

# =============================================================================
## @class  CSpline_pdf
#  The special spline for non-negative monothonic convex/concave function
#  @see Ostap::Models::ConvexSpline 
#  @see Ostap::Math::ConvexSpline 
#  @see http://en.wikipedia.org/wiki/I-spline
#  @see http://en.wikipedia.org/wiki/M-spline
#  @see http://en.wikipedia.org/wiki/B-spline
#  @see http://link.springer.com/chapter/10.1007%2F978-3-0348-7692-6_6
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class CSpline_pdf(PolyBase) :
    """A positive monothonic convex/concave spline

    >>> mass   = ... ## the variable
    >>> order  = 3   ## spline order
    
    ## create uniform spline 
    >>> inner  = 3   ## number of inner knots between min and max 
    >>> spline = Gaudi.Math.ConvexSpline( mass.xmin() , mass.xmax() , inner , order , True , True )

    ## create non-uniform spline with
    >>> knots = std.vector('double)()
    >>> knots.push_back ( mass.xmin() )
    >>> knots.push_back ( mass.xmax() )
    >>> knots.push_back ( ... )
    >>> spline = Gaudi.Math.ConvexSpline( knots , order, True  , True )

    >>> bkg = MSpline_pdf ( 'Spline' , mass , spline ) 
    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   mass             ,   ## the variable
                   spline           ) : ## the spline object Ostap::Math::ConvexSpline
        #
        PolyBase   .__init__ ( self , name  , spline.npars() , mass )
        #
        self.spline = spline 
        
        # make PDF 
        self.pdf  = Ostap.Models.ConvexSpline (
            'is_%s'            % name ,
            'ConvexSpline(%s)' % name ,
            self.mass                ,
            self.spline              , 
            self.phi_list            )

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
    >>> spline = Gaudi.Math.ConvexOnlySpline( mass.xmin() , mass.xmax() , inner , order , True )

    ## create non-uniform spline with
    >>> knots = std.vector('double)()
    >>> knots.push_back ( mass.xmin() )
    >>> knots.push_back ( mass.xmax() )
    >>> knots.push_back ( ... )
    >>> spline = Gaudi.Math.ConvexonlySpline( knots , order, True )

    >>> bkg = CPSpline_pdf ( 'Spline' , mass , spline ) 
    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   mass             ,   ## the variable
                   spline           ) : ## the spline object Ostap::Math::ConvexOnlySpline
        #
        PolyBase.__init__ ( self , name , spline.npars () , mass )
        #
        self.spline = spline 
        
        # make PDF 
        self.pdf  = Ostap.Models.ConvexOnlySpline (
            'is_%s'                % name ,
            'ConvexOnlySpline(%s)' % name ,
            self.mass                ,
            self.spline              , 
            self.phi_list            )

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
                   mass             ,   ## the varibale
                   m1               ,   ## the first  mass (constant)
                   m2               ) : ## the second mass (constant)
        #
        ## initialize the base 
        mass = makeVar ( mass , 'bmass' , 'background-mass' ) 
        PDF.__init__ ( self , name , mass )
        #

        if self.mass.getMax() < abs ( m1 ) + abs ( m2 ) :
            logger.error('PS2_pdf(%s): senseless setting of edges/threshold' % self.name ) 

        self.pdf  = Ostap.Models.PhaseSpace2 (
            'ps2_%s'          % name ,
            'PhaseSpace2(%s)' % name ,
            self.mass            ,
            m1  , m2 )
        
    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar
    

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
                   mass             ,   ## the variable
                   N                ,   ## N 
                   left   = None    ) : 
        #
        ## initialize the base 
        mass = makeVar ( mass , 'bmass' , 'background-mass' ) 
        PDF.__init__ ( self , name , mass )
        #
        self.__left = makeVar ( left                ,
                              'left_%s'    % name ,
                              'm_left(%s)' % name ,
                              None , mass.getMin() , mass.getMax() )

        if self.left.getMin() >= self.mass.getMax() :
            logger.error('PSLeft_pdf(%s): invalid setting!' % name )
            
        self.pdf  = Ostap.Models.PhaseSpaceLeft (
            'psl_%s'             % name ,
            'PhaseSpaceLeft(%s)' % name ,
            self.mass ,
            self.left ,
            N         ) 

    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar
    
    @property
    def left( self ) :
        """(Left) threshold for N-body phase space"""
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
                   mass             ,   ## the variable
                   L                ,   ## L  
                   N                ,   ## N
                   right   = None   ) : 
        #
        ## initialize the base 
        mass = makeVar ( mass , 'bmass' , 'background-mass' ) 
        PDF.__init__ ( self , name , mass )
        #
        self.__right = makeVar ( right ,
                                 'right_%s'      % name ,
                                 'm_{right}(%s)' % name ,
                                 None , mass.getMin() , mass.getMax() )
        
        if self.right.getMax() <= self.mass.getMax() :
            logger.error('PSRight_pdf(%s): invalid setting!' % name )
            
        self.pdf  = Ostap.Models.PhaseSpaceRight (
            'psr_%s'              % name ,
            'PhaseSpaceRight(%s)' % name ,
            self.mass  ,
            self.right ,
            L          , 
            N          )
        
    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar
    
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
                   mass             ,   ## the variable 
                   L                ,   ## L  
                   N                ,   ## N
                   left  = None     , 
                   right = None     ) : 
        ##
        #
        ## initialize the base 
        mass  = makeVar( mass , 'bmass' , 'background-mass' ) 
        PDF.__init__ ( self , name , mass )
        #
        #
        mmin = mass.getMin()
        mmax = mass.getMax()
        #
        self.__left  = makeVar ( left ,
                               'left_%s'        % name ,
                               'm_{left}(%s)'   % name , left  , 
                               0.9 * mmin + 0.1 * mmax ,
                               mmin ,
                               mmax ) 
        
        self.__right = makeVar ( right ,
                               'right_%s'       % name ,
                               'm_{right}(%s)'  % name , right , 
                               0.1 * mmin + 0.9 * mmax ,
                               mmin ,
                               mmax ) 
        
        if self.left.getMin()  >= self.mass.getMax() :
            logger.error('PSNL_pdf(%s): invalid setting!' % name )
            
        if self.right.getMax() <= self.mass.getMax() :
            logger.error('PSNL_pdf(%):  invalid setting!' % name )

        self.pdf  = Ostap.Models.PhaseSpaceNL (
            'psnl_%s'          % name ,
            'PhaseSpaceNL(%s)' % name ,
            self.mass  ,
            self.left  ,
            self.right ,
            L          , 
            N          )
        
    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar
    
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
                   mass             ,   ## the variable
                   m1               ,   ## mass the first particle  (const)
                   m2               ,   ## mass the second particle (const)
                   m3               ,   ## mass the third particle  (const)
                   m                ,   ## mass of the whole system (const)
                   L                ,   ## orbital momenutm between (1,2) and 3
                   l  = 0           ) : ## orbital momentum between 1 and 2
        #
        ## initialize the base 
        mass = makeVar( mass , 'bmass' , 'background-mass' ) 
        PDF.__init__ ( self , name , mass )
        #
        self.pdf  = Ostap.Models.PhaseSpace23L (
            'ps23l_%s'          % name ,
            'PhaseSpace23L(%s)' % name ,
            self.mass                  ,
            m1 , m2 , m3 , m , L , l )
        
    @property 
    def mass ( self ) :
        """``mass''-variable for the fit (the same as ``x'' or ``xvar'')"""
        return self.xvar
    

models.append ( PS23L_pdf ) 



# =============================================================================
if '__main__' == __name__ : 

    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )

# =============================================================================
# The END 
# =============================================================================
