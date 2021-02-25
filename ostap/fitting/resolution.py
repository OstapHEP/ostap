#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/resolution.py
#  Set of useful resolution models:
#  - single Gaussian                     (gaussian   tails)
#  - double Gaussian                     (gaussian   tails)
#  - symmetric Apollonios                (exponenial tails)
#  - Sech/hyperbolic  secant             (exponenial tails)
#  - Bukin                               (exponential or gaussian tails)
#  - symmetric double-sided Crystal Ball (power-law  tails)
#  - symmetric Student-T                 (power-law  tails)
#  - symmetric Sinh-Asinh model          (tails can be heavy or light)
#  - symmetric JohnsonSU  model          (tails can be heavy or light)
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-07-13
# =============================================================================
"""Set of useful resolution models:
- single Gaussian                     (gaussian    tails)
- double Gaussian                     (gaussian    tails)
- symmetric Apollonios                (exponential tails)
- Sech/hyperbolic  secant             (exponential tails)
- Logistic/Sech-squared               (exponential tails) 
- symmetric Bukin                     (exponential or gaussian tails)
- Symmetric double-sided Crystal Ball (power-law  tails)
- Student-T                           (power-law  tails)
- symmetric Sinh-Asinh model          (tails can be heavy or light)
- symmetric JohnsonSU  model          (tails can be heavy or light)
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'ResoGauss'     , ## simple single-Gaussian resolution model,
    'ResoGauss2'    , ## double-Gaussian resolutin model,
    'ResoApo2'      , ## symmetric Apollonios resolution model,
    'ResoCB2'       , ## symmetric double-sided Crystal Ball resolution model,
    'ResoStudentT'  , ## Student-T resolution model,
    'ResoSech'      , ## Sech/hyperbolic secant  resolution model
    'ResoLogistic'  , ## Logistic ("sech-squared") resoltuion model
    'ResoBukin'     , ## symmetric Bukin resolution model
    'ResoJohnsonSU' , ## symmetric Jonnson's SU resolution model 
    'ResoSinhAsinh' , ## symmetric Sinh-Asinh resolution model 
    )
# =============================================================================
import ROOT
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.resolution' )
else                       : logger = getLogger ( __name__                   )
# =============================================================================
from ostap.fitting.basic import Ostap, RESOLUTION
# =============================================================================    
models = set() 
# =============================================================================
## sigle gaussian model for resolution
# =============================================================================
## @class ResoGauss
#  Trivial single gaussian resolution model
class ResoGauss(RESOLUTION) :
    """Trivial single gaussian resolution model
    """
    def __init__ ( self         ,
                   name         ,   ## the  name 
                   xvar         ,   ## the variable 
                   sigma        ,   ## the first sigma
                   fudge = 1    ,   ## fudge-factor 
                   mean  = None ) : ## mean-value
        
        ## initialize the base
        super(ResoGauss,self).__init__( name  = name  ,
                                        xvar  = xvar  ,
                                        sigma = sigma ,                                     
                                        mean  = mean  ,
                                        fudge = fudge )

        #
        ## build gaussian resolution model
        #
        # self.gauss = ROOT.RooGaussModel(
        self.gauss = ROOT.RooGaussian (
            'ResoGauss_%s'  + name ,
            'ResoGauss(%s)' % name ,
            self.xvar              ,
            self.mean              , 
            self.sigma_corr        ) ## ATTENTION!
        
        self.pdf = self.gauss

        ##  save   the configuration
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mean'  : self.mean  ,
            'sigma' : self.sigma ,
            'fudge' : self.fudge ,
            }
        

models.add ( ResoGauss ) 
# =============================================================================
## @class ResoGauss2
#  Double Gaussian model for  resolution.
#  Parameters: 
#  - sigma of core Gaussian
#  - ratio of wide/core widths
#  - fraction of core(narrow) component
class ResoGauss2(RESOLUTION) :
    """Double-gaussian resolution model
    - sigma of core Gaussian
    - ratio of wide/core widths
    - fraction of core component
    """        
    def __init__ ( self            ,
                   name            ,   ## the name 
                   xvar            ,   ## the variable 
                   sigma           ,   ## the core sigma
                   scale    = 1.2  ,   ## sigma2/sigma1 ratio 
                   fraction = 0.5  ,   ## fraction of
                   fudge    = 1    ,   ## fudge-factor 
                   mean     = None ) : ## the mean value
        
        ## initialize the base 
        super(ResoGauss2,self). __init__ ( name  = name  ,
                                           xvar  = xvar  ,
                                           sigma = sigma ,
                                           mean  = mean  ,
                                           fudge = fudge )
        ## fraction of sigma1-component 
        self.__fraction = self.make_var (
            fraction                   , 
            'CoreFraction_'     + name ,
            'CoreFraction(%s)'  % name , fraction , 0 ,  1 ) 
        
        ## sigma2/sigma1 width ratio;
        self.__scale = self.make_var (
            scale ,
            'SigmaScale_'       + name ,
            'SigmaScale(%s)'    % name , scale    , 1 , 10 ) 

        #
        ## build resolution model
        # 
        self.pdf = Ostap.Models.DoubleGauss (
            "Reso2Gauss_"       + name ,
            "Reso2Gauss(%s)"    % name ,
            self.xvar       ,
            self.sigma_corr , ## ATTENTION! 
            self.fraction   ,
            self.scale      ,
            self.mean    
            )
        
        ##  save   the configuration
        self.config = {
            'name'     : self.name     ,
            'xvar'     : self.xvar     ,
            'mean'     : self.mean     ,
            'sigma'    : self.sigma    ,
            'scale'    : self.scale    ,
            'fraction' : self.fraction ,
            'fudge'    : self.fudge    ,
            }

    @property
    def fraction ( self  ) :
        """``fraction'' parameter for double Gaussian resolution function
        """
        return self.__fraction
    @fraction.setter
    def fraction ( self , value ) :
        value = float ( value )
        assert 0<= value <= 1, "``Fraction'' must be in  (0,1) range"
        self.__fraction.setVal ( value )

    @property
    def scale ( self  ) :
        """``scale'' parameter for double Gaussian resolution function
        """
        return self.__scale
    @scale.setter
    def scale ( self , value ) :
        value = float ( value )
        assert 1 < value, "``scale''-parameter must be >1"
        self.__scale.setVal ( value )
  
            
models.add ( ResoGauss2 ) 
# =============================================================================
## @class ResoApo2
#  Symmetrical  Apollonios  model for resolution
#   - Gaussian core 
#   - exponential tails
#  @see Ostap::Models::Apollonios2 
#  @see Ostap::Math::Apollonios2 
class ResoApo2(RESOLUTION) :
    """Symmetric variant of Apollonios model for the resolution function
    - Gaussian core 
    - exponential tails
    see Ostap.Models.Apollonios2 
    see Ostap.Math.Apollonios2 
    """
    def __init__ ( self         ,
                   name         ,   ## the  name 
                   xvar         ,   ## the variable 
                   sigma        ,   ## the sigma
                   beta  = 1    ,   ## beta parameter 
                   fudge = 1    ,   ## fudge-factor 
                   mean  = None ) : ## the mean value 

        ##  initlialize the base 
        super(ResoApo2,self).__init__ ( name  = name  ,
                                        xvar  = xvar  ,
                                        sigma = sigma ,
                                        mean  = mean  ,
                                        fudge = fudge )
        
        self.__beta    = self.make_var (
            beta ,
            'ResoBeta_%s'  % name  ,
            'ResoBeta(%s)' % name  , beta , 0.0001 , 10000 )

        #
        ## build resolution model
        #
        self.apo2  = Ostap.Models.Apollonios2 (
            "ResoApollonios_"    + name ,
            "ResoApollonios(%s)" % name ,
            self.xvar       ,
            self.mean       ,
            self.sigma_corr ,
            self.sigma      ,
            self.beta       ) 

        self.pdf = self.apo2

        ##  save   the configuration
        self.config = {
            'name'     : self.name     ,
            'xvar'     : self.xvar     ,
            'mean'     : self.mean     ,
            'sigma'    : self.sigma    ,
            'beta'     : self.beta     ,
            'fudge'    : self.fudge    ,
            }
        
    @property
    def beta ( self  ) :
        """``beta'' parameter for symmetric Apollonios resolution function"""
        return self.__beta
    @beta.setter
    def beta ( self , value ) :
        value = float ( value )
        assert 0< value , "``beta''-parameter must be positive"
        self.__beta.setVal ( value )

    
models.add ( ResoApo2 ) 
# =============================================================================
## @class ResoCB2
#  Symmetrical double-sided Crystal Ball model for resolution
#   - Gaussian core 
#   - power-law tails
#  @see Ostap::Math::CrystalBallDS
#  @see Ostap::Models::CrystalBallDS
class ResoCB2(RESOLUTION) :
    """Symmetric double-sided Crystal Ball model for resolution
    - Gaussian core 
    - power-law tails
    see Ostap.Math.CrystalBallDS
    see Ostap.Models.CrystalBallDS
    """
    def __init__ ( self         , 
                   name         ,   ## the  name 
                   xvar         ,   ## the  variable 
                   sigma        ,   ## core r esolution
                   alpha = 1.5  ,   ## alpha  
                   n     = 5    ,   ## power-law exponent
                   fudge = 1    ,   ## fudge-factor 
                   mean  = None ) : ## the mean value

        ## initialize the base 
        super(ResoCB2,self).__init__ ( name  = name  ,
                                       xvar  = xvar  ,
                                       sigma = sigma ,
                                       mean  = mean  ,
                                       fudge = fudge )
        
        self.__alpha = self.make_var (
            alpha                  ,
            'ResoAlpha_'    + name ,
            'ResoAlpha(%s)' % name , alpha , 0.1   , 10 )
        
        self.__n     = self.make_var (
            n                  ,
            'ResoN_'        + name ,
            'ResoN(%s)'     % name , n     , 1.e-6 , 50 )
        
        ## gaussian 
        self.cb2 = Ostap.Models.CrystalBallDS (
            'ResoCB2_'   + name ,
            'ResoCB2(%s' % name ,
            self.xvar           ,
            self.mean           , 
            self.sigma_corr     , ## ATTENTION!
            self.alpha          ,
            self.n              ,
            self.alpha          ,
            self.n              )
        
        ## the final PDF 
        self.pdf = self.cb2

        ##  save   the configuration
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'mean'     : self.mean  ,
            'sigma'    : self.sigma ,
            'alpha'    : self.alpha ,
            'n'        : self.n     ,
            'fudge'    : self.fudge ,
            }

    @property
    def alpha ( self  ) :
        """``alpha'' parameter for double-sided symmetric resolution function
        """
        return self.__alpha
    @alpha.setter
    def alpha ( self , value ) :
        value = float ( value )
        assert 0.2<= value<=5 , "``alpha''-parameter must be in 0.2,10 interval"
        self.__alpha.setVal ( value )

    @property
    def n ( self  ) :
        """``n'' parameter for double-sided symmetric resolution function"""
        return self.__n
    @n.setter
    def n ( self , value ) :
        value = float ( value )
        assert 1.e-4 <= value <= 40,  "``n'' must be in [1.e-4,40] interval"
        self.__n.setVal ( value )

models.add ( ResoCB2 )

# =============================================================================
## @class ResoStudentT
#  (symmetric) Student-T model for the resolution
#   - power-law tails 
#  @see Ostap::Models::StudentT
#  @see Ostap::Math::StudentT
#  @see http://en.wikipedia.org/wiki/Student%27s_t-distribution
class ResoStudentT(RESOLUTION) :
    """Student-T model for the resolution
    - power-law tails 
    - see http://en.wikipedia.org/wiki/Student%27s_t-distribution
    see Ostap.Models.StudentT
    see Ostap.Math.StudentT    
    """
    def __init__ ( self         ,
                   name         , ## the name 
                   xvar         , ## the variable
                   sigma        , ## the sigma
                   n            , ## N-parameter
                   fudge = 1    , ## fudge parameter 
                   mean  = None ) :
        
        ## initialize the base 
        super(ResoStudentT,self).__init__ ( name  = name  ,
                                            xvar  = xvar  ,
                                            sigma = sigma ,
                                            mean  = mean  ,
                                            fudge = fudge )
        
        self.__n     = self.make_var (
            n                      ,
            'ResoN_'        + name ,
            'ResoN(%s)'     % name , n , 1.e-6 , 100 )
        
        # 
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.StudentT (
            "ResoStT_"    + name ,
            "ResoStT(%s)" % name ,
            self.xvar       , 
            self.mean       ,
            self.sigma_corr , ## ATTENTION!
            self.n          )
        
        ##  save   the configuration
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'sigma'    : self.sigma ,
            'n'        : self.n     ,
            'mean'     : self.mean  ,
            'fudge'    : self.fudge ,
            }
        
    @property
    def n ( self  ) :
        """``n'' parameter for symmetric Student-T resolution function"""
        return self.__n
    
    @n.setter
    def n ( self , value ) :
        value = float ( value )
        assert 1.e-5 <= value , "``n'' must be positive!"
        self.__n.setVal ( value )

models.add ( ResoStudentT )


# =============================================================================
## @class ResoSech
#  Sech/hyperbolic  secant model for the resolution: exponential tails, leptokurtic
#   - Gaussian-like core
#   - exponential tails 
#  @see Ostap::Models::Sech 
#  @see Ostap::Math::Sech 
#  @see https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
class ResoSech(RESOLUTION) :
    """Sech/hyperbolic secant  for the resolution:
    - Gaussian-like core
    - exponential tails 
    - leptokurtic 
    - see https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
    see Ostap.Models.Sech 
    see Ostap.Math.Sech 
    """
    def __init__ ( self         ,
                   name         , ## the name 
                   xvar         , ## the variable
                   sigma        , ## the sigma
                   fudge = 1    , ## fudge-factor 
                   mean  = None ) :
        
        ## initialize the base 
        super(ResoSech,self).__init__ ( name  = name  ,
                                        xvar  = xvar  ,
                                        sigma = sigma ,
                                        mean  = mean  ,
                                        fudge = fudge )

        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Sech (
            "ResoSech_"    + name ,
            "ResoSech(%s)" % name ,
            self.xvar       ,
            self.mean       ,
            self.sigma_corr )  ## ATTENTION!
        
        ##  save   the configuration
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'sigma'    : self.sigma ,
            'mean'     : self.mean  ,
            'fudge'    : self.fudge
            }
        
models.add ( ResoSech )

# =============================================================================
## @class ResoBukin
#  Resolution function, described as symmetric Bukin's function
#   - Gaussian-like core
#   - exponential or Gaussian tails 
#  @see Ostap::Models::Bukin
#  @see Ostap::Math::Bukin
class ResoBukin (RESOLUTION) :
    """Resolution function, described as symmetric Bukin's function
    - Gaussian-like core
    - exponential or Gaussian tails 
    see Ostap::Models::Bukin
    see Ostap::Math::Bukin
    """
    def __init__ ( self         ,
                   name         , ## the name 
                   xvar         , ## the variable
                   sigma        , ## the sigma
                   rho   = 0    , ## the rho-parameter 
                   fudge = 1    , ## fudge-factor 
                   mean  = None ) :
        
        ## initialize the base 
        super(ResoBukin,self).__init__ ( name  = name  ,
                                         xvar  = xvar  ,
                                         sigma = sigma ,
                                         mean  = mean  ,
                                         fudge = fudge )
        
        ## parameter xi is zero! 
        self.__xi = ROOT.RooRealConstant.value ( 0 ) 
        
        ## rho 
        self.__rho = self.make_var ( rho               ,
                                     "rho_%s"   % name ,
                                     "#rho(%s)" % name , rho , 0 , 0 , 15 )        
        
        # 
        ## create PDF
        # 
        self.pdf = Ostap.Models.Bukin (
            "ResoBukin_"    + name ,
            "ResoBukin(%s)" % name ,
            self.xvar       ,
            self.mean       ,
            self.sigma_corr , ## ATTENTION!
            self.xi         ,
            self.rho        ,
            self.rho        )

        ##  save   the configuration
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'sigma'    : self.sigma ,
            'mean'     : self.mean  ,
            'rho'      : self.rho   ,
            'fudge'    : self.fudge
            }

    @property
    def xi ( self ) :
        """``xi''-parameter (asymmetry) for Bukin function"""
        return self.__xi    
    @property
    def rho ( self ) :
        """``rho''-parameter (tail) for Bukin function"""
        return self.__rho


# =============================================================================
## @class ResoJohnsonSU
#
#  Resolution model based on symmetric form of Johnson's SU distribution
# 
#  Johnson, N. L. (1949) 
#  "Systems of frequency curves generated by methods of translation"
#  Biometrika 36: 149–176 JSTOR 2332539
#  @see https://en.wikipedia.org/wiki/Johnson_SU_distribution
#
#  When variable \f$x\f$ follows Johnson-SU distribution, 
#  the variable 
#  \f$ z = \gamma + \delta \sinh^{-1}\dfrac{ x - \xi}{\lambda} \f$
#  follows normal distribtion with mean 0 and sigma 1.
#
#  @see Ostap::Math::JohnsonSU
#  @see Ostap::Models::JohnsonSU
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-07-019
class ResoJohnsonSU(RESOLUTION) :
    """
    Resolution model based on symmetric form of Johnson's SU distribution
    
    Johnson, N. L. (1949) 
    Systems of frequency curves generated by methods of translation
    Biometrika 36: 149–176 JSTOR 2332539
    
    https://en.wikipedia.org/wiki/Johnson_SU_distribution
        
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   lambd     = None ,   ## related to sigma 
                   delta     = 1    ,
                   fudge     = 1    ,
                   mean      = None ) :
        
        #
        ## initialize the base
        #
        super ( ResoJohnsonSU, self ) .__init__ ( name        = name                 ,
                                                  xvar        = xvar                 ,
                                                  sigma       = lambd                ,
                                                  mean        = mean                 ,
                                                  fudge       = fudge                ,
                                                  sigma_name  = 'lambda_%s'   % name ,
                                                  sigma_title = '#lambda(%s)' % name )
        self.__xi    = self.mean
        self.__lambd = self.sigma
        self.__gamma = ROOT.RooRealConstant.value ( 0 )
        
        self.lambd.setMax ( self.lambd.getMax() * 100 ) ## adjust it! 
    
        self.__delta = self.make_var ( delta                 ,
                                       'delta_%s'     % name ,
                                       '#delta(%s)'   % name , delta ,
                                       1 , 1.e-6 , 1000   )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.JohnsonSU (
            "ResoJSU_"          + name ,
            "ResoJohnsonSU(%s)" % name ,
            self.xvar       ,
            self.xi         ,
            ## self.lambd      ,
            self.sigma_corr , ## ATTENTION! as lambda ....
            self.delta      ,
            self.gamma      )

        ## save the configuration
        self.config = {
            'name'      : self.name    ,
            'xvar'      : self.xvar    ,
            'mean'      : self.mean    ,
            'lambd'     : self.lambd   ,
            'delta'     : self.delta   ,
            'fudge'     : self.fudge   , 
            }

    @property
    def delta ( self ) :
        """``delta''-parameter for Johnson-SU function"""
        return self.__delta
    @delta.setter
    def delta ( self, value ) :
        value = float ( value )
        assert   0 < value, "``delta''-parameter must be positive"
        self.__delta.setVal ( value ) 

    @property
    def gamma ( self ) :
        """``gamma''-parameter for Johnson-SU function"""
        return self.__gamma
        
    @property
    def xi ( self ) :
        """``xi''-parameter (location) for Johnson-SU function (the   same as ``mean'')"""
        return self.__xi
    @xi.setter
    def xi ( self, value ) :
        value = float ( value )
        self.__xi.setVal ( value ) 

    @property
    def lambd ( self ) :
        """``lambda''-parameter (scale) for Johnson-SU function (the  same  as ``sigma'')"""
        return self.__lambd
    @lambd.setter
    def lambd ( self, value ) :
        value = float ( value )
        assert   0 < value, "``lambda''-parameter must be positive"        
        self.__lambd.setVal ( value ) 

    
# =============================================================================
## @class ResoSinhAsinh
#
#  Resoltuion model based on symmetric form of "SinhAsinh" distribution
# 
#  @see Ostap::Math::SinhAsinh
#  @see Ostap::Models::SinhAsinh
#  Jones, M. C.; Pewsey, A. (2009). 
#  "Sinh-arcsinh distributions". Biometrika 96 (4): 761. 
#   doi:10.1093/biomet/asp053
#   http://oro.open.ac.uk/22510
#
#   Location & scale  parameters are the usual representation of the family of 
#   distributions:
#    - \f$\epsilon\f$ parameter control the skewness 
#    - \f$\delta\f$   parameter control the kurtosis 
#   Normal distribution reappears as \f$\epsilon=0\f$ and \f$\delta=1\f$ 
#  The heavy tails correspond to \f$\delta<1\f$, 
#  light tails correpond to \f$\delta>1\f$
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2020-12-10
class ResoSinhAsinh(RESOLUTION) :
    """ Resolution model based on symmetric form of `SinhAsinh`-function:
    
    see Jones, M. C.; Pewsey, A. (2009).
    ``Sinh-arcsinh distributions''. Biometrika 96 (4): 761. 
    doi:10.1093/biomet/asp053
    http://oro.open.ac.uk/22510
    
    Normal distribution reappears as delta=1 
    The heavy tails correspond to delta<1, light tails correpond to delta>1
    
    Parameters 
    - mean     : location
    - sigma    : scale 
    - delta>0  : parameter to control the tails/kurtosis 
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   sigma     = None ,
                   delta     = 1    ,
                   fudge     = 1    ,   ## fudge factor 
                   mean      = None ) : ## mu 

        ## initialize the base 
        super(ResoSinhAsinh,self).__init__ ( name  = name  ,
                                             xvar  = xvar  ,
                                             sigma = sigma ,
                                             mean  = mean  ,
                                             fudge = fudge )
        
        ## parameter epsilon is  fixed to zero! 
        self.__epsilon = ROOT.RooRealConstant.value ( 0 )
        self.__delta   = self.make_var ( delta ,
                                         'delta_%s'   % name ,
                                         '#delta(%s)' % name , delta ,
                                         1 , 1.e-6 , 1000   )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.SinhAsinh (
            "ResoSinhAsinhT_"   + name ,
            "ResoSinhAsinh(%s)" % name ,
            self.xvar       ,
            self.mean       ,
            self.sigma_corr , ## ATTENTION! 
            self.epsilon    ,
            self.delta      )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            'delta'     : self.delta ,
            'fudge'     : self.fudge ,
            }

    @property
    def epsilon( self ) :
        """``epsilon''-parameter for Sinh-Asinh function"""
        return self.__epsilon
    @property
    def delta ( self ) :
        """``delta-parameter'' for Sinh-Asinh function"""
        return self.__delta
    @delta.setter
    def delta ( self, value ) :
        value = float ( value )
        assert 0  < value, "``delta''-parameter must be positive"         
        self.__delta.setVal ( value ) 


# =============================================================================
## @class ResoLogistic
#  Logistic, aka "sech-square" PDF
#  \f$ f(x;\mu;s) = \dfrac{1}{4s}sech^2\left(\dfrac{x-\mu}{2s}\right)\f$, 
#   where
#   \f$  s = \sigma \dfrac{\sqrt{3}}{\pi}\f$
#  @see https://en.wikipedia.org/wiki/Logistic_distribution
#  @see Ostap::Math::Logistic
#  @see Ostap::Models::Logistic
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-06-14
class ResoLogistic(RESOLUTION) :
    """ Logistic, aka ``sech-square'' PDF
     \f$ f(x;\mu;s) = \dfrac{1}{4s}sech^2\left(\dfrac{x-\mu}{2s}\right)\f$, 
     where
     \f$  s = \sigma \dfrac{\sqrt{3}}{\pi}\f$
     - see https://en.wikipedia.org/wiki/Logistic_distribution
     - see Ostap::Math::Logistic
     - see Ostap::Models::Logistic
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   sigma     = None ,   ## related to sigma
                   fudge     = 1    , 
                   mean      = None ) : ## related to mean 
        
        
        ## initialize the base 
        super(ResoLogistic,self).__init__ ( name  = name  ,
                                            xvar  = xvar  ,
                                            sigma = sigma ,
                                            mean  = mean  ,
                                            fudge = fudge )
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Logistic (
            "logistic_"     + name ,
            "Logistic(%s)"  % name ,
            self.xvar       ,
            self.mean       ,
            self.sigma_corr ) 
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            'fudge'     : self.fudge ,
            }

        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
