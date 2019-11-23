#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/resolution.py
#  Set of useful resolution models:
#  - single Gaussian                     (gaussian   tails)
#  - double Gaussian                     (gaussian   tails)
#  - symmetric Apolonious                (exponenial tails)
#  - Sech/hyperbolic  secant             (exponenial tails)
#  - symmetric double-sided Crystal Ball (power-law  tails)
#  - symmetric Student-T                 (power-law  tails)
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2017-07-13
# =============================================================================
"""Set of useful resolution models:
- single Gaussian                     (gaussian   tails)
- double Gaussian                     (gaussian   tails)
- symmetric Apolonious                (exponenial tails)
- Sech/hyperbolic  secant             (exponenial tails)
- symmetric double-sided Crystal Ball (power-law  tails)
- Student-T                           (power-law  tails)
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'ResoGauss'     , ## simple single-Gaussian resolution model,
    'ResoGauss2'    , ## double-Gaussian resolutin model,
    'ResoApo2'      , ## symmetric Apolonios resolution model,
    'ResoCB2'       , ## symmetric double-sided Crystal Ball resolution model,
    'ResoStudentT'  , ## Student-T resolution model,
    'ResoSech'      , ## Sech/hyperbolic secant  resolution model 
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
        
        if mean is None : mean = ROOT.RooConstVar(
            'mean_%s'  % name ,
            'mean(%s)' % name , 0 )
        
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
#  Double Gaussian model for  resoltuion
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
        
        if mean is None : mean = ROOT.RooConstVar(
            'mean_%s'  % name ,
            'mean(%s)' % name , 0 )                 
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
#  Symmetrical  Apolonios  model for resolution
class ResoApo2(RESOLUTION) :
    """Symmetric variant of Apolonios model for the resolution function
    """
    def __init__ ( self         ,
                   name         ,   ## the  name 
                   xvar         ,   ## the variable 
                   sigma        ,   ## the sigma
                   beta  = 1    ,   ## beta parameter 
                   fudge = 1    ,   ## fudge-factor 
                   mean  = None ) : ## the mean value 

        if mean is None :  mean = ROOT.RooConstVar(
            'mean_%s'  % name ,
            'mean(%s)' % name , 0 )                 
            
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
        self.apo2  = Ostap.Models.Apolonios2 (
            "ResoApolonios_"    + name ,
            "ResoApolonios(%s)" % name ,
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
        """``beta'' parameter for Apolonious resolution function"""
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
class ResoCB2(RESOLUTION) :
    """Symmetric double-sided Crystal Ball model for resolution
    """
    def __init__ ( self         , 
                   name         ,   ## the  name 
                   xvar         ,   ## the  variable 
                   sigma        ,   ## core r esolution
                   alpha = 1.5  ,   ## alpha  
                   n     = 5    ,   ## power-law exponent
                   fudge = 1    ,   ## fudge-factor 
                   mean  = None ) : ## the mean value

        if mean is None : mean = ROOT.RooConstVar(
            'mean_%s'  % name ,
            'mean(%s)' % name , 0 )                 
        
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
#  @see Ostap::Models::StudentT
#  @see Ostap::Math::StudentT
#  @see http://en.wikipedia.org/wiki/Student%27s_t-distribution
class ResoStudentT(RESOLUTION) :
    """Student-T model for the resolution
    - see http://en.wikipedia.org/wiki/Student%27s_t-distribution
    """
    def __init__ ( self         ,
                   name         , ## the name 
                   xvar         , ## the variable
                   sigma        , ## the sigma
                   n            , ## N-parameter
                   fudge = 1    , ## fudge parameter 
                   mean  = None ) :
        
        if mean is None : mean = ROOT.RooConstVar(
            'mean_%s'  % name ,
            'mean(%s)' % name , 0 )                 
        
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
#  @see Ostap::Models::Sech 
#  @see Ostap::Math::Sech 
#  @see https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
class ResoSech(RESOLUTION) :
    """Sech/hyperbolic secant  for the resolution:
    - exponential tails
    - leptokurtic 
    - see https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
    """
    def __init__ ( self         ,
                   name         , ## the name 
                   xvar         , ## the variable
                   sigma        , ## the sigma
                   fudge = 1    , ## fudge-factor 
                   mean  = None ) :
        
        if mean is None : mean = ROOT.RooConstVar(
            'mean_%s'  % name ,
            'mean(%s)' % name , 0 )                 
        
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
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
# The END 
# =============================================================================
