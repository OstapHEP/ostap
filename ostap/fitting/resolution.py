#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/resolution.py
#  Set of useful resolution models:
#  - single Gaussian                     (gaussian   tails)
#  - double Gaussian                     (gaussian   tails)
#  - Apollonios-2                        (exponenial tails)
#  - Sech/hyperbolic  secant             (exponenial tails)
#  - Bukin                               (exponential or gaussian tails)
#  - double-sided Crystal Ball           (power-law  tails)
#  - Student-T                           (power-law  tails)
#  - Sinh-Asinh model                    (tails can be heavy or light)
#  - JohnsonSU  model                    (tails can be heavy or light)
#  - Hyperbolic model                    (tails are exponential)
#  - generalized Hyperbolic model        (tails are exponential or heavier)
#  - Hypatia model                       (tails are exponential or heavier)
#  - Das model                           (gaussian with exponential tails)
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-07-13
# =============================================================================
"""Set of useful resolution models:
- single Gaussian                     (gaussian    tails)
- double Gaussian                     (gaussian    tails)
- Apollonios-2                        (exponential tails)
- Sech/hyperbolic  secant             (exponential tails)
- Logistic/Sech-squared               (exponential tails) 
- Bukin                               (exponential or gaussian tails)
- double-sided Crystal Ball           (power-law  tails)
- Student-T                           (power-law  tails)
- Sinh-Asinh model                    (tails can be heavy or light)
- JohnsonSU  model                    (tails can be heavy or light)
- Hyperbolic                          (tails are exponential)
- generalized Hyperbolic              (tails are exponential or heavier)
- Hypatia model                       (tails are exponential or heavier)
- Das model                           (gaussian with exponential tails)
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    ##
    'ResoGauss'         , ## single-Gaussian resolution model,
    'ResoGauss2'        , ## double-Gaussian resolutin model,
    'ResoApo2'          , ## Apollonios-2 resolution model,
    'ResoCB2'           , ## double-sided Crystal Ball resolution model,
    'ResoStudentT'      , ## Student-T resolution model,
    'ResoSech'          , ## Sech/hyperbolic secant  resolution model
    'ResoLogistic'      , ## Logistic ("sech-squared") resoltuion model
    'ResoBukin'         , ## Bukin resolution model
    'ResoJohnsonSU'     , ## Jonnson's SU resolution model 
    'ResoSinhAsinh'     , ## Sinh-Asinh resolution model
    'ResoHyperbolic'    , ## Hyperbolic resolution model
    'ResoGenHyperbolic' , ## Generalised Hyperbolic resolution model
    'ResoHypatia'       , ## Hypatia resoltuion model
    'ResoDas'           , ## Das resolution model     
    )
# =============================================================================
import ROOT
# =============================================================================
from ostap.fitting.basic import Ostap, RESOLUTION, CheckMean 
# =============================================================================    
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.resolution' )
else                       : logger = getLogger ( __name__                   )
# =============================================================================
## @var ZERO : zero constant 
ZERO   = ROOT.RooRealConstant.value ( 0 ) 
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
                   mean  = None ,   ## mean-value
                   kappa = None ) : ## asymmetry
        
        ## initialize the base
        super(ResoGauss,self).__init__( name  = name  ,
                                        xvar  = xvar  ,
                                        sigma = sigma ,                                     
                                        mean  = mean  ,
                                        fudge = fudge )

        ## asymmetry parameter 
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa , 
                                       'kappa_%s'   % self.name ,
                                       '#kappa(%s)' % self.name ,
                                       ZERO if kappa is None else kappa , 0 , -1 , +1 ) 

        if kappa is None :

            self.__sigmaL = self.sigma_corr
            self.__sigmaR = self.sigma_corr

        else :
            
            self.__sigmaL , self.__sigmaR = self.vars_from_asymmetry (
                self.sigma_corr                                   , ## mean/average sigma
                self.kappa                                        , ## asymmetry parametet
                v1name  =  self.roo_name ( 'sigmaL' , self.name ) ,
                v2name  =  self.roo_name ( 'sigmaR' , self.name ) ,
                v1title = '#sigma_L: #sigma #times (1+#kappa)'    , 
                v2title = '#sigma_R: #sigma #times (1-#kappa)'    )
                #
        ## build gaussian resolution model
        #
        # self.gauss = ROOT.RooGaussModel(
        if kappa is None :
            
            self.gauss = ROOT.RooGaussian (
                self.roo_name ( 'rgauss_' )       ,
                "Resolution Gauss %s" % self.name ,
                self.xvar              ,
                self.mean              , 
                self.sigma_corr        ) ## ATTENTION!

        else :
            
            self.gauss = Ostap.Models.BifurcatedGauss ( 
                self.roo_name ( 'rbfgauss_' )                , 
                "Resolution Bifurcated Gauss %s" % self.name ,
                self.xvar   ,
                self.mean   , 
                self.sigmaL ,
                self.sigmaR )   
                
        self.pdf = self.gauss

        ##  save   the configuration
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mean'  : self.mean  ,
            'sigma' : self.sigma ,
            'fudge' : self.fudge ,
            'kappa' : kappa if kappa is None else self.kappa 
            }

    @property
    def kappa ( self ) :
        """``kappa'' : asymmetry parameter"""
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.set_value ( self.__kappa , value )

    @property
    def sigmaL ( self ) :
        """``sigmaL'': left sigma"""
        return self.__sigmaL
    @property
    def sigmaR ( self ) :
        """``sigmaR'': left sigma"""
        return self.__sigmaR
        
        
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
        ## fraction of sigma-1-component 
        self.__fraction = self.make_var (
            fraction                   , 
            'CoreFraction_'     + name ,
            'CoreFraction(%s)'  % name , fraction , 0 ,  1 ) 
        
        ## sigma-2/sigma-1 width ratio;
        self.__scale = self.make_var (
            scale ,
            'SigmaScale_'       + name ,
            'SigmaScale(%s)'    % name , scale    , 1 , 10 ) 

        #
        ## build the resolution model
        # 
        self.pdf = Ostap.Models.DoubleGauss (
            self.roo_name ( 'rgauss2_' )       ,
            "Resolution double Gauss %s" % self.name ,
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
                   mean  = None ,   ## the mean value 
                   kappa = None ) : ## asymmetry parameter
        
        ##  initlialize the base 
        super(ResoApo2,self).__init__ ( name  = name  ,
                                        xvar  = xvar  ,
                                        sigma = sigma ,
                                        mean  = mean  ,
                                        fudge = fudge )

        if kappa is None :
            
            self.__kappa  = ZERO
            self.__sigmaL = self.sigma_corr
            self.__sigmaR = self.sigma_corr
            
        else :

            self.__kappa = self.make_var ( kappa        ,
                                           'kappa_%s'   % self.name ,
                                           '#kappa(%s)' % self.name , kappa , 0 , -1 , +1 ) 
            
            self.__sigmaL , self.__sigmaR = self.vars_from_asymmetry (
                self.sigma_corr                                   , ## mean/average sigma
                self.kappa                                        , ## asymmetry parametet
                v1name  =  self.roo_name ( 'sigmaL' , self.name ) ,
                v2name  =  self.roo_name ( 'sigmaR' , self.name ) ,
                v1title = '#sigma_L: #sigma #times (1+#kappa)'    , 
                v2title = '#sigma_R: #sigma #times (1-#kappa)'    )
                
        self.__beta    = self.make_var (
            beta ,
            'beta_%s'   % name  ,
            '#beta(%s)' % name  , beta , 0.0001 , 10000 )

        #
        ## build the resolution model
        #
        self.apo2  = Ostap.Models.Apollonios2 (
            self.roo_name ( 'rapo2_' )       ,
            "Resolution Apollonios2 %s" % self.name ,
            self.xvar       ,
            self.mean       ,
            self.sigmaL     ,
            self.sigmaR     ,
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
            'kappa'    : None if kappa is None else self.kappa 
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

    @property
    def kappa ( self ) :
        """``kappa'' : asymmetry parameter
        """
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.set_value ( self.__kappa , value )
        
    @property
    def sigmaL ( self )  :
        """``sigmaL'' : left sigma-parameter"""
        return self.__sigmaL
    @property
    def sigmaR ( self )  :
        """``sigmaR'' : right  sigma-parameter"""
        return self.__sigmaR
    
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
    def __init__ ( self          , 
                   name          ,   ## the  name 
                   xvar          ,   ## the  variable 
                   sigma         ,   ## core resolution
                   alpha  = 1.5  ,   ## alpha  
                   n      = 5    ,   ## power-law exponent
                   fudge  = 1    ,   ## fudge-factor 
                   mean   = None ,   ## the mean value
                   kappaN = None ,   ## asymmetru for N
                   kappaA = None ) : ## asymmetru for alpha
                   

        ## initialize the base 
        super(ResoCB2,self).__init__ ( name  = name  ,
                                       xvar  = xvar  ,
                                       sigma = sigma ,
                                       mean  = mean  ,
                                       fudge = fudge )
        
        self.__alpha = self.make_var (
            alpha                  ,
            'alpha_'     + name ,
            '#alpha(%s)' % name , alpha , 0.05   , 10 )
        
        self.__n     = self.make_var (
            n                  ,
            'n_'         + name ,
            'n(%s)'      % name , n     , 1.e-6 , 50 )
        
        if kappaN is None :

            self.__kappaN = ZERO
            self.__nL     = self.n
            self.__nR     = self.n
            
        else :
            
            self.__kappaN = self.make_var ( kappaN                     ,
                                            'kappa_n%s'    % self.name ,
                                            '#kappa_{n}(%s)' % self.name , kappaN , 0 , -1 , +1 ) 
            
            self.__nL , self.__nR = self.vars_from_asymmetry (
                self.n                                        , ## mean/average n
                self.kappaN                                   , ## asymmetry parametet
                v1name  =  self.roo_name ( 'nL' , self.name ) ,
                v2name  =  self.roo_name ( 'nR' , self.name ) ,
                v1title = 'n_{L}: n #times (1+#kappa_{n})'    , 
                v2title = 'n_{R}: n #times (1-#kappa_{n})'    )

        if kappaA is None :

            self.__kappaA = ZERO
            self.__alphaL = self.alpha
            self.__alphaR = self.alpha

        else :
            
            self.__kappaA = self.make_var ( kappaA                            ,
                                            'kappa_a%s'           % self.name ,
                                            '#kappa_{#alpha}(%s)' % self.name , kappaA , 0 , -1 , +1 ) 
            
            
            self.__alphaL , self.__alphaR = self.vars_from_asymmetry (
                self.alpha                                                , ## mean/average alpha 
                self.kappaA                                               , ## asymmetry parametet
                v1name  =  self.roo_name ( 'alphaL' , self.name )         ,
                v2name  =  self.roo_name ( 'alphaR' , self.name )         ,
                v1title = '#alpha_{L}: #alpha #times (1+#kappa_{#alpha})' , 
                v2title = '#alpha_{R}: #alpha #times (1-#kappa_{#alpha})' )
        
        ## actual PDF 
        self.cb2 = Ostap.Models.CrystalBallDS (
            self.roo_name ( 'rcb2_' )       ,
            "Resolution double-sided Crystal Ball %s" % self.name ,
            self.xvar           ,
            self.mean           , 
            self.sigma_corr     , ## ATTENTION!
            self.alphaL         ,
            self.nL             ,
            self.alphaR         ,
            self.nR             )
        
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
            'kappaN'   : kappaN if kappaN is None else self.kappaN ,
            'kappaA'   : kappaA if kappaA is None else self.kappaA ,            
            }

    @property
    def alpha ( self  ) :
        """``alpha'' parameter for double-sided symmetric resolution function
        """
        return self.__alpha
    @alpha.setter
    def alpha ( self , value ) :
        value = float ( value )
        assert 0.1<= value<=6 , "``alpha''-parameter must be in [0.1,6] interval"
        self.set_value ( self.__alpha , value  ) 

    @property
    def n ( self  ) :
        """``n'' parameter for double-sided symmetric resolution function"""
        return self.__n
    @n.setter
    def n ( self , value ) :
        value = float ( value )
        assert 1.e-6 <= value <= 100,  "``n'' must be in [1.e-6,100] interval"
        self.set_value ( self.__n , value  ) 

    @property
    def kappaN ( self ) :
        """``kappaN'' : asymmetry for parameter ``n''
        """
        return self.__kappaN
    @kappaN.setter
    def kappaN ( self , value ) :
        self.set_value ( self.__kappaN , value )

    @property
    def kappaA ( self ) :
        """``kappaA'' : asymmetry for parameter ``alpha''
        """
        return self.__kappaA
    @kappaA.setter
    def kappaA ( self , value ) :
        self.set_value ( self.__kappaA , value )

    @property
    def nL ( self ) :
        """``nL'' : parameter ``n'' for left tail
        """
        return self.__nL
    @property
    def nR ( self ) :
        """``nR'' : parameter ``n'' for right tail
        """
        return self.__nR

    @property
    def alphaL ( self ) :
        """``alphaL'' : parameter ``alpha'' for left tail
        """
        return self.__alphaL

    @property
    def alphaR ( self ) :
        """``alphaR'' : parameter ``alpha'' for right tail
        """
        return self.__alphaR

    

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
    - when asymmetry is activates use `BifurcatedStudentT`
    see Ostap.Models.BigurcatedStudentT
    see Ostap.Math.BofurcatedStudentT    
    """
    def __init__ ( self           ,
                   name           ,   ## the name 
                   xvar           ,   ## the variable
                   sigma          ,   ## the sigma
                   n              ,   ## N-parameter
                   fudge  = 1     ,   ## fudge parameter 
                   mean   = None  ,   ## mean 
                   kappaS = None  ,   ## asymmetry for sigma 
                   kappaN = None  ) : ## asymmetry for N 
        
        ## initialize the base 
        super(ResoStudentT,self).__init__ ( name  = name  ,
                                            xvar  = xvar  ,
                                            sigma = sigma ,
                                            mean  = mean  ,
                                            fudge = fudge )
        
        self.__n     = self.make_var ( n                      ,
                                       'ResoN_'        + name ,
                                       'ResoN(%s)'     % name , n , 1 , 1.e-6 , 100 )

        if kappaN is None :

            self.__kappaN = ZERO 
            self.__nL     = self.n
            slef.__nR     = self.n

        else :
            
            self.__kappaN = self.make_var ( ZERO if kappaN is None else kappaN ,
                                            "kappaN_%s"      % name ,
                                            "#kappa_{n}(%s)" % name ,
                                            ZERO if kappaN is None else kappaN ,
                                            0 , -1 , 1 )
            
            self.__nL , self.__nR = self.vars_from_asymmetry (
                self.n                                        , ## mean/average n
                self.kappaN                                   , ## asymmetry parametet
                v1name  =  self.roo_name ( 'nL' , self.name ) ,
                v2name  =  self.roo_name ( 'nR' , self.name ) ,
                v1title = 'n_{L}: n #times (1+#kappa_{n})'    , 
                v2title = 'n_{R}: n #times (1-#kappa_{n})'    )

        if kappaS is None :

            self.__kappaS = ZERO 
            self.__sigmaL = self.sigma_corr
            self.__sigmaR = self.sigma_corr 
            
        else:
            
            self.__sigmaL , self.__sigmaR = self.vars_from_asymmetry (
                self.sigma_corr                                   , ## mean/average sigma
                self.kappa                                        , ## asymmetry parametet
                v1name  =  self.roo_name ( 'sigmaL' , self.name ) ,
                v2name  =  self.roo_name ( 'sigmaR' , self.name ) ,
                v1title = '#sigma_L: #sigma #times (1+#kappa)'    , 
                v2title = '#sigma_R: #sigma #times (1-#kappa)'    )
                
                 
        # 
        ## finally build pdf
        #
        if kappaL is None and kappaS is None :
            
            self.pdf = Ostap.Models.StudentT (
                self.roo_name ( 'rstt_' )       ,
                "Resolution Student's t %s" % self.name ,
                self.xvar       , 
                self.mean       ,
                self.sigma_corr , ## ATTENTION!
                self.n          )
            
        else :
            
            self.pdf = Ostap.Models.BifurcatedStudentT (
                self.roo_name ( 'rbfstt_' )       ,
                "Resolution Bifurcated Student's t %s" % self.name ,
                self.xvar       , 
                self.mean       ,
                self.sigmaL     , 
                self.sigmaR     , 
                self.nL         ,
                self.nR         )

        
        ##  save   the configuration
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'sigma'    : self.sigma ,
            'n'        : self.n     ,
            'mean'     : self.mean  ,
            'fudge'    : self.fudge ,
            'kappaS'   : kappaS if kappaS is None else self.kappaS ,  
            'kappaN'   : kappaN if kappaN is None else self.kappaN ,
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

    @property
    def kappaN ( self  ):
        """``kappaN'' : asymmetry for ``n''-parameter"""
        return self.__kappaN
    @kappaN.setter
    def kappaN ( self , value ) :
        self.set_value ( self.__kappaN , value )

    @property
    def kappaS ( self  ):
        """``kappaS'' : asymmetry for ``sigma''-parameter"""
        return self.__kappaS
    @kappaS.setter
    def kappaS ( self , value ) :
        self.set_value ( self.__kappaS , value )

    @property
    def nL ( self ) :
        """``nL'' : ``n''parameter for left  part"""
        return self.__nL        
    @property
    def nR ( self ) :
        """``nR'' : ``n''parameter for right part"""
        return self.__nR

    @property
    def sigmaL ( self ) :
        """``sigmaL'' : ``sigma''parameter for left  part"""
        return self.__sigmaL        
    @property
    def sigmaR ( self ) :
        """``sigmaR'' : ``sigma''parameter for right part"""
        return self.__sigmaR        
    
    
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
            self.roo_name ( 'rsech_' )       ,
            "Resolution Sech %s" % self.name ,
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
    def __init__ ( self            ,
                   name            ,   ## the name 
                   xvar            ,   ## the variable
                   sigma           ,   ## the sigma
                   rho      = 0    ,   ## the rho-parameter 
                   fudge    = 1    ,   ## fudge-factor 
                   mean     = None ,   ## mean 
                   xi       = None ,   ## core asymmetry parameter 
                   kappa    = None ) : ## tail asymmetry parameter 
        
        ## initialize the base 
        super(ResoBukin,self).__init__ ( name  = name  ,
                                         xvar  = xvar  ,
                                         sigma = sigma ,
                                         mean  = mean  ,
                                         fudge = fudge )
        
        ## parameter xi is zero! 
        self.__xi = self.make_var ( ZERO if xi is None else xi ,
                                    "xi_%s"   % name ,
                                    "#xi(%s)" % name ,
                                    ZERO if xi is None else xi , 0 , -10 , +10 )
        
        ## rho 
        self.__rho = self.make_var   ( rho               ,
                                       "rho_%s"   % name ,
                                       "#rho(%s)" % name , rho , 0 , 0 , 25 )        
        
        ## parameter 
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa ,
                                       "kappa_%s"          % name ,
                                       "#kappa_{#rho}(%s)" % name ,
                                       ZERO if kappa is None else kappa , 0 , -1 , +1 )
        
        if kappa is None :

            self.__rhoL = self.__rho
            self.__rhoR = self.__rho

        else :
            
            self.__rhoL , self.__rhoR = self.vars_from_asymmetry (
                self.rho                                            , ## mean/average rho 
                self.kappa                                          , ## asymmetry parameter
                v1name  =  self.roo_name ( 'rhoL' , self.name )     ,
                v2name  =  self.roo_name ( 'rhoR' , self.name )     ,
                v1title = '#rho_{L}: #rho #times (1+#kappa_{#rho})' , 
                v2title = '#rho_{R}: #rho #times (1-#kappa_{$rho})' )
                
        # 
        ## create PDF
        # 
        self.pdf = Ostap.Models.Bukin (
            self.roo_name ( 'rbukin_' )       ,
            "Resolution Bukin %s" % self.name ,
            self.xvar       ,
            self.mean       ,
            self.sigma_corr , ## ATTENTION!
            self.xi         ,
            self.rhoL       ,
            self.rhoR       )

        ##  save   the configuration
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'sigma'    : self.sigma ,
            'mean'     : self.mean  ,
            'rho'      : self.rho   ,
            'fudge'    : self.fudge ,
            'xi'       : None if xi    is None else self.xi     , 
            'kappa'    : None if kappa is None else self.kappa  ,         
            }

    @property
    def xi ( self ) :
        """``xi''-parameter (asymmetry) for Bukin function"""
        return self.__xi
    @xi.setter
    def xi ( self , value ) :
        self.set_value ( self.__xi , value ) 
    
    @property
    def rho ( self ) :
        """``rho''-parameter (tail) for Bukin function"""
        return self.__rho
    @rho.setter
    def rho ( self , value ) :
        self.set_value ( self.__rho , value ) 

    @property
    def kappa ( self ) :
        """``kappa''-parameter (tail asymmetry)  for Bukin function"""
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.set_value ( self.__kappa , value ) 

    @property
    def rhoL ( self ) :
        """``rhoL''-parameter (left tail) for Bukin function"""
        return self.__rhoL
    @property
    def rhoR ( self ) :
        """``rhoR''-parameter (right tail) for Bukin function"""
        return self.__rhoR


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
                   mean      = None ,   ## related to mean 
                   gamma     = None ) : ## related to asymmetry 
        
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

        self.lambd.setMax ( self.lambd.getMax() * 100 ) ## adjust it! 
    
        ## asymmetry parameter 
        self.__gamma = self.make_var ( ZERO if gamma is None else gamma , 
                                       'gamma_%s'     % name ,
                                       '#gamma(%s)'   % name ,
                                       ZERO if gamma is None else gamma ,                                        
                                       0 , -1000 , 1000 ) 
        
        self.__delta = self.make_var ( delta                 ,
                                       'delta_%s'     % name ,
                                       '#delta(%s)'   % name , delta ,
                                       1 , 1.e-6 , 1000   )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.JohnsonSU (
            self.roo_name ( 'rjsu_' )       ,
            "Resolution Johnson's SU %s" % self.name ,
            self.xvar       ,
            self.xi         ,
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
            'gamma'     : None if gamma is None else self.gamma 
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
        """``gamma''-parameter for Johnson-SU function - related to asymmetry"""
        return self.__gamma
    @gamma.setter
    def gamma ( self , value ) :
        self.set_value ( self.__gamma , value )
        
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
    - epsilon  : asymmetry parameter 
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   sigma     = None ,
                   delta     = 1    ,
                   fudge     = 1    ,   ## fudge factor 
                   mean      = None ,   ## mean value 
                   epsilon   = None ) : ## asymmety parameter 

        ## initialize the base 
        super(ResoSinhAsinh,self).__init__ ( name  = name  ,
                                             xvar  = xvar  ,
                                             sigma = sigma ,
                                             mean  = mean  ,
                                             fudge = fudge )
        
        ## parameter epsilon:
        self.__epsilon =  self.make_var ( ZERO if epsilon is None else epsilon ,
                                          'epsilon_%s'  % name ,
                                          '#epsilon(%s)' % name ,
                                          ZERO if epsilon is None else epsilon , 0 , -1000 , 1000   )

        ## parameter delta 
        self.__delta  = self.make_var ( delta ,
                                        'delta_%s'   % name ,
                                        '#delta(%s)' % name , delta ,
                                        1 , 1.e-6 , 1000   )
    
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.SinhAsinh (
            self.roo_name ( 'rsash_' )       ,
            "Resolution SinhAsinh %s" % self.name ,
            self.xvar       ,
            self.mean       ,
            self.sigma_corr , ## ATTENTION! 
            self.epsilon    ,
            self.delta      )
        
        ## save the configuration
        self.config = {
            'name'      : self.name   ,
            'xvar'      : self.xvar   ,
            'mean'      : self.mean   ,
            'sigma'     : self.sigma  ,
            'delta'     : self.delta  ,
            'fudge'     : self.fudge  ,
            'epsilon'   : epsilon if epsilon is None else self.epsilon , 
            }

    @property
    def epsilon( self ) :
        """``epsilon''-parameter for Sinh-Asinh function"""
        return self.__epsilon
    @epsilon.setter
    def epsilon ( self , value ) :
        self.set_value ( self.__epsilon , value )
        
    @property
    def delta ( self ) :
        """``delta-parameter'' for Sinh-Asinh function"""
        return self.__delta
    @delta.setter
    def delta ( self, value ) :
        value = float ( value )
        assert 0  < value, "``delta''-parameter must be positive"         
        self.set_value ( self.__delta , value  ) 
        
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
    r""" Logistic, aka ``sech-square'' PDF
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
            self.roo_name ( 'rlog_' )       ,
            "Resolution Logistic %s" % self.name ,
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
## @class ResoHyperbolic
#  @see Ostap::Math::Hyperbolic
#  @see Ostap::Models::Hyperbolic
class ResoHyperbolic(RESOLUTION) :
    """Symmetric Hyperbolic distribution
    - see Ostap::Math::Hyperbolic
    - see Ostap::Models::Hyperbolic
    - see Hyperbolic_pdf 
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   sigma     = None ,   ## related to sigma
                   zeta      = 1    ,   ## shape parameter 
                   fudge     = 1    ,   ## fudge parameter 
                   mu        = None ,   ## related to mean 
                   kappa     = None ) : ## asymmetry
        
        ## initialize the base 
        super(ResoHyperbolic,self).__init__ ( name       = name  ,
                                              xvar       = xvar  ,
                                              sigma      = sigma ,
                                              mean       = mu    ,
                                              fudge      = fudge ,
                                              mean_name  = 'mu_%s'   % name ,
                                              mean_title = '#mu(%s)' % name )
        
        ## Zeta
        self.__zeta  = self.make_var ( zeta                ,
                                       'zeta_%s'    % name ,
                                       '#zeta(%s)'  % name , zeta  , 1 , 1.e-10  , +100 )
        
        ## mu 
        self.__mu    = self.mean 

        ## parameter kappa - asymmetry
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa  ,
                                       'kappa_%s'    % name ,
                                       '#kappa(%s)'  % name ,
                                       ZERO if kappa is None else kappa  , 0 , -100 , +100 ) 
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Hyperbolic (
            self.roo_name ( 'rhyp_' ) , 
            "Resolution Hyperbolic %s" % self.name ,
            self.xvar       ,
            self.mu         ,
            self.sigma_corr ,
            self.zeta       ,
            self.kappa      )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'sigma'     : self.sigma ,
            'zeta'      : self.zeta  ,
            'fudge'     : self.fudge ,
            'mu'        : mu    if mu    is None else self.mu    ,  
            'kappa'     : kappa if kappa is None else self.kappa 
            }

    @property
    def mu ( self ) :
        """``mu'' : location parameter, same as ``mean'')"""
        return self.__mu

    @property 
    def zeta  ( self ) :
        """``zeta'' : dimensioneless parameter, related to kurtosis"""
        return self.__zeta
    @zeta.setter  
    def zeta ( self , value ) :
        self.set_value ( self.__zeta , value )
    
    @property
    def kappa ( self ) :
        """``kappa'' : dimensionless parameter, related to asymmetry"""
        return self.__kappa
    @kappa.setter  
    def kappa ( self , value ) :
        self.set_value ( self.__kappa, value )


    @property
    def alpha ( self ) :
        """``alpha'' : value of canonical parameter ``alpha''"""
        self.pdf.setPars ()
        return self.pdf.function().alpha ()

    @property
    def beta ( self ) :
        """``beta'' : value of canonical parameter ``beta''"""
        self.pdf.setPars ()
        return self.pdf.function().beta ()

    @property
    def gamma ( self ) :
        """``gamma'' : value of canonical parameter ``gamma''"""
        self.pdf.setPars ()
        return self.pdf.function().gamma ()
    
    @property
    def delta ( self ) :
        """``delta'' : value of canonical parameter ``delta''"""
        self.pdf.setPars ()
        return self.pdf.function().delta ()

    @property
    def nominal_mean ( self ) :
        """``nominal_mean'' : actual mean of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().delta ()

    @property
    def nominal_mode ( self ) :
        """``nominal_mode'' : actual mode of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().mode ()
    
    @property
    def nominal_variance ( self ) :
        """``nominal_variance'' : actual variance of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().variance ()

    @property
    def nominal_rms ( self ) :
        """``nominal_rms'' : actual RMS of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().rms ()


# =============================================================================
## @class ResoGenHyperbolic
#  @see Ostap::Math::GenHyperbolic
#  @see Ostap::Models::GenHyperbolic
class ResoGenHyperbolic(RESOLUTION) :
    """Symmetric generalised Hyperbolic distribution
    - see Ostap::Math::GenHyperbolic
    - see Ostap::Models::GenHyperbolic
    - see GenHyperbolic_pdf 
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   sigma     = None ,   ## related to sigma
                   zeta      = 1    ,   ## related to shape 
                   lambd     = 1    ,   ## related to shape 
                   fudge     = 1    ,   ## fudge-parameter 
                   mu        = None ,   ## related to mean 
                   kappa     = None ) : ## asymmetry 
        
        ## initialize the base 
        super(ResoGenHyperbolic,self).__init__ ( name       = name             ,
                                                 xvar       = xvar             ,
                                                 sigma      = sigma            ,
                                                 mean       = mu               ,
                                                 fudge      = fudge            ,
                                                 mean_name  = 'mu_%s'   % name ,
                                                 mean_title = '#mu(%s)' % name )
        
        
        ## Zeta
        self.__zeta  = self.make_var ( zeta                ,
                                       'zeta_%s'    % name ,
                                       '#zeta(%s)'  % name , zeta  , 1 , 1.e-10 , +100 ) 

        ## parameter kappa
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa ,
                                       'kappa_%s'    % name ,
                                       '#kappa(%s)'  % name ,
                                       ZERO if kappa is None else kappa , 0 , -100 , +100 ) 
        
        ## lambda 
        self.__lambda = self.make_var ( lambd               ,
                                        'lambda_%s'   % name ,
                                        '#lambda(%s)' % name , lambd , -2 , -100 , 100 )
        
        ## mu 
        self.__mu    = self.mean 
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.GenHyperbolic (
            self.roo_name ( 'rghyp_' ) , 
            "Resolution GenHyperbolic %s" % self.name ,
            self.xvar       ,
            self.mu         ,
            self.sigma_corr ,
            self.zeta       ,
            self.kappa      ,
            self.lambd      )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'sigma'     : self.sigma ,
            'zeta'      : self.zeta  ,
            'lambd'     : self.lambd ,
            'fudge'     : self.fudge ,
            'mu'        : None  if mu    is None else self.mu    ,
            'kappa'     : kappa if kappa is None else self.kappa , 
            }

    @property
    def mu ( self ) :
        """``mu'' : location parameter, same as ``mean'')"""
        return self.__mu

    @property 
    def zeta  ( self ) :
        """``zeta'' : dimensioneless parameter, related to shape"""
        return self.__zeta
    @zeta.setter  
    def zeta ( self , value ) :
        self.set_value ( self.__zeta , value )
    
    @property
    def kappa ( self ) :
        """``kappa'' : dimensionless parameter, related to asymmetry"""
        return self.__kappa
    @kappa.setter  
    def kappa ( self , value ) :
        self.set_value ( self.__kappa, value )
    
    @property
    def lambd ( self ) :
        """``lambd'' : dimensionless parameter, related to shape """
        return self.__lambda
    @lambd.setter
    def lambd ( self , value ) :    
        self.set_value ( self.__lambd , value )
        
    @property
    def alpha ( self ) :
        """``alpha'' : value of canonical parameter ``alpha''"""
        self.pdf.setPars ()
        return self.pdf.function().alpha ()

    @property
    def beta ( self ) :
        """``beta'' : value of canonical parameter ``beta''"""
        self.pdf.setPars ()
        return self.pdf.function().beta ()

    @property
    def gamma ( self ) :
        """``gamma'' : value of canonical parameter ``gamma''"""
        self.pdf.setPars ()
        return self.pdf.function().gamma ()
    
    @property
    def delta ( self ) :
        """``delta'' : value of canonical parameter ``delta''"""
        self.pdf.setPars ()
        return self.pdf.function().delta ()

    @property
    def nominal_mean ( self ) :
        """``nominal_mean'' : actual mean of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().delta ()

    @property
    def nominal_mode ( self ) :
        """``nominal_mode'' : actual mode of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().mode ()
    
    @property
    def nominal_variance ( self ) :
        """``nominal_variance'' : actual variance of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().variance ()

    @property
    def nominal_rms ( self ) :
        """``nominal_rms'' : actual RMS of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().rms ()
     
# =============================================================================
## @class ResoHypatia
# Variant of Hypatia pdf
# @see D. Martinez Santos, F. Duipertois,
#      "Mass distributions marginalized over per-event errors",
#       Nucl.Instrum.Meth.A 764 (2014) 150,
#       arXiv:1312.5000 [hep-ex]
# @see https://doi.org/10.1016/j.nima.2014.06.081
# @see https://arxiv.org/abs/1312.5000
# Actually this function corresponds to Hypatia function with
# \f$ a\rigaharrow +\infty, n=1\f$ 
#
# Convolution of Generalized Hyperbolic distrobution with "offset"
# Gaussian distribution
# @see GenHyperbolic_pdf
# @see Ostap::Math::GenHyperbolic
# @see Ostap::Models::GenHyperbolic
class ResoHypatia(RESOLUTION):
    r""" Variant of Hypatia pdf
    Convolution of Generalized Hyperbolic distrobution with ``offset''
    Gaussian distribution
    
    - see D. Martinez Santos, F. Duipertois,
    ``Mass distributions marginalized over per-event errors'',
    Nucl.Instrum.Meth.A 764 (2014) 150,
    arXiv:1312.5000 [hep-ex]
    - see https://doi.org/10.1016/j.nima.2014.06.081
    - see https://arxiv.org/abs/1312.5000

    Actually this function corresponds to Hypatia function with
    a -> +infinity, n=0
    
    - see Hypatia_pdf
    - see GenHyperbolic_pdf
    - see Ostap.Math.GenHyperbolic
    - see Ostap.Models.GenHyperbolic
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   sigma     =  1   ,   ## relatd  to width  
                   zeta      =  1   ,   ## related to shape 
                   lambd     = -2   ,   ## related to shape 
                   sigma0    = None ,   ## width of the ``offset'' Gaussian
                   fudge     =  1   ,   ## fudge factor 
                   mu        = None ,   ## related to mean
                   kappa     = None ,   ## related to asymmetry
                   cnvpars   = {}   ) : ## convolution parameters
        
        ## initialize the base 
        super(ResoHypatia,self).__init__ ( name       = name  ,
                                           xvar       = xvar  ,
                                           sigma      = sigma ,
                                           mean       = mu    ,
                                           fudge      = fudge ,
                                           mean_name  = 'mu_%s'   % name ,
                                           mean_title = '#mu(%s)' % name )
        
        self.__mu    = self.mean 

        ## Zeta
        self.__zeta   = self.make_var ( zeta                 ,
                                        'zeta_%s'     % name ,
                                        '#zeta(%s)'   % name , zeta  ,  1 , 1.e-10 , 1.e+5 ) 
        ## kappa  
        self.__kappa  = self.make_var ( ZERO if kappa is None else kappa , 
                                        'kappa_%s'    % name ,
                                        '#kappa(%s)'  % name ,
                                        ZERO if kappa is None else kappa , 0 ,  -50  ,  50 ) 
        
        ## lambda 
        self.__lambda = self.make_var ( lambd               ,
                                        'lambda_%s'   % name ,
                                        '#lambda(%s)' % name , lambd ,  -2 , -100   , 100 ) 
        

        ## create a generalized hyperbolic PDF 
        hname  = self.generate_name ( prefix = self.name , suffix = 'GHD' )
        from   ostap.fitting.signals import GenHyperbolic_pdf 
        with CheckMean ( False ) :
            self.__genhyp = GenHyperbolic_pdf ( name  = hname           , 
                                                xvar  = self.xvar       ,
                                                mu    = self.mu         , 
                                                sigma = self.sigma_corr , ## ATTENTION HERE!
                                                zeta  = zeta            ,
                                                kappa = kappa           ,
                                                lambd = lambd           )
        
        ## prepare FFT convolution
        from ostap.fitting.resolution import ResoGauss 
        gname = self.generate_name ( prefix = self.name , suffix = 'offset' ) 
        self.__resolution = ResoGauss     ( name  = gname     , 
                                            xvar  = self.xvar ,
                                            sigma = sigma0    )
        
        self.__cnvpars = {}
        self.__cnvpars.update ( cnvpars ) 
        
        cname = self.generate_name ( prefix = self.name , suffix = 'cnv' )  
        from ostap.fitting.convolution import Convolution_pdf 
        self.__convolved = Convolution_pdf ( name       = cname             , 
                                             pdf        = self.genhyp       ,
                                             xvar       = self.xvar         ,
                                             resolution = self.__resolution ,
                                             **self.cnvpars                 ) 
        
        ## final 
        self.pdf = self.convolved.pdf
        
        self.config = {
            'name'    : self.name    ,
            'xvar'    : self.xvar    ,
            'sigma'   : self.sigma   ,
            'zeta'    : self.zeta    ,
            'lambd'   : self.lambd   ,
            'sigma0'  : self.sigma0  ,
            'cnvpars' : self.cnvpars ,
            'fudge'   : self.fudge   , 
            'mu'      : None  if mu    is None else self.mu    ,
            'kappa'   : kappa if kappa is None else self.kappa }
        
    @property
    def genhyp ( self ) :
        """``genhyp'': get underlying generalized hyperbilis PDF"""
        return self.__genhyp
    
    @property
    def convolved ( self ) :
        """``convolved'' : get PDF as convolution"""
        return self.__convolved
    
    @property
    def sigma0    ( self ) :
        """``sigma0'' : width for the ``offset'' Gaussian"""
        return self.__resolution.sigma
    @sigma0.setter
    def sigma0    ( self , value ) :
        self.__resolution.sigma = value

    
    @property
    def cnvpars ( self ) :
        """``cnvpars'' : parameters for convolution"""
        return self.__cnvpars 

    @property
    def mu ( self ) :
        """``mu'' : location parameter (same as ``mean'')"""
        return self.__mu
    @mu.setter
    def mu ( self , value ) :    
        self.set_value ( self.__mu , value )
   

    @property 
    def zeta  ( self ) :
        """``zeta'' : dimensioneless parameter, related to shape """
        return self.__zeta
    @zeta.setter  
    def zeta ( self , value ) :
        self.set_value ( self.__zeta , value )
    
    @property
    def kappa ( self ) :
        """``kappa'' : dimensionless parameter, related to asymmetry"""
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :    
        self.set_value ( self.__kappa , value )

    @property
    def lambd ( self ) :
        """``lambd'' : dimensionless parameter, related to shape """
        return self.__lambda
    @lambd.setter
    def lambd ( self , value ) :    
        self.set_value ( self.__lambd , value )

    @property
    def alpha ( self ) :
        """``alpha'' : value of canonical parameter ``alpha''"""
        return self.genhyp.alpha 

    @property
    def beta ( self ) :
        """``beta'' : value of canonical parameter ``beta''"""
        return self.genhyp.beta

    @property
    def gamma ( self ) :
        """``gamma'' : value of canonical parameter ``gamma''"""
        return self.genhyp.gamma
    
    @property
    def delta ( self ) :
        """``delta'' : value of canonical parameter ``delta''"""
        return self.genhyp.delta 


# =============================================================================
## @class ResoDas
#  @see Ostap::Math::Das
#  @see Ostap::Models::Das
class ResoDas(RESOLUTION) :
    """Das resoltuoonmodel: gaussian with exponential tails 
    - see Ostap::Math::Das
    - see Ostap::Models::Das
    - see Das_pdf 
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   sigma     = None ,   ## related to sigma
                   k         = None ,   ## tail parameter  
                   fudge     = 1    ,   ## fudge-parameter 
                   mean      = None ,   ## related to mean 
                   kappa     = None ) : ## asymmetry 
        
        ## initialize the base 
        super(ResoDas,self).__init__ ( name  = name  ,
                                       xvar  = xvar  ,
                                       sigma = sigma ,
                                       mean  = mean  ,
                                       fudge = fudge )
        
        self.__k     = self.make_var ( k                    ,
                                       'k_%s'      % name   ,  
                                       'k(%s)'     % name   , k , 1 , 1.e-6 , 1000 )
        
        ## parameter kappa  
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa ,
                                       'kappa_%s'    % name ,
                                       '#kappa(%s)'  % name ,
                                       ZERO if kappa is None else kappa , 0 , -1 , +1 ) 
        
        if kappa  is None :

            self.__kL = self.k
            self.__kR = self.k

        else :

            self.__kL , self.__kR = self.vars_from_asymmetry (
                self.k                                        , ## mean/average k
                self.kappa                                    , ## asymmetry parametet
                v1name  =  self.roo_name ( 'kL' , self.name ) ,
                v2name  =  self.roo_name ( 'kR' , self.name ) ,
                v1title = 'k_{L}: k #times (1+#kappa)'        , 
                v2title = 'k_{R}: k #times (1-#kappa)'        )
            
        ## mu 
        self.__mu    = self.mean 
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Das (
            self.roo_name ( 'rdas_' ) , 
            "Resolution Das %s" % self.name ,
            self.xvar       ,
            self.mu         ,
            self.sigma_corr ,
            self.kL         ,
            self.kR         )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            'k'         : self.k     ,
            'kappa'     : None if kappa is None else self.kappa , 
            'fudge'     : self.fudge }

    @property
    def mu ( self ) :
        """``mu'' : location parameter, same as ``mean'')"""
        return self.__mu
    @mu.setter
    def mu ( self , value ) :
        self.mean = value

    @property
    def k ( self ) :
        """``k'' : tail parameter (exponential slope)
        """
        return self.__k
    @k.setter
    def k ( self , value ) :
        self.set_value ( self.__k , value )
        
    @property
    def kappa ( self ) :
        """``kappa'' : dimensionless parameter, related to asymmetry"""
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.set_value ( self.__kappa , value )

    @property
    def kL  ( self ) :
        """``kL'' : left tail parameter
        """
        return self.__kL

    @property
    def kR  ( self ) :
        """``kR'' : right tail parameter
        """
        return self.__kR
     
# =============================================================================



# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
