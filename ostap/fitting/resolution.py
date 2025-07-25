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
#  - Needham's function                  (variant of Crystal Ball(
#  - Student-T                           (power-law  tails)
#  - Sinh-Asinh model                    (tails can be heavy or light)
#  - JohnsonSU  model                    (tails can be heavy or light)
#  - Hyperbolic model                    (tails are exponential)
#  - generalized Hyperbolic model        (tails are exponential or heavier)
#  - Hypatia model                       (tails are exponential or heavier)
#  - Das model                           (gaussian with exponential tails)
#  - Normal Laplace model                (gaussian with exponential tails)
#  - PearsonIV model                     (power-law tails+asymmetry) 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-07-13
# =============================================================================
"""Set of useful resolution models:
- single Gaussian                     (gaussian    tails)
- double Gaussian                     (gaussian    tails)
- Apollonios-2                        (exponential tails)
- Sech/hyperbolic  secant             (exponential tails)
- Logistic/Sech-squared               (exponential tails) 
- GenLogisticIV                       (exponential tails+asymmetry) 
- Bukin                               (exponential or gaussian tails)
- double-sided Crystal Ball           (power-law  tails)
- Needham's function                  (variant of Crystal Ball(
- Student-T                           (power-law  tails)
- Sinh-Asinh model                    (tails can be heavy or light)
- JohnsonSU  model                    (tails can be heavy or light)
- Hyperbolic                          (tails are exponential)
- generalized Hyperbolic              (tails are exponential or heavier)
- Hypatia model                       (tails are exponential or heavier)
- Normal Laplace model                (gaussian with exponential tails) 
- Buin2  model                        (gaussian with exponential tails) 
- Das model                           (gaussian with exponential tails)
- Generalized Gaussian v1             (family that included Gaussian, Laplace, uniform etc...)
- PearsonIV model                     (power-law tails+asymmetry)
- Bates-shape models                  (spline-like finite model)
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
    'ResoCB2a'          , ## double-sided Crystal Ball resolution model,
    'ResoNeedham'       , ## variant of Crystal Ball funxtion
    'ResoStudentT'      , ## Student-T resolution model,
    'ResoPearsonIV'     , ## Pearson Tyep IV resolution model
    'ResoSkewGenT'      , ## Skewed Generalized t-distribution 
    'ResoSkewGenError'  , ## Skewed Generalized Error-distribution 
    'ResoSech'          , ## Sech/hyperbolic secant  resolution model
    'ResoLogistic'      , ## Logistic ("sech-squared") resolution model
    'ResoGenLogisticIV' , ## Generalized Logistic Type IV resolution model
    'ResoBukin'         , ## Bukin resolution model
    'ResoJohnsonSU'     , ## Jonnson's SU resolution model 
    'ResoSinhAsinh'     , ## Sinh-Asinh resolution model
    'ResoHyperbolic'    , ## Hyperbolic resolution model
    'ResoGenHyperbolic' , ## Generalised Hyperbolic resolution model
    'ResoHypatia'       , ## Hypatia resoltuion model
    'ResoDas'           , ## Das resolution model
    'ResoBukin2'        , ## Bukin2 resolution model
    'ResoNormalLaplace' , ## Normal Laplace resolution model
    'ResoGenGaussV1'    , ## Generalized  Gaussian v1
    'ResoBatesShape'    , ## Bates-shape models (spline-like finite model)

    )
# =============================================================================
from   ostap.core.core          import Ostap
from   ostap.fitting.pdfbasic   import Generic1D_pdf
from   ostap.fitting.fithelpers import ZERO
from   ostap.fitting.fit1d      import RESOLUTION, CheckMean 
import ROOT, math 
# =============================================================================    
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.resolution' )
else                       : logger = getLogger ( __name__                   )
# =============================================================================
models = set()
# =============================================================================
## single gaussian model for resolution
# =============================================================================
## @class ResoGauss
#  Trivial single gaussian resolution model
class ResoGauss(RESOLUTION) :
    """ Trivial single gaussian resolution model
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
                                       None , 0 , -1.0+1.e-9 , +1.0-1.e-9 )
        
        if kappa is None :            
            self.__AV_SIGMA = self.asymmetry_vars ( 'sigma'                   ,
                                                    var1    = self.sigma_corr ,
                                                    var2    = self.sigma_corr )
        else :            
            self.__AV_SIGMA = self.asymmetry_vars ( 'sigma'                   ,
                                                    halfsum = self.sigma_corr ,
                                                    kappa   = self.__kappa    )
            
        # self.gauss = ROOT.RooGaussModel(
        if kappa is None or self.__kappa is ZERO : 
            
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
        """'kappa' : asymmetry parameter"""
        return self.__AV_SIGMA.kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.__AV_SIGMA.kappa = value 

    @property
    def sigmaL ( self ) :
        """'sigmaL': left sigma"""
        return self.__AV_SIGMA.var1
    @property
    def sigmaR ( self ) :
        """'sigmaR': left sigma"""
        return self.__AV_SIGMA.var2 
        
        
models.add ( ResoGauss ) 
# =============================================================================
## @class ResoGauss2
#  Double Gaussian model for  resolution.
#  Parameters: 
#  - sigma of core Gaussian
#  - ratio of wide/core widths
#  - fraction of core(narrow) component
class ResoGauss2(RESOLUTION) :
    """ Double-Gaussian resolution model
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
        self.__fraction = self.make_var ( fraction                   , 
                                          'CoreFraction_'     + name ,
                                          'CoreFraction(%s)'  % name ,
                                          None , 0.75 , 0.01 ,  1 ) 
        
        ## sigma-2/sigma-1 width ratio;
        self.__scale = self.make_var   ( scale ,
                                         'SigmaScale_'       + name ,
                                         'SigmaScale(%s)'    % name ,
                                         None , 1.5  , 1.01 , 10 ) 
        
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
        """'fraction' : fraction parameter for the double Gaussian resolution function
        """
        return self.__fraction
    @fraction.setter
    def fraction ( self , value ) :
        value = float ( value )
        assert 0 <= value <= 1, "'Fraction' must be in  (0,1) range!"
        self.set_value ( self.__fraction , value ) 

    @property
    def scale ( self  ) :
        """'scale' : scale  parameter for double Gaussian resolution function
        """
        return self.__scale
    @scale.setter
    def scale ( self , value ) :
        value = float ( value )
        assert 1 < value, "'scale'-parameter must be larger than 1"
        self.set_value ( self.__scale , value )
  
models.add ( ResoGauss2 ) 
# =============================================================================
## @class ResoApo2
#  (A)Symmetrical  Apollonios  model for resolution
#   - (asymmetrical)Gaussian core 
#   - exponential tails
#  @see Ostap::Models::Apollonios2 
#  @see Ostap::Math::Apollonios2 
class ResoApo2(RESOLUTION) :
    """ (A)Symmetric variant of Apollonios model for the resolution function
    - (asymmetrical) Gaussian core 
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
        
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa   ,
                                       'kappa_%s'   % self.name           ,
                                       '#kappa(%s)' % self.name           ,
                                       None , 0 , -1.0+1.e-9 , +1.0-1.e-9 ) 
        
        if kappa is None or self.__kappa is ZERO :
        
            self.__AV_SIGMA = self.asymmetry_vars ( 'sigma' ,
                                                    var1    = self.sigma_corr ,
                                                    var2    = self.sigma_corr )
        else :
            
            self.__AV_SIGMA = self.asymmetry_vars ( 'sigma' ,
                                                    halfsum = self.sigma_corr ,
                                                    kappa   = self.__kappa    )
            
        self.__beta    = self.make_var ( beta ,
                                         'beta_%s'   % name  ,
                                         '#beta(%s)' % name  ,
                                         None , 0.0001 , 10000 )

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
        """'beta' parameter for symmetric Apollonios resolution function"""
        return self.__beta
    @beta.setter
    def beta ( self , value ) :
        value = float ( value )
        assert 0< value , "'beta'-parameter must be positive!"
        self.set_value ( self.__beta , value )

    @property
    def kappa ( self ) :
        """'kappa' : asymmetry parameter
        """
        return self.__AV_SIGMA.kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.__AV_SIGMA_kappa = value
        
    @property
    def sigmaL ( self )  :
        """'sigmaL' : left sigma-parameter"""
        return self.__AV_SIGMA.var1
    @property
    def sigmaR ( self )  :
        """'sigmaR' : right  sigma-parameter"""
        return self.__AV_SIGMA.var2 
    
models.add ( ResoApo2 )

# =============================================================================
## @class ResoCB2a
#  (A)Symmetrical double-sided Crystal Ball model for resolution
#   - (symmetric) Gaussian core 
#   - power-law tails
#  @see Ostap::Math::CrystalBallDS
#  @see Ostap::Models::CrystalBallDS
class ResoCB2a(RESOLUTION) :
    """ (A)Symmetric double-sided Crystal Ball model for resolution
    - Gaussian core 
    - power-law tails
    see Ostap.Math.CrystalBallDS
    see Ostap.Models.CrystalBallDS
    """
    def __init__ ( self          , 
                   name          ,   ## the   name 
                   xvar          ,   ## the   variable 
                   sigma         ,   ## core  resolution
                   alphaL = 1.5  ,   ## left  alpha
                   alphaR = None ,   ## right alpha (the same as left alpha if None)
                   nL     = 5    ,   ## left  N 
                   nR     = None ,   ## right N 
                   fudge  = 1    ,   ## fudge-factor 
                   mean   = None ) : ## the mean value

        ## initialize the base 
        super(ResoCB2a,self).__init__ ( name  = name  ,
                                        xvar  = xvar  ,
                                        sigma = sigma ,
                                        mean  = mean  ,
                                        fudge = fudge )
        
        alpha_pars = None , 2.0 , 0.1 , 4.0        
        self.__alphaL = self.make_var ( alphaL               ,
                                        'alphaL_'     + name ,
                                        '#alphaL(%s)' % name , *alpha_pars )
        
        self.__alphaR = self.make_var ( self.alphaL if alphaR is None else alphaR , 
                                        'alphaR_'     + name ,
                                        '#alphaR(%s)' % name , *alpha_pars )
        
        n_pars = None , 5 , 1.e-8 , 200 
        self.__nL     = self.make_var ( nL                   ,
                                        'nL_'         + name ,
                                        'nL(%s)'      % name , *n_pars)        
        self.__nR     = self.make_var ( self.nL if nR is None else nR ,
                                        'nR_'    + name   ,
                                        'nR(%s)' % name   , *n_pars ) 
        
        self.__AV_ALPHA = self.asymmetry_vars ( 'alpha' , var1 = self.alphaL , var2 = self.alphaR )
        self.__AV_N     = self.asymmetry_vars ( 'n'     , var1 = self.nL     , var2 = self.nR     )
        
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
            'name'     : self.name                       ,
            'xvar'     : self.xvar                       ,
            'mean'     : self.mean                       ,
            'sigma'    : self.sigma                      ,
            'alphaL'   : self.alphaL                     ,
            'nL'       : self.nL                         ,
            'alphaR'   : self.alphaR if alphaR else None ,
            'nR'       : self.nR     if nR     else None ,
            'fudge'    : self.fudge                      ,
        }

    @property
    def alpha ( self  ) :
        """'alpha' parameter for double-sided (a)symmetric resolution function
        """
        return self.__AV_ALPHA.halfsum
    
    @property
    def n    ( self  ) :
        """'n' parameter for double-sided (a)symmetric resolution function
        """
        return self.__AV_N.halfsum 
    
    @property
    def alphaL ( self  ) :
        """'alphaL' parameter for double-sided (a)symmetric resolution function
        """
        return self.__alphaL
    @alphaL.setter
    def alphaL ( self , value ) :
        self.set_value ( self.__alphaL , value  ) 

    @property
    def alphaR ( self  ) :
        """'alphaR' parameter for double-sided (a)symmetric resolution function
        """
        return self.__alphaR
    @alphaR.setter
    def alphaR ( self , value ) :
        self.set_value ( self.__alphaR , value  ) 

    @property
    def aL ( self  ) :
        """'aL' parameter for double-sided (a)symmetric resolution function
        """
        return self.__alphaL
    @aL.setter
    def aL ( self , value ) :
        self.set_value ( self.__alphaL , value  ) 

    @property
    def aR ( self  ) :
        """'aR' parameter for double-sided (a)symmetric resolution function
        """
        return self.__alphaR
    @aR.setter
    def aR ( self , value ) :
        self.set_value ( self.__alphaR , value  ) 
        
    @property
    def nL ( self  ) :
        """'nL' parameter for double-sided (a)symmetric resolution function"""
        return self.__nL 
    @nL.setter
    def nL ( self , value ) :
        self.set_value ( self.__nL , value  )
        
    @property
    def nR ( self  ) :
        """'nR' parameter for double-sided (a)symmetric resolution function"""
        return self.__nR 
    @nR.setter
    def nR ( self , value ) :
        self.set_value ( self.__nR , value  ) 

    @property
    def kappaA ( self ) :
        """`kappaA` : asymmetry for `alpha`"""
        return self.__AV_ALPHA.kappa
        
    @property
    def kappaN ( self ) :
        """`kappaN` : asymmetry for `n`"""
        return self.__AV_N.kappa

    @property
    def psiA ( self ) :
        """`psiA` : skew for `alpha`"""
        return self.__AV_ALPHA.psi
        
    @property
    def psiN ( self ) :
        """`psiN` : skew for `n`"""
        return self.__AV_N.psi
        
models.add ( ResoCB2a )
# ===============================================================================
## @class ResoNeedham
#  Needham's functiobn:
#  - variant of Crystall Ball function with alpha = alpha(sigma)
#
# - alpha is parameterized as function of sigma 
#  \f$ \alpha(\sigma) = c_0\frac{ (\sigma/c_1)^{c_2}}{ 1 + (\sigma/c_1)^{c_2} }\f$ 
#
#  @attention For majority of physics cases <code>n</code> 
#             can be fixed <code>n=0</code> (corresponds to <code>N=1</code>
#
#  @attention parameter \f$ c_1 \f$ is inverse with respect to the original 
#             Matt's code
#
#  Reasonable values:
#  - for \f$ c_0 \f$ :  \f$ 1.7 \le c_0 \le 3.5 \f$ 
#  - for \f$ c_1 \f$ :  \f$ c_2 \approx O(\sigma) \f$ 
#  - for \f$ c_2 \f$ :  \f$ c_2 \approx O(10) \f$
# 
#  @see Ostap::Math::Needham
#  @see Ostap::Models::Needham
#  @see Nededham_pdf
class ResoNeedham(RESOLUTION) :
    """ Nedham's function
    - variant of Crystal Ball function with alpha = alpha(sigma)     
    """
    def __init__ ( self           ,  
                   name           ,   ## the  name 
                   xvar           ,   ## the  variable 
                   sigma          ,   ## core resolution                                 
                   c0             ,   ## c0: between 1.8 and 3.5
                   c1             ,   ## c1: close to sigma
                   c2             ,   ## c2: close to 10
                   n       = ROOT.RooFit.RooConst ( 0 ) , 
                   fudge   = 1    ,   ## fudge-factor
                   mean    = None ) : ## the mean value
        
        super(ResoNeedham,self).__init__ ( name  = name  ,
                                           xvar  = xvar  ,
                                           sigma = sigma ,
                                           mean  = mean  ,
                                           fudge = fudge )
        
        self.__c0 = self.make_var ( c0                  ,
                                    "c0_%s"     % name  ,
                                    "c_{0}(%s)" % name  ,
                                    True , 2.5 , 1.5    , 4.0 )
        
        s_minmax = self.sigma.minmax ()         
        if s_minmax :
            smin, smax = s_minmax
            c1limits = 0.1 * smin , 10 * smax
        else :
            c1limits = () 
        
        self.__c1 = self.make_var ( c1                  ,
                                    "c1_%s"     % name  ,
                                    "c_{1}(%s)" % name  ,
                                    True , *c1limits    )  
        
        self.__c2 = self.make_var ( c2                  ,
                                    "c2_%s"     % name  ,
                                    "c_{2}(%s)" % name  ,
                                    True , 10 , 1 , 50  )
        
        ## effective parameteter 
        self.__n     = self.make_var ( n   ,
                                       'n_%s'            % name ,
                                       'n_{CB}(%s)'      % name ,
                                       None , 0 , -1 , 100 )
        
        ## true parameter 
        self.__N = Ostap.MoreRooFit.TailN ( 'N_%s' % name , self.__n )
        
        self.needham = Ostap.Models.Needham (
            self.roo_name ( 'needham_' ) , 
            'Needham function %s' % self.name ,
            self.xvar  ,
            self.mean  ,
            self.sigma_corr , ## ATTENTION HERE
            self.c0    ,
            self.c1    ,
            self.c2    ,
            self.n     , 
            )
        
        self.pdf = self.needham
        
        ## save the configuration
        self.config = {
            'name'   : self.name  ,
            'xvar'   : self.xvar  ,
            'mean'   : self.mean  ,
            'sigma'  : self.sigma ,
            'c0'     : self.c0    ,
            'c1'     : self.c1    ,
            'c2'     : self.c2    ,
            'n'      : self.n     ,
            'fudge'  : self.fudge ,
        }

    @property
    def c0 ( self ) :
        """'c0'-parameter for Needham' function"""
        return self.__c0
    @c0.setter
    def c0 ( self, value ) :
        self.set_value ( self.__c0 , value ) 

    @property
    def c1 ( self ) :
        """'c1'-parameter for Needham' function, INVERSE with respect ot originam Matt's code!"""
        return self.__c1
    @c1.setter
    def c1 ( self, value ) :
        self.set_value ( self.__c1 , value ) 

    @property
    def c2 ( self ) :
        """'c2'-parameter for Needham' function"""
        return self.__c2
    @c2.setter
    def c2 ( self, value ) :
        self.set_value ( self.__c2 , value )
        
    @property
    def n ( self ) :
        """n-parameter for Crystal Ball tail"""
        return self.__n
    @n.setter
    def n ( self, value ) :
        self.set_value ( self.__n , value ) 

    @property
    def N ( self ) :
        """`N` : actual N-parameter used for Crystal Ball """
        return self.__N
  

models.add ( ResoNeedham )

# ===============================================================================
## @class ResoCB2
#  (A)Symmetrical double-sided Crystal Ball model for resolution
#   - Gaussian core 
#   - power-law tails
#  @see Ostap::Math::CrystalBallDS
#  @see Ostap::Models::CrystalBallDS
class ResoCB2(RESOLUTION) :
    """ (A)Symmetric double-sided Crystal Ball model for resolution
    - (symmetric) Gaussian core 
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
                   kappaN = None ,   ## asymmetry for N
                   kappaA = None ) : ## asymmetry for alpha

        ## initialize the base 
        super(ResoCB2,self).__init__ ( name  = name  ,
                                       xvar  = xvar  ,
                                       sigma = sigma ,
                                       mean  = mean  ,
                                       fudge = fudge )

        self.__alpha = self.make_var ( alpha               ,
                                       'alpha_'     + name ,
                                       '#alpha(%s)' % name ,
                                       None , 2.0 , 0.05   , 5 )
        
        self.__n     = self.make_var ( n                  ,
                                       'n_'         + name ,
                                       'n(%s)'      % name ,
                                       None , 5 , 1.e-6 , 200 )
        
        self.__kappaN = self.make_var ( ZERO if kappaN is None else kappaN , 
                                        'kappa_n%s'      % self.name       ,
                                        '#kappa_{n}(%s)' % self.name       ,
                                        None , 0 , -1 , +1 )
        self.__kappaA = self.make_var ( ZERO if kappaA is None else kappaA , 
                                        'kappa_a%s'      % self.name       ,
                                        '#kappa_{a}(%s)' % self.name       ,
                                        None , 0 , -1 , +1 )
        
        self.__AV_ALPHA = self.asymmetry_vars ( 'alpha' , halfsum = self.alpha , kappa = self.kappaA )
        self.__AV_N     = self.asymmetry_vars ( 'n'     , halfsum = self.n     , kappa = self.kappaN )
        
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
        """'alpha' parameter for double-sided symmetric resolution function
        """
        return self.__alpha
    @alpha.setter
    def alpha ( self , value ) :
        self.set_value ( self.__alpha , value  ) 

    @property
    def n ( self  ) :
        """'n' parameter for double-sided symmetric resolution function"""
        return self.__n
    @n.setter
    def n ( self , value ) :
        self.set_value ( self.__n , value  ) 

    @property
    def kappaN ( self ) :
        """'kappaN' : asymmetry for parameter 'n'
        """
        return self.__kappaN
    @kappaN.setter
    def kappaN ( self , value ) :
        self.set_value ( self.__kappaN , value )

    @property
    def kappaA ( self ) :
        """'kappaA' : asymmetry for parameter 'alpha'
        """
        return self.__kappaA
    @kappaA.setter
    def kappaA ( self , value ) :
        self.set_value ( self.__kappaA , value )

    @property
    def nL ( self ) :
        """'nL' : parameter 'n' for left tail
        """
        return self.__AV_N.var1
    
    @property
    def nR ( self ) :
        """'nR' : parameter 'n' for right tail
        """
        return self.__AV_N.var2 

    @property
    def alphaL ( self ) :
        """'alphaL' : parameter 'alpha' for left tail
        """
        return self.__AV_ALPHA.var1 

    @property
    def alphaR ( self ) :
        """'alphaR' : parameter 'alpha' for right tail
        """
        return self.__AV_ALPHA.var2 

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
    - when asymmetry is activated use `BifurcatedStudentT`
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
        
        self.__n      = self.make_var ( n                      ,
                                        'ResoN_'        + name ,
                                        'ResoN(%s)'     % name ,
                                        None , 1 , 1.e-6 , 200 )
        
        self.__kappaN = self.make_var ( ZERO if kappaN is None else kappaN ,
                                        "kappaN_%s"      % name ,
                                        "#kappa_{n}(%s)" % name ,
                                        None , 0 , -1.0 + 1.e-9 , +1.0 - 1.e-9 )
        
        self.__kappaS = self.make_var ( ZERO if kappaS is None else kappaS ,
                                        "kappaS_%s"           % name ,
                                        "#kappa_{#sigma}(%s)" % name ,
                                        None , 0 , -1.0+1.e-9 , 1.0-1.e-9 )
        
        ## n,nL,nR    
        if kappaN is None or self.__kappaN is ZERO :            
            self.__AV_N = self.asymmetry_vars ( 'n'  , var1    = self.__n , var2  = self.__n  )
        else :            
            self.__AV_N = self.asymmetry_vars ( 'n'  , halfsum = self.__n , kappa = self.__kappaN )

        if kappaS is None or self.__kappaS is ZERO :
            self.__AV_SIGMA = self.asymmetry_vars ( 'sigma' , var1    = self.sigma_corr , var2  = self.sigma_corr )
        else :            
            self.__AV_SIGMA = self.asymmetry_vars ( 'sigma' , halfsum = self.sigma_corr , kappa = self.__kappaS   )
                 
        # 
        ## finally build pdf
        #
        if ( kappaN is None or self.__kappaN is ZERO ) and \
           ( kappaS is None or self.__kappaS is ZERO ) : 
            
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
            'kappaS'   : None if kappaS is None else self.kappaS ,  
            'kappaN'   : None if kappaN is None else self.kappaN ,
            }
        
    @property
    def n ( self  ) :
        """'n' parameter for symmetric Student-T resolution function"""
        return self.__AV_N.halfsum 
    @n.setter
    def n ( self , value ) :
        self.__AV_N.halfsum = value 

    @property
    def kappaN ( self  ):
        """'kappaN' : asymmetry for 'n'-parameter"""
        return self.__kappaN
    @kappaN.setter
    def kappaN ( self , value ) :
        self.set_value ( self.__kappaN , value )

    @property
    def kappaS ( self  ):
        """'kappaS' : asymmetry for 'sigma'-parameter"""
        return self.__kappaS
    @kappaS.setter
    def kappaS ( self , value ) :
        self.set_value ( self.__kappaS , value )

    @property
    def nL ( self ) :
        """'nL' : 'n'-parameter for left  part"""
        return self.__AV_N.var1 
    @property
    def nR ( self ) :
        """'nR' : 'n'-parameter for right part"""
        return self.__AV_N_var2

    @property
    def sigmaL ( self ) :
        """'sigmaL' : 'sigma'-parameter for left  part"""
        return self.__AV_SIGMA.var1      
    @property
    def sigmaR ( self ) :
        """'sigmaR' : 'sigma'-parameter for right part"""
        return self.__AV_SIGMA.var2    
        
models.add ( ResoStudentT )

# =============================================================================
## @class ResoPearsonIV
#  (asymmetric) Pearson Type IV model for the resolution
#   - power-law tails
#   - asymmetry
#   Pearson Type IV distribution  
#   \f$ f(x;\mu, n, \kappa) = 
#   C \left( 1 + y^{2}\right)^{-(\frac{1}{2}+n)}
#   \mathrm{e}^{ -\kappa \atan y }}\f$, where 
#   - \f$  y = \frac{x-\mu}{\sigma}\f$,
#   - \f$ 0 < n \f$  
#  @see https://en.wikipedia.org/wiki/Pearson_distribution
#  For $\kappa=0\f$ one gets Student's t-distribution
#  @see J. Heinrich, "A guide to the Pearson Type IV distribution", 
#       CDF/MEMO/STATISTICS/PUBLIC/6820, 2004 
#  @see http://www-cdf.fnal.gov/physics/statistics/notes/cdf6820_pearson4.pdf
#  @see Ostap::Models::PearsonIV 
#  @see Ostap::Math::PearsonIV 
class ResoPearsonIV(RESOLUTION) :
    """ (asymmetric) Pearson Type IV model for the resolution
    - power-law tails
    - asymmetry
    """
    def __init__ ( self           ,
                   name           ,   ## the name 
                   xvar           ,   ## the variable
                   varsigma       ,   ## the width-related parameter 
                   n              ,   ## N-parameter
                   fudge  = 1     ,   ## fudge parameter 
                   mean   = None  ,   ## mean 
                   kappa  = None  ) : ## asymmetry 
        
        ## initialize the base 
        super(ResoPearsonIV,self).__init__ ( name  = name     ,
                                             xvar  = xvar     ,
                                             sigma = varsigma ,
                                             mean  = mean     ,
                                             fudge = fudge    )
        
        self.__n     = self.make_var ( n                      ,
                                       'ResoN_'        + name ,
                                       'ResoN(%s)'     % name ,
                                       False , 2 , 1.e-6 , 200 )

        
        ## asymmetry parameter 
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa                      , 
                                       'kappa_%s'         % name ,
                                       '#kappa_{PIV}(%s)' % name ,
                                       False , 0 , -100 , 100 )
        
        ## finally build PDF 
        self.pdf = Ostap.Models.PearsonIV (
            self.roo_name ( 'p4_' )       ,
            "Resolution Pearson Type IV %s" % self.name ,
            self.xvar       , 
            self.mean       ,
            self.sigma_corr , ## ATTENTION!
            self.n          ,
            self.kappa      )
        
        ##  save   the configuration
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'varsigma' : self.sigma ,
            'n'        : self.n     ,
            'kappa'    : self.kappa ,
            ## 
            'mean'     : self.mean  ,
            'fudge'    : self.fudge ,
            }
        
    @property
    def n ( self  ) :
        """'n' parameter for Pearson Type IV resolution function"""
        return self.__n
    @n.setter
    def n ( self , value ) :
        self.set_value ( self.__n , value ) 

    @property
    def kappa ( self  ):
        """'kappa' : asymmetry parameter for Pearson Type IV function"""
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.set_value ( self.__kappa , value )

    @property
    def varsigma ( self ) :
        """'varsigma' : 'sigma'-related parameter"""
        return self.sigma        
    @varsigma.setter 
    def varsigma ( self , value ) :
        self.sigma = value
    
models.add ( ResoPearsonIV )

# =============================================================================
## @class ResoSkewGenT
#  Skewed generalised t-distribution
#  @see https://en.wikipedia.org/wiki/Skewed_generalized_t_distribution
#  Original function is parameterised in terms of parameters 
#  - \f$ \mu \$ related to locartion 
#  - \f$ \sigma \$ related to width/scale 
#  - \f$ -1 < \lambda < 1 \f$ related to asymmetry/skewness  
#  - \f$ 0<p, 0<q \f$ related to kutsosis
#
#  Mean value is defined if \f$ 1 < pq \f$ 
#  RMS si defined for \f$ 2 < pq \f$
# 
#  In this view here we adopt sligth reparameterisation in terms of 
#  - \f$ 0 < r \f$, such as  \f$  r = \frac{1}{p} 
#  - \f$ 0< \zeta \f$, such as \f$ pq = \zeta + 4 \f$
#  - \f$ -\infty < \xi < +\infty \f$, such as \f$ \lambda  = \tanh \xi \f$   
#
#  Usage of \f$ \zeta\f$ ensures the existance of the  mean, RMS, sewness & kurtosis
# 
#  Special limitnig cases:
#  - \f$ q\rigtharrow +\infty (\zeta \rightarrow +\infty) \f$ 
#     Generalized Error Distribution 
#  - \f$ \lambda=0 (\xi = 0)  \f$ Generalized t-distribution 
#  - \f$ p=2(r=\frac{1}{2}) \f$  Skewed t-distribution 
#  - \f$ p=1(r=1), q\rigtharrow +\infty (\zeta\rightarrow+\infty) \f$
#     Skewed Laplace distribution 
#  - \f$ \lambda=0, q\rigtharrow +\infty (\zeta\rightarrow+\infty) \f$
#     Generalized Error Distribution 
#  - \f$ p=2(r=\frac{1}{2}), q\rigtharrow +\infty (\zeta\rightarrow+\infty) \f$
#     Skewed Normal distribution 
#  - \f$ \sigma=1, \lambda=0,p=2(r=\frac{1}{2},  q=\frac{n+2}{2} (\alpha=n) \f$
#     Student's t-distribution 
#  - \f$ \lambda=0, p=1(r=1), q\rigtharrow +\infty (\zeta\rightarrow+\infty) \f$
#     Laplace distribution 
#  - \f$ \lambda=0, p=2(r=\frac{1}{2}, q\rigtharrow +\infty (\zeta\rightarrow+\infty) \f$
#     Skewed Normal distribution 
#  @see Ostap::Math::SkewGenT 
#  @see Ostap::Models::SkewGenT 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2022-12-01
class ResoSkewGenT(RESOLUTION) :
    """ Skewed Generilized t-distribution
    - power-law, exponential and gauisian tails    
    - asymmetry
    - see https://en.wikipedia.org/wiki/Skewed_generalized_t_distribution
    - see `Ostap.Math.SkewGenT`
    - see `Ostap.Mdoels.SkewGenT`
    """
    def __init__ ( self           ,
                   name           ,   ## the name 
                   xvar           ,   ## the variable
                   sigma          ,   ## the width-related parameter 
                   zeta           ,   ## zeta-parameter
                   r              ,   ## r-parameter                                      
                   fudge  = 1     ,   ## fudge parameter 
                   mean   = None  ,   ## mean 
                   kappa  = None  ) : ## asymmetry (same as xi)
        
        ## initialize the base 
        super(ResoSkewGenT,self).__init__ ( name  = name  ,
                                            xvar  = xvar  ,
                                            sigma = sigma ,
                                            mean  = mean  ,
                                            fudge = fudge )
        
        
        ## r parameter (shape)
        self.__r      = self.make_var ( r                         ,
                                        'r_%s'            % name  ,
                                        'r_{SGT}(%s)'     % name  ,
                                        False , 1 , 0.0001 , 1000 ) 
        ## zeta parameter (shape)
        self.__zeta   = self.make_var ( zeta                     ,
                                        'zeta_%s'         % name ,
                                        '#zeta_{SGT}(%s)' % name ,
                                        False , 1 , 0    , 1000  ) 
        
        ## asymmetry parameter 
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa                      , 
                                       'xi_%s'         % name ,
                                       '#xi_{SGT}(%s)' % name ,
                                       False , 0 , -100 , 100 )
            
        ## finally build PDF 
        self.pdf = Ostap.Models.SkewGenT (
            self.roo_name ( 'sgt_' )       ,
            "Resolution Skewed Generalised t:  %s" % self.name ,
            self.xvar       , 
            self.mean       ,
            self.sigma_corr , ## ATTENTION!
            self.xi         , ## same as kappa 
            self.r          ,
            self.zeta       )
        
        ##  save   the configuration
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'sigma'    : self.sigma ,
            'zeta'     : self.zeta  ,
            'r'        : self.r     ,
            'kappa'    : self.kappa ,
            ## 
            'mean'     : self.mean  ,
            'fudge'    : self.fudge ,
            }
        
    @property
    def kappa ( self  ):
        """'kappa' : asymmetry parameter for SkewGenT function (same as 'xi')"""
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.set_value ( self.__kappa , value )

    @property
    def xi ( self  ):
        """'xi' : asymmetry parameter for SkewGenT function (same as 'kappa')"""
        return self.__kappa
    @xi.setter
    def xi ( self , value ) :
        self.set_value ( self.__kappa , value )

    @property
    def r ( self  ):
        """'r' : shape parameter for SkewGenT function"""
        return self.__r
    @r.setter
    def r ( self , value ) :
        self.set_value ( self.__r , value )

    @property
    def zeta ( self  ):
        """'zeta' : shape parameter for SkewGenT function"""
        return self.__zeta
    @zeta.setter
    def kappa ( self , value ) :
        self.set_value ( self.__zeta , value )

models.add ( ResoSkewGenT )

# =============================================================================
## @class ResoSkewGenError
#  Skewed gheneralised error districbution 
#  @see https://en.wikipedia.org/wiki/Skewed_generalized_t_distribution#Skewed_generalized_error_distribution
#
#  The Special  case of Skewwed Generaliaed T-distribution 
#  @see Ostap::Math::SkewGenT 
# 
#  Original function is parameterised in terms of parameters 
#  - \f$ \mu \$ related to location  
#  - \f$ \sigma \$ related to width/scale 
#  - \f$ -1 < \lambda < 1 \f$ related to asymmetry/skewness  
#  - \f$ 0<p \f$ shape parameters 
#
#  \f[ f(x;\mu,\sigma,\lambda,p) = 
#    \frac{p}{2v\sigma\Gamma(1/p)} \mathrm{e}^{ - \Delta^{p}},  
#   \f]
#  where 
#   - \f$ v = \sqrt{ \frac{ \pi \Gamma(1/p)}{  \pi(1+3\lambda^2)\Gamma(3/p) 
#            -16^{1/p} \lambda^2 \Gamma(1/2+1/p)^2\Gamma(1/p) }  }\f$,
#   - \f$ \Delta = \frac{\left| \delta x \right|}{v\sigma ( 1+ \lambda \sign \delta x )} \f$
#   - \f$ \delta x = x - \mu + m \f$
#   - \f$ m =  2^{2/p} v \sigma \Gamma( 1/2+ 1/p)/\sqrt{\pi}\f$ 
#
#  Here we adopt sligth reparameterisation in terms of 
#  - \f$ -\infty < \xi < +\infty \f$, such as \f$ \lambda  = \tanh \xi \f$   
# 
#  special cases: 
#  - \f$ \xi=0 (\lambda=0), p=2\$ corresponds to Gaussian function 
#  - \f$ \xi=0 (\lambda=0), p=1\$ corresponds to Laplace case 
#
#  @see Ostap::Math::SkewGenError
#  @see Ostap::Math::SkewGenT 
#  @see Ostap::Models::SkewGenT 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2022-12-01
class ResoSkewGenError(RESOLUTION) :
    """ Skewed Generilized Error-distribution
    - exponential and gauisian tails    
    - asymmetry
    - see https://en.wikipedia.org/wiki/Skewed_generalized_t_distribution#Skewed_generalized_error_distribution
    - see `Ostap.Math.SkewGenError`
    - see `Ostap.Models.SkewGenError`
    - see `Ostap.Math.SkewGenError`
    - see `Ostap.Models.SkewGenError`
    """
    def __init__ ( self           ,
                   name           ,   ## the name 
                   xvar           ,   ## the variable
                   sigma          ,   ## the width-related parameter 
                   p              ,   ## r-parameter                                      
                   fudge  = 1     ,   ## fudge parameter 
                   mean   = None  ,   ## mean 
                   kappa  = None  ) : ## asymmetry (same as xi)
        
        ## initialize the base 
        super(ResoSkewGenError,self).__init__ ( name  = name  ,
                                                xvar  = xvar  ,
                                                sigma = sigma ,
                                                mean  = mean  ,
                                                fudge = fudge )
        
        
        ## p-parameter (shape)
        self.__p      = self.make_var ( p                         ,
                                        'p_%s'            % name  ,
                                        'p_{SGE}(%s)'     % name  ,
                                        False , 2 , 0.01 , 100 ) 
        
        ## asymmetry parameter 
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa , 
                                       'xi_%s'         % name ,
                                       '#xi_{SGE}(%s)' % name ,
                                       False , 0 , -100 , 100 )
            
        ## finally build PDF 
        self.pdf = Ostap.Models.SkewGenError  (
            self.roo_name ( 'sge_' )       ,
            "Resolution Skewed Generalised Error:  %s" % self.name ,
            self.xvar       , 
            self.mean       ,
            self.sigma_corr , ## ATTENTION!
            self.xi         , ## same as kappa 
            self.p          )
        
        ##  save   the configuration
        self.config = {
            'name'     : self.name  ,
            'xvar'     : self.xvar  ,
            'sigma'    : self.sigma ,
            'p'        : self.p     ,
            'kappa'    : self.kappa ,
            ## 
            'mean'     : self.mean  ,
            'fudge'    : self.fudge ,
            }
        
    @property
    def kappa ( self  ):
        """'kappa' : asymmetry parameter for Skewed Generalized Error shape (same as 'xi')"""
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.set_value ( self.__kappa , value )

    @property
    def xi ( self  ):
        """'xi' : asymmetry parameter for Skewed Generalized Error shape  (same as 'kappa')"""
        return self.__kappa
    @xi.setter
    def xi ( self , value ) :
        self.set_value ( self.__kappa , value )

    @property
    def p ( self  ):
        """'p' : shape parameter for Skewed Generalized Error shape"""
        return self.__p
    @p.setter
    def p ( self , value ) :
        self.set_value ( self.__p , value )

models.add ( ResoSkewGenError )

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
                                    None , 0 , -10 , +10 )
        
        ## rho 
        self.__rho = self.make_var   ( rho               ,
                                       "rho_%s"   % name ,
                                       "#rho(%s)" % name ,
                                       None , 0 , -1 , 25 )        
        
        ## parameter 
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa ,
                                       "kappa_%s"          % name ,
                                       "#kappa_{#rho}(%s)" % name ,
                                       None , 0 , -1 , +1 )
        
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
                v2title = '#rho_{R}: #rho #times (1-#kappa_{#rho})' )
                
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
        """'xi'-parameter (asymmetry) for Bukin function"""
        return self.__xi
    @xi.setter
    def xi ( self , value ) :
        self.set_value ( self.__xi , value ) 
    
    @property
    def rho ( self ) :
        """'rho'-parameter (tail) for Bukin function"""
        return self.__rho
    @rho.setter
    def rho ( self , value ) :
        self.set_value ( self.__rho , value ) 

    @property
    def kappa ( self ) :
        """'kappa'-parameter (tail asymmetry)  for Bukin function"""
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.set_value ( self.__kappa , value ) 

    @property
    def rhoL ( self ) :
        """'rhoL'-parameter (left tail) for Bukin function"""
        return self.__rhoL
    @property
    def rhoR ( self ) :
        """'rhoR'-parameter (right tail) for Bukin function"""
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
    """ Resolution model based on symmetric form of Johnson's SU distribution
    
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
                                       None , 0 , -1000 , 1000 ) 
        
        self.__delta = self.make_var ( delta                 ,
                                       'delta_%s'     % name ,
                                       '#delta(%s)'   % name ,
                                       None , 1 , 1.e-6 , 1000 )
        
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
        """'delta'-parameter for Johnson-SU function"""
        return self.__delta
    @delta.setter
    def delta ( self, value ) :
        value = float ( value )
        self.set_value ( self.__delta , value )

    @property
    def gamma ( self ) :
        """'gamma'-parameter for Johnson-SU function - related to asymmetry"""
        return self.__gamma
    @gamma.setter
    def gamma ( self , value ) :
        self.set_value ( self.__gamma , value )
        
    @property
    def xi ( self ) :
        """'xi'-parameter (location) for Johnson-SU function (the   same as 'mean')"""
        return self.__xi
    @xi.setter
    def xi ( self, value ) :
        self.set_value ( self.__xi , value )

    @property
    def lambd ( self ) :
        """'lambda'-parameter (scale) for Johnson-SU function (the  same  as 'sigma')"""
        return self.__lambd
    @lambd.setter
    def lambd ( self, value ) :
        self.set_value ( self.__lambd , value )
    
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
    'Sinh-arcsinh distributions'. Biometrika 96 (4): 761. 
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
                                          None , 0 , -1000 , 1000   )
        
        ## parameter delta 
        self.__delta  = self.make_var ( delta ,
                                        'delta_%s'   % name ,
                                        '#delta(%s)' % name ,
                                        None , 1 , 1.e-6 , 1000 )
    
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
        """'epsilon'-parameter for Sinh-Asinh function"""
        return self.__epsilon
    @epsilon.setter
    def epsilon ( self , value ) :
        self.set_value ( self.__epsilon , value )
        
    @property
    def delta ( self ) :
        """'delta'-parameter for Sinh-Asinh function"""
        return self.__delta
    @delta.setter
    def delta ( self, value ) :
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
    r""" Logistic, aka 'sech-square' PDF
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
## @class ResoGenLogisticIV
#  Generalized Logistic Type IV 
#  Type I   : beta  = 1 
#  Type II  : alpha = 1 
#  Type III : alpha = beta         
#  @see https://en.wikipedia.org/wiki/Generalized_logistic_distribution
#  @see Ostap::Math::GenLogisticIV
#  @see Ostap::Models::GenLogisticIV
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-06-14
class ResoGenLogisticIV(RESOLUTION) :
    """ Generalized Logistic Type IV
    - Type I   : beta  = 1 
    - Type II  : alpha = 1 
    - Type III : alpha = beta
    - see https://en.wikipedia.org/wiki/Generalized_logistic_distribution
    - see Ostap::Math::GenLogisticIV
    - see Ostap::Models::GenLogisticIV
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   sigma     = None ,   ## related to sigma
                   gamma     = 1    ,   ## "mean"-tail: 0.5*(alpha+gamma)
                   fudge     = 1    , 
                   mean      = None ,   ## related to mean
                   kappa     = None ) : ## related to asymmetry 0.5*(alpha-gamma)
        
        
        ## initialize the base 
        super(ResoGenLogisticIV,self).__init__ ( name  = name  ,
                                                 xvar  = xvar  ,
                                                 sigma = sigma ,
                                                 mean  = mean  ,
                                                 fudge = fudge )
        
        ## asymmetry parameter 
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa , 
                                       'kappa_%s'   % self.name ,
                                       '#kappa(%s)' % self.name ,
                                       None  , 0 , -1 , +1 ) 
        ## tail 
        self.__gamma = self.make_var ( gamma               ,
                                       'gamma_%s'   % name ,
                                       '#gamma(%s)' % name ,
                                       None , 1 , 1.e-5 , 100 ) 
        
        if kappa is None :
            
            self.__alpha = self.gamma
            self.__beta  = self.gamma
            
        else :
            
            self.__alpha , self.__beta = self.vars_from_asymmetry (
                self.gamma                                       , ## mean/average tail
                self.kappa                                       , ## asymmetry parametet
                v1name  =  self.roo_name ( 'alpha' , self.name ) ,
                v2name  =  self.roo_name ( 'beta'  , self.name ) ,
                v1title = '#alpha: #gamma #times (1+#kappa)'     , 
                v2title = '#beta : #gamma #times (1-#kappa)'     )
                      
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.GenLogisticIV (
            self.roo_name ( 'rgl4_' ) , 
            "Resolution GenLogisticIV %s" % self.name ,
            self.xvar      ,
            self.mean      ,
            self.sigma     , 
            self.alpha     , 
            self.beta      ) 

        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            'gamma'     : self.gamma ,
            'kappa'     : self.kappa ,
            'fudge'     : self.fudge ,
            }

    @property
    def mu ( self ) :
        """'mu' : location parameter (same as 'mean')"""
        return self.__mu
    @mu.setter
    def mu ( self , value ) :
        self.mean = value
    
    @property
    def kappa ( self ) :
        """'kappa' : dimensionless parameter, related to asymmetry"""
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.set_value ( self.__kappa , value )

    @property
    def gamma ( self ) :
        """'gamma' : mean tail-parameter: 0.5 * alpha + beta )"""
        return self.__gamma
    @gamma.setter
    def gamma ( self , value ) :
        self.set_value ( self.__gamma , value )

    @property
    def alpha ( self ) :
        """`alpha'- parameter for Generalized Logistic Type IV distribution
        """
        return self.__alpha

    @property
    def beta ( self ) :
        """`beta'- parameter for Generalized Logistic Type IV distribution
        """
        return self.__beta

# =============================================================================
## @class ResoHyperbolic
#  @see Ostap::Math::Hyperbolic
#  @see Ostap::Models::Hyperbolic
class ResoHyperbolic(RESOLUTION) :
    """ Symmetric Hyperbolic distribution
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
                                       '#zeta(%s)'  % name ,
                                       None , 1 , 1.e-10  , +100 )
        
        ## mu 
        self.__mu    = self.mean 

        ## parameter kappa - asymmetry
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa  ,
                                       'kappa_%s'    % name ,
                                       '#kappa(%s)'  % name ,
                                       None , 0 , -200 , +200 ) 
        
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
        """'mu' : location parameter, same as 'mean')"""
        return self.__mu

    @property 
    def zeta  ( self ) :
        """'zeta' : dimensioneless parameter, related to kurtosis"""
        return self.__zeta
    @zeta.setter  
    def zeta ( self , value ) :
        self.set_value ( self.__zeta , value )
    
    @property
    def kappa ( self ) :
        """'kappa' : dimensionless parameter, related to asymmetry"""
        return self.__kappa
    @kappa.setter  
    def kappa ( self , value ) :
        self.set_value ( self.__kappa, value )


    @property
    def alpha ( self ) :
        """'alpha' : value of canonical parameter 'alpha' """
        self.pdf.setPars ()
        return self.pdf.function().alpha ()

    @property
    def beta ( self ) :
        """'beta' : value of canonical parameter 'beta'"""
        self.pdf.setPars ()
        return self.pdf.function().beta ()

    @property
    def gamma ( self ) :
        """'gamma' : value of canonical parameter 'gamma'"""
        self.pdf.setPars ()
        return self.pdf.function().gamma ()
    
    @property
    def delta ( self ) :
        """'delta' : value of canonical parameter 'delta'"""
        self.pdf.setPars ()
        return self.pdf.function().delta ()

    @property
    def nominal_mean ( self ) :
        """'nominal_mean' : the actual mean of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().mean  ()

    @property
    def nominal_mode ( self ) :
        """'nominal_mode' : the actual mode of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().mode ()
    
    @property
    def nominal_variance ( self ) :
        """'nominal_variance' : actual variance of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().variance ()

    @property
    def nominal_rms ( self ) :
        """'nominal_rms' : actual RMS of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().rms ()


# =============================================================================
## @class ResoGenHyperbolic
#  @see Ostap::Math::GenHyperbolic
#  @see Ostap::Models::GenHyperbolic
class ResoGenHyperbolic(RESOLUTION) :
    """ Symmetric generalised Hyperbolic distribution
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
                                       '#zeta(%s)'  % name ,
                                       None , 1 , 1.e-10 , +100 ) 
        
        ## parameter kappa
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa ,
                                       'kappa_%s'    % name ,
                                       '#kappa(%s)'  % name ,
                                       None , 0 , -200 , +200 ) 
        
        ## lambda 
        self.__lambda = self.make_var ( lambd               ,
                                        'lambda_%s'   % name ,
                                        '#lambda(%s)' % name ,
                                        True ,  -2 , -200 , 200 )
        
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
        """'mu' : location parameter (same as 'mean')"""
        return self.__mu

    @property 
    def zeta  ( self ) :
        """'zeta' : dimensioneless parameter, related to shape"""
        return self.__zeta
    @zeta.setter  
    def zeta ( self , value ) :
        self.set_value ( self.__zeta , value )
    
    @property
    def kappa ( self ) :
        """'kappa' : dimensionless parameter, related to asymmetry"""
        return self.__kappa
    @kappa.setter  
    def kappa ( self , value ) :
        self.set_value ( self.__kappa, value )
    
    @property
    def lambd ( self ) :
        """'lambd' : dimensionless parameter, related to shape """
        return self.__lambda
    @lambd.setter
    def lambd ( self , value ) :    
        self.set_value ( self.__lambd , value )
        
    @property
    def alpha ( self ) :
        """'alpha' : value of canonical parameter 'alpha'"""
        self.pdf.setPars ()
        return self.pdf.function().alpha ()

    @property
    def beta ( self ) :
        """'beta' : value of canonical parameter 'beta'"""
        self.pdf.setPars ()
        return self.pdf.function().beta ()

    @property
    def gamma ( self ) :
        """'gamma' : value of canonical parameter 'gamma'"""
        self.pdf.setPars ()
        return self.pdf.function().gamma ()
    
    @property
    def delta ( self ) :
        """'delta' : value of canonical parameter 'delta'"""
        self.pdf.setPars ()
        return self.pdf.function().delta ()

    @property
    def nominal_mean ( self ) :
        """'nominal_mean' : the actual mean of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().mean  ()

    @property
    def nominal_mode ( self ) :
        """'nominal_mode' : actual mode of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().mode ()
    
    @property
    def nominal_variance ( self ) :
        """'nominal_variance' : actual variance of distribution"""
        self.pdf.setPars ()
        return self.pdf.function().variance ()

    @property
    def nominal_rms ( self ) :
        """'nominal_rms' : actual RMS of distribution"""
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
class ResoHypatia(ResoGenHyperbolic) : 
    r""" Variant of Hypatia pdf
    Convolution of Generalized Hyperbolic distrobution with 'offset'
    Gaussian distribution
    
    - see D. Martinez Santos, F. Duipertois,
    'Mass distributions marginalized over per-event errors',
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
                   sigma0    = None ,   ## width of the 'offset" Gaussian
                   fudge     =  1   ,   ## fudge factor 
                   mu        = None ,   ## related to mean
                   kappa     = None ,   ## related to asymmetry
                   cnvpars   = {}   ) : ## convolution parameters

        ## initialize the base 
        super(ResoHypatia,self).__init__ ( name       = name  ,
                                           xvar       = xvar  ,
                                           sigma      = sigma ,
                                           zeta       = zeta  ,
                                           fudge      = fudge ,
                                           mu         = mu    ,
                                           kappa      = kappa )

        
        ## safe created generic PDF 
        self.__genhyp = Generic1D_pdf ( self.pdf , xvar = self.xvar , name = self.new_name ( 'GenHyp' ) )
        
        ## prepare FFT convolution
        from ostap.fitting.resolution import ResoGauss 
        gname = self.generate_name ( 'gauss' , suffix = 'offset' ) 
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
        """'genhyp': get underlying generalized hyperbolic PDF"""
        return self.__genhyp
    
    @property
    def convolved ( self ) :
        """'convolved' : get PDF as convolution"""
        return self.__convolved
    
    @property
    def sigma0    ( self ) :
        """'sigma0' : width for the 'offset' Gaussian function"""
        return self.__resolution.sigma
    @sigma0.setter
    def sigma0    ( self , value ) :
        self.__resolution.sigma = value
    
    @property
    def cnvpars ( self ) :
        """'cnvpars' : parameters for convolution"""
        return self.__cnvpars 

# =============================================================================
## @class ResoDas
#  @see Ostap::Math::Das
#  @see Ostap::Models::Das
class ResoDas(RESOLUTION) :
    """ Das resoltuon model: gaussian with exponential tails 
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
                                       'k(%s)'     % name   ,
                                       None , 1 , 1.e-6 , 1000 )
        
        ## parameter kappa  
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa ,
                                       'kappa_%s'    % name ,
                                       '#kappa(%s)'  % name ,
                                       None , 0 , -1 , +1 ) 
        
        if kappa  is None or self.__kappa is ZERO :
            self.__AV_K = self.asymmetry_vars ( 'k' , var1    = self.__k , var2 = self.__k )
        else :
            self.__AV_K = self.asymmetry_vars ( 'k' , halfsum = self.__k , kappa = self.__kappa  )
            
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
        """'mu' : location parameter (same as 'mean')"""
        return self.__mu
    @mu.setter
    def mu ( self , value ) :
        self.mean = value

    @property
    def k ( self ) :
        """'k' : tail parameter (exponential slope)"""
        return self.__k
    @k.setter
    def k ( self , value ) :
        self.set_value ( self.__k , value ) 
        
    @property
    def kappa ( self ) :
        """'kappa' : dimensionless parameter, related to asymmetry"""
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.set_value ( self.__kappa , value )

    @property
    def kL  ( self ) :
        """'kL' : left tail parameter"""
        return self.__AV_K.var1 

    @property
    def kR  ( self ) :
        """'kR' : right tail parameter"""
        return self.__AV_K.var2 

# =============================================================================
## @class ResoGenGaussV1
#  Simple class that implements the generalized normal distribution v1
#  @see http://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1
#
#  Known also as the exponential power distribution, or the generalized error distribution,
#  this is a parametric family of symmetric distributions.
#  It includes all normal and Laplace distributions, and as limiting cases it includes all
#  continuous uniform distributions on bounded intervals of the real line.
#
#  This family includes the normal distribution when beta=2
#  (with mean mu and variance alpha^2/2)
#  and it includes the Laplace distribution when beta=1
#  As beta->inf, the density converges pointwise to a uniform density on (mu-alpha,mu+alpha)
# 
#  This family allows for tails that are either heavier than normal (when beta<2)
#  or lighter than normal (when beta>2).
#  It is a useful way to parametrize a continuum of symmetric, platykurtic densities
#  spanning from the normal (beta=2) to the uniform density (beta=inf),
#  and a continuum of symmetric, leptokurtic densities spanning from the Laplace
#  (beta=1) to the normal density (beta=2).
#    
#  Parameters:
#  - mu         : location/mean  
#  - alpha > 0  : scale 
#  - beta  > 0  : shape   (beta=2 corresponds to Gaussian distribution)
#
#  @see Ostap::Models::GenGaussV1 
#  @see Ostap::Math::GenGaussV1
#  @see GenGaussV1_pdf 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
#  @see M. T. Subbotin, “On the Law of Frequency of Error”, Mat. Sb., 31:2 (1923), 296–301
#  @see http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=sm&paperid=6854&option_lang=eng
#  @see Nadarajah, Saralees (September 2005). "A generalized normal distribution".
#       Journal of Applied Statistics. 32 (7): 685–694. doi:10.1080/02664760500079464.
#  @see https://doi.org/10.1080%2F02664760500079464
class ResoGenGaussV1(RESOLUTION) :
    """ Generalized Normal distribution v1
    see http://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_1

    Known also as the exponential power distribution, or the generalized error distribution,
    this is a parametric family of symmetric distributions.
    It includes all normal and Laplace distributions, and as limiting cases it includes all
    continuous uniform distributions on bounded intervals of the real line.

    This family includes the normal distribution when beta=2
    (with mean mu and variance alpha^2/2)
    and it includes the Laplace distribution when beta=1
    As beta->inf, the density converges pointwise to a uniform density on (mu-alpha,mu+alpha)
    
    This family allows for tails that are either heavier than normal (when beta<2)
    or lighter than normal (when beta>2).
    It is a useful way to parametrize a continuum of symmetric, platykurtic densities
    spanning from the normal (beta=2) to the uniform density (beta=inf),
    and a continuum of symmetric, leptokurtic densities spanning from the Laplace
    (beta=1) to the normal density (beta=2).
    
    - see M. T. Subbotin, 'On the Law of Frequency of Error', Mat. Sb., 31:2 (1923), 296–301
    - see http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=sm&paperid=6854&option_lang=eng
    - see Nadarajah, Saralees (September 2005). 'A generalized normal distribution'.
    Journal of Applied Statistics. 32 (7): 685–694. doi:10.1080/02664760500079464.
    - see https://doi.org/10.1080%2F02664760500079464
    
    Parameters:
    - mu         : location/mean  
    - alpha > 0  : scale 
    - beta  > 0  : shape   (beta=2 corresponds to Gaussian distribution)
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   alpha     = None ,
                   beta      = 2    ,  ## beta=2 is gaussian distribution
                   fudge     = 1    , 
                   mean      = None ) :

        ## initialize the base
        super(ResoGenGaussV1,self).__init__( name        = name  ,
                                             xvar        = xvar  ,
                                             sigma       = alpha ,                                     
                                             mean        = mean  ,
                                             fudge       = fudge ,
                                             sigma_name  = 'alpha_%s'   % name , 
                                             sigma_title = '#alpha(%s)' % name ) 
        
        self.__alpha = self.sigma         
        self.__beta  = self.make_var ( beta ,
                                       'beta_%s'        % name ,
                                       '#beta_{v1}(%s)' % name ,  
                                       None , 2 , 0.5 , 20     ) 
        #
        ## finally build PDF
        #
        self.pdf = Ostap.Models.GenGaussV1 (
            self.roo_name ( 'gaussgv1_' ) , 
            "Resolution generalized Gauss-V1 %s" % self.name ,
            self.xvar       ,
            self.mean       ,
            self.alpha_corr , ## attention 
            self.beta       )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'alpha'     : self.alpha ,            
            'beta'      : self.beta  ,
            'fudge'     : self.fudge 
            }


    @property
    def alpha ( self ) :
        """'alpha'-parameter for Generalized V1 Gaussian (the same as 'sigma')"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        self.set_value ( self.__alpha , value ) 

    @property
    def alpha_corr ( self ) :
        """Corrected 'alpha'-parameter for Generalized V1 Gaussian (the same as 'sigma_corr')"""
        return self.sigma_corr

    @property
    def beta ( self ) :
        """'beta'-parameter for Generalized  V1 Gaussian"""
        return self.__beta
    @beta.setter
    def beta ( self, value ) :
        self.set_value ( self.__beta , value )

# =============================================================================
## @class ResoNormalLaplace
#  Distribution for a sum of Gaussian and (asymmertric) Laplace variables 
#  It behaves like core Gaussian with exponential tails 
#  @see Wiliam J. Reed, "The Normal-Laplace Distribution Relatives", 
#  October, 2004
#  @see https://www.math.uvic.ca/faculty/reed/NL.draft.1.pdf
#  @see Reed, W.J, "The Normal-Laplace Distribution and Its Relatives". 
#        In: Balakrishnan, N., Sarabia, J.M., Castillo, E. (eds) 
#       "Advances in Distribution Theory, Order Statistics, and Inference. 
#       Statistics for Industry and Technology". Birkhäuser Boston. 
#  @see https://doi.org/10.1007/0-8176-4487-3_4
#
#  \f$ f(x; \mu, \sigma, k_L , k_R ) = 
#   \frac{1}{\sigma ( k_L + k_R) } 
#   \phi ( z ) \left( R ( \frac{1}{k_R} - z ) + 
#                     R ( \frac{1}{k_L} + z ) \right) 
#   \f$, where
#   - \f$ k_L,k_R \ge 0 \f$ 
#   - \f$ z = \frac{x-\mu}{\sigma} \f$ 
#   - \f$ \phi(z) \f$ is Gaussian PDF  
#   - \f$  R(x)   \f$ is Mill's ratio 
#    
class ResoNormalLaplace(RESOLUTION) :
    """ Distribution for a sum of Gaussian and (asymmertric) Laplace variables 
    It behaves like core Gaussian with exponential tails 
    - see Wiliam J. Reed, "The Normal-Laplace Distribution Relatives", 
    October, 2004
    -  see https://www.math.uvic.ca/faculty/reed/NL.draft.1.pdf
    - see Reed, W.J, 'The Normal-Laplace Distribution and Its Relatives'. 
    In: Balakrishnan, N., Sarabia, J.M., Castillo, E. (eds) 
    'Advances in Distribution Theory, Order Statistics, and Inference. 
    Statistics for Industry and Technology'. Birkhäuser Boston. 
    - see https://doi.org/10.1007/0-8176-4487-3_4
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   varsigma  = None ,
                   k         = 0    ,  ## k = 0 gives gaussian Resolution
                   fudge     = 1    , 
                   mean      = None ,
                   kappa     = None ) : ## tail symmetry parameter 

        ## initialize the base
        super(ResoNormalLaplace,self).__init__( name        = name     ,
                                                xvar        = xvar     ,
                                                sigma       = varsigma ,                                     
                                                mean        = mean     ,
                                                fudge       = fudge    ,
                                                sigma_name  = 'varsigma_%s'   % name ,                                               
                                                sigma_title = '#varsigma(%s)' % name )
        
        self.__varsigma = self.sigma

        ## k-parameter 
        self.__k     = self.make_var ( k ,
                                       'k_%s'  % self.name ,
                                       'k(%s)' % self.name ,
                                       None  , 1 , 0 , +100 ) 
        
        ## asymmetry parameter 
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa , 
                                       'kappa_%s'   % self.name ,
                                       '#kappa(%s)' % self.name ,
                                       None , 0 , -1 , +1 ) 
                
        if kappa is None or self.__kappa is ZERO : 
            self.__AV_K = self.asymmetry_vars ( 'k' , var1    = self.__k , var2  = self.__k     )
        else :
            self.__AV_K = self.asymmetry_vars ( 'k' , halfsum = self.__k , kappa = self.__kappa )
  
        #
        ## finally build PDF
        #
        self.pdf = Ostap.Models.NormalLaplace (
            self.roo_name ( 'normlapl_' ) , 
            "Resolution Normal Laplace %s" % self.name ,
            self.xvar       ,
            self.mean       ,
            self.sigma_corr , ## attention 
            self.kL         ,
            self.kR         )
        
        ## save the configuration
        self.config = {
            'name'      : self.name     ,
            'xvar'      : self.xvar     ,
            'mean'      : self.mean     ,
            'varsigma'  : self.varsigma ,
            'k'         : self.k        ,
            'kappa'     : self.kappa    ,
            'fudge'     : self.fudge 
            }
        
    @property
    def varsigma    ( self ) :
        """'varsigma' : varsigma parameter for Normal Laplace function (same as 'sigma')"""
        return self.sigma
    @varsigma.setter
    def varsigma    ( self , value ) :
        self.set_value ( self.__varsigma , value )

    @property 
    def k  ( self ) :
        """'k' :  (dimensioneless) k-parameter"""
        return self.__k
    @k.setter  
    def k ( self , value ) :
        self.set_value ( self.__k , value )

    @property 
    def kL ( self ) :
        """'kL' :  (dimensioneless) kL-parameter"""
        return self.__AV_K.var1 
        
    @property 
    def kR ( self ) :
        """'kR' :  (dimensioneless) kR-parameter"""
        return self.__AV_K.var2 

    @property
    def kappa ( self ) :
        """'kappa' : tail asymmetry parameter 0.5(kL-kR)/(kL+_kR)
        """
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.set_value ( self.__kappa , value )


# =============================================================================
## @class ResoBukin2
#  Resolution function based on Buni2_pdf: sum of two
#  exponental modified Gaussina functions
class ResoBukin2(RESOLUTION) :
    """ Resolution function based on Buni2_pdf: sum of two
    exponental modified Gaussina functions
    """
    def __init__ ( self            ,
                   name            ,
                   xvar            ,
                   varsigma = None ,
                   kA       = 0    ,  ## k = 0 gives gaussian Resolution
                   kB       = None ,  ## k = 0 gives gaussian Resolution
                   phi      = None , 
                   fudge    = 1    , 
                   mean     = None ,
                   kappa    = None ) : ## sigma symmetry parameter 
        
        ## initialize the base
        super(ResoBukin2,self).__init__( name        = name     ,
                                         xvar        = xvar     ,
                                         sigma       = varsigma ,                                     
                                         mean        = mean     ,
                                         fudge       = fudge    ,
                                         sigma_name  = 'varsigma_%s'   % name ,                                         
                                         sigma_title = '#varsigma(%s)' % name )
        
        ## asymmetry parameter 
        self.__kappa = self.make_var ( ZERO if kappa is None else kappa , 
                                       'kappa_%s'   % self.name ,
                                       '#kappa(%s)' % self.name ,
                                       None , 0 , -1 , +1 )

        
        if kappa is None or self.__kappa is ZERO  :
            self.__AV_S = self.asymmetry_vars ( 'varsigma' ,
                                                var1 = self.sigma_corr    ,
                                                var2 = self.sigma_corr    , left_right = "AB" )
        else :
            self.__AV_S = self.asymmetry_vars ( 'varsigma' ,
                                                halfsum = self.sigma_corr ,
                                                kappa   = self.__kappa    , left_right = "AB" )
            
        ## k-parameters 
        if   kA is None and kB is None :
                
            self.__kA , self.__kB = ZERO, ZERO 
            
        elif kA is None :
            
            self.__kB  = self.make_var      ( kB,
                                              'kB_%s'  % self.name ,
                                              'kB(%s)' % self.name ,
                                              None  , 0 , -1000 , +1000 ) 
            self.__kA  = self.vars_multiply ( self.kB , -1.0 ,
                                              'kA_%s'  % self.name ,
                                              'kA(%s)' % self.name )
        elif kB is None :
            
            self.__kA  = self.make_var      ( kA,
                                              'kA_%s'  % self.name ,
                                              'kA(%s)' % self.name ,
                                              None  , 0 , -1000 , +1000 ) 
            self.__kB  = self.vars_multiply ( self.kA , -1.0 ,
                                              'kB_%s'  % self.name ,
                                              'kB(%s)' % self.name )
            
        else :
            
            self.__kA  = self.make_var     ( kA,
                                             'kA_%s'  % self.name ,
                                             'kA(%s)' % self.name ,
                                             None  , -1.e-5 , -1000 , 0     ) 
            self.__kB  = self.make_var     ( kB,
                                             'kB_%s'  % self.name ,
                                             'kB(%s)' % self.name ,
                                             None  , +1.e-5 , 0     , +1000 ) 
            
        ## phi parameter 
        self.__phi = self.make_var ( ZERO if phi is None else phi , 
                                     'phi_%s'   % self.name ,
                                     '#phi(%s)' % self.name ,
                                     None , 0 , -4 * math.pi , +4 * math.pi )
        
        #
        ## finally build PDF
        #
        self.pdf = Ostap.Models.Bukin2 (
            self.roo_name ( 'bukin2_' ) , 
            "Resolution Bukin2 %s" % self.name ,
            self.xvar       ,
            self.mean       ,            
            self.varsigmaA  , ## attention
            self.varsigmaB  , ## attention
            self.kA         ,
            self.kB         ,
            self.phi        )
        
        ## save the configuration
        self.config = {
            'name'      : self.name     ,
            'xvar'      : self.xvar     ,
            'mean'      : self.mean     ,
            'varsigma'  : self.sigma    ,
            'kappa'     : self.kappa    ,
            'kA'        : self.kA       ,
            'kB'        : self.kB       ,
            'phi'       : self.phi      ,
            'fudge'     : self.fudge 
            }
        
    @property
    def varsigma    ( self ) :
        """'varsigma' : varsigma parameter for Bukin2 function (same as 'sigma')"""
        return self.sigma
    @varsigma.setter
    def varsigma    ( self , value ) :
        self.set_value ( self.__varsigma , value )

    @property
    def kappa ( self ) :
        """'kappa' : sigma asymmetry parameter
        """
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.set_value ( self.__kappa , value )

    @property
    def varsigmaA   ( self ) :
        """'varsigmaA' : varsigmaA parameter for Bukin2 function"""
        return self.__AV_S.var1

    @property
    def varsigmaB   ( self ) :
        """'varsigmaB' : varsigmaB parameter for Bukin2 function"""
        return self.__AV_S.var2  

    @property 
    def kA  ( self ) :
        """'kA' :  (dimensioneless) kA-parameter"""
        return self.__kA
    @kA.setter  
    def kA ( self , value ) :
        self.set_value ( self.__kA , value )

    @property 
    def kB  ( self ) :
        """'kB' :  (dimensioneless) kB-parameter"""
        return self.__kB
    @kB.setter  
    def kB ( self , value ) :
        self.set_value ( self.__kB , value )

    @property
    def phi ( self ) :
        """'phi' : phi parameter
        """
        return self.__phi
    @phi.setter
    def phi ( self , value ) :
        self.set_value ( self.__phi , value )

# =============================================================================
## @class ResoBateShape 
#  @see Ostap::Math::BatesShape 
#  @see Ostap::Math::Bates
#  @see Ostap::Math::IrwinHall
#  @see Ostap::Models::BatesShape
#  For good fits parameter n shoul be rather large. e.g.  n>10 
class ResoBatesShape(RESOLUTION) :
    """ BatesShape resolution model (spline-like finite model) 
    - see Ostap::Math::BatesShape
    - see Ostap::Math::Bates
    - see Ostap::Math::IrwinHall
    - see Ostap::Models::BatesShape 
    - see BatesShape_pdf 
    For good fits parameter n should be large. e.g.  n>10 
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   sigma     = None ,   ## related to sigma
                   n         = 2    ,   ## shape 
                   fudge     = 1    ,   ## fudge-parameter 
                   mean      = None ) : ## related to mean 
        
        ## initialize the base 
        super(ResoBatesShape,self).__init__ ( name  = name  ,
                                              xvar  = xvar  ,
                                              sigma = sigma ,
                                              mean  = mean  ,
                                              fudge = fudge )
        
        assert isinstance ( n , int ) and 1<= n, "Invalid parameter n!"
        self.__n      = int ( n )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.BatesShape (
            self.roo_name ( 'rbs_' ) , 
            "Resolution Bates-shape %s" % self.name ,
            self.xvar       ,
            self.mean       ,
            self.sigma_corr ,
            self.n          )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            'n'         : self.n     ,
            'fudge'     : self.fudge }

    @property
    def n ( self ) :
        """'n' : shape parameter"""
        return self.__n
        
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
