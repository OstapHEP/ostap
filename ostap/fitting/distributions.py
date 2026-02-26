#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/distributions.py
#  A set of various smooth shapes and PDFs 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
# =============================================================================
""" A set of various smooth shapes and PDFs

- GammaDist_pdf        : Gamma-distributuon in shape/scale parameterization
- GenGammaDist_pdf     : Generalized Gamma-distribution
- Amoroso_pdf          : Another view of generalized Gamma distribution
- LogGammaDist_pdf     : Gamma-distributuon in shape/scale parameterization
- Log10GammaDist_pdf   : Gamma-distributuon in shape/scale parameterization

- LogGamma_pdf         : Log-Gamma distribution

- Landau_pdf           : Landau distribution 
- Slash_pdf'           : Slash: symmetric bump with very heavy tails

- Argus_pdf            : ARGUS distribution 
- GenArgus_pdf         : Generalized ARGUS distribution 

- Beta_pdf             : Beta distribution
- BetaPrime_pdf        : Beta-prime distribution 
- GenBeta_pdf          : Generalized Bet distribution 
- GenBetaPrime_pdf     : generalized Beta-prime distribution 

- TwoExpos_pdf         : Difference of two exponents

- Gumbel_pdf           : Gumbel distributions
- Rice_pdf             : Rice distribution
- GenInvGauss_pdf      : Generalized Inverse Gaussian distribution
- Weibull_pdf          : Weibull distributions

- GenPareto_pdf        : Generalised Pareto distribution
- ExGenPareto_pdf      : Exponentiated Generalised Pareto distribution

- GEV_pdf              : Generalised Extreme Value distribution
- Benini_pdf           : Benini distribution
- MPERT_pdf            : Modified PERT distribution
- BirnbaumSaunders_pdf : Birnbaum-Saunders distribution
- Frechet_pdf          : Frechet distribution
- Dagum_pdf            : Dagum distribution
- BenktanderI_pdf      : Benktander Type I distribution
- BenktanderII_pdf     : Benktander Type 2 distribution
- LogNormal_pdf        : Log-normal distribution
- ExpoLog_pdf          : Expo-Log   distribution
- Davis_pdf            : Davis distribution
- Kumaraswami_pdf      : Kumaraswami distribution
- InverseGamma_pdf     : Inverse-Gamma distribution
- Burr_pdf             : Burr Type XII distribution

- Tsallis_pdf          : Tsallis PDF 
- QGSM_pdf             : QGSM PDF 
- Hagedorn_pdf         : Hagedorn PDF 
- Tsallis2_pdf         : 2D Tsallis PDF 
"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    #
    'GammaDist_pdf'        , ## Gamma-distributuon in shape/scale parameterization
    'GenGammaDist_pdf'     , ## Generalized Gamma-distribution
    'Amoroso_pdf'          , ## another view of generalized Gamma distribution
    'LogGammaDist_pdf'     , ## Gamma-distributuon in shape/scale parameterization
    'Log10GammaDist_pdf'   , ## Gamma-distributuon in shape/scale parameterization
    #
    'LogGamma_pdf'         , ## Log-Gamma distribution
    #
    'Landau_pdf'           , ## Landau distribution 
    'Slash_pdf'            , ## Slash: symmetric bump with very heavy tails 
    #
    'Argus_pdf'            , ## ARGUS distribution 
    'GenArgus_pdf'         , ## Generalized ARGUS distribution
    #
    'Beta_pdf'             , ## Beta distribution
    'BetaPrime_pdf'        , ## Beta-prime distribution 
    'GenBeta_pdf'          , ## generalized Beta distribution
    'GenBetaPrime_pdf'     , ## generalized Beta-prime distribution
    #
    'TwoExpos_pdf'         , ## difference of two exponents
    #
    'Gumbel_pdf'           , ## Gumbel distributions
    'Rice_pdf'             , ## Rice distribution 
    'GenInvGauss_pdf'      , ## Generalized Inverse Gaussian distribution
    'Weibull_pdf'          , ## Weibull distributions
    #
    'GenPareto_pdf'        , ## Generalised Pareto distribution
    'ExGenPareto_pdf'      , ## Exponentiated Generalised Pareto distribution
    #
    'GEV_pdf'              , ## Generalised Extreme Value distribution
    'Benini_pdf'           , ## Benini distribution
    'MPERT_pdf'            , ## Modified PERT distribution
    #
    'BirnbaumSaunders_pdf' , ## Birnbaum-Saunders distribution
    'Frechet_pdf'          , ## Frechet distribution
    'Dagum_pdf'            , ## Dagum distribution
    'BenktanderI_pdf'      , ## Benktander Type I distribution
    'BenktanderII_pdf'     , ## Benktander Type 2 distribution
    'LogNormal_pdf'        , ## Log-normal distribution
    'ExpoLog_pdf'          , ## Expo-Log   distribution
    'Davis_pdf'            , ## Davis distribution
    'Kumaraswami_pdf'      , ## Kumaraswami distribution
    'InverseGamma_pdf'     , ## Inverse-Gamma distribution
    'Burr_pdf'             , ## Burr Type XII distribution
    # 
    'Tsallis_pdf'          , ## Tsallis PDF 
    'QGSM_pdf'             , ## QGSM PDF 
    'Hagedorn_pdf'         , ## Hagedorn PDF 
    'Tsallis2_pdf'         , ## 2D Tsallis PDF
    #
    )
# =============================================================================
from   ostap.core.core          import Ostap, VE 
from   ostap.math.base          import isfinite, pos_infinity, neg_infinity 
from   ostap.fitting.pdfbasic   import PDF1, PDF2
from   ostap.fitting.fithelpers import ( ShiftAndScale , Shift , Scale ,
                                         AlphaAndBeta  , Alpha , Beta  ,
                                         P , Q , R , PQ ) 
import ROOT, math
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.distributions' )
else                       : logger = getLogger ( __name__                      )
# =============================================================================
models  = []
spectra = []
# =============================================================================
## @class GammaDist_pdf
#  Gamma-distribution with shape/scale parameters
#  http://en.wikipedia.org/wiki/Gamma_distribution
#  It suits nicely for fits of multiplicity and/or chi2 distributions
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::GammaDist 
#  @see Ostap::Math::GammaDist 
class GammaDist_pdf(PDF1,ShiftAndScale) :
    """ Gamma-distribution with shape/scale parameters
    http://en.wikipedia.org/wiki/Gamma_distribution
    It suits nicely for fits of multiplicity and/or, especially chi2 distributions
    
    In probability theory and statistics, the gamma distribution is a two-parameter
    family of continuous probability distributions.
    The common exponential distribution and chi-squared distribution are special
    cases of the gamma distribution.
    Here a shape parameter k and a scale parameter theta are used for parameterization
    
    The parameterization with k and theta appears to be more common in econometrics
    and certain other applied fields, where e.g. the gamma distribution is frequently
    used to model waiting times. For instance, in life testing, the waiting time
    until death is a random variable that is frequently modeled with a gamma distribution.

    If k is an integer, then the distribution represents an Erlang distribution;
    i.e., the sum of k independent exponentially distributed random variables,
    each of which has a mean of theta (which is equivalent to a rate parameter of 1/theta).

    If X ~ Gamma(nu/2, 2), then X is identical to chi2(nu),
    the chi-squared distribution with nu degrees of freedom.
    Conversely, if Q ~ chi2(nu) and c is a positive constant,
    then cQ ~ Gamma(nu/2, 2c).
    
    Parameters
    - k>0     : shape
    - theta>0 : scale 
    """
    ## constructor
    def __init__ ( self        , * , 
                   xvar        , ## the variable
                   name  = ''  , ## the name 
                   k     = 2   , ## k-parameter
                   theta = 1   , ## theta-parameter (scale)
                   shift = ROOT.RooFit.RooConst ( 0 ) ) : ## shift
        
        ## initiailze the 1st base 
        PDF1.__init__ ( self , name = name , xvar = xvar )
        ## Initialize the 2nd base
        ShiftAndScale.__init__ ( self ,
                                 shift       = shift        ,
                                 scale       = theta        ,
                                 scale_name  = 'theta_%s'   % self.name ,
                                 scale_title = '#theta(%s)' % self.name )
        
        ## shape 
        self.__k     = self.make_var ( k     ,
                                       'k_%s'                % self.name ,
                                       'k_{#Gamma}(%s)'      % self.name ,
                                       None , 2 , 1.e-3 , 100 )
        
        self.pdf  = Ostap.Models.GammaDist (
            self.roo_name ( 'gamma_' ) ,
            'Gamma distribution %s' % self.name , 
            self.x                 ,
            self.k                 ,
            self.theta             )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'k'     : self.k     ,
            'theta' : self.theta ,            
            'shift' : self.shift ,            
            }

    @property
    def k ( self ) :
        """`k'-parameter of Gamma distribution   (k>0)"""
        return self.__k
    @k.setter 
    def k ( self , value ) :
        self.set_value ( self.__k , value ) 

    ## ALIAS
    theta = ShiftAndScale.scale
    low   = ShiftAndScale.shift

    
models.append ( GammaDist_pdf )

# =============================================================================
## @class GenGammaDist_pdf 
#  Generalized Gamma-distribution with additional shift parameter 
#  http://en.wikipedia.org/wiki/Generalized_gamma_distribution
#  Special cases : 
#  - p == 1      : Gamma  distribution
#  - p == k      : Weibull distribution
#  - p == k == 1 : Exponential distribution
#  - p == k == 2 : Rayleigh    distribution
#  @see Ostap::Math::GenGammaDist
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Math::GenGammaDist 
#  @see Ostap::Models::GenGammaDist 
class GenGammaDist_pdf(GammaDist_pdf,P) :
    """ Generalized Gamma-distribution with additional shift parameter 
    http://en.wikipedia.org/wiki/Generalized_gamma_distribution
    Special cases : 
    - p == 1      : Gamma       distribution
    - p == k      : Weibull     distribution
    - p == k == 1 : Exponential distribution
    - p == k == 2 : Rayleigh    distribution
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       , ## the variable
                   name  = '' , ## the name 
                   k     = 2  , ## k-parameter
                   theta = 1  , ## theta-parameter, scale 
                   logp  = 0  , ## p-parameter
                   shift = ROOT.RooFit.RooConst ( 0 ) ) : ## shift-parameter
        
        ## Initialize the 1st base 
        GammaDist_pdf.__init__ ( self ,
                                 name  = name  ,
                                 xvar  = xvar  ,
                                 k     = k     ,
                                 theta = theta ,
                                 shift = shift )
        ## initiailze the second base
        P.__init__ ( self  , logp = logp )
        ## 
        self.pdf  = Ostap.Models.GenGammaDist (
            self.roo_name ( 'ggamma_' ) ,
            'Generalized Gamma distribution %s' % self.name , 
            self.x         ,
            self.k         ,
            self.theta     ,
            self.logp      , 
            self.shift     )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'k'     : self.k     ,
            'theta' : self.theta ,            
            'logp'  : self.logp  ,            
            'shift' : self.shift ,            
            }

models.append ( GenGammaDist_pdf ) 
# =============================================================================
## @class Amoroso_pdf
#  Another view on generalized gamma distribution
#  http://arxiv.org/pdf/1005.3274
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Math::Amoroso
#  @see Ostap::Models::Amoroso
class Amoroso_pdf(PDF1,ShiftAndScale,AlphaAndBeta) :
    """ Another view on generalized gamma distribution
    http://arxiv.org/pdf/1005.3274
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name 
                   theta = 1  ,   ## theta-parameter/ scale 
                   alpha = 1  ,   ## alpha-parameter
                   beta  = 1  ,   ## beta-parameter
                   a     = 0  ) : ## a-parameter/shift 
        
        ## The 1st base 
        PDF1         .__init__ ( self , name = name , xvar = xvar )
        ## The 2nd base
        ShiftAndScale.__init__ ( self ,
                                 scale       = theta ,
                                 shift       = a     ,
                                 scale_name  = 'theta_%s'             % self.name ,
                                 scale_title = '#theta_{Amoroso}(%s)' % self.name ,
                                 shift_name  = 'a_%s'                 % self.name ,
                                 shift_title = 'a_{Amoroso}(%s)'      % self.name )
        ## The 3rd base
        AlphaAndBeta .__init__  ( self          ,
                                  alpha = alpha ,
                                  beta  = beta  )
                                  
        self.pdf  = Ostap.Models.Amoroso (
            self.roo_name ( 'amo_' ) ,
            'Amoroso %s' % self.name , 
            self.x         ,
            self.theta     ,
            self.alpha     ,
            self.beta      ,
            self.a         )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'theta' : self.theta ,            
            'alpha' : self.alpha ,
            'beta'  : self.beta  ,            
            'a'     : self.a     ,            
            }
        
    ## ALIAS
    theta = ShiftAndScale.scale 
    a     = ShiftAndScale.shift


models.append ( Amoroso_pdf ) 
# =============================================================================
## @class LogGammaDist_pdf
#  Distribution for log(x), where x follows Gamma distribution
#  It suits nicely for fits of log(multiplicity) and/or log(chi2) distributions
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::LogGammaDist 
#  @see Ostap::Math::LogGammaDist 
#  @see Ostap::Models::GammaDist 
#  @see Ostap::Math::GammaDist 
class LogGammaDist_pdf(GammaDist_pdf) :
    """ Distribution for log(x), where x follows Gamma distribution
    It suits nicely for fits of log(multiplicity) and/or log(chi2) distributions
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       , ## the variable
                   name  = '' , ## the name 
                   k     = 2  , ## k-parameter
                   theta = 1  , ## theta-parameter
                   shift = ROOT.RooFit.RooConst ( 0 ) ) : ## shift-parameter
        #
        GammaDist_pdf.__init__ ( self  ,
                                 name  = name  ,
                                 xvar  = xvar  ,
                                 k     = k     ,
                                 theta = theta ,
                                 shift = shift )
        #
        self.pdf  = Ostap.Models.LogGammaDist (
            self.roo_name ( 'lgamma_' ) ,
            'Log-Gamma %s' % self.name , 
            self.x                 ,
            self.k                 ,
            self.theta             )

        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'k'     : self.k     ,            
            'theta' : self.theta ,            
            'shift' : self.shift ,            
            }

models.append ( LogGammaDist_pdf ) 
# =============================================================================
## @class Log10GammaDist_pdf
#  Distribution for log10(x), where x follows Gamma distribution
#  It suits nicely for fits of log10(multiplicity) and/or log10(chi2) distributions
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::Log10GammaDist 
#  @see Ostap::Math::Log10GammaDist 
#  @see Ostap::Models::LogGammaDist 
#  @see Ostap::Math::LogGammaDist 
#  @see Ostap::Models::GammaDist 
#  @see Ostap::Math::GammaDist 
class Log10GammaDist_pdf(LogGammaDist_pdf) :
    """ Distribution for log10(x), where x follows Gamma distribution
    It suits nicely for fits of log10(multiplicity) and/or log10(chi2) distributions
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       , 
                   name  = '' , ## the name 
                   k     = 2  , ## k-parameter
                   theta = 1  , ## theta-parameter
                   shift = ROOT.RooFit.RooConst ( 0 ) ) : ## shift-parameter
        #
        LogGammaDist_pdf.__init__ ( self ,
                                    name  = name  ,
                                    xvar  = xvar  ,
                                    k     = k     ,
                                    theta = theta ,
                                    shift = shift )
        ## 
        self.pdf  = Ostap.Models.Log10GammaDist (
            self.roo_name ( 'l10gamma_' ) ,
            'Log10 Gamma %s' % self.name , 
            self.x                 ,
            self.k                 ,
            self.theta             )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'k'     : self.k     ,            
            'theta' : self.theta ,            
            'shift' : self.shift ,            
            }
    
models.append ( Log10GammaDist_pdf ) 
# =============================================================================
## @class LogGamma_pdf
#  - http://arxiv.org/pdf/1005.3274
#  - Prentice, R. L. (1974). A log gamma model and its maximum likelihood
#                            estimation. Biometrika 61, 539
#  - Johnson, N. L., Kotz, S., and Balakrishnan, N. (1995). Continuous
#            univariate distributions, 2nd ed. Vol. 2. Wiley, New York.
#  - Bartlett, M. S. and G., K. M. (1946). The statistical analysis of
#                  variance-heterogeneity and the logarithmic transformation. 
#                 J. Roy. Statist. Soc. Suppl. 8, 1, 128.
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::LogGamma
#  @see Ostap::Math::LogGamma
class LogGamma_pdf(PDF1,ShiftAndScale) :
    """ Log-Gamma distribution
    - http://arxiv.org/pdf/1005.3274
    - Prentice, R. L. (1974). A log gamma model and its maximum likelihood
    estimation. Biometrika 61, 539
    - Johnson, N. L., Kotz, S., and Balakrishnan, N. (1995). Continuous
    univariate distributions, 2nd ed. Vol. 2. Wiley, New York.
    - Bartlett, M. S. and G., K. M. (1946). The statistical analysis of
    variance-heterogeneity and the logarithmic transformation. 
    J. Roy. Statist. Soc. Suppl. 8, 1, 128.
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name 
                   nu    = 0  ,   ## nu-parameter (shift) 
                   lambd = 1  ,   ## lambda-parameter (scale) 
                   alpha = 1  ) : ## alpha-parameter  (shape) 
        
        ## The first base 
        PDF1         .__init__ ( self , name = name , xvar = xvar )
        ## The second base 
        ShiftAndScale.__init__ ( self          ,
                                 shift = nu    ,
                                 scale = lambd ,
                                 shift_name  = 'nu_%s'                    % self.name ,
                                 shift_title = '#nu_{#log#Gamma}(%s)'     % self.name ,
                                 scale_name  = 'lambda_%s'                % self.name ,
                                 scale_title = '#lambda_{#log#Gamma}(%s)' % self.name )
        ## 
        self.__alpha  = self.make_var ( alpha    ,
                                        'alpha_%s'                 % name ,
                                        '#alpha_{#log#Gamma}(%s)'  % name ,
                                        None , 1 , 1.e-3 , 1000 )
        
            
        self.pdf  = Ostap.Models.LogGamma (
            self.roo_name ( 'loggamma_' ) ,
            'Log-Gamma %s' % self.name    , 
            self.x     ,
            self.nu    ,
            self.lambd ,
            self.alpha )

        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'nu'    : self.nu    ,            
            'lambd' : self.lambd ,            
            'alpha' : self.alpha ,            
            }

    ## ALIASES 
    nu     = ShiftAndScale.shift
    lambd  = ShiftAndScale.scale 
    Lambda = ShiftAndScale.scale 

    @property
    def alpha  ( self ) :
        """`alpha'-parameter of log-Gamma distribution   (alpha>0)"""
        return self.__alpha
    @alpha.setter 
    def alpha ( self , value ) :
        self.set_value ( self.__alpha , value )

models.append ( LogGamma_pdf )

# =============================================================================
## @class Landau_pdf
#  http://en.wikipedia.org/wiki/Landau_distribution
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::Landau
#  @see Ostap::Math::Landau
class Landau_pdf(PDF1,ShiftAndScale) :
    """ Landau distribution 
    - http://en.wikipedia.org/wiki/Landau_distribution
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name 
                   scale = 1  ,   ## scale-parameter 
                   shift = 0  ) : ## shift-parameter 
        #
        PDF1          .__init__ ( self  , name = name , xvar = xvar )
        ShiftAndScale .__init__ ( self  ,
                                  scale = scale ,
                                  shift = shift ,
                                  shift_name  = 'mu_%s'            % self.name ,
                                  shift_title = '#mu_{Landau}(%s)' % self.name ,
                                  scale_name  = 'c_%s'             % self.name ,
                                  scale_title = 'c_{Landau}(%s)'   % self.name )
        
        self.pdf  = Ostap.Models.Landau (
            self.roo_name ( 'landau_' )   ,
            'Landau %s' % self.name  , 
            self.x     ,
            self.scale ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'scale' : self.scale ,            
            'shift' : self.shift ,            
            }
        
    ## ALIASES
    mu = ShiftAndScale.shift
    c  = ShiftAndScale.scale 
    
models.append ( Landau_pdf )

# =============================================================================
## @class Slash_pdf
#  Symmetric function with very heavy tails.
#  @see https://en.wikipedia.org/wiki/Slash_distribution
#  The tails  are so heavy that moments does not exists
#  @see Ostap::Math::Slash
#  @see Ostap::Models::Slash
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-02-004
class Slash_pdf(PDF1,ShiftAndScale) :        
    """ Symmetric function with very heavy tails.
    - see https://en.wikipedia.org/wiki/Slash_distribution
    The tails  are so heavy that moments does not exists
    - see Ostap::Math::Slash
    - see Ostap::Models::Slash
    """
    def __init__ ( self         , * , 
                   xvar         ,   ## variable 
                   name  = ''   ,   ## name 
                   mu    = 0    ,   ## related to location 
                   scale = 1    ) : ## related to scale
        
        ## initialize the first base
        PDF1.__init__ ( self , name = name , xvar = xvar )
        ## initialize the second base
        ShiftAndScale.__init__ ( self ,
                                 shift       = mu                  , 
                                 scale       = scale               ,
                                 shift_name  = 'mu_%s'      % self.name ,
                                 shift_title = '#mu(%s)'    % self.name )
        
        ## finally build pdf
        self.pdf = Ostap.Models.Slash (
            self.roo_name ( 'slash_' ) , 
            "Slash %s" % self.name ,
            self.xvar      ,
            self.mu        ,
            self.scale     )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mu'        : self.mu    ,
            'scale'     : self.scale ,
            }

    ## ALIAS
    mu = ShiftAndScale.shift
    
models.append ( Slash_pdf )      

# =============================================================================
## @class Argus_pdf
#  http://en.wikipedia.org/wiki/ARGUS_distribution
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::Argus
#  @see Ostap::Math::Argus
class Argus_pdf(PDF1,ShiftAndScale) :
    """ Argus distribution
    - http://en.wikipedia.org/wiki/ARGUS_distribution
    - support:   mu-c < x < mu
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       , ## the variable
                   name  = '' , ## the name 
                   chi   = 1  , ## parameter chi 
                   c     = 1  , ## (m
                   mu    = 1  ) :
        
        ## Initiailze the 1st base 
        PDF1         .__init__ ( self , name = name , xvar = xvar )
        ShiftAndScale.__init__ ( self ,
                                 shift       = mu ,
                                 scale       = c  , 
                                 scale_name  = 'c_%s'            % self.name ,
                                 scale_title = 'c_{Argus}(%s)'   % self.name ,
                                 shift_name  = 'mu_%s'           % self.name ,
                                 shift_title = '#mu_{Argus}(%s)' % self.name )

        ## shape 
        self.__chi  = self.make_var ( chi     ,
                                      'chi_%s'           % name ,
                                      '#chi_{Argus}(%s)' % name ,
                                      None , 1 , 1.e-6 , 20  )

        ## create PDF 
        self.pdf  = Ostap.Models.Argus (
            self.roo_name ( 'argus_' )   ,
            'ARGUS %s' % self.name  , 
            self.x     ,
            self.mu    ,
            self.c     ,
            self.chi   )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mu'    : self.mu    ,            
            'c'     : self.c     ,            
            'chi'   : self.chi   ,            
            }

    ## ALIAS
    mu = ShiftAndScale.shift
    c  = ShiftAndScale.scale 
    
    @property
    def chi ( self ) :
        """`chi'-parameter of Argus distribution"""
        return self.__chi
    @chi.setter 
    def chi ( self , value ) :
        self.set_value ( self.__chi , value )

models.append ( Argus_pdf ) 

# =============================================================================
## @class GenArgus_pdf
#  http://en.wikipedia.org/wiki/ARGUS_distribution
#  Generalised Argus 
#  - support: \f$ \mu-c < x < \mu \f$
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::GenArgus
#  @see Ostap::Math::GenArgus
class GenArgus_pdf(Argus_pdf) :
    """ Generalized Argus distribution
    - http://en.wikipedia.org/wiki/ARGUS_distribution
    - support m mu-c < x < mu 
    """
    ## constructor
    def __init__ ( self        , * , 
                   xvar        , ## the variable
                   name = ''   , ## the name 
                   chi  = 1    ,
                   c    = 1.0  ,
                   dp   = 1.5  , ## corresponds to regular Argus 
                   mu   = 1    ) :
        
        ## initialize the base 
        Argus_pdf.__init__ ( self        ,
                             name = name ,
                             xvar = xvar ,
                             c    = c    ,
                             chi  = chi  ,
                             mu   = mu   )         
        #
        self.__dp  = self.make_var ( dp     ,
                                     'dp_%s'               % self.name ,
                                     '#deltap_{Argus}(%s)' % self.name ,
                                     None , 1.5 , 1.e-5 , 20 )
        
        ## create PDF 
        self.pdf  = Ostap.Models.GenArgus (
            self.roo_name ( 'genargus_' )   ,
            'Generalized ARGUS %s' % self.name  , 
            self.x     ,
            self.mu    ,
            self.c     ,
            self.chi   ,
            self.dp    )
                    
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mu'    : self.mu    ,            
            'c'     : self.c     ,            
            'chi'   : self.chi   ,            
            'dp'    : self.dp    ,            
            }

    @property
    def dp ( self ) :
        """`dp'-parameter of Generalized Argus distribution (p=dp-1)"""
        return self.__dp
    @dp.setter 
    def dp ( self , value ) :
        self.set_value ( self.__dp , value )

models.append ( GenArgus_pdf ) 

# =============================================================================
## @class Beta_pdf
#  Beta distribution 
#  @see https://en.wikipedia.org/wiki/Beta_distribution
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2025-08-10
#  @see Ostap::Models::Beta
#  @see Ostap::Math::Beta
class Beta_pdf(PDF1,ShiftAndScale,PQ) :
    """ Beta distribution 
    - see https://en.wikipedia.org/wiki/Beta_distribution
    with scale and shift 
    """
    ## constructor
    def __init__ ( self        , * , 
                   xvar        , ## the variable
                   name   = '' , ## the name 
                   logp   = 0  ,
                   logq   = 0  ,
                   scale  = 1  ,
                   shift  = 0  ) :
        
        # =====================================================================
        # Initiailze the base
        PDF1         .__init__ ( self , name  = name  , xvar  = xvar  )
        ShiftAndScale.__init__ ( self , scale = scale , shift = shift )
        PQ           .__init__ ( self , logp  = logp  , logq  = logq  ) 
        # =====================================================================
        ## create PDF 
        self.pdf  = Ostap.Models.Beta (
            self.roo_name ( 'beta_' ) ,
            'Beta %s' % self.name     , 
            self.x     ,            
            self.logp  ,
            self.logq  ,
            self.scale ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'logp'  : self.logp  ,            
            'logq'  : self.logq  ,            
            'scale' : self.scale ,            
            'shift' : self.shift ,            
            }

models.append ( Beta_pdf ) 

# =============================================================================
## @class BetaPrime_pdf
#  http://en.wikipedia.org/wiki/Beta_prime_distribution
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::BetaPrime
#  @see Ostap::Math::BetaPrime
class BetaPrime_pdf(Beta_pdf) :
    """ Beta-prime disribution 
    - http://en.wikipedia.org/wiki/Beta_prime_distribution
    """
    ## constructor
    def __init__ ( self        , * , 
                   xvar        , ## the variable
                   name   = '' , ## the name 
                   logp   = 0  ,
                   logq   = 0  ,
                   scale  = 1  ,
                   shift  = 0  ) :
        ## 
        Beta_pdf.__init__ ( self  ,
                            name  = name  ,
                            xvar  = xvar  ,
                            logp  = logp  ,
                            logq  = logq  ,
                            scale = scale ,
                            shift = shift )
        
        # 
        self.pdf  = Ostap.Models.BetaPrime (
            self.roo_name ( 'betap_' ) ,
            "Beta' %s" % self.name     , 
            self.x     ,
            self.logp  ,
            self.logq  ,
            self.scale ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'logp'  : self.logp  ,            
            'logq'  : self.logq  ,            
            'scale' : self.scale ,            
            'shift' : self.shift ,            
            }

models.append ( BetaPrime_pdf ) 

# =============================================================================
## @class GenBetaPrime_pdf
#  Generalized beta'-distribution 
#  http://en.wikipedia.org/wiki/Beta_prime_distribution
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::GenBetaPrime
#  @see Ostap::Math::GenBetaPrime
#  @see Ostap::Models::BetaPrime
#  @see Ostap::Math::BetaPrime
class GenBetaPrime_pdf(Beta_pdf) :
    """ Generalized Beta-prime disribution 
    - http://en.wikipedia.org/wiki/Beta_prime_distribution
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name                   
                   a     = 1  ,   ## alpha-parameter                   
                   logp  = 0  ,   ## p-parameter, p>0 
                   logq  = 0  ,   ## q-parameter, q>0 
                   scale = 1  ,   ## scale-parameter 
                   shift = 0  ) : ## shift-parameter
        
        ## inialize the base 
        Beta_pdf.__init__ ( self  ,
                            xvar  = xvar  , 
                            name  = name  ,
                            logp  = logp  ,
                            logq  = logq  ,
                            scale = scale ,
                            shift = shift ) 
        
        ##  a: It also could  be negative!  
        self.__a = self.make_var ( a,
                                   'a_%s'          % self.name ,
                                   'a_{#beta}(%s)' % self.name ,
                                   None  , 1 , 1.e-3 , 100 )
        
        ## make PDF 
        self.pdf  = Ostap.Models.GenBetaPrime (
            self.roo_name ( 'genbetap_' )   ,
            "gen-Beta' %s"% self.name , 
            self.x     ,
            self.a     ,
            self.logp  ,
            self.logq  ,
            self.scale ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'a'     : self.a     ,            
            'logp'  : self.logp  ,            
            'logq'  : self.logq  ,            
            'scale' : self.scale ,            
            'shift' : self.shift ,            
            }
        
    @property
    def a ( self ) :
        """`e` : parameter a for gen-Beta(') distribution
        """
        return self.__a
    @a.setter
    def a ( self , value ) :
        self.set_value ( self.__a , value )
        
models.append ( GenBetaPrime_pdf ) 

# =============================================================================
## @class GenBeta1_pdf
#  Generalized beta-distribution
#  Generalized Beta distribution of 1st kind 
#  @see https://en.wikipedia.org/wiki/Generalized_beta_distribution
#  @see Ostap::Models::GenBeta1 
#  @see Ostap::Math::GenBeta1 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class GenBeta1_pdf(GenBetaPrime_pdf) :
    """ Generalized Beta-prime disribution of 1st kind 
    - http://en.wikipedia.org/wiki/Beta_prime_distribution
    """
    ## constructor
    def __init__ ( self        , * , 
                   xvar        ,   ## the variable
                   name  = ''  ,   ## the name 
                   a     = 1   ,   ## alpha-parameter (shape) 
                   logp  = 0   ,   ## log(p)-parameter
                   logq  = 0   ,   ## log(q_-parameter
                   scale = 1   ,   ## scale-parameter
                   shift = 0   ) : ## shift-parameter
        
        ## inialize the base 
        GenBetaPrime_pdf.__init__ ( self ,
                                    xvar  = xvar , 
                                    name  = name ,
                                    a     = a    , 
                                    logp  = logp ,
                                    logq  = logq ,
                                    scale = scale ,
                                    shift = shift )
        
        self.pdf  = Ostap.Models.GenBeta1 (
            self.roo_name ( 'genbeta1_' )   ,
            "gen-Beta-1 %s"% self.name , 
            self.x     ,
            self.a     ,
            self.logp  ,
            self.logq  ,
            self.scale ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'a'     : self.a     ,            
            'logp'  : self.logp  ,            
            'logq'  : self.logq  ,            
            'scale' : self.scale ,            
            'shift' : self.shift ,            
            }

models.append ( GenBeta1_pdf ) 

# =============================================================================
## @class GenBeta2_pdf
#  Generalized beta-distribution
#  Generalized Beta distribution of 2nd kind 
#  @see https://en.wikipedia.org/wiki/Generalized_beta_distribution
#  @see Ostap::Models::GenBeta2 
#  @see Ostap::Math::GenBeta2 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class GenBeta2_pdf(GenBeta1_pdf) :
    """ Generalized Beta-prime disribution of 2nd kind 
    - http://en.wikipedia.org/wiki/Beta_prime_distribution
    """
    ## constructor
    def __init__ ( self        , * , 
                   xvar        ,   ## the variable
                   name  = ''  ,   ## the name 
                   a     = 1   ,   ## alpha-parameter (shape) 
                   logp  = 0   ,   ## log(p)-parameter
                   logq  = 0   ,   ## log(q_-parameter
                   scale = 1   ,   ## scale-parameter
                   shift = 0   ) : ## shift-parameter
        
        ## inialize the base 
        GenBeta1_pdf.__init__ ( self ,
                                xvar  = xvar , 
                                name  = name ,
                                a     = a    , 
                                logp  = logp ,
                                logq  = logq ,
                                scale = scale ,
                                shift = shift )
        
        self.pdf  = Ostap.Models.GenBeta2 (
            self.roo_name ( 'genbeta2_' )   ,
            "gen-Beta-2 %s"% self.name , 
            self.x     ,
            self.a     ,
            self.logp  ,
            self.logq  ,
            self.scale ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'a'     : self.a     ,            
            'logp'  : self.logp  ,            
            'logq'  : self.logq  ,            
            'scale' : self.scale ,            
            'shift' : self.shift ,            
            }

models.append ( GenBeta2_pdf ) 

# =============================================================================
## @class GenBeta_pdf
#  Generalized beta-distribution
#  Generalized Beta distribution (the most general form)
#  - \f$ c \equiv \sin^2 \frac{\pi\gamma}{2} \f$ to ensure \f$ 0 \le c \le 1 \f$ 
#  @see https://en.wikipedia.org/wiki/Generalized_beta_distribution
#  @see Ostap::Models::GenBeta 
#  @see Ostap::Math::GenBeta 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
class GenBeta_pdf(GenBeta2_pdf,R) :
    """ Generalized Beta-prime disribution 
    - http://en.wikipedia.org/wiki/Beta_prime_distribution
    """
    ## constructor
    def __init__ ( self        , * , 
                   xvar        ,   ## the variable
                   name  = ''  ,   ## the name 
                   a     = 1   ,   ## alpha-parameter (shape) 
                   logr  = 0   ,   ## log(r)
                   logp  = 0   ,   ## log(p)-parameter
                   logq  = 0   ,   ## log(q_-parameter
                   scale = 1   ,   ## scale-parameter
                   shift = 0   ) : ## shift-parameter
        
        ## inialize the base 
        GenBeta2_pdf.__init__ ( self ,
                                xvar  = xvar , 
                                name  = name ,
                                a     = a    , 
                                logp  = logp ,
                                logq  = logq ,
                                scale = scale ,
                                shift = shift )
        ## the second base 
        R.__init__ ( self , logr = logr )
        
        self.pdf  = Ostap.Models.GenBeta (
            self.roo_name ( 'genbeta_' )   ,
            "gen-Beta %s"% self.name , 
            self.x     ,
            self.a     ,
            self.logr  ,
            self.logp  ,
            self.logq  ,
            self.scale ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'a'     : self.a     ,            
            'logr'  : self.logr  ,            
            'logp'  : self.logp  ,            
            'logq'  : self.logq  ,            
            'scale' : self.scale ,            
            'shift' : self.shift ,            
            }
        
    @property
    def logr ( self ) :
        """`logr'-parameter: logarithm of R-parameter """
        return self.__logr
    @logr.setter 
    def logr ( self , value ) :
        self.set_value ( self.__logr , value )

models.append ( GenBeta_pdf ) 

# =============================================================================
## @class TwoExpos_pdf
#  simple difference of two exponents
#  \f$ f \propto 
#        \mathrm{e}^{-a_1    x}       -\mathrm{e}^{-a_2 x} = 
#        \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f$
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @see Ostap::Models::TwoExpos
#  @see Ostap::Math::TwoExpos
class TwoExpos_pdf(PDF1,Shift) :
    r""" Simple difference of two exponents:    
    \f$ f \propto 
    \mathrm{e}^{-a_1    x}       -\mathrm{e}^{-a_2 x} = 
    \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f$    
    """
    ## constructor
    def __init__ ( self        , * , 
                   xvar        ,   ## the variable
                   name  = ''  ,   ## the name 
                   alpha = 1   ,   ## shape-parameter 
                   delta = 1   ,   ## high-parameter 
                   x0    = 0   ) : ## low-parameter 
        #
        PDF1  .__init__ ( self , name = name , xvar = xvar )
        Shift .__init__ ( self ,
                          shift = x0 ,  
                          shift_name  = 'x0_%s'       % self.name ,
                          shift_title = 'x0_{2e}(%s)' % self.name ) 
        #
        self.__alpha  = self.make_var ( alpha      ,
                                        'alpha_%s'        % self.name ,
                                        '#alpha_{2e}(%s)' % self.name , 
                                        None , 1 , 1.e-4 , 100 )
        self.__delta  = self.make_var ( delta     ,
                                        'delta_%s'        % self.name ,
                                        '#delta_{2e}(%s)' % self.name , 
                                        None , 1 , 1.e-4 , 100 )
        
        self.pdf  = Ostap.Models.TwoExpos (
            self.roo_name ( 'exp2_' )   ,
            'Two exponentials %s' % self.name  , 
            self.x     ,
            self.alpha ,
            self.delta ,
            self.x0    )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'alpha' : self.alpha ,            
            'delta' : self.delta ,            
            'x0'    : self.x0    ,            
            }

    ## ALIAS
    x0 = Shift.shift
    
    @property
    def alpha ( self ) :
        """`alpha'-parameter (slope of leading exponent) of the TwoExpo function"""
        return self.__alpha
    @alpha.setter
    def alpha ( self , value ) :
        self.set_value ( self.__alpha , value )
    
    @property
    def delta ( self ) :
        """`delta'-parameter (second exponent slope is `alpha+delta') of the TwoExpo function"""
        return self.__delta
    @delta.setter
    def delta ( self , value ) :
        self.set_value ( self.__delta , value )
    
models.append ( TwoExpos_pdf )

# =============================================================================
## @class Gumbel_pdf
#  Gumbel distribution
#  @see https://en.wikipedia.org/wiki/Gumbel_distribution
#  \f$  f(x,\mu,\beta) = \frac{1}{\left|\beta\right|} e^{-e^{-z}} \f$,
#  where \f$ z = \frac{x-\mu}{\beta}\f$
#  Very useful and important case: 
#  if \f$ g(x) \propto exp(-\tau x ) \f$ and \f$ z = \log(x) \f$,
#  than \f$ F(z) = g(x) = f(z; -log(\tau) ,  1) \f$
#  It means that sum exponential components will be represented by a set  of
#  peak-like shifted Gumbel' structures with \f$ \beta=1, \mu=-log(\tau) \f$ 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2017-09-02
#  @see Ostap::Models::Gumbel
class Gumbel_pdf(PDF1,ShiftAndScale) :
    r""" Gumbel distribution
    - see https://en.wikipedia.org/wiki/Gumbel_distribution
    \f$  f(x,\mu,\beta) = \frac{1}{\left|\beta\right|} e^{-e^{-z}} \f$,
    where \f$ z = \frac{x-\mu}{\beta}\f$
    - Very useful and important case: 
    if \f$ g(x) \propto exp(-\tau x ) \f$ and \f$ z = \log(x) \f$,
    than \f$ F(z) = g(x) = f(z; -log(\tau) ,  1) \f$
    It means that sum of exponential components will be represented by a set of
    peak-like Gumbel' structures with \f$ \beta=1, \mu=-log(\tau) \f$ 
    """
    ## constructor
    def __init__ ( self      , * , 
                   xvar      ,   ## the variable 
                   name = '' ,   ## the name 
                   mu   = 0  ,   ## shift parameter/mode
                   beta = 1  ) : ## scale parameter 
        #
        PDF1         .__init__ ( self , name = name , xvar = xvar )
        ShiftAndScale.__init__ ( self ,
                                 shift       = mu   , 
                                 scale       = beta ,
                                 scale_name  = 'beta_%s'       % self.name , 
                                 scale_title = '#beta_{G}(%s)' % self.name , 
                                 shift_name  = 'mu_%s'         % self.name , 
                                 shift_title = '#mu_{G}(%s)'   % self.name )
        
        self.pdf  = Ostap.Models.Gumbel (
            self.roo_name ( 'gumbel_' )   ,
            'Gumbel %s' % self.name  , 
            self.x              ,
            self.mu             ,
            self.beta           )

        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mu'    : self.mu    ,            
            'beta'  : self.beta  ,            
            }

    mu   = ShiftAndScale.shift
    beta = ShiftAndScale.scale 

models.append ( Gumbel_pdf ) 

# =============================================================================
## @class Rice_pdf
#  Rice distribution
#  @see Ostap::Math::Rice
#  Rice distribution 
#  \f$ f(x; \nu , \varsigma) = 
#  \frac{\delta x}{\varsigma^2} \mathrm{e}^{-\frac{ \delta x^2+\nu^2}{2\varsigma^2} } 
#   I_0 (\frac{\delta x\nu}{\varsigma^2}) \f$, 
#   where  \f$ \delta x = x - \x_0\f$ and 
#    - \f$ x\ge x_0\f$  
#    - \f$ \nu \ge 0 \f$ 
#    - \f$ \varsigma \ge 0 \f$ 
#  @see https://en.wikipedia.org/wiki/Rice_distribution
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @see Ostap::Math::Rice
#  @see Ostap::Models::Rice 
class Rice_pdf(PDF1, ShiftAndScale) :
    """ Rice distribution
    - see https://en.wikipedia.org/wiki/Rice_distribution
    """
    ## constructor
    def __init__ ( self          , * , 
                   xvar          , ## the variable 
                   name     = '' , ## the name 
                   nu       = 0  , ## parameter nu 
                   varsigma = 1  , ## parameter varsigma
                   shift    = ROOT.RooFit.RooConst ( 0.0 )  ) : ## shift parameter
        #
        PDF1         .__init__ ( self , name = name , xvar = xvar )
        ShiftAndScale.__init__ ( self ,
                                 shift = shift    ,
                                 scale = varsigma , 
                                 scale_name  = 'varsigma_%s'       % self.name , 
                                 scale_title = '#varsigma_{R}(%s)' % self.name  )
        
        self.__nu       = self.make_var ( nu        ,
                                          'nu_%s'                  % self.name ,
                                          '#nu_{Rice}(%s)'         % self.name ,
                                          None , 0 , 0.0    , 1.e+6 )
        self.pdf  = Ostap.Models.Rice (
            self.roo_name ( 'rice_' )   ,
            'Rice %s' % self.name  , 
            self.x        ,
            self.nu       ,
            self.varsigma ,
            self.shift    )
        
        ## save the configuration:
        self.config = {
            'name'     : self.name     ,
            'xvar'     : self.xvar     ,
            'nu'       : self.nu       ,            
            'varsigma' : self.varsigma ,            
            'shift'    : self.shift    ,            
            }
        
    @property
    def nu ( self ) :
        """`nu'-parameter of Rice function"""
        return self.__nu
    @nu.setter
    def nu ( self , value ) :
        self.set_value ( self.__nu , value )

    ## ALIAS 
    varsigma = ShiftAndScale.scale
    
models.append ( Rice_pdf ) 

# =============================================================================
## @class GenInvGauss_pdf
#  Generalized Inverse Gaussian distribution
#  @see Ostap::Math::GenInvGauss
#  Generalised Inverse Gaussian distribution using
#  \f$ (\theta,\eta) \f$ parameterisation  
#  - |f$ \theta = \sqrt{ab}\$ 
#  - |f$ \eta   = \sqrt{\frac{b}{a}}\$ 
#   @see https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution

#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @see Ostap::Math::GenInvGauss
#  @see Ostap::Models::GenInvGauss
class GenInvGauss_pdf(PDF1,Shift) :
    """ Generalized Inverse Gaussian distribution
    Generalised Inverse Gaussian distribution using (theta,eta) parameterisation
    - see https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution
    - see Ostap::Math::GenInvGauss
    """
    ## constructor
    def __init__ ( self        , * , 
                   xvar        , ## the variable 
                   name  = ''  , ## the name 
                   theta = 1   , ## parameter theta
                   eta   = 1   , ## parameter eta
                   p     = 0   , ## parameter p
                   shift = ROOT.RooFit.RooConst ( 0.0 )  ) : ## shift parameter
        #
        PDF1 .__init__ ( self , name  = name  , xvar = xvar )
        Shift.__init__ ( self , shift = shift )
        #
        self.__theta    = self.make_var ( theta             ,
                                          'theta_%s'          % self.name ,
                                          '#theta_{GIG}(%s)'  % self.name ,
                                          None , 1.0 , 1.e-8 , 100 )
        self.__eta      = self.make_var ( eta               ,
                                          'eta_%s'            % self.name ,
                                          '#eta_{GIG}(%s)'    % self.name ,
                                          None ,  1.0 , 1.e-8 , 100 )
        self.__p        = self.make_var ( p                 ,
                                          'p_%s'              % self.name ,
                                          'p_{GIG}(%s)'       % self.name ,
                                          None , 0   , -100  , 100 )        
        
        self.pdf  = Ostap.Models.GenInvGauss (
            self.roo_name ( 'gig_' )   ,
            'GenInvGauss %s' % self.name  , 
            self.x      ,
            self.theta  ,
            self.eta    ,
            self.p      ,
            self.shift  )
        
        ## save the configuration:
        self.config = {
            'name'     : self.name     ,
            'xvar'     : self.xvar     ,
            'theta'    : self.theta    ,
            'eta'      : self.eta      ,
            'p'        : self.p        ,
            'shift'    : self.shift    ,            
            }
        
    @property
    def theta ( self ) :
        """`theta'-parameter of Generalized Inverse Gaussian  function"""
        return self.__theta
    @theta.setter
    def theta ( self , value ) :
        self.set_value ( self.__theta , value )

    @property
    def eta ( self ) :
        """`eta'-parameter of Generalized Inverse Gaussian  function"""
        return self.__eta
    @eta.setter
    def eta ( self , value ) :
        self.set_value ( self.__eta , value )

    @property
    def p   ( self ) :
        """`p'-parameter of Generalized Inverse Gaussian  function"""
        return self.__p
    @p.setter
    def p ( self , value ) :
        self.set_value ( self.__p , value )

models.append ( GenInvGauss_pdf ) 

# =============================================================================
## @class Weibull_pdf 
#  3-parameter  Weibull distribution 
#  \f$ f(x,\lambda,k,x_0) = \frac{k}{\lambda}  y^{k-1} e^{-y^k}\f$, where 
#  \f$ y \equiv \frac{x-x_0}{\lambda}\f$
#  @see https://en.wikipedia.org/wiki/Weibull_distribution
#  Shape parameter:
#  - for k>2 is has peak-like, 'signal'-like shape
#  - for 1<k<2 is has smooth shape, 'threshold'-like
#  - for k<1 is is smooth decreasing shape
#  @see Ostap::Models::Weibull 
#  @see Ostap::Math::Weibull 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-02-27
class Weibull_pdf(PDF1,ShiftAndScale) :
    r""" 3-parameter  Weibull distribution 
    \f$ f(x,\lambda,k,x_0) = \frac{k}{\lambda}  y^{k-1} e^{-y^k}\f$, where 
    \f$ y \equiv \frac{x-x_0}{\lambda}\f$
    
    Shape parameter:
    
    - for k>2 is has peak-like, signal-like shape
    - for 1<k<2 is has smooth shape,  'threshold'-like
    - for k<1 is is smooth decreasing shape
    
    - see https://en.wikipedia.org/wiki/Weibull_distribution
    - see Ostap::Models::Weibull 
    - see Ostap::Math::Weibull 
    """
    def __init__ ( self       , * , 
                   xvar       ,
                   name  = '' ,
                   scale = 1  ,   ## scale/lambda 
                   shape = 1  ,   ## shape/k 
                   shift = 0  ) : ## shift/x0      
        #
        ## initialize the base
        #        
        PDF1         .__init__  ( self , name  = name   , xvar = xvar  )
        ShiftAndScale.__init__  ( self , scale = scale , shift = shift ) 
        ## shape 
        self.__shape = self.make_var ( shape                   ,
                                       'shape_%s'              % self.name ,
                                       'shape_{Weibull}(%s)'   % self.name ,
                                       None , 1 , 1.e-8 , 1.e+6 )
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Weibull (
            self.roo_name ( 'weibull_' )   ,
            'Weibull %s' % self.name  , 
            self.xvar  ,
            self.scale ,
            self.shape ,
            self.shift )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'scale'     : self.scale ,
            'shape'     : self.shape ,
            'shift'     : self.shift ,
            }

    @property
    def shape ( self ) :
        """`shape'-parameter, shape>0
        - for shape>2    : peak-like, signal-like shape
        - for 1<shape<2  : 'threshold'-like
        - for shape<1    : smooth decreasing         
        """
        return self.__shape    
    @shape.setter
    def shape ( self, value ) :
        self.set_value ( self.__shape , value )

models.append ( Weibull_pdf )      

# =============================================================================
## @class GenPareto_pdf
#  Generalized Pareto distribution
#  @see https://en.wikipedia.org/wiki/Generalized_Pareto_distribution
#  @see Ostap::Models::GenPareto
#  @see Ostap::Math::GenPareto
class GenPareto_pdf(PDF1,ShiftAndScale) :
    """ Generalized Pareto Distirbutution
    - see https://en.wikipedia.org/wiki/Generalized_Pareto_distribution
    - see `Ostap.Models.GenPareto`
    - see `Ostap.Math.GenPareto`
    """
    ## constructor
    def __init__ ( self         , * , 
                   xvar         ,   ## the variable
                   name  = ''   ,   ## the name 
                   mu    = 0    ,   ## location parameter
                   scale = 1    ,   ## scale parameter
                   shape = 1    ) : ## shape parameter
        #
        PDF1         .__init__ ( self , name , xvar )
        ShiftAndScale.__init__ ( self ,
                                 scale = scale ,
                                 shift = mu    , 
                                 shift_name  = 'mu_%s'        % self.name , 
                                 shift_title = '#mu_{GP}(%s)' % self.name  )
        #
        self.__shape = self.make_var ( shape ,
                                       'shape_%s'    % self.name ,
                                       'shape(%s)'   % self.name ,
                                       None , -1+1.e-6 , 100 )
        
        self.pdf  = Ostap.Models.GenPareto (
            self.roo_name ( 'gpd_' ) ,
            'GenPareto %s' % self.name , 
            self.x     ,
            self.mu    ,
            self.scale ,
            self.shape )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mu'    : self.mu    ,            
            'scale' : self.scale ,            
            'shape' : self.shape ,            
            }

    ## ALIAS 
    mu = ShiftAndScale.shift 

    @property
    def shape ( self ) :
        """'shape'- shape parameter"""
        return self.__shape 
    @shape.setter 
    def shape ( self , value ) :
        self.set_value ( self.__shape , value )

        
models.append ( GenPareto_pdf ) 

# =============================================================================
## @class ExGenPareto_pdf
#  Exponentiated Generalized Pareto distribution
#  @see https://en.wikipedia.org/wiki/Generalized_Pareto_distribution
#  @see Ostap::Models::ExGenPareto
#  @see Ostap::Math::ExGenPareto
class ExGenPareto_pdf(GenPareto_pdf) :
    """ Exponentiated Generalized Pareto Distirbutution
    - see https://en.wikipedia.org/wiki/Generalized_Pareto_distribution
    - see `Ostap.Models.ExGenPareto`
    - see `Ostap.Math.ExGenPareto`
    """
    ## constructor
    def __init__ ( self        , * , 
                   xvar        ,   ## the variable
                   name  = ''  ,   ## the name 
                   mu    = 0   ,   ## location parameter
                   scale = 1   ,   ## scale parameter
                   shape = 1   ) : ## shape parameter
        #
        GenPareto_pdf.__init__ ( self ,
                                 name  = name  ,
                                 xvar  = xvar  ,
                                 mu    = mu    ,
                                 scale = scale ,
                                 shape = shape )
        
        self.pdf  = Ostap.Models.ExGenPareto (
            self.roo_name ( 'egpd_' ) ,
            'ExGenPareto %s' % self.name , 
            self.xvar  ,
            self.mu    ,
            self.scale ,
            self.shape )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mu'    : self.mu    ,            
            'scale' : self.scale ,            
            'shape' : self.shape ,            
            }
        
models.append ( ExGenPareto_pdf ) 

# =============================================================================
## @class GEV_pdf
#  Generalized Extreme Value distribution
#  @see https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution
#  @see Ostap::Models::GEV
#  @see Ostap::Math::GEV
class GEV_pdf(GenPareto_pdf) :
    """ Generalized Extreme Value distribution
    - see https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution
    - see `Ostap.Models.GEV`
    - see `Ostap.Math.GEV`
    """
    ## constructor
    def __init__ ( self        , * , 
                   xvar        ,   ## the variable
                   name  = ''  ,   ## the name 
                   mu    = 0   ,   ## location parameter
                   scale = 1   ,   ## scale parameter
                   shape = 1   ) : ## shape parameter
        #
        GenPareto_pdf.__init__ ( self ,
                                 name  = name  ,
                                 xvar  = xvar  ,
                                 mu    = mu    ,
                                 scale = scale ,
                                 shape = shape )
        
        self.pdf  = Ostap.Models.GEV (
            self.roo_name ( 'gev_' ) ,
            'GEV %s' % self.name , 
            self.x     ,
            self.mu    ,
            self.scale ,
            self.shape )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mu'    : self.mu    ,            
            'scale' : self.scale ,            
            'shape' : self.shape ,            
            }
        
models.append ( GEV_pdf ) 

# =============================================================================
## @class Benini_pdf
# Modified version of Benini distribution 
# @see https://en.wikipedia.org/wiki/Benini_distribution
# Parameters 
#  - \f$ 0 < \alpha  \f$ linear in log 
#  - \f$ 0 < \beta   \f$ quadratic in log 
#  - \f$ 0 < \gamma  \f$ cubic in log 
#  - \f$ 0 < \delta  \f$ 4th order  in log 
#  - shift  \f$ \mu  \f$ 
#  - scale  \f$ s \f$ 
#
#  For standard Benini one has \f$ (\frac{x}{\sigma}\f$,
#  here one has \f$ (\frac{x-\delta}{s}\f$
#
#  Standard Benini distribution: 
#  - \f$ \mu=0    \f$ 
#  - \f$ \gamma=0 \f$ 
#  - \f$ \delta=0 \f$ 
# 
#  @see Ostap::Models::Benini
#  @see Ostap::Math::Benini
class Benini_pdf(PDF1, ShiftAndScale) :
    """ Modified version of Benini distribution 
    - see https://en.wikipedia.org/wiki/Benini_distribution
    Parameters 
    - 0 <  alpha  : linear     in log 
    - 0 <= beta   : quadractic in log 
    - 0 <= gamma  : cubic      in log 
    - 0 <= delta  : cubic     in log 
    - shift mu  
    - scale s 
    
    For standard Benini pne has (x/sigma) and 
    here one has (x-delta)/s 
    
    Standard Benini distribution: 
    - mu=0
    - gamma=0
    - delta=0
    
    - see `Ostap.Models.Benini`
    - see `Ostap.Math.Benini`
    """
    ## constructor
    def __init__ ( self                    , * , 
                   xvar                    ,   ## the variable
                   name  = ''              ,   ## the name 
                   shape =  ( 1 , 1 , 1 )  ,   ## shape  parameters
                   scale = 1               ,   ## scale  parameter
                   shift = 0               ) : ## shift  parameter
        #
        PDF1         .__init__ ( self , name , xvar )
        ShiftAndScale.__init__ ( self , scale = scale , shift = shift ) 
        #
        ##
        xmnmx = self.xminmax()

        self.__shape = ROOT.RooArgList()
        for  i , p in enumerate ( shape ) :  
            pp = self.make_var ( p ,
                                 'p_%s_%s'  % ( i , self.name ) ,
                                 'p_%s(%s)' % ( i , self.name ) ,
                                 None , 0 , 1000 )
            self.__shape.add  ( pp )

        assert 2 <= len ( self.__shape ) , 'Invalid number of shape parameters!' 
        
        self.pdf  = Ostap.Models.Benini (
            self.roo_name ( 'benini_' ) ,
            'Benini %s' % self.name , 
            self.x     ,
            self.shape ,
            self.scale ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'shape' : self.shape ,            
            'scale' : self.scale ,            
            'shift' : self.shift ,            
            }
        
    @property
    def alpha ( self ) :
        """'alpha'- shape parameter for Benini-distribution"""
        return self.__shape[0]
    @alpha.setter 
    def alpha ( self , value ) :
        self.set_value ( self.__shape [ 0 ] , value )
        
    @property
    def beta ( self ) :
        """'beta'- shape parameter for Benini-distribution"""
        assert 1 < len ( self.__shape ) , 'Invalid parameter!'
        return self.__shape[1] 
    @beta.setter 
    def beta ( self , value ) :
        assert 1 < len ( self.__shape ) , 'Invalid parameter!'
        self.set_value ( self.__shape [ 1 ] , value )

    @property
    def gamma ( self ) :
        """'gamma'- shape parameter for Benini-distribution"""
        assert 2 < len ( self.__shape ) , 'Invalid parameter!'        
        return self.__shape [ 2 ] 
    @gamma.setter 
    def gamma ( self , value ) :
        assert 2 < len ( self.__shape ) , 'Invalid parameter!'        
        self.set_value ( self.__shape [ 2 ]  , value )

    @property
    def delta ( self ) :
        """'delta'- shape parameter for Benini-distribution"""
        assert 3 < len ( self.__shape ) , 'Invalid parameter!'        
        return self.__shape [ 3 ] 
    @delta.setter 
    def delta ( self , value ) :
        assert 3 < len ( self.__shape ) , 'Invalid parameter!'
        self.set_value ( self.__shape [ 3 ]  , value )

    @property
    def shape  ( self ) :
        """'shape' : get all shape parameters"""
        return self.__shape
    @shape.setter
    def shape ( self , values ) :
        self.component_setter ( self.__shape , values )
        
models.append ( Benini_pdf ) 


# =============================================================================
## @class MPERT_pdf
#  Modified PERT distribution 
#  @see https://en.wikipedia.org/wiki/PERT_distribution
#  @see https://www.vosesoftware.com/riskwiki/ModifiedPERTdistribution.php
#  @see Ostap::Models::MPERT
#  @see Ostap::Math::MPERT
class MPERT_pdf(PDF1) :
    """ Modified PERT distribution 
    - see https://en.wikipedia.org/wiki/PERT_distribution
    - see https://www.vosesoftware.com/riskwiki/ModifiedPERTdistribution.php
    - see `Ostap.Models.MPERT`
    - see `Ostap.Math.MPERT`
    """
    ## constructor
    def __init__ ( self                  , * , 
                   xvar                  ,   ## the variable
                   name   = ''           ,   ## the name 
                   xi     = 0.5          ,   ## mode parameter 
                   gamma  = 4            ,   ## gamma/shape 
                   Xmin   = pos_infinity ,   ## xmin 
                   Xmax   = neg_infinity ) : ## xmax
        #
        PDF1.__init__ ( self , name = name , xvar = xvar )
        #
        xmnmx = self.xminmax()
        if xmnmx :
            xmn , xmx = xmnmx
            if isfinite ( Xmin ) and xmn <= Xmin <= xmx : xmn = Xmin
            if isfinite ( Xmax ) and xmn <= Xmax <= xmx : xmx = Xmax
            assert xmn < xmx , 'Invalid setting xmn/xmx/Xmin/Xmax: %s/%s/%s/%s' % ( xmn , xmx , Xmin , Xmax ) 
            xmid        = 0.5 * ( xmn + xmx )
            limits      = xmid , xmn , xmx
            Xmin , Xmax = xmn  , xmx 
        elif isfinite ( Xmin ) and isfinite ( Xmax ) and Xmin < Xmax :
            xmid        = 0.5 * ( Xmin + Xmax )
            limits      = xmid , Xmin , Xmax        
        else :
            raise TypeError ( "Cannot deduce Xmin/Xmax parameters!" )
                
        ## the mode 
        self.__xi    = self.make_var ( xi                   ,
                                       'xi_%s'       % self.name ,
                                       '#xi(%s)'     % self.name ,
                                       None , *limits )
        
        self.__gamma = self.make_var ( gamma                   ,
                                       'gamma_%s'       % self.name ,
                                       '#gamma(%s)'     % self.name ,
                                       None , 0 , 100 )
        
        self.__Xmin = Xmin
        self.__Xmax = Xmax

        self.pdf  = Ostap.Models.MPERT (
            self.roo_name ( 'mpert_' ) ,
            'MPERT %s' % self.name     , 
            self.x                     ,
            self.xi                    ,
            self.gamma                 , 
            self.Xmin                  ,
            self.Xmax                  )
        
        ## save the configuration:
        self.config = {
            'name'   : self.name  ,
            'xvar'   : self.xvar  ,
            'xi'     : self.xi    , 
            'gamma'  : self.gamma ,            
            'Xmin'   : self.Xmin  ,
            'Xmax'   : self.Xmax  ,         
            }
        
    @property
    def xi ( self ) :
        """'xi'- mode parameter"""
        return self.__xi
    @xi.setter 
    def xi ( self , value ) :
        self.set_value ( self.__xi , value )

    @property
    def gamma ( self ) :
        """gamma'- shape parameter"""
        return self.__gamma
    @gamma.setter 
    def gamma ( self , value ) :
        self.set_value ( self.__gamma , value )

    @property
    def Xmin ( self ) :
        """'xmin'- parameter"""
        return self.__Xmin
    
    @property
    def Xmax ( self ) :
        """'xmin'- parameter"""
        return self.__Xmax 
        
models.append ( MPERT_pdf ) 

# =============================================================================
## @class BirnbaumSaunders_pdf
#  Birnbaum-Saunders distribution 
#  @see https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution
#  \f[ f(x;\mu, \beta,\gamma) = 
#   \frac{ z + z^{-1}}{2\gamma(x-\mu)}\phi( \frac{1}{\gamma}(z-z^{-1}) \f]
#  where
#   - \f$ z=\frac{x-\mu}{\beta}\f$
#   - \f$ \phi\f$ is Gaussian PDF 
#  @see Ostap::Models::BirnbaumSaunders
#  @see Ostap::Math::BirnbaumSaunders
class BirnbaumSaunders_pdf(PDF1,ShiftAndScale) :
    """ Birnbaum-Saunders distribution 
    - see https://en.wikipedia.org/wiki/Birnbaum%E2%80%93Saunders_distribution
    - see `Ostap.Models.BirnbaumSaunders`
    - see `Ostap.Math.BirnbaumSaunders`
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name 
                   gamma = 1  ,   ## shape  
                   mu    = 0  ,   ## shift/location
                   beta  = 1  ) : ## scale parameter         
        #
        PDF1         .__init__ ( self , name = name , xvar = xvar )
        ShiftAndScale.__init__ ( self ,
                                 shift       = gamma , 
                                 scale       = beta  , 
                                 scale_name  = 'beta_%s'        % self.name , 
                                 scale_title = '#beta_{BS}(%s)' % self.name  , 
                                 shift_name  = 'mu_%s'          % self.name , 
                                 shift_title = '#mu_{BS}(%s)'   % self.name  )
        
        self.__gamma  = self.make_var ( gamma    ,
                                        'gamma_%s'              % self.name ,
                                        '#gamma_{BS}(%s)'       % self.name ,
                                        None , 1 , 0.001 , 1000 )
        

        ## create PDF 
        self.pdf  = Ostap.Models.BirnbaumSaunders (
            self.roo_name ( 'bs_' ) ,
            'Birnbaum-Saunders %s' % self.name     , 
            self.x      , 
            self.mu     , 
            self.beta   , 
            self.gamma  ) 
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'gamma' : self.gamma ,            
            'mu'    : self.mu         ,            
            'beta'  : self.beta       ,            
            }
        
    beta = ShiftAndScale.scale
    mu   = ShiftAndScale.shift 
            
    @property
    def gamma ( self ) :
        """`gamma`- shape parameter for Birnbaum-Saunders distribution"""
        return self.__gamma 
    @gamma.setter 
    def gamma ( self , value ) :
        self.set_value ( self.__gamma , value )
        
models.append ( BirnbaumSaunders_pdf ) 

# =============================================================================
## @class Frechet_pdf
#  Frechet distribution
#  @see https://en.wikipedia.org/wiki/Fr%C3%A9chet_distribution
#  @see Ostap::Models::Frechet
#  @see Ostap::Math::Frechet
class Frechet_pdf(PDF1,ShiftAndScale,Alpha) :
    """ Frechet distribution
    - see https://en.wikipedia.org/wiki/Fr%C3%A9chet_distribution
    - see Ostap::Models::Frechet
    - see Ostap::Math::Frechet
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name 
                   alpha = 1  ,   ## shape 
                   scale = 1  ,   ## scale 
                   shift = 0  ) : ## shift  
        
        #
        PDF1         .__init__ ( self  , name = name , xvar = xvar )
        ShiftAndScale.__init__ ( self  ,
                                 scale = scale ,
                                 shift = shift )
        Alpha        .__init__ ( self , alpha = alpha )
        ##
        ## create PDF 
        self.pdf  = Ostap.Models.Frechet (
            self.roo_name ( 'fr_' )  ,
            'Frechet %s' % self.name , 
            self.x     ,
            self.alpha ,
            self.scale ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'alpha' : self.alpha ,            
            'scale' : self.scale ,            
            'shift' : self.shift }
        
models.append ( Frechet_pdf ) 

# =============================================================================
## @class Dagum_pdf
#  Dagum distribution
#  @see https://en.wikipedia.org/wiki/Dagum_distribution
#  @see Ostap::Models::Dagum
#  @see Ostap::Math::Dagum 
class Dagum_pdf(PDF1,ShiftAndScale) :
    """ Dagum distribution
    - see https://en.wikipedia.org/wiki/Dagum_distribution
    - see Ostap::Models::Dagum
    - see Ostap::Math::Dagum 
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name 
                   p     = 1  ,
                   a     = 1  ,   
                   b     = 1  ,   ## b-parametes is a scale 
                   shift = 0  ) : ## shift  
        
        #
        PDF1         .__init__ ( self  , name = name , xvar = xvar )
        ShiftAndScale.__init__ ( self  ,
                                 scale       = b     ,
                                 shift       = shift , 
                                 scale_name  = 'b_%s'      % self.name , 
                                 scale_title = 'b_{D}(%s)' % self.name )
        
        #
        self.__a  = self.make_var ( a    ,
                                    'a_%s'       % self.name ,
                                    'a_{D}(%s)'  % self.name ,
                                    None , a , 1.e-6 , 100 )
        self.__p  = self.make_var ( p    ,
                                    'p_%s'       % self.name ,
                                    'p_{D}(%s)'  % self.name ,
                                    None , p , 1.e-6 , 100 )
        
        
        ## create PDF 
        self.pdf  = Ostap.Models.Dagum (
            self.roo_name ( 'dagum_' )  ,
            'Dagum %s' % self.name , 
            self.x     ,
            self.p     ,
            self.a     ,
            self.b     ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'p'     : self.p     ,            
            'a'     : self.a     ,            
            'b'     : self.b     ,            
            'shift' : self.shift }

    ## ALIAS
    b  = ShiftAndScale.scale
    
    @property
    def a ( self ) :
        """`a` : a-parameter for Dagum distribution"""
        return self.__a
    @a.setter 
    def a ( self , value ) :
        self.set_value ( self.__a , value )

    @property
    def p ( self ) :
        """`p` : p-parameter for Dagum distribution"""
        return self.__p
    @p.setter 
    def p ( self , value ) :
        self.set_value ( self.__p , value )
        
models.append ( Dagum_pdf ) 

# =============================================================================
## @class BenktanderI_pdf
#  Variant of Benktander type 1 distribution
#  @see https://en.wikipedia.org/wiki/Benktander_type_I_distribution
#  - \f$ z = \frac{x-x_0}{\sigma} + 1 \f$ for \f$ 0\le x \f$
#  - \f$ b = p \frac{a(a+1)}{2}\f$ , where \f$ 0 < p \le 1 \f$
#  - \f$ p = 1 / \sqrt{ r^2 + 1 } \f$
#  @see Ostap::Models::BenktanderI
#  @see Ostap::Math::BenktanderI
class BenktanderI_pdf(PDF1,ShiftAndScale) :
    r"""Variant of Benktander type 1 distribution
    - see https://en.wikipedia.org/wiki/Benktander_type_I_distribution
    - \f$ z = \frac{x-x_0}{\sigma} + 1 \f$ for \f$ 0\le x \f$
    - \f$ b = p \frac{a(a+1)}{2}\f$ , where \f$ 0 < p \le 1 \f$
    - \f$ p = 1 / \sqrt{ r^2 + 1 } \f$
    see Ostap::Models::BenktanderI
    see Ostap::Math::BenktanderI
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name 
                   a     = 1  ,
                   r     = 1  ,   
                   scale = 1  ,   ## scale 
                   shift = 0  ) : ## shift  
        
        #
        PDF1         .__init__ ( self  , name = name , xvar = xvar )
        ShiftAndScale.__init__ ( self  ,
                                 scale = scale ,
                                 shift = shift )
        #
        self.__a  = self.make_var ( a    ,
                                    'a_%s'       % self.name ,
                                    'a_{B}(%s)'  % self.name ,
                                    None , a , 1.e-3 , 1000 )
        self.__r  = self.make_var ( r    ,
                                    'r_%s'       % self.name ,
                                    'r_{D}(%s)'  % self.name ,
                                    None , r  )
        
        ## create PDF 
        self.pdf  = Ostap.Models.BenktanderI (
            self.roo_name ( 'bkI_' )  ,
            'BenktanderI %s' % self.name , 
            self.x     ,
            self.a     ,
            self.r     ,
            self.scale ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'a'     : self.a     ,            
            'r'     : self.r     ,            
            'scale' : self.scale ,            
            'shift' : self.shift }
        
    @property
    def a ( self ) :
        """`a` : a-parameter for Benktander I&II distributions"""
        return self.__a
    @a.setter 
    def a ( self , value ) :
        self.set_value ( self.__a , value )

    @property
    def r ( self ) :
        """`r` : r-parameter for Benktander I&&II distributions"""
        return self.__r
    @r.setter 
    def r ( self , value ) :
        self.set_value ( self.__r , value )

        
models.append ( BenktanderI_pdf )

# =============================================================================
## @class BenktanderII_pdf
#  Variant of Benktander type II distribution
#  @see https://en.wikipedia.org/wiki/Benktander_type_II_distribution
#  - \f$ z = \frac{x-x_0}{\sigma} + 1 \f$ for \f$ 0\le x \f$
#  - \f$ b = 1 / \sqrt{ r^1 + 1 }  \f$ 
#  @see Ostap::Models::BenktanderII
#  @see Ostap::Math::BenktanderII
class BenktanderII_pdf(BenktanderI_pdf) :
    r"""Variant of Benktander type II distribution
    - see https://en.wikipedia.org/wiki/Benktander_type_II_distribution
    - \f$ z = \frac{x-x_0}{\sigma} + 1 \f$ for \f$ 0\le x \f$
    - \f$ b = 1 / \sqrt{ r^1 + 1 }  \f$ 
    - see Ostap::Models::BenktanderII
    - see Ostap::Math::BenktanderII
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name 
                   a     = 1  ,
                   r     = 1  ,   
                   scale = 1  ,   ## scale 
                   shift = 0  ) : ## shift  
        
        #
        BenktanderI_pdf.__init__ ( self  ,
                                   name  = name  ,
                                   xvar  = xvar  ,
                                   a     = a     ,
                                   r     = r     ,
                                   scale = scale ,
                                   shift = shift ) 
        
        ## create PDF 
        self.pdf  = Ostap.Models.BenktanderII (
            self.roo_name ( 'bkII_' )  ,
            'BenktanderII %s' % self.name , 
            self.x     ,
            self.a     ,
            self.r     ,
            self.scale ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'a'     : self.a     ,            
            'r'     : self.r     ,            
            'scale' : self.scale ,            
            'shift' : self.shift }


models.append ( BenktanderII_pdf )

# =============================================================================
## @class LogNormal_pdf
#  Log-normal distribution
#  @see https://en.wikipedia.org/wiki/Log-normal_distribution
#  - We add here an shift parameter
#  - and use "mu = log(scale)"
#  - and we use "sigma" as shape parameter
#  @see Ostap::Models::LogNormal 
#  @see Ostap::Math::LogNormal 
class LogNormal_pdf(PDF1,ShiftAndScale) :
    """ Log-normal distribution
    - see https://en.wikipedia.org/wiki/Log-normal_distribution
    - We add here an shift parameter
    - and use "mu = log(scale)"
    - and we use "sigma" as shape parameter
    - see Ostap::Models::LogNormal 
    - see Ostap::Math::LogNormal 
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name 
                   shape = 1  ,   
                   scale = 1  ,   ## scale 
                   shift = 0  ) : ## shift  
        
        #
        PDF1         .__init__ ( self  , name = name , xvar = xvar )
        ShiftAndScale.__init__ ( self  ,
                                 scale       = scale ,
                                 shift       = shift )
        #
        self.__shape  = self.make_var ( shape ,
                                        'shape_%s'         % self.name ,
                                        '#sigma_{LN}(%s)'  % self.name ,
                                        None , shape , 1.e-3 , 1000 )

        ## create PDF 
        self.pdf  = Ostap.Models.LogNormal (
            self.roo_name ( 'ln_' )  ,
            'LogNormal %s' % self.name , 
            self.x     ,
            self.shape ,
            self.scale ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'shape' : self.shape ,            
            'scale' : self.scale ,            
            'shift' : self.shift }
        
    @property
    def shape ( self ) :
        """`shape` : shape/sigma-parameter for LogNormal  distributions"""
        return self.__shape
    @shape.setter 
    def shape ( self , value ) :
        self.set_value ( self.__shape , value )

models.append ( LogNormal_pdf )

# =============================================================================
## @class ExpoLog_pdf
#  Exponential-logarithmic distribution
#  @see https://en.wikipedia.org/wiki/Exponential-logarithmic_distribution
# - We have added a shift parameter
# - to ensure \f$ 0 < p < 1 \f$  we use
# \f$ p = \frac{1}{2}\left[ 1 + \tanh \psi \right] \f$ 
# see Ostap::Math::ExpoLog 
class ExpoLog_pdf(PDF1,Shift,Beta) :
    """ Log-normal distribution
    - see https://en.wikipedia.org/wiki/Log-normal_distribution
    - We add here an shift parameter
    - and use "mu = log(scale)"
    - and we use "sigma" as shape parameter
    - see Ostap::Models::LogNormal 
    - see Ostap::Math::LogNormal 
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name 
                   beta  = 1  ,   
                   psi   = 0  ,   
                   shift = 0  ) : ## shift  
        
        ##
        PDF1  .__init__ ( self , name  = name  , xvar = xvar  )
        Shift .__init__ ( self , shift = shift )
        Beta  .__init__ ( self , beta  = beta  ) 
        ## 
        self.__psi  = self.make_var ( psi ,
                                      'psi_%s'         % self.name ,
                                      '#psi_{EL}(%s)'  % self.name ,
                                      None , psi , -20 , 20  ) 
        
        
        ## create PDF 
        self.pdf  = Ostap.Models.ExpoLog (
            self.roo_name ( 'expolog_' )  ,
            'Expo-Log%s' % self.name , 
            self.x     ,
            self.beta  ,
            self.psi   ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'beta'  : self.beta  ,            
            'psi'   : self.psi   ,            
            'shift' : self.shift }
        
    @property
    def psi ( self ) :
        """`psi` : psi-parameter for Expo-Log distribution: p = [ 1+ tanh ( psi ) ] / 2 """
        return self.__psi
    @psi.setter 
    def psi ( self , value ) :
        self.set_value ( self.__psi , value )
        
models.append ( ExpoLog_pdf )

# =============================================================================
## @class Davis_pdf
#  Davis distribution 
#  @see https://en.wikipedia.org/wiki/Davis_distribution
#  @see Ostap::Math::Davis 
class Davis_pdf(PDF1,ShiftAndScale) :
    """ Davis distribution 
    - see https://en.wikipedia.org/wiki/Davis_distribution
    - see Ostap::Math::Davis 
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name 
                   n     = 10 ,   ## shape parametr
                   b     = 1  ,   ## scale parameter
                   mu    = 0  ) : ## shift/location parameer
        ## 
        PDF1         .__init__ ( self  , name = name , xvar = xvar  )
        ShiftAndScale.__init__ ( self  ,
                                 scale       = b  ,
                                 shift       = mu , 
                                 scale_name  = 'b_%s'        % self.name , 
                                 scale_title = 'b_{D}(%s)'   % self.name , 
                                 shift_name  = 'mu_%s'       % self.name , 
                                 shift_title = '#mu_{D}(%s)' % self.name )
        
        #
        self.__n   = self.make_var ( n ,
                                     'n_%s'       % self.name ,
                                     'n_{D}(%s)'  % self.name ,
                                     None , n , 1.e-3 , 100 )

        ## create PDF 
        self.pdf  = Ostap.Models.Davis (
            self.roo_name ( 'davis_' )  ,
            'Davis %s' % self.name , 
            self.x     ,
            self.b     ,
            self.n     ,
            self.mu    )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'n'     : self.n     ,            
            'b'     : self.b     ,            
            'mu'    : self.mu    }

    ## aliases
    b  = ShiftAndScale.scale
    mu = ShiftAndScale.shift
        
    @property
    def n ( self ) :
        """`n` : n/shape-parameter for Davis distribution"""
        return self.__n
    @n.setter 
    def n ( self , value ) :
        self.set_value ( self.__n , value )

models.append ( Davis_pdf )

# =============================================================================
## @class Kumaraswami_pdf
#  Kumaraswami distribution with scale and shift
#  @see https://en.wikipedia.org/wiki/Kumaraswamy_distribution
#  @see Ostap::Models::Kumaraswami 
#  @see Ostap::Math::Kumaraswami 
class Kumaraswami_pdf(PDF1,ShiftAndScale,AlphaAndBeta) :
    """ Kumaraswami distribution with scale and shift
    - see https://en.wikipedia.org/wiki/Kumaraswamy_distribution
    - see Ostap::Models::Kumaraswami 
    - see Ostap::Math::Kumaraswami 
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,  ## the variable
                   name  = '' ,  ## the name 
                   alpha = 2  ,
                   beta  = 2  , 
                   scale = 1  ,
                   shift = 0  ) :
        ## 
        PDF1         .__init__ ( self , name  = name  , xvar  = xvar  )
        ShiftAndScale.__init__ ( self , scale = scale , shift = shift )
        AlphaAndBeta .__init__ ( self , alpha = alpha , beta  = beta  ) 
        ## 
        ## create PDF 
        self.pdf  = Ostap.Models.Kumaraswami(
            self.roo_name ( 'kumaraswami_' )  ,
            'Kumaraswami %s' % self.name , 
            self.x     ,
            self.alpha ,
            self.beta  ,
            self.scale ,
            self.shift )
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'alpha' : self.alpha ,            
            'beta'  : self.beta  ,            
            'scale' : self.scale ,            
            'shift' : self.shift }
        
models.append ( Kumaraswami_pdf )

# =============================================================================
## @class InverseGamma_pdf
#  Inverse Gamma distribution (with shift)
#  @see https://en.wikipedia.org/wiki/Inverse-gamma_distribution
#  @see Ostap::Math::InverseGamma
class InverseGamma_pdf(PDF1,ShiftAndScale,AlphaAndBeta) :
    """ Inverse Gamma distribution (with shift)
    - see https://en.wikipedia.org/wiki/Inverse-gamma_distribution
    - see Ostap::Math::InverseGamma
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,   ## the variable
                   name  = '' ,   ## the name 
                   alpha = 8  ,   ## shape parameter 
                   beta  = 1  ,   ## scale parameter 
                   shift = 0  ) : ## shift parameter 
        ## 
        PDF1         .__init__ ( self , name = name , xvar = xvar  )
        ShiftAndScale.__init__ ( self ,
                                 scale = beta  , 
                                 shift = shift , 
                                 scale_name  = 'beta_%s'        % self.name , 
                                 scale_title = '#beta_{IG}(%s)' % self.name )
        AlphaAndBeta.__init__  ( self , alpha = alpha , beta = self.scale )
        ## 
        ## create PDF 
        self.pdf  = Ostap.Models.InverseGamma (
            self.roo_name ( 'invgamma_' )  ,
        'Inverse-Gamma %s' % self.name , 
            self.x     ,
            self.alpha ,
            self.beta  ,
            self.shift )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'alpha' : self.alpha ,            
            'beta'  : self.beta  ,            
            'shift' : self.shift }
        
    ## ALIAS 
    beta = ShiftAndScale.scale 
                                
models.append ( InverseGamma_pdf )

# =============================================================================
## @class Burr_pdf
#  Type XII Burr distribution
#  @see https://en.wikipedia.org/wiki/Burr_distribution
#
#  We have added two parameters: 
#  - scale
#  - shift 
#  @see Ostap::Models::Burr 
#  @see Ostap::Math::Burr
class Burr_pdf(PDF1,ShiftAndScale) :
    """ Burr Type XII distribution with scale and shift
    - see https://en.wikipedia.org/wiki/Burr_distribution
    - see Ostap::Models::Burr 
    - see Ostap::Math::Burr 
    """
    ## constructor
    def __init__ ( self       , * , 
                   xvar       ,  ## the variable
                   name  = '' ,  ## the name 
                   c     = 1  ,
                   k     = 8  , 
                   scale = 1  ,
                   shift = 0  ) :
        ## 
        PDF1         .__init__ ( self  , name  = name  , xvar  = xvar )
        ShiftAndScale.__init__ ( self  , scale = scale , shift = shift )
        #
        self.__c   = self.make_var ( c ,
                                     'c_%s'         % self.name ,
                                     'c_{Burr}(%s)' % self.name ,
                                     None , c , 1.e-5 , 1000  )
        self.__k   = self.make_var ( k ,
                                     'k_%s'         % self.name ,
                                     'k_{Burr}(%s)' % self.name ,
                                     None , k , 1.e-5 , 1000  )
        
        ## create PDF 
        self.pdf  = Ostap.Models.Burr (
            self.roo_name ( 'burr_' )  ,
            'Burr Type XII  %s' % self.name , 
            self.x     ,
            self.c     ,
            self.k     ,
            self.scale ,
            self.shift )
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'c'     : self.c     ,            
            'k'     : self.k     ,            
            'scale' : self.scale ,            
            'shift' : self.shift }
        
    @property
    def c ( self ) :
        """`c` : a-parameter for Burr Type XII distribution"""
        return self.__c
    @c.setter 
    def c ( self , value ) :
        self.set_value ( self.__c , value )
        
    @property
    def k ( self ) :
        """`k` : b-parameter for Burr Type XII distribution"""
        return self.__k
    @k.setter 
    def k ( self , value ) :
        self.set_value ( self.__k , value )

models.append ( Burr_pdf )

# =============================================================================
## @class Tsallis_pdf
#  Useful function to describe pT-spectra of particles 
#
#  - C. Tsallis, 
#  "Possible generalization of Boltzmann-Gibbs statistics,
#  J. Statist. Phys. 52 (1988) 479.
#  - C. Tsallis, 
#  Nonextensive statistics: theoretical, experimental and computational 
#  evidences and connections, Braz. J. Phys. 29 (1999) 1.
# 
#  \f[ \frac{d\sigma}{dp_T} \propto  
#    p_T\times \left( 1 + \frac{E_{kin}}{Tn}\right)^{-n}\f],
#  where \f$E_{kin} = \sqrt{p_T^2-M^2}-M\f$ 
#  is transverse kinetic energy 
#
#  @see Ostap::Models::Tsallis
#  @see Ostap::Math::Tsallis
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Tsallis_pdf(PDF1) :
    r""" Useful function to describe pT-spectra of particles 
    
    - C. Tsallis, 
    Possible generalization of Boltzmann-Gibbs statistics,
    J. Statist. Phys. 52 (1988) 479.
    - C. Tsallis, 
    Nonextensive statistics: theoretical, experimental and computational 
    evidences and connections, Braz. J. Phys. 29 (1999) 1.
    
    \f[ \frac{d\sigma}{dp_T} \propto  
    p_T\times \left( 1 + \frac{E_{kin}}{Tn}\right)^{-n}\f],
    
    where \f$E_{kin} = \sqrt{p_T^2-M^2}-M\f$
    
    is transverse kinetic energy 
    """
    def __init__ ( self                   , * , 
                   xvar                   ,   ## pT-variable (for fitting) 
                   name      = ''         , 
                   m0        = 0.135      ,   ## particle mass (may be fixed)
                   n         = None       ,   ## shape parameter
                   T         = None       ) : ## temperature parameter                   

        ## initialize the base 
        PDF1.__init__  ( self , name = name , xvar = xvar )
        
        self.__m0   = self.make_var ( m0                   ,
                                      'm0_%s'  % self.name , 
                                      'm0(%s)' % self.name ,
                                      True , 1.e-10   , 1e+6 )
        
        self.__n    = self.make_var ( n               ,
                                      'n_%s'   % self.name , 
                                      'n(%s) ' % self.name ,
                                      False , 1 , 0.01  , 1000 )  
        
        self.__T    = self.make_var ( T               ,
                                      'T_%s'   % self.name , 
                                      'T(%s) ' % self.name ,
                                      False , 1 , 1.e-4 , 1e+6 )
        
        self.pdf  = Ostap.Models.Tsallis (
            self.roo_name ( 'tsallis_' )   ,
            'Tsallis %s' % self.name  , 
            self.pt               ,
            self.n                ,
            self.T                ,
            self.m0               )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'n'     : self.n     ,            
            'T'     : self.T     ,            
            'm0'    : self.m0    ,            
            }
    
    @property
    def pt ( self ) :
        """`pt'-variable for Tsallis distribution (the same as ``x'')"""
        return self.xvar
    
    @property
    def m0 ( self ) :
        """`m0'-parameter of Tsallis' function"""
        return self.__m0
    @m0.setter
    def m0 ( self , value ) :
        self.set_value ( self.__m0 , value )

    @property
    def n ( self ) :
        """`n'-parameter of Tsallis' function"""
        return self.__n
    @n.setter
    def n ( self , value ) :
        self.set_value ( self.__n , value )

    @property
    def T ( self ) :
        """`T'-parameter of Tsallis' function"""
        return self.__T
    @T.setter
    def T ( self , value ) :
        self.set_value ( self.__T , value )
 
        
models  .append ( Tsallis_pdf )
spectra .append ( Tsallis_pdf )

# =============================================================================
## @class QGSM_pdf
#  Useful function to describe pT-spectra of particles 
#
# - A. B. Kaidalov and O. I. Piskunova, Z. Phys. C 30 (1986) 145.
# - O. I. Piskounova, arXiv:1301.6539 [hep-ph]; 
# - O. I. Piskounova, arXiv:1405.4398 [hep-ph].
# - A. A. Bylinkin and O. I. Piskounova, 
#  "Transverse momentum distributions of baryons at LHC energies",
#  arXiv:1501.07706.
#
#  \f[ \frac{d\sigma}{dp_T} \propto 
#  p_T \times \mathrm{e}^{ -b_0 (m_T-m)} \f], 
#  where transverse mass is defined as \f$m_T = \sqrt{p_T^2+m^2}\f$
# 
#  @see Ostap::Models::QGSM
#  @see Ostap::Math::QGSM
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class QGSM_pdf(PDF1) :
    r""" Useful function to describe pT-spectra of particles 
    
    - A. B. Kaidalov and O. I. Piskunova, Z. Phys. C 30 (1986) 145.
    - O. I. Piskounova, arXiv:1301.6539 [hep-ph]; 
    - O. I. Piskounova, arXiv:1405.4398 [hep-ph].
    - A. A. Bylinkin and O. I. Piskounova, 
    'Transverse momentum distributions of baryons at LHC energies',
    arXiv:1501.07706.

    \f[ \frac{d\sigma}{dp_T} \propto p_T \times \mathrm{e}^{ -b_0 (m_T-m)} \f],
    
    where transverse mass is defined as \f$m_T = \sqrt{p_T^2+m^2}\f$
    """
    def __init__ ( self             ,
                   name             , 
                   xvar             ,   ## pT-variable (for fitting) 
                   m0        = 0    ,   ## particle mass (may be fixed)
                   b         = None ) : ## slope parameter
        
        ## initialize the base 
        PDF1.__init__  ( self , name = name , xvar = xvar )

        
        self.__m0   = self.make_var ( m0              ,
                                      'm0_%s'  % self.name , 
                                      'm0(%s)' % self.name ,
                                      True , 0  , 1e+6 )
        
        self.__b    = self.make_var ( b               ,
                                      'b_%s'   % self.name , 
                                      'b(%s) ' % self.name ,
                                      False , 0. , 1e+6 )  
        
        self.pdf  = Ostap.Models.QGSM (
            self.roo_name ( 'qgsm_' ) ,
            'QGSM %s' % self.name , 
            self.pt               ,
            self.b                ,
            self.m0               )

        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'b'     : self.b     ,            
            'm0'    : self.m0    ,            
            }
    
    @property
    def pt ( self ) :
        """`pt'-variable for QGSM distribution (the same as ``x'')"""
        return self.xvar
    
    @property
    def m0 ( self ) :
        """`m0'-parameter of QGSM function"""
        return self.__m0
    @m0.setter
    def m0 ( self , value ) :
        self.set_value ( self.__m0 , value )

    @property
    def b ( self ) :
        """`b'-parameter of QGSM function"""
        return self.__b
    @b.setter
    def b ( self , value ) :
        self.set_value ( self.__b , value )

models  .append ( QGSM_pdf )
spectra .append ( QGSM_pdf )

# =============================================================================
## @class Hagedorn_pdf
#  Useful function to describe pT-spectra of particles 
#  @see R.Hagedorn, "Multiplicities, p_T distributions and the 
#       expected hadron \to Quark - Gluon Phase Transition", 
#       Riv.Nuovo Cim. 6N10 (1983) 1-50
#  @see https://doi.org/10.1007/BF02740917 
#  @see https://inspirehep.net/literature/193590
#  
#  \f[ f(p_T; m, T) \propto 
#   p_T \sqrt{p^2_T + m^2} K_1( \beta \sqrt{ p^2_T+m^2} ) \f] 
#
#  where \f$ \beta \f$ is inverse temporature 
#  \f$ \beta = \frac{1}{T} f$ 
#
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @date 2022-12-06
#  @see Ostap::Models::Hagedorn
#  @see Ostap::Math::Hagedorn
class Hagedorn_pdf(PDF1) :
    r"""Useful function to describe pT-spectra of particles 
    
    - see R.Hagedorn, 'Multiplicities, p_T distributions and the 
    expected hadron \to Quark - Gluon Phase Transition', 
    Riv.Nuovo Cim. 6N10 (1983) 1-50
    - see https://doi.org/10.1007/BF02740917 
    - see https://inspirehep.net/literature/193590
    
    \f[ f(p_T; m, T) \propto 
    p_T \sqrt{p^2_T + m^2} K_1( \beta \sqrt{ p^2_T+m^2} ) \f] 
    
    where \f$ \beta \f$ is inverse temporature 
    \f$ \beta = \frac{1}{T} f$ 
    """
    def __init__ ( self             ,
                   name             , 
                   xvar             ,   ## pT-variable (for fitting) 
                   m0        = 0    ,   ## particle mass (may be fixed)
                   beta      = None ) : ## inverse temperature
        
        ## initialize the base 
        PDF1.__init__  ( self , name  = name  , xvar = xvar )

        
        self.__m0   = self.make_var ( m0              ,
                                      'm0_%s'  % self.name , 
                                      'm0(%s)' % self.name ,
                                      True , 0  , 1e+6 )
        
        self.__beta = self.make_var ( beta                 ,
                                      'beta_%s'    % self.name  , 
                                      '#beta(%s) ' % self.name  ,
                                      False , 1.e-6 , 1e+6 )  
        
        self.pdf  = Ostap.Models.Hagedorn (
            self.roo_name ( 'hage_' ) ,
            'Hagedorn %s' % self.name , 
            self.pt               ,
            self.beta             ,
            self.m0               )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'beta'  : self.beta  ,             
            'm0'    : self.m0    ,            
            }
    
    @property
    def pt ( self ) :
        """'pt'-variable for Hagedorn distribution (the same as ``x'')"""
        return self.xvar
    
    @property
    def m0 ( self ) :
        """'m0'-parameter of Hagedorn function"""
        return self.__m0
    @m0.setter
    def m0 ( self , value ) :
        self.set_value ( self.__m0 , value )

    @property
    def beta ( self ) :
        """'beta'-parameter (inverse temperature) of Hagedorn  function"""
        return self.__beta
    @beta.setter
    def beta ( self , value ) :
        self.set_value ( self.__beta , value )
       
models  .append ( Hagedorn_pdf )
spectra .append ( Hagedorn_pdf )

# =============================================================================
## @class Tsallis2_pdf
#  Useful function to describe pT and saidity -spectra of particles 
#
#  2D particle density distribution as function of pt and rapidity 
#  @see L. Marques, J. Cleymans, A. Deppman, 
#       "Description of High-Energy pp Collisions 
#        Using Tsallis Thermodynamics: 
#        Transverse Momentum and Rapidity Distributions", 
#        Phys. Rev. D 91, 054025, 	arXiv:1501.00953 
#  @see https://arxiv.org/abs/1501.00953
#  @see https://doi.org/10.1103/PhysRevD.91.054025
#  @see Ostap::Models::Tsallis2
#  @see Ostap::Models::Tsallis
#  @see Ostap::Math::Tsallis2
#  @see Ostap::Math::Tsallis
class Tsallis2_pdf(PDF2) :
    """ 2D particle density distribution as function of pt and rapidity 
    @see L. Marques, J. Cleymans, A. Deppman, 
    ``Description of High-Energy pp Collisions 
    Using Tsallis Thermodynamics: 
    Transverse Momentum and Rapidity Distributions'', 
    Phys. Rev. D 91, 054025, 	arXiv:1501.00953 
    - see https://arxiv.org/abs/1501.00953
    - see https://doi.org/10.1103/PhysRevD.91.054025
    - see `Ostap.Models.Tsallis2`
    - see `Ostap.Models.Tsallis`
    - see `Ostap.Math.Tsallis2`
    - see `Ostap.Math.Tsallis`
    """
    def __init__ ( self                   ,
                   xvar                   ,   ## pT-observable )
                   yvar                   ,   ## rapidity observable
                   m0        = 1          ,   ## partile mass (presumably constant) 
                   T         = None       ,   ## temperature parameter
                   q         = 1.1        ,   ## q-parameter
                   mu        = 0          ,   ## chemical potential                    
                   name      = ''         ) :
        
        ## initialize the base 
        PDF2.__init__  ( self , name = name , xvar = xvar , yvar =  yvar )
        
        self.__m0   = self.make_var ( m0              ,
                                      'm0_%s'  % self.name , 
                                      'm0(%s)' % self.name ,
                                      True , 0     , 1e+6 )
        
        self.__q    = self.make_var ( q               ,
                                      'q_%s'   % self.name , 
                                      'q(%s) ' % self.name ,
                                      None , 1.1 , 1.e-6, 10 )  
        
        self.__T    = self.make_var ( T               ,
                                      'T_%s'   % self.name , 
                                      'T(%s) ' % self.name ,
                                      None , 1 , 1.e-4 , 1e+6 )
        
        self.__mu   = self.make_var ( mu               ,
                                      'mu_%s'   % self.name , 
                                      '#mu(%s) '% self.name ,
                                      True  , 0 , 0 , 100 )
        
        self.pdf  = Ostap.Models.Tsallis2 (
            self.roo_name ( 'tsallis2_' )   ,
            'Tsallis2 %s' % self.name  , 
            self.pt               ,
            self.y                ,
            self.m0               ,
            self.T                ,
            self.q                ,
            self.mu               )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'yvar'  : self.yvar  ,
            'm0'    : self.m0    ,            
            'q'     : self.q     ,            
            'T'     : self.T     ,            
            'mu'    : self.mu    ,            
            }
        
    @property
    def pt ( self ) :
        """'pt'-observable (transverse momentum) for Tsallis distribution (the same as 'xvar')"""
        return self.xvar
    
    @property
    def y ( self ) :
        """'y'-observable (rapidity) for Tsallis distribution (the same as 'yvar')"""
        return self.yvar

    @property
    def m0 ( self ) :
        """'m0'-parameter (particle mass) of Tsallis' function"""
        return self.__m0
    @m0.setter
    def m0 ( self , value ) :
        self.set_value ( self.__m0 , value )

    @property
    def q ( self ) :
        """'q'-parameter (shape) of Tsallis' function"""
        return self.__q
    @q.setter
    def q ( self , value ) :
        self.set_value ( self.__q , value )

    @property
    def T ( self ) :
        """'T'-parameter (temperature) of Tsallis' function"""
        return self.__T
    @T.setter
    def T ( self , value ) :
        self.set_value ( self.__T , value )

    @property
    def mu ( self ) :
        """'mu'-parameter (chemical potential) of Tsallis' function"""
        return self.__mu
    @mu.setter
    def mu ( self , value ) :
        self.set_value ( self.__mu , value )
 
        
models  .append ( Tsallis2_pdf ) 
spectra .append ( Tsallis2_pdf )

# =============================================================================
if '__main__' == __name__ : 
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )

# =============================================================================
##                                                                      The END 
# =============================================================================
