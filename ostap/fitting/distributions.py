#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/distributions.py
#  A set of various smooth shapes and PDFs 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
# =============================================================================
"""A set of various smooth shapes and PDFs

- GammaDist_pdf      : Gamma-distributuon in shape/scale parameterization
- GenGammaDist_pdf   : Generalized Gamma-distribution
- Amoroso_pdf        : another view of generalized Gamma distribution
- LogGammaDist_pdf   : Gamma-distributuon in shape/scale parameterization
- Log10GammaDist_pdf : Gamma-distributuon in shape/scale parameterization
- LogGamma_pdf       : Log-Gamma distribution  
- BetaPrime_pdf      : Beta-prime distribution 
- Landau_pdf         : Landau distribution 
- Argus_pdf          : ARGUS distribution 
- GenArgus_pdf       : Generalized ARGUS distribution 
- TwoExpos_pdf       : Difference of two exponents
- Gumbel_pdf         : Gumbel distributions
- Rice_pdf           : Rice distribution
- GenInvGauss_pdf    : Generalized Inverse Gaussian distribution
- Weibull_pdf        : Weibull distributions
- Tsallis_pdf        : Tsallis PDF 
- QGSM_pdf           : QGSM PDF 
- Tsallis2_pdf       : 2D Tsallis PDF 

"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
__all__     = (
    'GammaDist_pdf'      , ## Gamma-distributuon in shape/scale parameterization
    'GenGammaDist_pdf'   , ## Generalized Gamma-distribution
    'Amoroso_pdf'        , ## another view of generalized Gamma distribution
    'LogGammaDist_pdf'   , ## Gamma-distributuon in shape/scale parameterization
    'Log10GammaDist_pdf' , ## Gamma-distributuon in shape/scale parameterization
    'LogGamma_pdf'       , ## Log-Gamma distribution 
    'BetaPrime_pdf'      , ## Beta-prime distribution 
    'Landau_pdf'         , ## Landau distribution 
    'Argus_pdf'          , ## ARGUS distribution 
    'GenArgus_pdf'       , ## Generalized ARGUS distribution 
    'TwoExpos_pdf'       , ## difference of two exponents
    'Gumbel_pdf'         , ## Gumbel distributions
    'Rice_pdf'           , ## Rice distribution 
    'GenInvGauss_pdf'    , ## Generalized Inverse Gaussian distribution
    'Weibull_pdf'        , ## Weibull distributions
    'Tsallis_pdf'        , ## Tsallis PDF 
    'QGSM_pdf'           , ## QGSM PDF 
    'Tsallis2_pdf'       , ## 2D Tsallis PDF 
    )
# =============================================================================
from   ostap.core.core        import Ostap, VE 
from   ostap.fitting.pdfbasic import PDF1, PDF2 
import ROOT, math
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.distributions' )
else                       : logger = getLogger ( __name__                      )
# =============================================================================
models = []
# =============================================================================
## @class GammaDist_pdf
#  Gamma-distribution with shape/scale parameters
#  http://en.wikipedia.org/wiki/Gamma_distribution
#  It suits nicely for fits of multiplicity and/or chi2 distributions
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::GammaDist 
#  @see Ostap::Math::GammaDist 
class GammaDist_pdf(PDF1) :
    """Gamma-distribution with shape/scale parameters
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
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   k     = None     ,   ## k-parameter
                   theta = None     ) : ## theta-parameter
        #
        PDF1.__init__ ( self , name , xvar )
        #
        self.__k     = self.make_var ( k       ,
                                       'k_%s'                % name ,
                                       'k_{#Gamma}(%s)'      % name ,
                                       None , 1 , 1.e-3 , 100 )
        self.__theta = self.make_var ( theta   ,
                                       'theta_%s'            % name ,
                                       '#theta_{#Gamma}(%s)' % name ,
                                       None , 1 , 1.e-3 , 100 )
        
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
            }

    @property
    def k ( self ) :
        """``k''-parameter of Gamma distribution   (k>0)"""
        return self.__k
    @k.setter 
    def k ( self , value ) :
        value = float ( value )
        assert 0 < value, '``k''-value must be positive'
        self.__k.setVal ( value ) 
        return self.__k.getVal() 

    @property
    def theta ( self ) :
        """``theta''-parameter of Gamma distribution   (theta>0)"""
        return self.__theta
    @theta.setter 
    def theta ( self , value ) :
        value = float ( value )
        assert 0 < value, "``theta''-value must be positive"
        self.__theta.setVal ( value ) 
        return self.__theta.getVal() 

   

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
class GenGammaDist_pdf(PDF1) :
    """Generalized Gamma-distribution with additional shift parameter 
    http://en.wikipedia.org/wiki/Generalized_gamma_distribution
    Special cases : 
    - p == 1      : Gamma  distribution
    - p == k      : Weibull distribution
    - p == k == 1 : Exponential distribution
    - p == k == 2 : Rayleigh    distribution
    """
    ## constructor
    def __init__ ( self                ,
                   name                ,   ## the name 
                   xvar                ,   ## the variable
                   k     = None  ,   ## k-parameter
                   theta = None  ,   ## theta-parameter
                   p     = None  ,   ## p-parameter
                   low   = None  ) : ## low-parameter
        #
        PDF1.__init__ ( self , name , xvar )
        #
        self.__k     = self.make_var ( k       ,
                                       'k_%s'                % name ,
                                       'k_{#Gamma}(%s)'      % name ,
                                       None , 1 , 1.e-3 , 100 )
        self.__theta = self.make_var ( theta   ,
                                       'theta_%s'            % name ,
                                       '#theta_{#Gamma}(%s)' % name , theta ,
                                       None , 1.e-3 , 100 )
        self.__p     = self.make_var ( p       ,
                                       'p_%s'                % name ,
                                       'p_{#Gamma}(%s)'      % name ,
                                       NNone , 1 , 1.e-3 ,   6 )

        limits_low = ()
        if   self.xminmax() :
            mn , mx   = self.xminmax()
            limits_low = mn , mn , mx
            
        self.__low   = self.make_var ( low      ,
                                       'low_%s'         % name ,
                                       'l_{#Gamma}(%s)' % name ,
                                       None , *limits_low )
        
        self.pdf  = Ostap.Models.GenGammaDist (
            self.roo_name ( 'ggamma_' ) ,
            'Generalized Gamma distribution %s' % self.name , 
            self.x         ,
            self.k         ,
            self.theta     ,
            self.p         , 
            self.low       )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'k'     : self.k     ,
            'theta' : self.theta ,            
            'p'     : self.p     ,            
            'low'   : self.low   ,            
            }

    @property
    def k ( self ) :
        """``k''-parameter of generalized Gamma distribution   (k>0)"""
        return self.__k
    @k.setter 
    def k ( self , value ) :
        value = float ( value )
        assert 0 < value, "``k''-value must be positive"
        self.__k.setVal ( value ) 
    
    @property
    def theta ( self ) :
        """``theta''-parameter of generalized Gamma distribution   (theta>0)"""
        return self.__theta
    @theta.setter 
    def theta ( self , value ) :
        value = float ( value )
        assert 0 < value, "``theta''-value must be positive"
        self.__theta.setVal ( value ) 

    @property
    def p ( self ) :
        """``p''-parameter of generalized Gamma distribution   (p>0)"""
        return self.__p
    @p.setter 
    def p ( self , value ) :
        value = float ( value )
        assert 0 < value, "``p''-value must be positive"
        self.__p.setVal ( value ) 

    @property
    def low ( self ) :
        """``x-low''-parameter of generalized Gamma distribution   (f(x)=0., for x < x_low)"""
        return self.__low
    @low.setter 
    def low ( self , value ) :
        value = float ( value )
        self.__low.setVal ( value ) 

models.append ( GenGammaDist_pdf ) 
# =============================================================================
## @class Amoroso_pdf
#  Another view on generalized gamma distribution
#  http://arxiv.org/pdf/1005.3274
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Math::Amoroso
#  @see Ostap::Models::Amoroso
class Amoroso_pdf(PDF1) :
    """Another view on generalized gamma distribution
    http://arxiv.org/pdf/1005.3274
    """
    ## constructor
    def __init__ ( self          ,
                   name          ,   ## the name 
                   xvar          ,   ## the variable
                   theta = None  ,   ## theta-parameter
                   alpha = None  ,   ## alpha-parameter
                   beta  = None  ,   ## beta-parameter
                   a     = None  ) : ## a-parameter
        
        #
        PDF1.__init__ ( self , name , xvar )
        #
        
        self.__theta = self.make_var ( theta   ,
                                       'theta_%s'             % name ,
                                       '#theta_{Amoroso}(%s)' % name ,
                                       None , 1 , 1.e-3 , 100 )
        self.__alpha = self.make_var ( alpha   ,
                                       'alpha_%s'             % name ,
                                       '#alpha_{Amoroso}(%s)' % name ,
                                       None , 1 , 1.e-3 , 100 )
        self.__beta  = self.make_var ( beta    ,
                                       'beta_%s'              % name ,
                                       '#beta_{Amoroso}(%s) ' % name ,
                                       None , 1 , 1.e-3 ,  10 )        
        self.__a     = self.make_var ( a       ,
                                       'a_%s'                 % name ,
                                       'a_{Amoroso}(%s)'      % name ,
                                       None , 1 , -10   ,  10  )
        
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
                                     
    @property
    def theta ( self ) :
        """``theta''-parameter of Amoroso function (theta>0)"""
        return self.__theta
    @theta.setter 
    def theta ( self , value ) :
        value = float ( value )
        assert 0 < value, "``theta''-value must be positive"
        self.__theta.setVal ( value ) 

    @property
    def alpha ( self ) :
        """``alpha''-parameter of Amoroso function (alpha>0)"""
        return self.__alpha
    @alpha.setter 
    def alpha ( self , value ) :
        value = float ( value )
        assert 0 < value, "``alpha''-value must be positive"
        self.__alpha.setVal ( value ) 

    @property
    def beta ( self ) :
        """``beta''-parameter of Amoroso function (beta>0)"""
        return self.__beta
    @beta.setter 
    def beta ( self , value ) :
        value = float ( value )
        assert 0 < value, "``beta''-value must be positive"
        self.__beta.setVal ( value ) 

    @property
    def a ( self ) :
        """``a''-parameter of Amoroso function"""
        return self.__a
    @a.setter 
    def a ( self , value ) :
        value = float ( value )
        assert 0 < value, "``a''-value must be positive"
        self.__a.setVal ( value ) 


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
class LogGammaDist_pdf(PDF1) :
    """Distribution for log(x), where x follows Gamma distribution
    It suits nicely for fits of log(multiplicity) and/or log(chi2) distributions
    """
    ## constructor
    def __init__ ( self         ,
                   name         ,   ## the name 
                   xvar         ,   ## the variable
                   k     = None ,   ## k-parameter
                   theta = None ) : ## theta-parameter
        #
        PDF1.__init__ ( self , name , xvar )
        #
        self.__k     = self.make_var ( k       ,
                                       'k_%s'                   % name ,
                                       'k_{log#Gamma}(%s)'      % name ,
                                       None , 1 , 1.e-5 , 1000 )
        self.__theta = self.make_var ( theta   ,
                                       'theta_%s'               % name ,
                                       '#theta_{log#Gamma}(%s)' % name ,
                                       None , 1 , 1.e-5 , 1000 )

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
            }

    @property
    def k ( self ) :
        """``k''-parameter of log(Gamma) distribution   (k>0)"""
        return self.__k
    @k.setter 
    def k ( self , value ) :
        value = float ( value )
        assert 0 < value, "``k''-value must be positive"
        self.__k.setVal ( value ) 

    @property
    def theta ( self ) :
        """``theta''-parameter of log(Gamma) distribution   (theta>0)"""
        return self.__theta
    @theta.setter 
    def theta ( self , value ) :
        value = float ( value )
        assert 0 < value, "``theta''-value must be positive"
        self.__theta.setVal ( value ) 

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
class Log10GammaDist_pdf(PDF1) :
    """Distribution for log10(x), where x follows Gamma distribution
    It suits nicely for fits of log10(multiplicity) and/or log10(chi2) distributions
    """
    ## constructor
    def __init__ ( self         ,
                   name         ,   ## the name 
                   xvar         ,   ## the variable
                   k     = None ,   ## k-parameter
                   theta = None ) :  ## theta-parameter
        #
        PDF1.__init__ ( self , name , xvar )
        #
        self.__k     = self.make_var ( k       ,
                                       'k_%s'                     % name ,
                                       'k_{log10#Gamma}(%s)'      % name ,
                                       None , 1 , 1.e-4 , 10000 )
        self.__theta = self.make_var ( theta   ,
                                       'theta_%s'                 % name ,
                                       '#theta_{log10#Gamma}(%s)' % name ,
                                       None , 1 , 1.e-4 , 10000 )

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
            }
    
    @property
    def k ( self ) :
        """``k''-parameter of log10(Gamma) distribution   (k>0)"""
        return self.__k
    @k.setter 
    def k ( self , value ) :
        value = float ( value )
        assert 0 < value, "``k''-value must be positive"
        self.__k.setVal ( value )

    @property
    def theta ( self ) :
        """``theta''-parameter of log10(Gamma) distribution   (theta>0)"""
        return self.__theta
    @theta.setter 
    def theta ( self , value ) :
        value = float ( value )
        assert 0 < value, "``theta''-value must be positive"
        self.__theta.setVal ( value ) 

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
class LogGamma_pdf(PDF1) :
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
    def __init__ ( self         ,
                   name         ,   ## the name 
                   xvar         ,   ## the variable
                   nu    = None ,   ## nu-parameter
                   lambd = None ,   ## lambda-parameter
                   alpha = None ) : ## nu-parameter
        #
        PDF1.__init__ ( self , name , xvar )
        #
        limits_nu = ()
        if   self.xminmax() :
            mn,mx = self.xminmax()
            dx = mx - mn
            xm = 0.5 * ( mn + mx ) , mn - 10* dx , mx + 10 *  dx 
            
        self.__nu     = self.make_var ( nu       ,
                                        'nu_%s'                    % name ,
                                        '#nu_{#log#Gamma}(%s)'     % name ,
                                        None , *limits_nu )
        
        self.__lambd  = self.make_var ( lambd      ,
                                        'lambda_%s'                % name ,
                                        '#lambda_{#log#Gamma}(%s)' % name ,
                                        None  , 2 , -1000 , 1000 )
        
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
    
    @property
    def nu ( self ) :
        """``nu''-parameter (location) of log-Gamma distribution"""
        return self.__nu
    @nu.setter 
    def nu ( self , value ) :
        value = float ( value )
        self.__nu.setVal ( value ) 
        return self.__nu.getVal() 
    
    @property
    def lambd ( self ) :
        """``lambda''-parameter (scale) of log-Gamma distribution"""
        return self.__lambd
    @lambd.setter 
    def lambd ( self , value ) :
        value = float ( value )
        self.__lambd.setVal ( value ) 
        return self.__lambd.getVal() 

    @property
    def alpha  ( self ) :
        """``alpha-parameter'' of log-Gamma distribution   (alpha>0)"""
        return self.__alpha
    @alpha.setter 
    def alpha ( self , value ) :
        value = float ( value )
        assert 0 < value, "``alpha''-value must be positive"
        self.__alpha.setVal ( value ) 
        return self.__alpha.getVal() 

models.append ( LogGamma_pdf ) 
# =============================================================================
## @class BetaPrime_pdf
#  http://en.wikipedia.org/wiki/Beta_prime_distribution
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::BetaPrime
#  @see Ostap::Math::BetaPrime
class BetaPrime_pdf(PDF1) :
    """Beta-prime disribution 
    - http://en.wikipedia.org/wiki/Beta_prime_distribution
    """
    ## constructor
    def __init__ ( self         ,
                   name         ,   ## the name 
                   xvar         ,   ## the variable
                   alpha = None ,   ## alpha-parameter
                   beta  = None ,   ## beta-parameter
                   scale = 1    ,   ## scale-parameter 
                   delta = 0    ) : ## shift-parameter 
        #
        PDF1.__init__ ( self , name , xvar )
        # 
        self.__alpha  = self.make_var ( alpha    ,
                                        'alpha_%s'                 % name ,
                                        '#alpha_{#beta#prime}(%s)' % name ,
                                        None , 1 , 1.e-3 , 1000 )
        self.__beta   = self.make_var ( beta     ,
                                        'beta_%s'                  % name ,
                                        '#beta_{#beta#prime}(%s)'  % name ,
                                        None , 1 , 1.e-3 , 1000 )        
        self.__scale  = self.make_var ( scale     ,
                                        'scale_%s'                 % name ,
                                        '#theta_{#beta#prime}(%s)' % name , 
                                        None , 1 , -1000 , 1000 )

        limits_delta = ()
        if self.xminmax() :
            mn, mx = self.xminmax()
            dx = mx - mn
            limits_delta = mn - 10 * dx , mx + 10 * dx
            
        self.__delta  = self.make_var ( delta     ,
                                        'delta_%s'                 % name ,
                                        '#delta_{#beta#prime}(%s)' % name , 
                                        None , 0 , *limits_delta )
        
        self.pdf  = Ostap.Models.BetaPrime (
            self.roo_name ( 'betap_' )   ,
            'Beta-prime %s' % self.name  , 
            self.x     ,
            self.alpha ,
            self.beta  ,
            self.scale ,
            self.delta )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'alpha' : self.alpha ,            
            'beta'  : self.beta  ,            
            'scale' : self.scale ,            
            'delta' : self.delta ,            
            }
    
    @property
    def alpha ( self ) :
        """``alpha''-parameter of Beta' distribution   (alpha>0)"""
        return self.__alpha
    @alpha.setter 
    def alpha ( self , value ) :
        value = float ( value )
        assert 0 < value, "``alpha''-value must be positive"
        self.__alpha.setVal ( value ) 

    @property
    def beta ( self ) :
        """``beta''-parameter of Beta' distribution   (alpha>0)"""
        return self.__beta
    @beta.setter 
    def beta ( self , value ) :
        value = float ( value )
        assert 0 < value, "``beta''-value must be positive"
        self.__beta.setVal ( value ) 

    @property
    def scale ( self ) :
        """``scale''-parameter of Beta' distribution"""
        return self.__scale
    @scale.setter 
    def scale ( self , value ) :
        value = float ( value )
        self.__scale.setVal ( value ) 

    @property
    def delta ( self ) :
        """``delta''-parameter of Beta' distribution"""
        return self.__delta
    @delta.setter 
    def delta ( self , value ) :
        value = float ( value )
        self.__delta.setVal ( value ) 

models.append ( BetaPrime_pdf ) 
# =============================================================================
## @class Landau_pdf
#  http://en.wikipedia.org/wiki/Landau_distribution
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::Landau
#  @see Ostap::Math::Landau
class Landau_pdf(PDF1) :
    """Landau distribution 
    - http://en.wikipedia.org/wiki/Landau_distribution
    """
    ## constructor
    def __init__ ( self      ,
                   name      ,   ## the name 
                   xvar      ,   ## the variable
                   scale = 1 ,   ## scale-parameter 
                   delta = 0 ) : ## shift-parameter 
        #
        PDF1.__init__ ( self , name , xvar ) 
        #
        self.__scale  = self.make_var ( scale     ,
                                        'scale_%s'            % name ,
                                        '#theta_{Landau}(%s)' % name ,
                                        None , 1 , -1000 , 1000 )
        
        
        limits_delta = ()
        if self.xminmax() :
            mn, mx = self.xminmax()
            dx = mx - mn
            limits_delta = mn - 10 * dx , mx + 10 * dx
            
        self.__delta  = self.make_var ( delta     ,
                                        'delta_%s'            % name ,
                                        '#delta_{Landau}(%s)' % name , 
                                        None , 0 , *limits_delta )
        self.pdf  = Ostap.Models.Landau (
            self.roo_name ( 'landau_' )   ,
            'Landau %s' % self.name  , 
            self.x     ,
            self.scale ,
            self.delta )
        
        ## save the configuration:
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'scale' : self.scale ,            
            'delta' : self.delta ,            
            }
        
    @property
    def scale ( self ) :
        """``scale''-parameter of Landau distribution"""
        return self.__scale
    @scale.setter 
    def scale ( self , value ) :
        value = float ( value )
        self.__scale.setVal ( value ) 

    @property
    def delta ( self ) :
        """``delta''-parameter of Landau distribution"""
        return self.__delta
    @delta.setter 
    def delta ( self , value ) :
        value = float ( value )
        self.__delta.setVal ( value ) 


models.append ( Landau_pdf )

# =============================================================================
## @class Argus_pdf
#  http://en.wikipedia.org/wiki/ARGUS_distribution
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::Argus
#  @see Ostap::Math::Argus
class Argus_pdf(PDF1) :
    """Argus distribution
    - http://en.wikipedia.org/wiki/ARGUS_distribution
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   chi              ,
                   c                , 
                   mu    = None     ) :
        #
        PDF1.__init__ ( self , name , xvar ) 
        #
        self.__c  = self.make_var ( c     ,
                                    'c_%s'          % name ,
                                    'c_{Argus}(%s)' % name ,
                                    None  , 1 , 1.e-6 , 200 )
        
        if mu is None :  self.__mu = self.c
        else          :
            lims = self.xminmax() 
            self.__mu = self.make_var ( mu     ,
                                        'mu_%s'           % name ,
                                        '#mu_{Argus}(%s)' % name ,
                                        None , *self.xminmax() )
            
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

    @property
    def mu ( self ) :
        """``mu''-parameter of Argus distribution"""
        return self.__mu
    @mu.setter 
    def mu ( self , value ) :
        self.set_value ( self.__mu , value )
        
    @property
    def c  ( self ) :
        """``c''-parameter of Argus distribution"""
        return self.__c
    @c.setter 
    def mu ( self , value ) :
        self.set_value ( self.__c , value )

    @property
    def chi ( self ) :
        """``chi''-parameter of Argus distribution"""
        return self.__chi
    @chi.setter 
    def mu ( self , value ) :
        self.set_value ( self.__chi , value )

models.append ( Argus_pdf ) 


# =============================================================================
## @class GenArgus_pdf
#  http://en.wikipedia.org/wiki/ARGUS_distribution
#  Generalised Argus 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::GenArgus
#  @see Ostap::Math::GenArgus
class GenArgus_pdf(Argus_pdf) :
    """Generalized Argus distribution
    - http://en.wikipedia.org/wiki/ARGUS_distribution
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   chi              ,
                   c                ,
                   dp    = 1.5      ,   ## corresponds to regular Argus 
                   mu    = None     ) :
        #
        Argus_pdf.__init__ ( self        ,
                             name = name ,
                             xvar = xvar ,
                             c    = c    ,
                             chi  = chi  ,
                             mu   = mu   ) 
                             
        #
        self.__dp  = self.make_var ( dp     ,
                                     'dp_%s'               % name ,
                                     '#deltap_{Argus}(%s)' % name ,
                                     None , 1.5 , 1.e-5 , 20 )
        
        ## create PDF 
        self.pdf  = Ostap.Models.GenArgus (
            self.roo_name ( 'gargus_' )   ,
            'Generlized ARGUS %s' % self.name  , 
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
        """``dp''-parameter of Generalized Argus distribution (p=dp-1)"""
        return self.__dp
    @dp.setter 
    def dp ( self , value ) :
        self.set_value ( self.__dp , value )

models.append ( GenArgus_pdf ) 


# =============================================================================
## @class TwoExpos_pdf
#  simple difference of two exponents
#  \f$ f \propto 
#        \mathrm{e}^{-a_1    x}       -\mathrm{e}^{-a_2 x} = 
#        \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f$
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @see Ostap::Models::TwoExpos
#  @see Ostap::Math::TwoExpos
class TwoExpos_pdf(PDF1) :
    r""" Simple difference of two exponents:    
    \f$ f \propto 
    \mathrm{e}^{-a_1    x}       -\mathrm{e}^{-a_2 x} = 
    \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f$    
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   xvar             ,   ## the variable
                   alpha = None     ,   ## shape-parameter 
                   delta = None     ,   ## high-parameter 
                   x0    = 0        ) : ## low-parameter 
        #
        PDF1.__init__ ( self , name , xvar ) 
        #
        self.__alpha  = self.make_var ( alpha      ,
                                        'alpha_%s'        % name ,
                                        '#alpha_{2e}(%s)' % name , 
                                        None , 1 , 1.e-4 , 100 )
        self.__delta  = self.make_var ( delta     ,
                                        'delta_%s'        % name ,
                                        '#delta_{2e}(%s)' % name , 
                                        None , 1 , 1.e-4 , 100 )
        
        limits_x0  = ()
        if  self.xminmax() :
            mn, mx = self.xminmax()
            dm  =  mx - mn
            limits_x0 = mn , mn - 0.2 *  dm , mx  + 0.1 *  dm

        self.__x0     = self.make_var ( x0    ,
                                        'x0_%s'       % name ,
                                        'x0_{2e}(%s)' % name ,
                                        None , *limits_x0 )
        
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

    @property
    def alpha ( self ) :
        """``alpha''-parameter (slope of leading exponnet) of the TwoExpo function"""
        return self.__alpha
    @alpha.setter
    def alpha ( self , value ) :
        value = float ( value )
        assert 0 <= value, "``alpha''-parameter must be non-negative"
        self.alpha.setVal ( value ) 
    
    @property
    def delta ( self ) :
        """``delta''-parameter (second exponent slope is ``alpha+delta'') of the TwoExpo function"""
        return self.__delta
    @delta.setter
    def delta ( self , value ) :
        value = float ( value )
        assert 0 <= value, "``delta''-parameter must be non-negative"
        self.delta.setVal ( value ) 
    
    @property
    def x0 ( self ) :
        """x0-parameter of the TwoExpo function  (f(x)=0 for x<x0)"""
        return self.__x0
    @x0.setter
    def x0 ( self , value ) :
        value = float ( value )
        self.x0.setVal ( value ) 
  
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
class Gumbel_pdf(PDF1) :
    r"""Gumbel distribution
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
    def __init__ ( self      ,   
                   name      ,   ## the name 
                   xvar      ,   ## the variable 
                   mu   = 0  ,   ## shift parameter/mode
                   beta = 1  ) : ## scale parameter 
        #
        PDF1.__init__ ( self , name , xvar )
        #
        self.__mu    = self.make_var ( mu      ,
                                       'mu_%s'                  % name ,
                                       'mu_{Gumbel}(%s)'        % name , None )
        self.__beta  = self.make_var ( beta        ,
                                       'beta_%s'                % name ,
                                       'beta_{Gumbel}(%s)'      % name , None ) 
        
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
        
    @property
    def mu ( self ) :
        """``mu''-parameter (shift) of the Gumbel function"""
        return self.__mu
    @mu.setter
    def mu ( self , value ) :
        value = float ( value )
        self.__mu.setVal ( value ) 

    @property
    def beta ( self ) :
        """``beta''-parameter (scale) of the Gumbel function"""
        return self.__beta
    @beta.setter
    def mu ( self , value ) :
        value = float ( value )
        self.__beta.setVal ( value ) 
  
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
class Rice_pdf(PDF1) :
    """Rice distribution
    - see https://en.wikipedia.org/wiki/Rice_distribution
    """
    ## constructor
    def __init__ ( self         ,   
                   name         ,                               ## the name 
                   xvar         ,                               ## the variable 
                   nu       = 0 ,                               ## parameter nu 
                   varsigma = 1 ,                               ## parameter varsigma
                   shift    = ROOT.RooFit.RooConst ( 0.0 )  ) : ## shift parameter
        #
        PDF1.__init__ ( self , name , xvar )
        #

        self.__nu       = self.make_var ( nu        ,
                                          'nu_%s'                  % name ,
                                          '#nu_{Rice}(%s)'         % name ,
                                          None , 0 , 0.0    , 1.e+6 )
        self.__varsigma = self.make_var ( varsigma ,
                                          'varsigma_%s'            % name ,
                                          '#varsigma_{Rice}(%s)'   % name ,
                                          None , 1 , 1.e-6  , 1.e+4 )
        self.__shift   = self.make_var ( shift        ,
                                         'shift_%s'                % name ,
                                         'shift_{Rice}(%s)'        % name ,
                                         None , 0 , -1.e+6 , 1.e+6 )
        
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
        """``nu''-parameter of Rice function"""
        return self.__nu
    @nu.setter
    def nu ( self , value ) :
        self.set_value ( self.__nu , value )

    @property
    def varsigma ( self ) :
        """``varsigma''-parameter of Rice function"""
        return self.__varsigma
    @varsigma.setter
    def varsigma ( self , value ) :
        self.set_value ( self.__varsigma , value )

    @property
    def shift ( self ) :
        """``shift''-parameter of Rice function"""
        return self.__shift
    @shift.setter
    def shift ( self , value ) :
        self.set_value ( self.__shift , value )

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
class GenInvGauss_pdf(PDF1) :
    """ Generalized Inverse Gaussian distribution
    Generalised Inverse Gaussian distribution using (theta,eta) parameterisation
    - see https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution
    - see Ostap::Math::GenInvGauss
    """
    ## constructor
    def __init__ ( self         ,   
                   name         ,                               ## the name 
                   xvar         ,                               ## the variable 
                   theta        ,                               ## parameter theta
                   eta          ,                               ## parameter eta
                   p            ,                               ## parameter p
                   shift = ROOT.RooFit.RooConst ( 0.0 )  ) : ## shift parameter
        #
        PDF1.__init__ ( self , name , xvar )
        #

        self.__theta    = self.make_var ( theta             ,
                                          'theta_%s'          % name ,
                                          '#theta_{GIG}(%s)'  % name ,
                                          None , 1.0 , 1.e-8 , 100 )
        self.__eta      = self.make_var ( eta               ,
                                          'eta_%s'            % name ,
                                          '#eta_{GIG}(%s)'    % name ,
                                          None ,  1.0 , 1.e-8 , 100 )
        self.__p        = self.make_var ( p                 ,
                                          'p_%s'              % name ,
                                          'p_{GIG}(%s)'       % name ,
                                          None , 0   , -100  , 100 )        
        self.__shift   = self.make_var ( shift        ,
                                         'shift_%s'           % name ,
                                         'shift_{GIG}(%s)'    % name ,
                                         None  , 0 , -1.e+6 , 1.e+6 )
        
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
        """``theta''-parameter of Generilized Inverse Gaussian  function"""
        return self.__theta
    @theta.setter
    def theta ( self , value ) :
        self.set_value ( self.__theta , value )

    @property
    def eta ( self ) :
        """``eta''-parameter of Generilized Inverse Gaussian  function"""
        return self.__eta
    @eta.setter
    def eta ( self , value ) :
        self.set_value ( self.__eta , value )

    @property
    def p   ( self ) :
        """``p''-parameter of Generilized Inverse Gaussian  function"""
        return self.__p
    @p.setter
    def p ( self , value ) :
        self.set_value ( self.__p , value )

    @property
    def shift ( self ) :
        """``shift''-parameter of Rice function"""
        return self.__shift
    @shift.setter
    def shift ( self , value ) :
        self.set_value ( self.__shift , value )

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
class Weibull_pdf(PDF1) :
    r"""3-parameter  Weibull distribution 
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
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   scale = 1        ,   ## scale/lambda 
                   shape = 1        ,   ## shape/k 
                   shift = None     ) : ## shift/x0      
        #
        ## initialize the base
        #
        
        PDF1.__init__  ( self , name , xvar )
        
        self.__scale = self.make_var ( scale                   ,
                                       'scale_%s'              % name ,
                                       'scale_{Weibull}(%s)'   % name ,
                                       None , 1 , 1.e-8 , 1.e+6 )
        self.__shape = self.make_var ( shape                   ,
                                       'shape_%s'              % name ,
                                       'shape_{Weibull}(%s)'   % name ,
                                       None , 1 , 1.e-8 , 1.e+6 )
        
        limits_shift = () 
        if self.xminmax() :
            mn , mx = self.xminmax()
            dx = ( mx -  mn )
            limits_shift = mn + 0.01 * dx , mn , mx
            
        self.__shift = self.make_var ( shift                   ,
                                       'shift_%s'              % name ,
                                       'shift_{Weibull}(%s)'   % name ,
                                       None , *limits_shift )
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
    def scale ( self ) :
        """``scale''-parameter, scale>0"""
        return self.__scale
    @scale.setter
    def scale ( self, value ) :
        value = float ( value ) 
        assert 0 < value , "``scale''-parameter must be positive"
        self.__scale.setVal ( value )
        
    @property
    def shape ( self ) :
        """``shape''-parameter, shape>0
        - for ``shape''>2   : peak-like, signal-like shape
        - for 1<``shape'<2  : 'threshold'-like
        - for ``shape''k<1  : smooth decreasing         
        """
        return self.__shape    
    @shape.setter
    def shape ( self, value ) :
        value = float ( value ) 
        assert 0 < value , "``shape''-parameter must be positive"
        self.__shape.setVal ( value )         

    @property
    def shift ( self ) :
        """``shift''-parameter, function is   zero for x<``shift''"""
        return self.__shift
    @shift.setter
    def scale ( self, value ) :
        value = float ( value )
        if self.xminmax() :
            mn , mx = self.xminmax() 
            assert value <= mx, "``shift''-parameter must be smaller ``xmax''"
        self.__shift.setVal ( value )
        
models.append ( Weibull_pdf )      

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
    r"""Useful function to describe pT-spectra of particles 
    
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
    def __init__ ( self                   ,
                   xvar                   ,   ## pT-variable (for fitting) 
                   m0        = 0          ,   ## particle mass (may be fixed)
                   n         = None       ,   ## shape parameter
                   T         = None       ,   ## temperature parameter                   
                   name      = ''         ) :

        ## initialize the base 
        PDF1.__init__  ( self , name , xvar )
        
        self.__m0   = self.make_var ( m0              ,
                                      'm0_%s'  % name , 
                                      'm0(%s)' % name ,
                                      None , 0     , 1e+6 )
        
        self.__n    = self.make_var ( n               ,
                                      'n_%s'   % name , 
                                      'n(%s) ' % name ,
                                      None , 1 , 0.01  , 1000 )  
        
        self.__T    = self.make_var ( T               ,
                                      'T_%s'   % name , 
                                      'T(%s) ' % name ,
                                      None , 1 , 1.e-4 , 1e+6 )
        
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
        """``pt''-variable for Tsallis distribution (the same as ``x'')"""
        return self.xvar
    
    @property
    def m0 ( self ) :
        """``m0''-parameter of Tsallis' function"""
        return self.__m0
    @m0.setter
    def m0 ( self , value ) :
        value = float ( value )
        self.__m0.setVal ( value ) 

    @property
    def n ( self ) :
        """``n''-parameter of Tsallis' function"""
        return self.__n
    @n.setter
    def n ( self , value ) :
        value = float ( value )
        self.__n.setVal ( value ) 

    @property
    def T ( self ) :
        """``T''-parameter of Tsallis' function"""
        return self.__T
    @T.setter
    def T ( self , value ) :
        value = float ( value )
        self.__T.setVal ( value ) 
 
        
models.append ( Tsallis_pdf ) 
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
    r"""Useful function to describe pT-spectra of particles 
    
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
                   xvar             ,   ## pT-variable (for fitting) 
                   m0        = 0    ,   ## particle mass (may be fixed)
                   b         = None ,   ## slope parameter
                   name      = ''   ) :
        
        ## initialize the base 
        PDF1.__init__  ( self , name , pt )

        
        self.__m0   = self.make_var ( m0              ,
                                      'm0_%s'  % name , 
                                      'm0(%s)' % name ,
                                      None , 0  , 1e+6 )
        
        self.__b    = self.make_var ( b               ,
                                      'b_%s'   % name , 
                                      'b(%s) ' % name ,
                                      None , 0. , 1e+6 )  
        
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
        """``pt''-variable for QGSM distribution (the same as ``x'')"""
        return self.xvar
    
    @property
    def m0 ( self ) :
        """``m0''-parameter of QGSM function"""
        return self.__m0
    @m0.setter
    def m0 ( self , value ) :
        value = float ( value )
        self.__b0.setVal ( value ) 

    @property
    def b ( self ) :
        """``b''-parameter of QGSM function"""
        return self.__b
    @b.setter
    def b ( self , value ) :
        value = float ( value )
        self.__b.setVal ( value ) 
       

models.append ( QGSM_pdf )



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
                                      'm0_%s'  % name , 
                                      'm0(%s)' % name ,
                                      True , 0     , 1e+6 )
        
        self.__q    = self.make_var ( q               ,
                                      'q_%s'   % name , 
                                      'q(%s) ' % name ,
                                      None , 1.1 , 1.e-6, 10 )  
        
        self.__T    = self.make_var ( T               ,
                                      'T_%s'   % name , 
                                      'T(%s) ' % name ,
                                      None , 1 , 1.e-4 , 1e+6 )
        
        self.__mu   = self.make_var ( mu               ,
                                      'mu_%s'   % name , 
                                      '#mu(%s) '% name ,
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
 
        
models.append ( Tsallis2_pdf ) 

# =============================================================================
if '__main__' == __name__ : 
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )

# =============================================================================
##                                                                      The END 
# =============================================================================
