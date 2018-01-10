#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file disributions.py
#  A set of various smooth shapes and PDFs 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
# =============================================================================
"""A set of various smooth shapes and PDFs"""
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
    'LogGamma_pdf'       , ## 
    'BetaPrime_pdf'      , ## Beta-prime distribution 
    'Landau_pdf'         , ## Landau distribution 
    'Argus_pdf'          , ## ARGUS distribution 
    'TwoExpos_pdf'       , ## difference of two exponents
    'Gumbel_pdf'         , ## Gumbel distributions
    'Tsallis_pdf'        , ## Tsallis PDF 
    'QGSM_pdf'           , ## QGSM PDF 
    )
# =============================================================================
import ROOT, math
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.distributions' )
else                       : logger = getLogger ( __name__                      )
# =============================================================================
from   ostap.core.core     import cpp, Ostap, VE 
from   ostap.fitting.basic import makeVar, PDF
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
class GammaDist_pdf(PDF) :
    """Gamma-distribution with shape/scale parameters
    http://en.wikipedia.org/wiki/Gamma_distribution
    It suits nicely for firs of multiplicity and/or, especially chi2 distributions
    
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
                   x                ,   ## the variable
                   k     = None     ,   ## k-parameter
                   theta = None     ) : ## theta-parameter
        #
        PDF.__init__ ( self , name , x )
        #
        self.__k     = makeVar ( k       ,
                                 'k_%s'                % name ,
                                 'k_{#Gamma}(%s)'      % name , k     , 1 , 1.e-6 , 1000 )
        self.__theta = makeVar ( theta   ,
                                 'theta_%s'            % name ,
                                 '#theta_{#Gamma}(%s)' % name , theta , 1 , 1.e-6 , 1000 )
        
        self.pdf  = Ostap.Models.GammaDist (
            'gd_%s'         % name ,
            'GammaDist(%s)' % name ,
            self.x                 ,
            self.k                 ,
            self.theta             )

    @property
    def k ( self ) :
        """k-parameter of Gamma distribution   (k>0)"""
        return self.__k
    @k.setter 
    def k ( self , value ) :
        value = float ( value )
        assert 0 < value, 'k-value must be positive!'
        self.__k.setVal ( value ) 
        return self.__k.getVal() 

    @property
    def theta ( self ) :
        """Theta-parameter of Gamma distribution   (theta>0)"""
        return self.__theta
    @theta.setter 
    def theta ( self , value ) :
        value = float ( value )
        assert 0 < value, 'theta-value must be positive!'
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
class GenGammaDist_pdf(PDF) :
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
                   x                   ,   ## the variable
                   k     = None  ,   ## k-parameter
                   theta = None  ,   ## theta-parameter
                   p     = None  ,   ## p-parameter
                   low   = None  ) : ## low-parameter
        #
        PDF.__init__ ( self , name , x  )
        #
        self.__k     = makeVar ( k       ,
                                 'k_%s'                % name ,
                                 'k_{#Gamma}(%s)'      % name , k     , 1 , 1.e-5 , 1000 )
        self.__theta = makeVar ( theta   ,
                                 'theta_%s'            % name ,
                                 '#theta_{#Gamma}(%s)' % name , theta , 1 , 1.e-5 , 1000 )
        self.__p     = makeVar ( p       ,
                                 'p_%s'                % name ,
                                 'p_{#Gamma}(%s)'      % name , p     , 1 , 1.e-5 ,   10 )
        
        self.__low   = makeVar ( low      ,
                                 'low_%s'              % name ,
                                 'l_{#Gamma}(%s)'      % name ,
                                 low      ,
                                 min ( 0 , x.getMin() ) , x.getMax() ) 
        
        self.pdf  = Ostap.Models.GenGammaDist (
            'ggd_%s'           % name ,
            'GenGammaDist(%s)' % name ,
            self.x         ,
            self.k         ,
            self.theta     ,
            self.p         , 
            self.low       )

    @property
    def k ( self ) :
        """k-parameter of generalized Gamma distribution   (k>0)"""
        return self.__k
    @k.setter 
    def k ( self , value ) :
        value = float ( value )
        assert 0 < value, 'k-value must be positive!'
        self.__k.setVal ( value ) 
        return self.__k.getVal() 

    @property
    def theta ( self ) :
        """Theta-parameter of generalized Gamma distribution   (theta>0)"""
        return self.__theta
    @theta.setter 
    def theta ( self , value ) :
        value = float ( value )
        assert 0 < value, 'theta-value must be positive!'
        self.__theta.setVal ( value ) 
        return self.__theta.getVal() 

    @property
    def p ( self ) :
        """p-parameter of generalized Gamma distribution   (p>0)"""
        return self.__p
    @p.setter 
    def theta ( self , value ) :
        value = float ( value )
        assert 0 < value, 'p-value must be positive!'
        self.__p.setVal ( value ) 
        return self.__p.getVal() 

    @property
    def low ( self ) :
        """``x-low''-parameter of generalized Gamma distribution   (f(x)=0., for x < x_low)"""
        return self.__low
    @low.setter 
    def theta ( self , value ) :
        value = float ( value )
        self.__low.setVal ( value ) 
        return self.__low.getVal() 


models.append ( GenGammaDist_pdf ) 
# =============================================================================
## @class Amoroso_pdf
#  Another view on generalized gamma distribution
#  http://arxiv.org/pdf/1005.3274
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Math::Amoroso
#  @see Ostap::Models::Amoroso
class Amoroso_pdf(PDF) :
    """ Another view on generalized gamma distribution
    http://arxiv.org/pdf/1005.3274
    """
    ## constructor
    def __init__ ( self          ,
                   name          ,   ## the name 
                   x             ,   ## the variable
                   theta = None  ,   ## theta-parameter
                   alpha = None  ,   ## alpha-parameter
                   beta  = None  ,   ## beta-parameter
                   a     = None  ) : ## s-parameter
        
        #
        PDF.__init__ ( self , name , x )
        #
        self.__theta = makeVar ( theta   ,
                                 'theta_%s'             % name ,
                                 '#theta_{Amoroso}(%s)' % name , theta , 1 , 1.e-6 , 1000 )
        self.__alpha = makeVar ( alpha   ,
                                 'alpha_%s'             % name ,
                                 '#alpha_{Amoroso}(%s)' % name , alpha , 1 , 1.e-6  , 1000 )
        self.__beta  = makeVar ( beta    ,
                                 'beta_%s'              % name ,
                                 '#beta_{Amoroso}(%s) ' % name , beta  , 1 , 1.e-6  ,  100 )        
        self.__a     = makeVar ( a       ,
                                 'a_%s'                 % name ,
                                 'a_{Amoroso}(%s)'      % name , a     , 1 , -100   ,  100 )
        
        self.pdf  = Ostap.Models.Amoroso (
            'amo_%s'      % name ,
            'Amoroso(%s)' % name ,
            self.x         ,
            self.theta     ,
            self.alpha     ,
            self.beta      ,
            self.a         )

    @property
    def theta ( self ) :
        """Theta-parameter of Amoroso function (theta>0)"""
        return self.__theta
    @theta.setter 
    def theta ( self , value ) :
        value = float ( value )
        assert 0 < value, 'theta-value must be positive!'
        self.__theta.setVal ( value ) 
        return self.__theta.getVal() 

    @property
    def alpha ( self ) :
        """Alpha-parameter of Amoroso function (alpha>0)"""
        return self.__alpha
    @alpha.setter 
    def alpha ( self , value ) :
        value = float ( value )
        assert 0 < value, 'theta-value must be positive!'
        self.__alpha.setVal ( value ) 
        return self.__alpha.getVal() 

    @property
    def beta ( self ) :
        """Beta-parameter of Amoroso function (beta>0)"""
        return self.__beta
    @beta.setter 
    def beta ( self , value ) :
        value = float ( value )
        assert 0 < value, 'beta-value must be positive!'
        self.__beta.setVal ( value ) 
        return self.__beta.getVal() 

    @property
    def a ( self ) :
        """A-parameter of Amoroso function"""
        return self.__a
    @a.setter 
    def a ( self , value ) :
        value = float ( value )
        assert 0 < value, 'a-value must be positive!'
        self.__a.setVal ( value ) 
        return self.__a.getVal() 


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
class LogGammaDist_pdf(PDF) :
    """Distribution for log(x), where x follows Gamma distribution
    It suits nicely for fits of log(multiplicity) and/or log(chi2) distributions
    """
    ## constructor
    def __init__ ( self         ,
                   name         ,   ## the name 
                   x            ,   ## the variable
                   k     = None ,   ## k-parameter
                   theta = None ) : ## theta-parameter
        #
        PDF.__init__ ( self , name , x  )
        #
        self.k     = makeVar ( k       ,
                               'k_%s'                   % name ,
                               'k_{log#Gamma}(%s)'      % name , k     , 1 , 1.e-6 , 10000 )
        self.theta = makeVar ( theta   ,
                               'theta_%s'               % name ,
                               '#theta_{log#Gamma}(%s)' % name , theta , 1 , 1.e-6 , 10000 )

        self.pdf  = Ostap.Models.LogGammaDist (
            'lgd_%s'           % name ,
            'LogGammaDist(%s)' % name ,
            self.x                 ,
            self.k                 ,
            self.theta             )

    @property
    def k ( self ) :
        """k-parameter of log(Gamma) distribution   (k>0)"""
        return self.__k
    @k.setter 
    def k ( self , value ) :
        value = float ( value )
        assert 0 < value, 'k-value must be positive!'
        self.__k.setVal ( value ) 
        return self.__k.getVal() 

    @property
    def theta ( self ) :
        """Theta-parameter of log(Gamma) distribution   (theta>0)"""
        return self.__theta
    @theta.setter 
    def theta ( self , value ) :
        value = float ( value )
        assert 0 < value, 'theta-value must be positive!'
        self.__theta.setVal ( value ) 
        return self.__theta.getVal() 


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
class Log10GammaDist_pdf(PDF) :
    """Distribution for log10(x), where x follows Gamma distribution
    It suits nicely for fits of log10(multiplicity) and/or log10(chi2) distributions
    """
    ## constructor
    def __init__ ( self         ,
                   name         ,   ## the name 
                   x            ,   ## the variable
                   k     = None ,   ## k-parameter
                   theta = None ) :  ## theta-parameter
        #
        PDF.__init__ ( self , name , x  )
        #
        self.__k     = makeVar ( k       ,
                                 'k_%s'                     % name ,
                                 'k_{log10#Gamma}(%s)'      % name , k     , 1 , 1.e-6 , 100000 )
        self.__theta = makeVar ( theta   ,
                                 'theta_%s'                 % name ,
                                 '#theta_{log10#Gamma}(%s)' % name , theta , 1 , 1.e-6 , 100000 )
        
        self.pdf  = Ostap.Models.Log10GammaDist (
            'lgd10_%s'           % name ,
            'Log10GammaDist(%s)' % name ,
            self.x                 ,
            self.k                 ,
            self.theta             )
        
    @property
    def k ( self ) :
        """k-parameter of log10(Gamma) distribution   (k>0)"""
        return self.__k
    @k.setter 
    def k ( self , value ) :
        value = float ( value )
        assert 0 < value, 'k-value must be positive!'
        self.__k.setVal ( value ) 
        return self.__k.getVal() 

    @property
    def theta ( self ) :
        """Theta-parameter of log10(Gamma) distribution   (theta>0)"""
        return self.__theta
    @theta.setter 
    def theta ( self , value ) :
        value = float ( value )
        assert 0 < value, 'theta-value must be positive!'
        self.__theta.setVal ( value ) 
        return self.__theta.getVal() 



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
class LogGamma_pdf(PDF) :
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
                   x            ,   ## the variable
                   nu    = None ,   ## nu-parameter
                   lam   = None ,   ## lambda-parameter
                   alpha = None ) : ## nu-parameter
        #
        PDF.__init__ ( self , name , x )
        #
        xmin = self.x.getMin()
        xmax = self.x.getMax()
        dx   = xmax - xmin
        # 
        self.__nu     = makeVar ( nu       ,
                                  'nu_%s'                    % name ,
                                  '#nu_{#log#Gamma}(%s)'     % name , nu    ,
                                  0.5 * ( xmin + xmax ) ,
                                  xmin - 10 * dx ,
                                  xmax + 10 * dx ) 
        
        self.__lambd  = makeVar ( lam      ,
                                  'lambda_%s'                % name ,
                                  '#lambda_{#log#Gamma}(%s)' % name , lam   , 2 , -10000 , 10000 )
        
        self.__alpha  = makeVar ( alpha    ,
                                  'alpha_%s'                 % name ,
                                  '#alpha_{#log#Gamma}(%s)'  % name , alpha , 1 , 1.e-6 , 10000 )
        
        self.pdf  = Ostap.Models.LogGamma (
            'lg_%s'        % name ,
            'LogGamma(%s)' % name ,
            self.x     ,
            self.nu    ,
            self.lambd ,
            self.alpha )

    @property
    def nu ( self ) :
        """nu-parameter (location) of log-Gamma distribution"""
        return self.__nu
    @nu.setter 
    def nu ( self , value ) :
        value = float ( value )
        self.__nu.setVal ( value ) 
        return self.__nu.getVal() 

    @property
    def lambd ( self ) :
        """lambda-parameter (scale) of log-Gamma distribution   (theta>0)"""
        return self.__lambd
    @lambd.setter 
    def lambd ( self , value ) :
        value = float ( value )
        self.__lambd.setVal ( value ) 
        return self.__lambd.getVal() 

    @property
    def alpha  ( self ) :
        """alpha-parameter of log-Gamma distribution   (alpha>0)"""
        return self.__alpha
    @alpha.setter 
    def alpha ( self , value ) :
        value = float ( value )
        assert 0 < value, 'alpha-value must be positive!'
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
class BetaPrime_pdf(PDF) :
    """Beta-prime disribution 
    - http://en.wikipedia.org/wiki/Beta_prime_distribution
    """
    ## constructor
    def __init__ ( self         ,
                   name         ,   ## the name 
                   x            ,   ## the variable
                   alpha = None ,   ## alpha-parameter
                   beta  = None ,   ## beta-parameter
                   scale = 1    ,   ## scale-parameter 
                   delta = 0    ) : ## shift-parameter 
        #
        PDF.__init__ ( self , name , x )
        # 
        self.__alpha  = makeVar ( alpha    ,
                                  'alpha_%s'                 % name ,
                                  '#alpha_{#beta#prime}(%s)' % name , alpha , 1 , 1.e-6 , 10000 )
        self.__beta   = makeVar ( beta     ,
                                  'beta_%s'                  % name ,
                                  '#beta_{#beta#prime}(%s)'  % name , beta  , 1 , 1.e-6 , 10000 )
        
        self.__scale  = makeVar ( scale     ,
                                  'scale_%s'                 % name ,
                                  '#theta_{#beta#prime}(%s)' % name , scale ,
                                  1 , -1000 , 1000 )

        _dm = self.x.getMax()  - self.x.getMin() 
        self.__delta  = makeVar ( delta     ,
                                  'delta_%s'                 % name ,
                                  '#delta_{#beta#prime}(%s)' % name , delta ,
                                  0 ,
                                  self.x.getMin() - 10 * _dm ,
                                  self.x.getMax() + 10 * _dm ) 
        
        self.pdf  = Ostap.Models.BetaPrime (
            'bp_%s'         % name ,
            'BetaPrime(%s)' % name ,
            self.x     ,
            self.alpha ,
            self.beta  ,
            self.scale ,
            self.delta )
        
    @property
    def alpha ( self ) :
        """alpha-parameter of Beta' distribution   (alpha>0)"""
        return self.__alpha
    @alpha.setter 
    def alpha ( self , value ) :
        value = float ( value )
        assert 0 < value, 'alpha-value must be positive!'
        self.__alpha.setVal ( value ) 
        return self.__alpha.getVal() 

    @property
    def beta ( self ) :
        """beta-parameter of Beta' distribution   (alpha>0)"""
        return self.__beta
    @beta.setter 
    def beta ( self , value ) :
        value = float ( value )
        assert 0 < value, 'beta-value must be positive!'
        self.__beta.setVal ( value ) 
        return self.__beta.getVal() 

    @property
    def scale ( self ) :
        """scale-parameter of Beta' distribution"""
        return self.__scale
    @scale.setter 
    def scale ( self , value ) :
        value = float ( value )
        self.__scale.setVal ( value ) 
        return self.__scale.getVal() 

    @property
    def delta ( self ) :
        """delta-parameter of Beta' distribution"""
        return self.__delta
    @delta.setter 
    def delta ( self , value ) :
        value = float ( value )
        self.__delta.setVal ( value ) 
        return self.__delta.getVal() 



models.append ( BetaPrime_pdf ) 
# =============================================================================
## @class Landau_pdf
#  http://en.wikipedia.org/wiki/Landau_distribution
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::Landau
#  @see Ostap::Math::Landau
class Landau_pdf(PDF) :
    """Landau distribution 
    - http://en.wikipedia.org/wiki/Landau_distribution
    """
    ## constructor
    def __init__ ( self      ,
                   name      ,   ## the name 
                   x         ,   ## the variable
                   scale = 1 ,   ## scale-parameter 
                   delta = 0 ) : ## shift-parameter 
        #
        PDF.__init__ ( self , name , x  ) 
        #
        self.__scale  = makeVar ( scale     ,
                                  'scale_%s'            % name ,
                                  '#theta_{Landau}(%s)' % name , scale ,
                                  1 , -1000 , 1000 )
        
        _dm = self.x.getMax()  - self.x.getMin() 
        self.__delta  = makeVar ( delta     ,
                                  'delta_%s'            % name ,
                                  '#delta_{Landau}(%s)' % name , delta ,
                                  0 ,
                                  self.x.getMin() - 10 * _dm ,
                                  self.x.getMax() + 10 * _dm ) 
        
        self.pdf  = Ostap.Models.Landau (
            'landau_%s'    % name ,
            'Landau(%s)' % name ,
            self.x     ,
            self.scale ,
            self.delta )

    @property
    def scale ( self ) :
        """scale-parameter of Landau distribution"""
        return self.__scale
    @scale.setter 
    def scale ( self , value ) :
        value = float ( value )
        self.__scale.setVal ( value ) 
        return self.__scale.getVal() 

    @property
    def delta ( self ) :
        """delta-parameter of Landau distribution"""
        return self.__delta
    @delta.setter 
    def delta ( self , value ) :
        value = float ( value )
        self.__delta.setVal ( value ) 
        return self.__delta.getVal() 


models.append ( Landau_pdf ) 
# =============================================================================
## @class Argus_pdf
#  http://en.wikipedia.org/wiki/ARGUS_distribution
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2013-05-11
#  @see Ostap::Models::Argus
#  @see Ostap::Math::Argus
class Argus_pdf(PDF) :
    """Argus distribution
    - http://en.wikipedia.org/wiki/ARGUS_distribution
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   x                ,   ## the variable
                   shape = None     ,   ## shape-parameter 
                   high  = None     ,   ## high-parameter 
                   low   = 0        ) : ## low-parameter 
        #
        PDF.__init__ ( self , name , x  ) 
        #
        self.__shape  = makeVar ( shape     ,
                                  'shape_%s'         % name ,
                                  '#chi_{Argus}(%s)' % name , shape ,
                                  1     ,
                                  1.e-4 , 20 )
        
        self.__high  = makeVar ( high      ,
                                 'high_%s'          % name ,
                                 'high_{Argus}(%s)' % name , high ,
                                 0.1 * self.x.getMin ()  + 0.9 * self.x.getMax () , 
                                 self.x.getMin () ,
                                 self.x.getMax () )
        
        _dm  = self.x.getMax()  - self.x.getMin()
        lmin = min ( 0 , self.x.getMin() - 10 * _dm )
        lmax =           self.x.getMax() + 10 * _dm 
        self.__low   = makeVar ( low      ,
                                 'low_%s'          % name ,
                                 'low_{Argus}(%s)' % name , low ,
                                 0.9 * self.x.getMin ()  + 0.1 * self.x.getMax () , 
                                 lmin , lmax ) 
        
        self.pdf  = Ostap.Models.Argus (
            'arg_%s'    % name ,
            'Argus(%s)' % name ,
            self.x     ,
            self.shape ,
            self.high  ,
            self.low   )

    @property
    def shape ( self ) :
        """shape-parameter of Argus distribution"""
        return self.__shape
    @shape.setter 
    def shape ( self , value ) :
        value = float ( value )
        self.__shape.setVal ( value ) 
        return self.__shape.getVal() 

    @property
    def high ( self ) :
        """``high''-parameter of Argus distribution"""
        return self.__high
    @high.setter 
    def shape ( self , value ) :
        value = float ( value )
        self.__high.setVal ( value ) 
        return self.__high.getVal() 

    @property
    def low ( self ) :
        """``low''-parameter of Argus distribution"""
        return self.__low
    @low.setter 
    def low ( self , value ) :
        value = float ( value )
        self.__low.setVal ( value ) 
        return self.__low.getVal() 



models.append ( Argus_pdf ) 


# =============================================================================
## @class TwoExpos_pdf
#  simple difference of two exponents
#  \f[ f \propto \mathrm{e}^{-a_1 x} -\mathrm{e}^{-a_2 x} = 
#        \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f]
#  @author Vanya Belyaev Ivan.Belyaev@itep.ru
#  @see Ostap::Models::TwoExpos
#  @see Ostap::Math::TwoExpos
class TwoExpos_pdf(PDF) :
    r"""Simple difference of two exponents
    \f$ f  \propto 
    \mathrm{e}^{-a_1    x}       -\mathrm{e}^{-a_2 x} = 
    \mathrm{e}^{-\alpha x}\left(1-\mathrm{e}^{-\delta x}\right) \f$
    """
    ## constructor
    def __init__ ( self             ,
                   name             ,   ## the name 
                   x                ,   ## the variable
                   alpha = None     ,   ## shape-parameter 
                   delta = None     ,   ## high-parameter 
                   x0    = 0        ) : ## low-parameter 
        #
        PDF.__init__ ( self , name , x  ) 
        #
        self.__alpha  = makeVar ( alpha      ,
                                  'alpha_%s'        % name ,
                                  '#alpha_{2e}(%s)' % name , alpha ,
                                  1     ,
                                  1.e-6 , 100 )
        self.__delta  = makeVar ( delta     ,
                                  'delta_%s'        % name ,
                                  '#delta_{2e}(%s)' % name , delta ,
                                  1     ,
                                  1.e-6 , 100 )
        
        xmin = self.x.getMin()
        xmax = self.x.getMin()
        dx   = xmax - xmin 
        self.__x0     = makeVar ( x0    ,
                                  'x0_%s'       % name ,
                                  'x0_{2e}(%s)' % name , x0 ,
                                  xmin ,
                                  xmin -0.1 * dx , xmax + 0.1 * dx )
        
        self.pdf  = Ostap.Models.TwoExpos (
            'exp2_%s'    % name ,
            '2Expos(%s)' % name ,
            self.x     ,
            self.alpha ,
            self.delta ,
            self.x0    )

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
#  @see Ostap::Math::Gumbel
class Gumbel_pdf(PDF) :
    """Gumbel distribution
    - see https://en.wikipedia.org/wiki/Gumbel_distribution
    f(x,\mu,\beta) = 1/abs(beta) e^{-e^{-z}}, where z = (x-mu)/beta 
    - Very useful and important case: 
    if  g(x) ~ exp(- tau x ) and  z = log(x), than F(z) = f(z; -log(tau),1)
    It means that sum of exponential components will be represented by a set of
    peak-like Gumbel' structures with beta=1, mu=-log(tau) 
    """
    ## constructor
    def __init__ ( self      ,   
                   name      ,   ## the name 
                   x         ,   ## the variable 
                   mu   = 0  ,   ## shift parameter/mode
                   beta = 1  ) : ## scale parameter 
        #
        PDF.__init__ ( self , name , x )
        #
        self.__mu    = makeVar ( mu      ,
                                 'mu_%s'                  % name ,
                                 'mu_{Gumbel}(%s)'        % name , mu   )
        self.__beta  = makeVar ( beta        ,
                                 'beta_%s'                % name ,
                                 'beta_{Gumbel}(%s)'      % name , beta ) 
        
        self.pdf  = Ostap.Models.Gumbel (
            'gumbel_%s'  % name ,
            'Gumbel(%s)' % name ,
            self.x              ,
            self.mu             ,
            self.beta           )

    @property
    def mu ( self ) :
        """mu-parameter (shift) of the Gumbel function"""
        return self.__mu
    @mu.setter
    def mu ( self , value ) :
        value = float ( value )
        self.__mu.setVal ( value ) 
        return self.__mu.getVal()

    @property
    def beta ( self ) :
        """beta-parameter (scale) of the Gumbel function"""
        return self.__beta
    @beta.setter
    def mu ( self , value ) :
        value = float ( value )
        self.__beta.setVal ( value ) 
        return self.__beta.getVal()
    
        
models.append ( Gumbel_pdf ) 

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
class Tsallis_pdf(PDF) :
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
                   pt                     ,   ## pT-variable (for fitting) 
                   mass      = 0          ,   ## particle mass (may be fixed)
                   n         = None       ,   ## shape parameter
                   T         = None       ,   ## temperature parameter                   
                   name      = ''         ) :

        ## initialize the base 
        PDF.__init__  ( self , name , pt )
        
        self.__m0   = makeVar ( mass            ,
                                'm0_%s'  % name , 
                                'm0(%s)' % name , mass , 0     , 1e+6 )
        
        self.__n    = makeVar ( n               ,
                                'n_%s'   % name , 
                                'n(%s) ' % name , n    , 0.01  , 1000 )  
        
        self.__T    = makeVar ( T               ,
                                'n_%s'   % name , 
                                'n(%s) ' % name , n    , 1.e-3 , 1e+6 )
        
        self.pdf  = Ostap.Models.Tsallis (
            'tsallis_'    + name  ,
            'Tsallis(%s)' % name  ,
            self.xvar             ,
            self.n                ,
            self.T                ,
            self.m0               )

    @property
    def pt ( self ) :
        """``pt''-variable for Tsallis distribution (the same as ``x'')"""
        return self.xvar
    
    @property
    def m0 ( self ) :
        """m0-parameter of Tsallis' function"""
        return self.__m0
    @m0.setter
    def m0 ( self , value ) :
        value = float ( value )
        self.__b0.setVal ( value ) 
        return self.__b0.getVal()

    @property
    def n ( self ) :
        """n-parameter of Tsallis' function"""
        return self.__n
    @n.setter
    def n ( self , value ) :
        value = float ( value )
        self.__n.setVal ( value ) 
        return self.__n.getVal()

    @property
    def T ( self ) :
        """T-parameter of Tsallis' function"""
        return self.__T
    @T.setter
    def T ( self , value ) :
        value = float ( value )
        self.__T.setVal ( value ) 
        return self.__T.getVal()
    
        
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
#  \f[ \frac{d\sigma}{dp_T} \propto p_T \times \mathrm{e}^{ -b_0 (m_T-m)} \f], 
#  where transverse mass is defined as \f$m_T = \sqrt{p_T^2+m^2}\f$
#
#  @see Ostap::Models::QGSM
#  @see Ostap::Math::QGSM
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class QGSM_pdf(PDF) :
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
                   pt               ,   ## pT-variable (for fitting) 
                   mass      = 0    ,   ## particle mass (may be fixed)
                   b         = None ,   ## slope parameter
                   name      = ''   ) :
        
        ## initialize the base 
        PDF.__init__  ( self , name , pt )

        self.__m0   = makeVar ( mass            ,
                                'm0_%s'  % name , 
                                'm0(%s)' % name , mass , 0     , 1e+6 )
        
        self.__b    = makeVar ( b               ,
                                'b_%s'   % name , 
                                'b(%s) ' % name , b    , 0.    , 1e+6 )  
        
        self.pdf  = Ostap.Models.QGSM (
            'qgsm_'    + name  ,
            'QGSM(%s)' % name  ,
            self.pt               ,
            self.b                ,
            self.m0               )
        
    @property
    def pt ( self ) :
        """``pt''-variable for QGSM distribution (the same as ``x'')"""
        return self.xvar
    
    @property
    def m0 ( self ) :
        """m0-parameter of QGSM function"""
        return self.__m0
    @m0.setter
    def m0 ( self , value ) :
        value = float ( value )
        self.__b0.setVal ( value ) 
        return self.__b0.getVal()

    @property
    def b ( self ) :
        """n-parameter of QGSM function"""
        return self.__b
    @b.setter
    def n ( self , value ) :
        value = float ( value )
        self.__b.setVal ( value ) 
        return self.__b.getVal()


models.append ( QGSM_pdf )

# =============================================================================
if '__main__' == __name__ : 
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )

# =============================================================================
# The END 
# =============================================================================
