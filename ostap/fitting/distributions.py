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
        PDF.__init__ ( self , name )
        #
        self.x     = x
        self.mass  = x  ## ditto
        #
        self.k     = makeVar ( k       ,
                               'k_%s'                % name ,
                               'k_{#Gamma}(%s)'      % name , k     , 1 , 1.e-3 , 100 )
        self.theta = makeVar ( theta   ,
                               'theta_%s'            % name ,
                               '#theta_{#Gamma}(%s)' % name , theta , 1 , 1.e-3 , 100 )
        
        if self.k.getMin() <= 0 :
            self.k.setMin ( 1.e-3 ) 
            logger.warning( 'GammaDist(%s): min(k)     is set %s ' % ( name , self.k.getMin() ) )
            
        if self.theta.getMin() <= 0 :
            theta.setMin ( 1.e-3 ) 
            logger.warning( 'GammaDist(%s): min(theta) is set %s ' % ( name , self.theta.getMin() ) )
            
        self.pdf  = Ostap.Models.GammaDist (
            'gd_%s'         % name ,
            'GammaDist(%s)' % name ,
            self.x                 ,
            self.k                 ,
            self.theta             )

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
        PDF.__init__ ( self , name )
        #
        self.x     = x
        self.mass  = x  ## ditto
        #
        self.k     = makeVar ( k       ,
                               'k_%s'                % name ,
                               'k_{#Gamma}(%s)'      % name , k     , 1 , 1.e-3 , 100 )
        self.theta = makeVar ( theta   ,
                               'theta_%s'            % name ,
                               '#theta_{#Gamma}(%s)' % name , theta , 1 , 1.e-3 , 100 )
        self.p     = makeVar ( p       ,
                               'p_%s'                % name ,
                               'p_{#Gamma}(%s)'      % name , p     , 1 , 1.e-3 ,   6 )
        
        self.low   = makeVar ( low      ,
                               'low_%s'              % name ,
                               'l_{#Gamma}(%s)'      % name ,
                               low      ,
                               min ( 0 , x.getMin() ) , x.getMax() ) 
        
        if self.k    .getMin() <= 0 :
            self.k   .setMin ( 1.e-3 ) 
            logger.warning( 'GenGammaDist(%s): min(k)     is set %s ' % ( name , self.k.getMin() ) ) 
            
        if self.theta.getMin() <= 0 :
            self.theta.setMin ( 1.e-3 ) 
            logger.warning( 'GenGammaDist(%s): min(theta) is set %s ' % ( name , self.theta.getMin() ) ) 
            
        if self.p    .getMin() <= 0 :
            self.p   .setMin ( 1.e-3 ) 
            logger.warning( 'GenGammaDist(%s): min(p)     is set %s ' % ( name , self.p.getMin() ) ) 
            
        self.pdf  = Ostap.Models.GenGammaDist (
            'ggd_%s'           % name ,
            'GenGammaDist(%s)' % name ,
            self.x         ,
            self.k         ,
            self.theta     ,
            self.p         , 
            self.low       )

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
        PDF.__init__ ( self , name )
        #
        self.x     = x
        self.mass  = x  ## ditto
        #
        self.theta = makeVar ( theta   ,
                               'theta_%s'             % name ,
                               '#theta_{Amoroso}(%s)' % name , theta , 1 , 1.e-3 , 100 )
        self.alpha = makeVar ( alpha   ,
                               'alpha_%s'             % name ,
                               '#alpha_{Amoroso}(%s)' % name , alpha , 1 , 1.e-3 , 100 )
        self.beta  = makeVar ( beta    ,
                               'beta_%s'              % name ,
                               '#beta_{Amoroso}(%s) ' % name , beta  , 1 , 1.e-3 ,  10 )        
        self.a     = makeVar ( a       ,
                               'a_%s'                 % name ,
                               'a_{Amoroso}(%s)'      % name , a     , 1 , -10   ,  10  )
           
        logger.debug ('Amoroso theta  %s' % self.theta  )
        logger.debug ('Amoroso alpha  %s' % self.alpha  )
        logger.debug ('Amoroso beta   %s' % self.beta   )
        logger.debug ('Amoroso a      %s' % self.a      )

        self.pdf  = Ostap.Models.Amoroso (
            'amo_%s'      % name ,
            'Amoroso(%s)' % name ,
            self.x         ,
            self.theta     ,
            self.alpha     ,
            self.beta      ,
            self.a         )

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
        PDF.__init__ ( self , name )
        #
        self.x     = x
        self.mass  = x  ## ditto
        #
        self.k     = makeVar ( k       ,
                               'k_%s'                   % name ,
                               'k_{log#Gamma}(%s)'      % name , k     , 1 , 1.e-5 , 1000 )
        self.theta = makeVar ( theta   ,
                               'theta_%s'               % name ,
                               '#theta_{log#Gamma}(%s)' % name , theta , 1 , 1.e-5 , 1000 )

        if self.k.getMin() <= 0 :
            self.k.setMin ( 1.e-3 ) 
            logger.warning( 'LogGammaDist(%s): min(k)     is set %s ' % ( name  , self.k.getMin() ) ) 
            
        if self.theta.getMin() <= 0 :
            theta.setMin ( 1.e-3 ) 
            logger.warning( 'LogGammaDist(%s): min(theta) is set %s ' % ( name , self.theta.getMin() ) )
            
        logger.debug ('LogGammaDist k      %s' % self.k      )
        logger.debug ('LogGammaDist theta  %s' % self.theta  )

        self.pdf  = Ostap.Models.LogGammaDist (
            'lgd_%s'           % name ,
            'LogGammaDist(%s)' % name ,
            self.x                 ,
            self.k                 ,
            self.theta             )

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
        PDF.__init__ ( self , name )
        #
        self.x     = x
        self.mass  = x  ## ditto
        #
        self.k     = makeVar ( k       ,
                               'k_%s'                     % name ,
                               'k_{log10#Gamma}(%s)'      % name , k     , 1 , 1.e-4 , 10000 )
        self.theta = makeVar ( theta   ,
                               'theta_%s'                 % name ,
                               '#theta_{log10#Gamma}(%s)' % name , theta , 1 , 1.e-4 , 10000 )
        
        if self.k.getMin() <= 0 :
            self.k.setMin ( 1.e-4 ) 
            logger.warning( 'Log10GammaDist(%s): min(k)     is set %s ' % ( name , self.k.getMin() ) ) 
            
            if self.theta.getMin() <= 0 :
                theta.setMin ( 1.e-4 ) 
            logger.warning( 'Log10GammaDist(%s): min(theta) is set %s ' % ( name , self.theta.getMin() ) ) 
            
        logger.debug ('Log10GammaDist k      %s' % self.k      )
        logger.debug ('Log10GammaDist theta  %s' % self.theta  )

        self.pdf  = Ostap.Models.Log10GammaDist (
            'lgd10_%s'           % name ,
            'Log10GammaDist(%s)' % name ,
            self.x                 ,
            self.k                 ,
            self.theta             )

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
        PDF.__init__ ( self , name )
        #
        self.x     = x
        self.mass  = x  ## ditto
        #
        xmin = self.x.getMin()
        xmax = self.x.getMax()
        dx   = xmax - xmin
        # 
        self.nu     = makeVar ( nu       ,
                                'nu_%s'                    % name ,
                                '#nu_{#log#Gamma}(%s)'     % name , nu    ,
                                0.5 * ( xmin + xmax ) ,
                                xmin - 10 * dx ,
                                xmax + 10 * dx ) 
        
        self.lam    = makeVar ( lam      ,
                                'lambda_%s'                % name ,
                                '#lambda_{#log#Gamma}(%s)' % name , lam   , 2 , -1000 , 1000 )
        
        self.alpha  = makeVar ( alpha    ,
                                'alpha_%s'                 % name ,
                                '#alpha_{#log#Gamma}(%s)'  % name , alpha , 1 , 1.e-3 , 1000 )
        
        if self.alpha.getMin() <= 0 :
            self.alpha.setMin ( 1.e-3 ) 
            logger.warning( 'LogGamma(%s): min(alpha) is set %s ' % ( name , self.alpha.getMin() ) ) 
            
        logger.debug ('LogGamma nu     %s' % self.nu     )
        logger.debug ('LogGamma lambda %s' % self.lam    )
        logger.debug ('LogGamma alpha  %s' % self.alpha  )

        self.pdf  = Ostap.Models.LogGamma (
            'lg_%s'        % name ,
            'LogGamma(%s)' % name ,
            self.x     ,
            self.nu    ,
            self.lam   ,
            self.alpha )

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
        PDF.__init__ ( self , name )
        # 
        self.x     = x
        self.mass  = x  ## ditto
        #
        self.alpha  = makeVar ( alpha    ,
                                'alpha_%s'                 % name ,
                                '#alpha_{#beta#prime}(%s)' % name , alpha , 1 , 1.e-3 , 1000 )
        self.beta   = makeVar ( beta     ,
                                'beta_%s'                  % name ,
                                '#beta_{#beta#prime}(%s)'  % name , beta  , 1 , 1.e-3 , 1000 )

        
        if self.alpha.getMin() <= 0 :
            self.alpha.setMin ( 1.e-3 ) 
            logger.warning( 'BetaPrime(%s): min(alpha) is set %s ' % ( name , self.alpha.getMin() ) ) 
            
        if self.beta .getMin() <= 0 :
            self.beta.setMin ( 1.e-3 ) 
            logger.warning( 'BetaPrime(%s): min(beta) is set %s ' %  ( name , self.beta.getMin  () ) ) 
    
        self.scale  = makeVar ( scale     ,
                                'scale_%s'                 % name ,
                                '#theta_{#beta#prime}(%s)' % name , scale ,
                                1 , -1000 , 1000 )

        _dm = self.x.getMax()  - self.x.getMin() 
        self.delta  = makeVar ( delta     ,
                                'delta_%s'                 % name ,
                                '#delta_{#beta#prime}(%s)' % name , delta ,
                                0 ,
                                self.x.getMin() - 10 * _dm ,
                                self.x.getMax() + 10 * _dm ) 

        logger.debug ("Beta' alpha     %s" % self.alpha     )
        logger.debug ("Beta' beta      %s" % self.beta      )
        logger.debug ("Beta' scale     %s" % self.scale     )
        logger.debug ("Beta' sdelta    %s" % self.delta     )

        self.pdf  = Ostap.Models.BetaPrime (
            'bp_%s'         % name ,
            'BetaPrime(%s)' % name ,
            self.x     ,
            self.alpha ,
            self.beta  ,
            self.scale ,
            self.delta )

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
        PDF.__init__ ( self , name ) 
        #
        self.x     = x
        self.mass  = x  ## ditto
        #
        self.scale  = makeVar ( scale     ,
                                'scale_%s'            % name ,
                                '#theta_{Landau}(%s)' % name , scale ,
                                1 , -1000 , 1000 )
        
        _dm = self.x.getMax()  - self.x.getMin() 
        self.delta  = makeVar ( delta     ,
                                'delta_%s'            % name ,
                                '#delta_{Landau}(%s)' % name , delta ,
                                0 ,
                                self.x.getMin() - 10 * _dm ,
                                self.x.getMax() + 10 * _dm ) 
                                
        logger.debug ('Landau scale  %s' % self.scale  )
        logger.debug ('Landau delta  %s' % self.delta  )
        
        self.pdf  = Ostap.Models.Landau (
            'land_%s'    % name ,
            'Landau(%s)' % name ,
            self.x     ,
            self.scale ,
            self.delta )

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
        PDF.__init__ ( self , name ) 
        #
        self.x     = x
        self.mass  = x  ## ditto
        #
        self.shape  = makeVar ( shape     ,
                                'shape_%s'         % name ,
                                '#chi_{Argus}(%s)' % name , shape ,
                                1     ,
                                1.e-4 , 20 )
        
        self.high  = makeVar ( high      ,
                               'high_%s'          % name ,
                               'high_{Argus}(%s)' % name , high ,
                               0.1 * self.x.getMin ()  + 0.9 * self.x.getMax () , 
                               self.x.getMin () ,
                               self.x.getMax () )
        
        _dm  = self.x.getMax()  - self.x.getMin()
        lmin = min ( 0 , self.x.getMin() - 10 * _dm )
        lmax =           self.x.getMax() + 10 * _dm 
        self.low   = makeVar ( low      ,
                               'low_%s'          % name ,
                               'low_{Argus}(%s)' % name , low ,
                               0.9 * self.x.getMin ()  + 0.1 * self.x.getMax () , 
                               lmin , lmax ) 
        
        logger.debug ('ARGUS shape  %s' % self.shape  )
        logger.debug ('ARGUS high   %s' % self.high   )
        logger.debug ('ARGUS low    %s' % self.low    )

        self.pdf  = Ostap.Models.Argus (
            'arg_%s'    % name ,
            'Argus(%s)' % name ,
            self.x     ,
            self.shape ,
            self.high  ,
            self.low   )

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
        PDF.__init__ ( self , name ) 
        #
        self.x     = x
        self.mass  = x  ## ditto
        #
        self.alpha  = makeVar ( alpha      ,
                                'alpha_%s'        % name ,
                                '#alpha_{2e}(%s)' % name , alpha ,
                                1     ,
                                1.e-4 , 50 )
        self.delta  = makeVar ( delta     ,
                                'delta_%s'        % name ,
                                '#delta_{2e}(%s)' % name , delta ,
                                1     ,
                                1.e-4 , 50 )
        
        xmin = self.x.getMin()
        xmax = self.x.getMin()
        dx   = xmax - xmin 
        self.x0     = makeVar ( x0    ,
                                'x0_%s'       % name ,
                                'x0_{2e}(%s)' % name , x0 ,
                                xmin ,
                                xmin -0.1 * dx , xmax + 0.1 * dx )

        logger.debug ('2expos alpha  %s' % self.alpha  )
        logger.debug ('2expos delta  %s' % self.delta  )
        logger.debug ('2expos x0     %s' % self.x0     )

        self.pdf  = Ostap.Models.TwoExpos (
            'exp2_%s'    % name ,
            '2Expos(%s)' % name ,
            self.x     ,
            self.alpha ,
            self.delta ,
            self.x0    )

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
        PDF.__init__ ( self , name )
        #
        self.x     = makeVar ( x ,
                               'x_%s'           % name ,
                               'x_{Gumbel}(%s)' % name , x )                       
        self.mass  = self.x  ## ditto
        #
        self.mu    = makeVar ( mu      ,
                               'mu_%s'                  % name ,
                               'mu_{Gumbel}(%s)'        % name , mu   )
        self.beta  = makeVar ( beta        ,
                               'beta_%s'                % name ,
                               'beta_{Gumbel}(%s)'      % name , beta ) 
        
        self.pdf  = Ostap.Models.Gumbel (
            'gumbel_%s'  % name ,
            'Gumbel(%s)' % name ,
            self.x              ,
            self.mu             ,
            self.beta           )
        
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
        PDF.__init__  ( self , name )
        if not isinstance ( pt , ROOT.RooAbsReal ) :
            raise AttributeError( "Tsallis(%s): invalid 'pt'-parameter %s" % ( name , pt ) )
        
        self.pt   = pt
        
        self.m    = self.pt
        self.mass = self.pt
        
        self.m0   = makeVar ( mass            ,
                              'm0_%s'  % name , 
                              'm0(%s)' % name , mass , 0     , 1e+6 )
        
        self.n    = makeVar ( n               ,
                              'n_%s'   % name , 
                              'n(%s) ' % name , n    , 0.01  , 1000 )  
        
        self.n    = makeVar ( T               ,
                              'n_%s'   % name , 
                              'n(%s) ' % name , n    , 1.e-3 , 1e+6 )
        
        self.pdf  = Ostap.Models.Tsallis (
            'tsallis_'    + name  ,
            'Tsallis(%s)' % name  ,
            self.pt               ,
            self.n                ,
            self.T                ,
            self.m0               )
        
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
        PDF.__init__  ( self , name )
        if not isinstance ( pt , ROOT.RooAbsReal ) :
            raise AttributeError( "QGSM(%s): invalid 'pt'-parameter %s" % ( name , pt ) )
        
        self.pt   = pt
        
        self.m    = self.pt
        self.mass = self.pt
        
        self.m0   = makeVar ( mass            ,
                              'm0_%s'  % name , 
                              'm0(%s)' % name , mass , 0     , 1e+6 )
        
        self.b    = makeVar ( b               ,
                              'b_%s'   % name , 
                              'b(%s) ' % name , b    , 0.    , 1e+6 )  
        
        self.pdf  = Ostap.Models.QGSM (
            'qgsm_'    + name  ,
            'QGSM(%s)' % name  ,
            self.pt               ,
            self.b                ,
            self.m0               )
        

models.append ( QGSM_pdf )

# =============================================================================
if '__main__' == __name__ : 
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )

# =============================================================================
# The END 
# =============================================================================
