#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# $Id$
# =============================================================================
## @file signals.py
#
#  Set of useful PDFs for various ``signal'' 1D and 2D fits
#  It includes
#  - soeme empricial PDFs to describe narrow peaks: Gauss, CrystalBall, ....
#  - some PDF to describe "wide" peaks: BreitWigner,LASS, Bugg, Flatter, ...
#  - some useful PDFs to describe smooth background: phase space ;
#    expo times polynomial; phase space times polynomial, ...
#  - set of smooth non-facrorizeable model for 2D fits 
#
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# 
#                    $Revision$
#  Last modification $Date$
#                 by $Author$
# =============================================================================
"""Set of useful PDFs for various ``signal'' 1D and 2D fits

It includes

Empricial PDFs to describe narrow peaks

  - Gauss 
  - Crystal Ball
  - right-side Crystal Ball
  - double-side Crystal Ball
  - Needham function for J/psi, psi' and Y peaks
  - Apolonios
  - Apolonios2 (bifurcated Apolonious)
  - bifurcated Gauissian
  - generalized normal v1 
  - generalized normal v2
  ## - skew Gaussian  (temporarily removed)
  - Bukin,
  - Student-T
  - bifurcated Student-T
  - SinhAsinh_pdf   
  - JohnsonSU_pdf   
  - Atlas_pdf   
  - Sech_pdf   
  - Logistic_pdf   
  
PDF to describe ``wide'' peaks

  - BreitWigner
  - LASS
  - Bugg
  - Flatte
  - Swanson's S=wave cusp 
  - ...

Special stuff:

  - Voigt & PseudoVoigt
  - BW23L

"""
# =============================================================================
__version__ = "$Revision:"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-07-25"
# =============================================================================
__all__ = (
    #
    ## empirical 1D signal models
    # 
    'Gauss_pdf'              , ## simple     Gauss
    'CrystalBall_pdf'        , ## Crystal-ball function
    'CrystalBallRS_pdf'      , ## right-side Crystal-ball function
    'CB2_pdf'                , ## double-sided Crystal Ball function    
    'Needham_pdf'            , ## Needham function for J/psi or Y fits 
    'Apolonios_pdf'          , ## Apolonios function         
    'Apolonios2_pdf'         , ## Apolonios function         
    'BifurcatedGauss_pdf'    , ## bifurcated Gauss
    'GenGaussV1_pdf'         , ## generalized normal v1  
    'GenGaussV2_pdf'         , ## generalized normal v2 
    ## 'SkewGauss_pdf'          , ## skewed gaussian (temporarily removed)
    'Bukin_pdf'              , ## generic Bukin PDF: skewed gaussian with exponential tails
    'StudentT_pdf'           , ## Student-T function 
    'BifurcatedStudentT_pdf' , ## bifurcated Student-T function
    'SinhAsinh_pdf'          , ## "Sinh-arcsinh distributions". Biometrika 96 (4): 761
    'JohnsonSU_pdf'          , ## JonhsonSU-distribution 
    'Atlas_pdf'              , ## modified gaussian with exponenital tails 
    'Sech_pdf'               , ## hyperboilic secant  (inverse-cosh) 
    'Logistic_pdf'           , ## Logistic aka "sech-squared"   
    #
    ## pdfs for "wide" peaks, to be used with care - phase space corrections are large!
    # 
    'BreitWigner_pdf'      , ## (relativistic) 2-body Breit-Wigner
    'Flatte_pdf'           , ## Flatte-function  (pipi)
    'Flatte2_pdf'          , ## Flatte-function  (KK) 
    'LASS_pdf'             , ## kappa-pole
    'Bugg_pdf'             , ## sigma-pole
    'Swanson_pdf'          , ## Swanson's S-wave cusp 
    ##
    'Voigt_pdf'            , ## Voigt-profile
    'PseudoVoigt_pdf'      , ## PseudoVoigt-profile
    'BW23L_pdf'            , ## BW23L
    #
    )
# =============================================================================
import ROOT, math
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.models_signal' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
from   ostap.core.core     import cpp , Ostap 
from   ostap.fitting.basic import makeVar, MASS 
# =============================================================================
models = [] 
# =============================================================================
## @class Gauss_pdf
#  simple wrapper over Gaussian-pdf
#  @see http://en.wikipedia.org/wiki/Normal_distribution
#  @see RooGaussian
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Gauss_pdf(MASS) :
    """Trivial Gaussian function:
    http://en.wikipedia.org/wiki/Normal_distribution
    """
    def __init__ ( self             ,
                   name             ,
                   mn        = None ,
                   mx        = None ,
                   mass      = None ,
                   mean      = None ,
                   sigma     = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name   ,
                         mn      , mx     , mass ,
                         mean    , sigma  )
        
        #
        ## build pdf
        # 
        self.pdf = ROOT.RooGaussian (
            'gauss_%s'  % name ,
            "Gauss(%s)" % name ,
            self.mass  ,
            self.mean  ,
            self.sigma )

models.append ( Gauss_pdf ) 
# =============================================================================
## @class CrystalBall_pdf
#  @see Ostap::Models::CrystalBall
#  @see Ostap::Math::CrystalBall
#  @see http://en.wikipedia.org/wiki/Crystal_Ball_function
#  - T. Skwarnicki, A study of the radiative CASCADE transitions between
#                   the Upsilon-Prime and Upsilon resonances,
#                   Ph.D Thesis, DESY F31-86-02(1986), Appendix E.
#                   @see http://inspirehep.net/record/230779/files/f31-86-02.pdf
#  - J. E. Gaiser,  Charmonium Spectroscopy from Radiative Decays of the J/Psi and Psi-Prime,
#                   Ph.D. Thesis, SLAC-R-255 (1982), Appendix F, p 178
#                   @see http://www.slac.stanford.edu/cgi-wrap/getdoc/slac-r-255.pdf
#  - M. J. Oreglia, A Study of the Reactions psi prime --> gamma gamma psi,
#                   Ph.D. Thesis, SLAC-R-236 (1980), Appendix D.
#                   @see http://www.slac.stanford.edu/pubs/slacreports/slac-r-236.html
#  @attention Unlike the original definition parameter <code>n</code> here is shifted by 1:
#                    \f$ n_{0} \leftarrow \left| n_{here} \right| + 1  \f$
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class CrystalBall_pdf(MASS) :
    """ Crystal Ball function
    http://en.wikipedia.org/wiki/Crystal_Ball_function
    
    - T. Skwarnicki,
    ``A study of the radiative cascade transitions between the Upsilon-Prime
    and Upsilon resonances'', Ph.D Thesis, DESY F31-86-02(1986), Appendix E.
    http://inspirehep.net/record/230779/files/f31-86-02.pdf
    - J. E. Gaiser,
    ``Charmonium Spectroscopy from Radiative Decays of the J/Psi and Psi-Prime'',
    Ph.D. Thesis, SLAC-R-255 (1982), Appendix F, p 178
    http://www.slac.stanford.edu/cgi-wrap/getdoc/slac-r-255.pdf
    - M. J. Oreglia,
    ``A Study of the Reactions psi prime --> gamma gamma psi'',
    Ph.D. Thesis, SLAC-R-236 (1980), Appendix D.
    http://www.slac.stanford.edu/pubs/slacreports/slac-r-236.html
    
    Note:
    - Unlike the original definition parameter 'n' here is shifted by 1: n <- |n| + 1
    - Typical value of parameter alpha for ``physical'' peaks is 1.5<alpha<3.0,
    - For large alpha (e.g. alpha>3), there is no sensitivity for n;
    similarly in the limit of large n, sensitivity for alpha is minimal
    """
    def __init__ ( self             ,
                   name             ,
                   mn       = None  ,
                   mx       = None  , 
                   mass     = None  ,
                   mean     = None  , 
                   sigma    = None  ,
                   alpha    = None  ,
                   n        = None  ) : 
                   
        #
        ## initialize the base
        #
        MASS.__init__ ( self    , name     ,
                        mn      , mx       , mass    ,
                        mean    , sigma    )
        
        self.alpha = makeVar ( alpha ,
                               'alpha_%s'        % name ,
                               '#alpha_{CB}(%s)' % name ,  alpha  ,
                               2.0 , 0  , 10 )
        
        self.n     = makeVar ( n   ,
                               'n_%s'            % name ,
                               'n_{CB}(%s)'      % name , n       ,
                               1.0 , 0  , 20 )
        
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.CrystalBall (
            'cb_%s'           % name ,
            'CrystalBall(%s)' % name ,
            self.mass  ,
            self.mean  ,
            self.sigma ,
            self.alpha ,
            self.n     )

models.append ( CrystalBall_pdf )    
# =============================================================================
## @class CrystalBallRS_pdf
#  Crystal Ball function with the right side tail.
#  to be rewritten
#  @see Ostap::Models::CrystalBallRS
#  @see Ostap::Math::CrystalBallRS
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class CrystalBallRS_pdf(MASS) :
    """Right-side CrystalBall
    """
    def __init__ ( self              ,
                   name              ,
                   mn       = None   ,
                   mx       = None   , 
                   mass     = None   ,
                   mean     = None   , 
                   sigma    = None   ,
                   alpha    = None   ,
                   n        = None   ) : 
                   
        
        MASS.__init__ ( self    , name     ,
                        mn      , mx       , mass    ,
                        mean    , sigma    ) 
        
        self.alpha = makeVar ( alpha ,
                               'alpha_%s'          % name ,
                               '#alpha_{CBRS}(%s)' % name , alpha , 
                               2.0 , 0  , 10      )
        
        self.n     = makeVar ( n     ,
                               'n_%s'              % name ,
                               'n_{CBRS}(%s)'      % name ,  n , 
                               1   ,  0 , 20      )
        
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.CrystalBallRS (
            'cbrs_%s'           % name ,
            'CrystalBallRS(%s)' % name ,
            self.mass  ,
            self.mean  ,
            self.sigma ,
            self.alpha ,
            self.n     )
        
models.append ( CrystalBallRS_pdf )    
# =============================================================================
## @class CB2_pdf
#  Double-sided Cristal Ball function
#  It appears to be very powerful and is used for many LHCb papers to describe
#   B-hadron mass signals, especially for \f$B \rightarrow J/\psi X\f$-final states 
#  @see Ostap::Models::CrystalBallDS
#  @see CrystalBall_pdf 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class CB2_pdf(MASS) :
    """Double sided Crystal Ball function with both left and rigth sides
    It appears to be very powerful and is used for many LHCb papers to describe
    B-hadron mass signals, especially for B->J/psi X final states 
    
    Note:
    - Similar to CrystalBall_pdf and unlike the original definition,
    the parameters 'n' here are shifted by 1: n <- |n| + 1
    - Typical value of parameters alpha for ``physical'' peaks is 1.5<alpha<3.0,
    - For large alpha (e.g. alpha>3), there is no sensitivity for n;
    similarly in the limit of large n, sensitivity for alpha is minimal
    """
    def __init__ ( self              ,
                   name              ,
                   mn        = None  ,
                   mx        = None  ,
                   mass      = None  , 
                   mean      = None  ,
                   sigma     = None  ,
                   alphaL    = None  ,
                   alphaR    = None  ,
                   nL        = None  ,
                   nR        = None  ) : 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name     ,
                         mn      , mx       , mass    ,
                         mean    , sigma    )
        #
        ## treat the specific parameters
        #
        self.aL    = makeVar ( alphaL                  ,
                               "aL_%s"          % name ,
                               "#alpha_{L}(%s)" % name , alphaL    , 2.0 , 0 , 10 )
        self.nL    = makeVar ( nL                      ,                     
                               "nL_%s"          % name ,
                               "n_{L}(%s)"      % name , nL        , 1   , 0 , 20 )
        self.aR    = makeVar ( alphaR ,
                               "aR_%s"          % name ,
                               "#alpha_{R}(%s)" % name , alphaR    , 2.0 , 0 , 10 )
        self.nR    = makeVar ( nR                      ,
                               "nR_%s"          % name ,
                               "n_{R}(%s)"      % name , nR        , 1   , 0 , 20 )
        
        self.pdf = Ostap.Models.CrystalBallDS(
            "cb2_"       + name ,
            "CB_{2}(%s)" % name ,
            self.mass    ,
            self.mean    ,
            self.sigma   ,
            self.aL      ,
            self.nL      ,
            self.aR      ,
            self.nR      )
        
models.append ( CB2_pdf )    
# =============================================================================
## @class Needham_pdf
#  Needham function: specific parameterisation of Crystal Ball function with
#   - \f$ n = 1 \f$  
#   - \f$ \alpha(\sigma) = a_0 + \sigma\times (a_1+\sigma \times a_2)
#  The function is very well sutable to fit
#  \f$J/\psi \rightarrow \mu^+\mu^-\f$,
#  \f$\psi^{\prime} \rightarrow \mu^+\mu^-\f$ and
#  \f$\Upsilon \rightarrow \mu^+\mu^-\f$ signals and
#  is has been used with great success for all LCHb papers
#  on quarkonia production in dimuon final states
#  @thank Matthew Needham
#  @see CrystalBall_pdf 
#  @see Ostap::Models::Needham 
#  @see Ostap::Math::Needham 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Needham_pdf(MASS) :
    """Needham function: specific parameterisation of Crystal Ball function with
    - n = 1 
    - alpha(sigma) = a_0 + sigma*(a_1+sigma*a_2)
    
    The function is well sutable to fit dimuon final states: 
    - J/psi   -> mu+ mu-
    - psi'    -> mu+ mu-
    - Upsilon -> mu+ mu-
    Is has been used with great success for all LCHb papers on
    quarkonia production in dimuon final states
    """
    def __init__ ( self                 ,
                   name                 ,
                   mn       = None      ,
                   mx       = None      , 
                   mass     = None      ,
                   mean     = 3.096     ,   ## GeV  
                   sigma    = 0.013     ,   ## GeV 
                   a0       = 1.975     ,
                   a1       = -0.0011   ,   ## GeV^-1
                   a2       = -0.00018  ) : ## GeV^-2  
        
        MASS.__init__ ( self    ,
                        name    ,
                        mn      , mx    ,
                        mass    ,
                        mean    , sigma )
        
        #
        unit = 1000
        #
        if   self.mass.getMin() <= 3.096 <= self.mass.getMax() : unit = 1000 
        elif self.mass.getMin() <=  3096 <= self.mass.getMax() : unit = 1
        elif self.mass.getMin() <= 9.460 <= self.mass.getMax() : unit = 1000 
        elif self.mass.getMin() <=  9460 <= self.mass.getMax() : unit = 1
        #
        self.a0 = makeVar ( a0                  ,
                            "a0_%s"     % name  ,
                            "a_{0}(%s)" % name  , a0 , 
                            1.975               ,   0           , 10           )
        self.a1 = makeVar ( a1                  ,
                            "a1_%s"     % name  ,
                            "a_{1}(%s)" % name  , a1 , 
                            -0.0011   * unit    , -10 * unit    , 10 * unit    )
        self.a2 = makeVar ( a2                  ,
                            "a2_%s"     % name  ,
                            "a_{2}(%s)" % name  , a2 , 
                            -0.00018  * unit**2 , -10 * unit**2 , 10 * unit**2 )
        #
        self.pdf = Ostap.Models.Needham (
            'needham_%s'  % name ,
            'needham(%s)' % name ,
            self.mass  ,
            self.mean  ,
            self.sigma ,
            self.a0    ,
            self.a1    ,
            self.a2
            )
        
models.append ( Needham_pdf )    
# =============================================================================
## @class Apolonios_pdf
#  simple wrapper over Apolonios PDF 
#  @see Ostap::Models::Apolonios 
#  The function is proposed by Diego Martinez Santos 
#  @see http://arxiv.org/abs/1312.5000
#  Here a bit modified version is used with redefined parameter <code>n</code>
#  to be coherent with local definitions of Crystal Ball
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Apolonios_pdf(MASS) :
    """Apolonios function
    http://arxiv.org/abs/1312.5000
    
    The function is proposed by Diego Martinez Santos 
    https://indico.itep.ru/getFile.py/access?contribId=2&resId=1&materialId=slides&confId=262633
    
    Note :
    - similar to CrystalBall case,  parameter n is redefined:   n <- |n|+1
    to be coherent with local definitions of Crystal Ball
    - unfortuately neither sigma nor b parameters allows easy interpretation
    - typical value of parameters alpha for ``physical'' peaks is 1.5<alpha<3.0,
    - for large alpha (e.g. alpha>3), there is no sensitivity for n;
    similarly in the limit of large n, sensitivity for alpha is minimal
    """
    def __init__ ( self                    ,
                   name                    ,
                   mn        = None        ,
                   mx        = None        ,
                   mass      = None        ,
                   mean      = None        ,
                   sigma     = None        ,
                   alpha     = None        ,
                   n         = None        ,
                   b         = None        ) : 
                   
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name  ,
                         mn      , mx    , mass    ,
                         mean    , sigma ) 
        
        self.alpha = makeVar ( alpha                     ,
                               'alpha_%s'         % name ,
                               '#alpha_{Apo}(%s)' % name , alpha , 
                               2.0   , 0 , 10 )
        
        self.n     = makeVar ( n                    ,
                               'n_%s'        % name ,
                               'n_{Apo}(%s)' % name , n ,
                               2.0   , 0 , 20 )
        
        self.b     = makeVar ( b                    ,
                               'b_%s'        % name ,
                               'b_{Apo}(%s)' % name ,  b  ,
                               1         , 0.01 , 10000 ) 
        
        #
        ## finally build PDF
        #
        self.pdf  = Ostap.Models.Apolonios (
            "apolo_"        + name ,
            "Apolonios(%s)" % name ,
            self.mass   ,
            self.mean   ,
            self.sigma  ,
            self.alpha  ,
            self.n      ,
            self.b      ) 


models.append ( Apolonios_pdf )    
# =============================================================================
## @class Apolonios2_pdf
#  "Bifurcated Apolonious"
#  Gaussian with exponential (asymmetrical) tails
#
#  A convinient reparameterization is applied to keep reduce 
#  the correlations between "sigma"s and "beta" 
#  \f[ f(x;\mu,\sigma_l,\sigma_r,\beta) \propto 
#         \mathrm{e}^{\left|\beta\right|( \left|\beta\right| -
#         \sqrt{ \beta^2+\left(\delta x\right)^2}} \f] 
#   where 
#  \f[ \delta x  = \left\{ \begin{array}{ccc}
#          \frac{x-\mu}{\sigma_l} & \text{for} & x \le \mu \\
#          \frac{x-\mu}{\sigma_r} & \text{for} & x \gt \mu \\
#          \end{array}
#          \right.\f]
#  Large betas corresponds to gaussian 
#      
#  @see Ostap::Models::Apolonios2 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-20
class Apolonios2_pdf(MASS) :
    """Bifurcated Apolonious:
    Gaussian with exponential (asymmetrical) tails
    
    f(x; mu, sigma_l, sigma_r, beta) ~ exp( |beta|(|\beta| - sqrt( beta^2+( delta x)^2 ))      
    with
    delta x  = (x-mu)/sigma_l for  x < mu
    delta x  = (x-mu)/sigma_r for  x > mu
    and
    sigma_{l,r} = sigma * ( 1 + kappa)
    with asymmetry parameter -1 < kappa < 1 
    
    The function is inspired by Appolonios function, but it has no power-law tails:
    instead the exponential tails are used.
    Reparameterization reduces the correlation between ``sigma'' and ``beta'',
    allowing their easy interpretation.
    Large ``beta'' and small ``asymmetry'' corresponds to gaussian
    
    """
    def __init__ ( self               ,
                   name               ,
                   mn        = None   ,
                   mx        = None   ,
                   mass      = None   ,
                   mean      = None   ,
                   sigma     = None   ,
                   asymmetry = None   ,
                   beta      = None   ) : 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name   ,
                         mn      , mx     , mass    ,
                         mean    , sigma  )
        
        self.asym = makeVar ( asymmetry                  ,
                              'asym_%s'           % name ,
                              '#kappa_{Apo2}(%s)' % name ,
                              asymmetry , -1 , 1  ) 
        
        self._lst_R = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.sigmaR = ROOT.RooFormulaVar (
            "sigmaR_%s"     % name   ,
            "sigma_{R}(%s)" % name   ,
            "%s*(1-%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self._lst_R   )
        
        self._lst_L = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.sigmaL = ROOT.RooFormulaVar (
            "sigmaL_%s"     % name   ,
            "sigma_{L}(%s)" % name   ,
            "%s*(1+%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self._lst_L   )
        
        self.beta    = makeVar ( beta ,
                                 'beta_%s'          % name  ,
                                 '#beta_{Apo2}(%s)' % name  ,
                                 beta , 0.01  , 10000 ) 
        #
        ## finally build PDF
        #
        self.pdf  = Ostap.Models.Apolonios2 (
            "apolo2_"        + name ,
            "Apolonios2(%s)" % name ,
            self.mass   ,
            self.mean   ,
            self.sigmaL ,
            self.sigmaR ,
            self.beta   ) 


models.append ( Apolonios2_pdf )    
# =============================================================================
## @class BifurcatedGauss_pdf
#  simple wrapper over bifurcated-gaussian
#  @see RooGaussian
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class BifurcatedGauss_pdf(MASS) :
    """Bifurcated Gauss :

    f(x; mu, sigma_l, sigma_r ) ~ exp ( -0.5 * dx^2 )
    
    with
    delta x  = (x-mu)/sigma_l for  x < mu
    delta x  = (x-mu)/sigma_r for  x > mu
    and
    sigma_{l,r} = sigma * ( 1 + kappa)
    with asymmetry parameter -1 < kappa < 1 
    
    """
    def __init__ ( self                  ,
                   name                  ,
                   mn        = None      ,
                   mx        = None      ,
                   mass      = None      ,
                   mean      = None      ,
                   sigma     = None      ,
                   asymmetry = None      ) : 
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name  ,
                         mn   , mx    , mass ,
                         mean , sigma ) 
        
        self.asym = makeVar ( asymmetry               ,
                              'asym_%s'        % name ,
                              '#xi_{asym}(%s)' % name ,
                              asymmetry , -1 , 1  ) 
        
        self._lst_R = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.sigmaR = ROOT.RooFormulaVar (
            "sigmaR_%s"     % name   ,
            "sigma_{R}(%s)" % name   ,
            "%s*(1+%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self._lst_R   )
        
        self._lst_L = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.sigmaL = ROOT.RooFormulaVar (
            "sigmaL_%s"     % name   ,
            "sigma_{L}(%s)" % name   ,
            "%s*(1-%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self._lst_L   )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.BifurcatedGauss (
            "fbgau_"         + name ,
            "BufurGauss(%s)" % name ,
            self.mass   ,
            self.mean   ,
            self.sigmaL ,
            self.sigmaR )

models.append ( BifurcatedGauss_pdf )    
# =============================================================================
## @class GenGaussV1_pdf
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
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class GenGaussV1_pdf(MASS) :
    """Generalized Normal distribution v1
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
    
    Parameters:
    - mu         : location/mean  
    - alpha > 0  : scale 
    - beta  > 0  : shape   (beta=2 corresponds to Gaussian distribution)
    """
    def __init__ ( self             ,
                   name             ,
                   mn        = None ,
                   mx        = None ,
                   mass      = None ,
                   mean      = None ,
                   alpha     = None ,
                   beta      = 2    ) :  ## beta=2 is gaussian distribution 
        
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self  , name  ,
                         mn    , mx    , mass  ,
                         mean  , alpha ) 
        
        #
        ## rename it!
        #
        self.alpha = self.sigma
        sname  = self.alpha.GetName  ()
        stitle = self.alpha.GetTitle ()
        gname  = sname .replace ( 'sigma' , 'alpha' )
        gtitle = stitle.replace ( 'sigma' , 'alpha' )
        self.alpha.SetName  ( gname  ) 
        self.alpha.SetTitle ( gtitle )
        
        self.beta  = makeVar ( beta ,
                               'beta_%s'        % name  ,
                               '#beta_{v1}(%s)' % name  , beta , 
                               2 , 1.e-4  , 1.e+6 ) 
        
        #
        ## finally build PDF
        #
        self.pdf = Ostap.Models.GenGaussV1 (
            "gengV1_"        + name ,
            "GenGaussV1(%s)" % name ,
            self.mass   ,
            self.mean   ,
            self.alpha  ,
            self.beta   )

models.append ( GenGaussV1_pdf )    
# =============================================================================
## @class GenGaussV2_pdf
#  Generalized Normal distribution v2
#  @see http://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_2
#  @see Ostap::Models::GenGaussV2
#  @see Ostap::Math::GenGaussV2 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class GenGaussV2_pdf(MASS) :
    """Generalized normal distribution v2
    see http://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_2

    This is a family of continuous probability distributions in which the shape
    parameter can be used to introduce skew.
    When the shape parameter is zero, the normal distribution results.
    Positive values of the shape parameter yield left-skewed distributions
    bounded to the right, and negative values of the shape parameter yield
    right-skewed distributions bounded to the left.

    Parameters:
    - mu      : location 
    - alpha>0 : scale
    - kappa   : shape   (kappa=0 corresponds to gaussian distribution)
    """
    def __init__ ( self               ,
                   name               ,
                   mn          = None ,
                   mx          = None ,
                   mass        = None ,
                   mean        = None ,
                   alpha       = None ,
                   kappa       = 0    ) : ## 0 corresponds to gaussian distribution 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name  ,
                         mn      , mx    , mass    ,
                         mean    , alpha ) 
        
        #
        ## rename it!
        #
        self.alpha = self.sigma        
        sname      = self.alpha.GetName  ()
        stitle     = self.alpha.GetTitle ()
        gname      = sname .replace ( 'sigma' , 'alpha' )
        gtitle     = stitle.replace ( 'sigma' , 'alpha' )
        self.alpha.SetName  ( gname  ) 
        self.alpha.SetTitle ( gtitle )
        
        self.xi   = self.mean 
        ## sname     = self.xi.GetName  ()
        ## stitle    = self.xi.GetTitle ()
        ## gname     = sname .replace   ( 'mean' , 'xi' )
        ## gtitle    = stitle.replace   ( 'mean' , 'xi' )
        ## self.xi.SetName              ( gname  ) 
        ## self.xi.SetTitle             ( gtitle )
        
        self.kappa = makeVar ( kappa ,
                               'kappa_%s'        % name  ,
                               '#kappa_{v2}(%s)' % name  , kappa , 
                               0 , -4  , 4 ) 
        
        #
        ## finally build PDF
        #
        self.pdf = Ostap.Models.GenGaussV2 (
            "gengV2_"        + name ,
            "GenGaussV2(%s)" % name ,
            self.mass   ,
            self.xi     ,
            self.alpha  ,
            self.kappa  )

models.append ( GenGaussV2_pdf )    
## # =============================================================================
## ## @class SkewGauss_pdf
## #  Simple class that implements the skew normal distribution
## #  @see http://en.wikipedia.org/wiki/Skew_normal_distribution
## #  @see Ostap::Models::SkewGauss 
## #  @see Ostap::Math::SkewGauss 
## #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
## #  @date 2011-07-25
## class SkewGauss_pdf(MASS) :
##     """Skew normal distribution
##     see http://en.wikipedia.org/wiki/Skew_normal_distribution
    
##     The skew normal distribution is a continuous probability distribution that
##     generalises the normal distribution to allow for non-zero skewness.

##     Parameters:
##     - location
##     - omega>0   : scale
##     - alpha     : shape   (alpha=0 corresponds to gaussian distribuition)
##     """
##     def __init__ ( self             ,
##                    name             ,
##                    mn        = None ,
##                    mx        = None ,
##                    mass      = None ,
##                    mean      = None ,
##                    omega     = None ,
##                    alpha     = 0    ) : ## alpha=0 correspond to gaussian 
##         #
##         ## initialize the base
##         # 
##         MASS.__init__  ( self    , name  ,
##                          mn      , mx    , mass    ,
##                          mean    , omega ) 
        
##         self.omega = self.sigma
##         sname  = self.omega.GetName  ()
##         stitle = self.omega.GetTitle ()
##         gname  = sname .replace ( 'sigma' , 'omega' )
##         gtitle = stitle.replace ( 'sigma' , 'omega' )
##         self.omega.SetName  ( gname  ) 
##         self.omega.SetTitle ( gtitle )

##         self.alpha = makeVar ( alpha                      ,
##                                'alpha_%s'          % name ,
##                                '#alpha_{Skew}(%s)' % name , alpha , 
##                                0 , -100 , 100  ) 
##         #
##         ## finally build pdf
##         # 
##         self.pdf = Ostap.Models.SkewGauss (
##             "skewg_"         + name ,
##             "SkewGauss(%s)" % name ,
##             self.mass   ,
##             self.mean   ,
##             self.omega  ,
##             self.alpha  )

## models.append ( SkewGauss_pdf )      
# =============================================================================
## @class Bukin_pdf
#  Bukin function, aka ``modified Novosibirsk function''
#  - asymmetrical gaussian-like core
#  - exponential (optionally gaussian) asymmetrical tails
#  @see http://journals.aps.org/prd/abstract/10.1103/PhysRevD.84.112007
#  @http://arxiv.org/abs/1107.5751
#  Here small reparameterization is applied to achieve more stable fits.
# 
#  It is very well suitable to describe high statistic charm meson peaks,
#  as well some asymmetric distributions, e.g. log( chi^2(IP) ) 
#
#  Parameters:
#  - mean      : peak position 
#  - sigma     : defines full width at half heigh 
#  - xi        : asymmetry
#  - rho_l>=0  : define gaussian component to left  tail
#  - rho_r>=0  : define gaussian component to right tail
#
#  Note: 
#  - for rho_{l,r}=0 left/right tails are exponential.
#  - for large asymmetry parameter function has weird shape
#
#  @see http://dx.doi.org/10.1007/JHEP06(2012)141     
#  @see Ostap::Math::Bukin
#  @see Analusis::Models::Bukin
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Bukin_pdf(MASS) :
    """ Bukin function, aka ``modified Novosibirsk function'':
    - asymmetrical gaussian-like core
    - exponential (optionally gaussian) asymmetrical tails
    see http://journals.aps.org/prd/abstract/10.1103/PhysRevD.84.112007
    see http://arxiv.org/abs/1107.5751
    see http://dx.doi.org/10.1007/JHEP06(2012)141     
    Here small reparameterization is applied to achieve more stable fits.
    
    It is very well suitable to describe high statistic charm meson peaks,
    as well some asymmetric distributions, e.g. log( chi^2(IP) ) 

    Parameters:
    - mean      : peak position 
    - sigma     : defines full width at half heigh 
    - xi        : asymmetry
    - rho_l>=0  : define gaussian component to left  tail
    - rho_r>=0  : define gaussian component to right tail
    Note: 
    - for rho_{l,r}=0 left/right tails are exponential
    - for large asymmetry parameter function has weird shape    
    """
    def __init__ ( self            ,
                   name            ,
                   mn       = None , ## low-edge   (not used if 'mass' is specified)
                   mx       = None , ## high-edge  (not used if 'mass' is specified) 
                   mass     = None , 
                   mean     = None ,
                   sigma    = None ,
                   xi       = None ,
                   rhol     = None ,
                   rhor     = None ) :

        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name  ,
                         mn      , mx    , mass    ,
                         mean    , sigma )
        #
        ## treat the specific parameters
        #
        ## asymmetry 
        self.xi    = makeVar ( xi                    ,
                               "xi_%s"        % name ,
                               "#xi(%s)"      % name , xi      , 0  , -1 , 1    )
        
        self.rhol  = makeVar ( rhol                  ,
                               "rhol_%s"      % name ,
                               "#rho_{L}(%s)" % name , rhol    , 0  ,  -10 , 10 )
        
        self.rhor  = makeVar ( rhor                  ,
                               "rhor_%s"      % name ,
                               "#rho_{R}(%s)" % name , rhor    , 0  ,  -10 , 10 )
        # 
        ## create PDF
        # 
        self.pdf = Ostap.Models.Bukin (
            "bkn_"      + name ,
            "Bukin(%s)" % name ,
            self.mass  ,
            self.mean  ,
            self.sigma ,
            self.xi    ,
            self.rhol  ,
            self.rhor  )

models.append ( Bukin_pdf )      
# =============================================================================
## @class StudentT_pdf
#  Student-T distribution
#  @see http://en.wikipedia.org/wiki/Student%27s_t-distribution
#
#  \f[  f(y) = \frac{1}{\sqrt{\pi n}} \frac { \Gamma( \frac{n+1}{2}) } { \Gamma( \frac{n}{2}  ) }
#  \left( 1 + \frac{y^2}{n} \right)^{ -\frac{n+1}{2}} \f], 
#  where \f$ y = \frac{x - \mu}{\sigma} \f$  
# 
#  @attention Unlike the original definition parameter 'n' here is shifted by 1: n <- |n| + 1
#  @see Ostap::Models::StudentT
#  @see Ostap::Math::StudentT
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class StudentT_pdf(MASS) :
    """Student-T distribution:
    http://en.wikipedia.org/wiki/Student%27s_t-distribution
    
    f(dx) ~ (1 + dx^2/n)^{-(n+1)/2 }
    with
    dx = ( x - mu )  / sigma
    
    Note:
    - unlike the original definition parameter 'n' here is shifted by 1: n <- |n| + 1
    
    Parameters 
    - mean    : location 
    - sigma   : scale 
    - n       : n-parameter, |n|+1 is used 
    """
    def __init__ ( self             ,
                   name             ,
                   mn        = None ,
                   mx        = None ,
                   mass      = None ,
                   mean      = None ,
                   sigma     = None ,
                   n         = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name  ,
                         mn      , mx    , mass    ,
                         mean    , sigma ) 
        
        # 
        self.n  = makeVar ( n                    ,
                            'n_%s'        % name ,
                            '#n_{ST}(%s)' % name , n , 
                            2 , 0 , 100  ) 
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.StudentT (
            "stT_"         + name ,
            "StudentT(%s)" % name ,
            self.mass   ,
            self.mean   ,
            self.sigma  ,
            self.n      )

models.append ( StudentT_pdf )      
# =============================================================================
## @class BifurcatedStudentT_pdf
#  bifurcated Student-T distribution
# 
#  @see Ostap::Models::BifurcatedStudentT
#  @see Ostap::Math::BifurcatedStudentT
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class BifurcatedStudentT_pdf(MASS) :
    """Bifurcated Student-T distribution:

    f(dx) ~ (1 + dx^2/n)^{-(n+1)/2 }
    where
    for x < mu   n=n_l and dx = ( x - mu )  / sigma_l
    for x > mu   n=n_l and dx = ( x - mu )  / sigma_r
    with 
    sigma_{l,r} = sigma * ( 1 + kappa)
    with asymmetry parameter -1 < kappa < 1 

    Parameters:
    - mean
    - sigma
    - -1<asymmetry<1 
    - n_L     
    - n_R
    """
    def __init__ ( self             ,
                   name             ,
                   mn        = None , 
                   mx        = None , 
                   mass      = None ,
                   mean      = None ,
                   sigma     = None ,
                   asymmetry = None ,
                   nL        = None , 
                   nR        = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name  ,
                         mn      , mx    , mass    ,
                         mean    , sigma )
        
        self.asym = makeVar ( asymmetry                  ,
                              'asym_%s'        % name ,
                              '#xi_{asym}(%s)' % name , asymmetry , 
                              0 , -1 , 1  ) 
        
        self._lst_R = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.sigmaR = ROOT.RooFormulaVar (
            "sigmaR_stt_%s"     % name   ,
            "sigma_{R}(%s)" % name   ,
            "%s*(1+%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self._lst_R   )
        
        self._lst_L = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.sigmaL = ROOT.RooFormulaVar (
            "sigmaL_stt_%s"     % name   ,
            "sigma_{L}(%s)" % name   ,
            "%s*(1-%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self._lst_L   )
        
        # 
        self.nL =  makeVar ( nL                     ,
                             'nL_%s'         % name ,
                             '#nL_{BST}(%s)' % name , nL , 
                             2  , 0 , 100  )
        
        self.nR =  makeVar ( nR                    ,
                             'nR_%s'         % name ,
                             '#nR_{BST}(%s)' % name , nR , 
                             2  , 0 , 100  ) 
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.BifurcatedStudentT (
            "bstT_"         + name ,
            "BStudentT(%s)" % name ,
            self.mass   ,
            self.mean   ,
            self.sigmaL ,
            self.sigmaR ,
            self.nL     ,
            self.nR     ) 
        
models.append ( BifurcatedStudentT_pdf )      
# =============================================================================
## @class SinhAsinh_pdf
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
#  @date   2014-08-02
class SinhAsinh_pdf(MASS) :
    """SinhAsinh-function: 
    see Jones, M. C.; Pewsey, A. (2009).
    ``Sinh-arcsinh distributions''. Biometrika 96 (4): 761. 
    doi:10.1093/biomet/asp053
    http://oro.open.ac.uk/22510
    
    Normal distribution reappears as epsilon=0 and delta=1 
    The heavy tails correspond to delta<1, light tails correpond to delta>1
    
    Parameters 
    - mean     : location
    - sigma    : scale 
    - epsilon  : parameter to control the skewness 
    - delta>0  : parameter to control the kurtosis 
    """
    def __init__ ( self             ,
                   name             ,
                   mn        = None ,
                   mx        = None , 
                   mass      = None ,
                   mean      = None , ## mu 
                   sigma     = None ,
                   epsilon   = 0    ,
                   delta     = 1    ) :

        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name  ,
                         mn      , mx    , mass    ,
                         mean    , sigma )

        self.mu      = self.mean
        self.epsilon = makeVar ( epsilon ,
                                 'epsilon_%s'   % name ,
                                 '#epsilon(%s)' % name , epsilon ,
                                 0 , -1000 , +1000 )
        self.delta   = makeVar ( delta ,
                                 'delta_%s'   % name ,
                                 '#delta(%s)' % name , delta ,
                                 1 , 1.e-6 , 1000   )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.SinhAsinh (
            "sinhaT_"        + name ,
            "SinhAsinh(%s)" % name ,
            self.mass      ,
            self.mean      ,
            self.sigma     ,
            self.epsilon   ,
            self.delta     ) 

models.append ( SinhAsinh_pdf )      


# =============================================================================
## @class JohnsonSU_pdf
#
#  Johnson, N. L. (1949) 
#  "Systems of frequency curves generated by methods of translation"
#  Biometrika 36: 149–176 JSTOR 2332539
#  @see https://en.wikipedia.org/wiki/Johnson_SU_distribution
#
#  When variable \f$x\f$ follows Johnson-SU distribution, 
#  the variable 
#  \f$ z = \gamma + \delta \sinh^{-1}\frac{ x - \xi}{\lambda} \f$
#  follows normal distribtion with mean 0 and sigma 1.
#
#  Note:
#  Symmetric case of JonhsonSU distribution is 
#  recovere by \f$\delta\rightarrow0\f$ for 
#  "sinh-asinh" distribution, see 
#  Jones, M. C.; Pewsey, A. (2009). 
#  "Sinh-arcsinh distributions". Biometrika 96 (4): 761. 
#  doi:10.1093/biomet/asp053
#  http://oro.open.ac.uk/22510
#  
#  @see Ostap::Math::JohnsonSU
#  @see Ostap::Models::JohnsonSU
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
#  @date   2015-07-019
class JohnsonSU_pdf(MASS) :
    """Johnson, N. L. (1949) 
    Systems of frequency curves generated by methods of translation
    Biometrika 36: 149–176 JSTOR 2332539
    
    https://en.wikipedia.org/wiki/Johnson_SU_distribution
    
    When variable x follows Johnson-SU distribution, 
    the variable  z = gamma +  delta  sinh^{-1} ( x - xi)/ lambda 
    follows normal distribtion with mean 0 and sigma 1.
    
    Note:
    Symmetric case of JonhsonSU distribution is 
    recovered by delta -> 0  for 
    'sinh-asinh' distribution, see 
    Jones, M. C.; Pewsey, A. (2009). 
    'Sinh-arcsinh distributions'. Biometrika 96 (4): 761. 
    doi:10.1093/biomet/asp053
    http://oro.open.ac.uk/22510
    
    """
    def __init__ ( self             ,
                   name             ,
                   mn        = None ,
                   mx        = None , 
                   mass      = None ,
                   xi        = None ,   ## related to mean 
                   lam       = None ,   ## related to sigma 
                   delta     = 1    ,
                   gamma     = 0    ) :  

        #
        ## initialize the base
        #
        
        MASS.__init__  ( self    , name   ,
                         mn      , mx     , mass    ,
                         xi      , lam    ) ## mean    , sigma  )

        self.xi = self.mean
        if None == xi : ## newly created 
            oname   = self.xi.GetName ()
            otitle  = self.xi.GetTitle()
            nname   = oname .replace ( 'mean' , 'xi' )
            ntitle  = otitle.replace ( 'mean' , 'xi' )
            self.xi.SetName  ( nname  )
            self.xi.SetTitle ( ntitle )
            
        self.lam = self.sigma
        if None == lam : ## newly created 
            oname    = self.lam.GetName ()
            otitle   = self.lam.GetTitle()
            nname    = oname .replace ( 'sigma' , 'lambda' )
            ntitle   = otitle.replace ( 'sigma' , 'lambda' )
            self.lam.SetName  ( nname  )
            self.lam.SetTitle ( ntitle )
            self.lam.setMax ( self.lam.getMax() * 10 ) ## adjust it! 


        ## provdie backup name 
        self.lambda_ = self.lam
        
        del self.mean
        del self.sigma
        
        self.delta   = makeVar ( delta                 ,
                                 'delta_%s'     % name ,
                                 '#delta(%s)'   % name , delta ,
                                 1 , 1.e-6 , 1000   )
        
        self.gamma   = makeVar ( gamma               ,
                                 'gamma_%s'   % name ,
                                 '#gamma(%s)' % name , gamma ,
                                 0 , -1000 , +1000 )
        
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.JohnsonSU (
            "jSU_"          + name ,
            "JohnsonSU(%s)" % name ,
            self.mass      ,
            self.xi        ,
            self.lam       ,
            self.delta     ,
            self.gamma     ) 

models.append ( JohnsonSU_pdf )      
# =============================================================================
## @class Atlas_pdf
#  Modified gaussian with exponential tails
#  \f$  f(x) \propto \exp( -frac{\delta x^{1+\frac{1}{1+\deltax/2}}}{2})\f$,
#  where \f$\delta x = \left| x - \mu \right|/\sigma\f$
#  Function is taken from http://arxiv.org/abs/arXiv:1507.07099
# 
#  @see http://arxiv.org/abs/hep-ex/0505008 
#  @see http://arxiv.org/abs/1206.3122
#  @see http://arxiv.org/abs/arXiv:1507.07099
#
#  @see Ostap::Math::Atlas 
#  @see Ostap::Models::Atlas
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-08-024
class Atlas_pdf(MASS) :
    """Modified gaussian with exponential tails
    \f$  f(x) \propto \exp( -frac{\delta x^{1+\frac{1}{1+\deltax/2}}}{2})\f$,
    where \f$\delta x = \left| x - \mu \right|/\sigma\f$
    Function is taken from http://arxiv.org/abs/arXiv:1507.07099    
    """
    def __init__ ( self             ,
                   name             ,
                   mn        = None ,
                   mx        = None , 
                   mass      = None ,
                   mean      = None ,   ## related to mean 
                   sigma     = None ) : ## related to sigma 

        #
        ## initialize the base
        #
        
        MASS.__init__  ( self    , name   ,
                         mn      , mx     , mass    ,
                         mean    , sigma  )

        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Atlas (
            "atlas_"    + name ,
            "ATLAS(%s)" % name ,
            self.mass      ,
            self.mean      ,
            self.sigma     ) 

models.append ( Atlas_pdf )      

# =============================================================================
## @class Sech_pdf
#  Hyperbolic secant distribution or "inverse-cosh" distribution
# 
#  The hyperbolic secant distribution shares many properties with the 
#  standard normal distribution: 
#  - it is symmetric with unit variance and zero mean, 
#    median and mode
#  -its pdf is proportional to its characteristic function. 
#
#  However, the hyperbolic secant distribution is leptokurtic; 
#  that is, it has a more acute peak near its mean, and heavier tails, 
#  compared with the standard normal distribution.
#
#  \f$ f(x,\mu,\sigma) \propto \frac{1}{2} \sech ( \frac{\pi}{2}\frac{x-\mu}{\sigma} )\f$ 
#   @see https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
#
#  @see Ostap::Math::Sech
#  @see Ostap::Models::Sech
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-04-025
class Sech_pdf(MASS) :
    """Hyperbolic secant distribution or ``inverse-cosh'' distribution
    
    The hyperbolic secant distribution shares many properties with the 
    standard normal distribution: 
    - it is symmetric with unit variance and zero mean, median and mode
    -its pdf is proportional to its characteristic function. 
     
    However, the hyperbolic secant distribution is leptokurtic; 
    that is, it has a more acute peak near its mean, and heavier tails, 
    compared with the standard normal distribution.
    
    f(x,\mu,\sigma) \propto \frac{1}{2} \sech ( \frac{\pi}{2}\frac{x-\mu}{\sigma} )
    
    @see https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
    """
    def __init__ ( self             ,
                   name             ,
                   mn        = None ,
                   mx        = None , 
                   mass      = None ,
                   mean      = None ,   ## related to mean 
                   sigma     = None ) : ## related to sigma 

        #
        ## initialize the base
        #
        
        MASS.__init__  ( self    , name   ,
                         mn      , mx     , mass    ,
                         mean    , sigma  )

        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Sech (
            "sech_"    + name ,
            "SECH(%s)" % name ,
            self.mass      ,
            self.mean      ,
            self.sigma     ) 
        
models.append ( Sech_pdf )      



# =============================================================================
## @class Logistic_pdf
#  Logistic, aka "sech-square" PDF
#  \f$ f(x;\mu;s) = \frac{1}{4s}sech^2\left(\frac{x-\mu}{2s}\right)\f$, 
#   where
#   \f$  s = \sigma \frac{\sqrt{3}}{\pi}\f$
#  @see https://en.wikipedia.org/wiki/Logistic_distribution
#  @see Ostap::Math::Logistic
#  @see Ostap::Models::Logistic
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-06-14
class Logistic_pdf(MASS) :
    """ Logistic, aka ``sech-square'' PDF
     \f$ f(x;\mu;s) = \frac{1}{4s}sech^2\left(\frac{x-\mu}{2s}\right)\f$, 
     where
     \f$  s = \sigma \frac{\sqrt{3}}{\pi}\f$
     - see https://en.wikipedia.org/wiki/Logistic_distribution
     - see Ostap::Math::Logistic
     - see Ostap::Models::Logistic
    """
    def __init__ ( self             ,
                   name             ,
                   mn        = None ,
                   mx        = None , 
                   mass      = None ,
                   mean      = None ,   ## related to mean 
                   sigma     = None ) : ## related to sigma 

        #
        ## initialize the base
        #
        
        MASS.__init__  ( self    , name   ,
                         mn      , mx     , mass    ,
                         mean    , sigma  )

        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Logistic (
            "logistic_"    + name ,
            "Logistic(%s)" % name ,
            self.mass      ,
            self.mean      ,
            self.sigma     ) 
        
models.append ( Logistic_pdf )      
# =============================================================================
## @class Voigt_pdf
#  Voigt-pdf distribution
#  @see Ostap::Models::Voigt
#  @see Ostap::Math::Voigt
#  The implementation relied on Faddeeva function 
#  @see http://en.wikipedia.org/wiki/Faddeeva_function
#  @see http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Voigt_pdf(MASS) :
    """Voigt function:
    Convolution of non-relativistic Breit-Wigner with gaussian resolution
    
    The implementation relied on Faddeeva function 
    http://en.wikipedia.org/wiki/Faddeeva_function
    http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package
    
    Parameters
    - mean  : location 
    - gamma : gamma for breight-wigner pole
    - sigma : resolution parameter for gaussian 
    """
    def __init__ ( self             ,
                   name             ,
                   mn        = None ,
                   mx        = None ,
                   mass      = None ,
                   mean      = None ,
                   sigma     = None ,
                   gamma     = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name , mn , mx ,
                         mass    ,
                         mean    , gamma ) 

        #
        ##  rename it
        #
        self.gamma = self.sigma
        sname  = self.gamma.GetName  ()
        stitle = self.gamma.GetTitle ()
        gname  = sname .replace ( 'sigma' , 'gamma' )
        gtitle = stitle.replace ( 'sigma' , 'Gamma' )
        self.gamma.SetName  ( gname  ) 
        self.gamma.SetTitle ( gtitle )

        self.sigma  = makeVar ( sigma               ,
                                'sigma_%s'   % name ,   
                                '#sigma(%s)' % name , sigma , 
                                0.0001 * ( mass.getMax()  - mass.getMin() ) , 
                                0    ,
                                0.3000 * ( mass.getMax()  - mass.getMin() ) )
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Voigt (
            "vgt_"       + name ,
            "Voigt(%s)" % name ,
            self.mass   ,
            self.mean   ,
            self.gamma  ,
            self.sigma  )

models.append ( Voigt_pdf )                          
# =============================================================================
## @class PseudoVoigt_pdf
#  PseudoVoigt-pdf distribution
#  @see Ostap::Models::PseudoVoigt
#  @see Ostap::Math::PseudoVoigt
#  CPU-efficient Approximation of Voight profile
#  @see T. Ida, M. Ando and H. Toraya, 
#       "Extended pseudo-Voigt function for approximating the Voigt profile"
#       J. Appl. Cryst. (2000). 33, 1311-1316
#  @see doi:10.1107/S0021889800010219
#  @see http://dx.doi.org/10.1107/S0021889800010219
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2016-06-15
class PseudoVoigt_pdf(MASS) :
    """Voigt function:
    Convolution of non-relativistic Breit-Wigner with gaussian resolution
    
    CPU-efficient Approximation of Voight profile
    -@see T. Ida, M. Ando and H. Toraya, 
    ``Extended pseudo-Voigt function for approximating the Voigt profile''
    J. Appl. Cryst. (2000). 33, 1311-1316
    - see doi:10.1107/S0021889800010219
    - see http://dx.doi.org/10.1107/S0021889800010219
    
    Parameters
    - mean  : location 
    - gamma : gamma for breight-wigner pole
    - sigma : resolution parameter for gaussian 
    """
    def __init__ ( self             ,
                   name             ,
                   mn        = None ,
                   mx        = None ,
                   mass      = None ,
                   mean      = None ,
                   sigma     = None ,
                   gamma     = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name , mn , mx ,
                         mass    ,
                         mean    , gamma ) 

        #
        ##  rename it
        #
        self.gamma = self.sigma
        sname  = self.gamma.GetName  ()
        stitle = self.gamma.GetTitle ()
        gname  = sname .replace ( 'sigma' , 'gamma' )
        gtitle = stitle.replace ( 'sigma' , 'Gamma' )
        self.gamma.SetName  ( gname  ) 
        self.gamma.SetTitle ( gtitle )

        self.sigma  = makeVar ( sigma               ,
                                'sigma_%s'   % name ,   
                                '#sigma(%s)' % name , sigma , 
                                0.0001 * ( mass.getMax()  - mass.getMin() ) , 
                                0    ,
                                0.3000 * ( mass.getMax()  - mass.getMin() ) )
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.PseudoVoigt (
            "pvgt_"           + name ,
            "PseudoVoigt(%s)" % name ,
            self.mass   ,
            self.mean   ,
            self.gamma  ,
            self.sigma  )

models.append ( PseudoVoigt_pdf )                          
# =============================================================================
## @class BreitWigner_pdf 
#  Relativistic Breit-Wigner function using Jackson's parameterization
#  J.D.Jackson,
#  "Remarks on the Phenomenological Analysis of Resonances",
#  In Nuovo Cimento, Vol. XXXIV, N.6
#  http://www.springerlink.com/content/q773737260425652/
#  Optional convolution with resolution function is possible a
#  @see Ostap::Models::BreitWigner 
#  @see Ostap::Math::BreitWigner
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class BreitWigner_pdf(MASS) :
    """Relativistic Breit-Wigner function using Jackson's parameterization
    J.D.Jackson, ``Remarks on the Phenomenological Analysis of Resonances'',
    In Nuovo Cimento, Vol. XXXIV, N.6
    http://www.springerlink.com/content/q773737260425652/

    >>> bw    = Ostap.Math.BreitWigner( m_X , g_X , m_Jpsi , m_pipi , 0 )
    >>> breit = Models.BreitWigner_pdf ( 'BW'          ,
    ...                                  bw            ,
    ...                                  mass  = mass  ,
    ...                                  mean  = m_X   ,
    ...                                  gamma = g_X   )
    
    Optional convolution with the resolution function is possible
    e.g. use gaussian resoltuion and fast-fourier convolution method:
    
    >>> breit = Models.BreitWigner_pdf ( 'BW'          ,
    ...                                  bw            ,
    ...                                  mass  = mass  ,
    ...                                  mean  = m_X   ,
    ...                                  gamma = g_X   ,
    ...                                  convolution = 1*MeV ,
    ...                                  useFFT      = True  ) 
    
    Other resolution functions can be used:
    
    >>> resolution_pdf = ...
    >>> breit = Models.BreitWigner_pdf ( 'BW'          ,
    ...                                  ...
    ...                                  convolution = resolution_pdf ,
    ...                                  useFFT      = True  ) 

    
    Parameters:
    - mean  : location
    - gamma : width of Breight-Wigner function
    
    """
    def __init__ ( self               ,
                   name               ,
                   breitwigner        , ## Ostap::Math::BreitWeigner object
                   mn          = None ,
                   mx          = None ,
                   mass        = None ,
                   mean        = None , 
                   gamma       = None ,
                   convolution = None ,
                   useFFT      = True ) :
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name     ,
                         mn      , mx       ,
                         mass    ,
                         mean    , gamma    )
        
        self.gamma = self.sigma
        sname  = self.gamma.GetName  ()
        stitle = self.gamma.GetTitle ()
        gname  = sname .replace ( 'sigma' , 'gamma' )
        gtitle = stitle.replace ( 'sigma' , 'Gamma' )
        self.gamma.SetName  ( gname  ) 
        self.gamma.SetTitle ( gtitle )
        del self.sigma
        
        #
        ## define the actual BW-shape using
        #      Ostap::Math::BreitWeigner object
        #
        self.breitwigner = breitwigner  ## Ostap::Math::BreitWeigner object
        self.bw          = breitwigner  
        
        ## create PDF 
        self.breit = Ostap.Models.BreitWigner ( 
            "rbw_"    + name ,
            "RBW(%s)" % name ,
            self.mass        ,
            self.mean        ,
            self.gamma       ,
            self.breitwigner )

        if  None is convolution : self.pdf = self.breit
        else :
            from ostap.fitting.basic import Convolution 
            self.conv = Convolution ( name        ,
                                      self.breit  , self.mass ,
                                      convolution , useFFT    ) 
            self.pdf  = self.conv.pdf

models.append ( BreitWigner_pdf )                          
# =============================================================================
## @class BW23L_pdf
#  The shape of Breit-Wigner resonace from 2-body decays, e.f.
#  @see Ostap::Models::BW23L 
#  @see Ostap::Math::BW23L
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-25
class BW23L_pdf(MASS) :
    """The shape of Breit-Wigner resonace from 3-body decays, e.f. X -> ( A B ) C
    In this case the phase space factors can be modified by the  orbital momentum
    between (AB) and C-systems
    """
    def __init__ ( self               ,
                   name               ,
                   breitwigner        , ## Ostap::Math::BW23L object
                   mn          = None ,
                   mx          = None ,
                   mass        = None ,
                   mean        = None , 
                   gamma       = None ) : 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name     ,
                         mn      , mx       ,
                         mass    ,
                         mean    , gamma    )
        
        self.gamma = self.sigma
        sname  = self.gamma.GetName  ()
        stitle = self.gamma.GetTitle ()
        gname  = sname .replace ( 'sigma' , 'gamma' )
        gtitle = stitle.replace ( 'sigma' , 'Gamma' )
        self.gamma.SetName  ( gname  ) 
        self.gamma.SetTitle ( gtitle )
        
        #
        ## define the actual BW-shape using
        #      Ostap::Math::BW23L object
        #
        self.breitwigner = breitwigner  ## Ostap::Math::BW23L object
        self.bw          = breitwigner  
        
        ## create PDF 
        self.pdf = Ostap.Models.BW23L ( 
            "rbw23_"    + name ,
            "RBW23(%s)" % name ,
            self.mass          ,
            self.mean          ,
            self.gamma         ,
            self.breitwigner   )

models.append ( BW23L_pdf )                          
# =============================================================================
## @class Flatte_pdf
#  Flatte function
#  S.M.Flatte, "Coupled-channel analysis of the \f$\pi\eta\f$ 
#    and \f$K\bar{K}\f$ systems near \f$K\bar{K}\f$ threshold  
#    Phys. Lett. B63, 224 (1976)
#  Well suitable for \f$\f_0(980)\rightarrow \pi^+ \pi^-\f$
#  @see http://www.sciencedirect.com/science/article/pii/0370269376906547
#  @see Ostap::Models::Flatte
#  @see Ostap::Math::Flatte
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-01-18
class Flatte_pdf(MASS) :
    """Flatte function:
    S.M.Flatte, ``Coupled-channel analysis of the (pi eta)
    and (KbarK) systems near (KbarK) threshold'' 
    Phys. Lett. B63, 224 (1976
    http://www.sciencedirect.com/science/article/pii/0370269376906547

    Typical case:    f0 -> pi+ pi- shape 
    """
    def __init__ ( self              ,
                   name              ,
                   flatte            , ## Ostap::Math::Flatte 
                   mn       = None   ,
                   mx       = None   ,
                   mass     = None   ,
                   m0_980   = None   ,    ## mass  of f0(980) resonance
                   m0g1     = 165000 ,    ## m0(f0(980))*gamma_1 
                   g2og1    = 4.21   ) :  ## gamma2/gamma1 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name     ,
                         mn      , mx       , mass    ,
                         m0_980  , None     )
        
        del self.sigma 

        self.flatte = flatte 

        self.m0_980 = self.mean 
        sname  = self.mean.GetName  ()
        stitle = self.mean.GetTitle ()
        gname  = sname .replace ( 'mean' , 'm0_980' )
        gtitle = stitle.replace ( 'mean' , 'm0_980' )
        self.m0_980.SetName  ( gname  ) 
        self.m0_980.SetTitle ( gtitle ) 
        
        self.m0g1 = makeVar  ( m0g1                          ,
                               'm0g1_%s'              % name ,
                               'm_{0}*\gamma_{1}(%s)' % name , m0g1 ,
                               165                           ,
                               1.e-5                         ,
                               1.e+5                         )
        
        self.g2og1 = makeVar ( g2og1    ,
                               'g2og1_%s'                  % name ,
                               '#gamma_{2}/#gamma_{1}(%s)' % name , g2og1 , 
                               4.21     , 
                               0.01     , 100 ) 


        #
        ## build the actual pdf
        #
    
        ## create PDF 
        self.pdf = Ostap.Models.Flatte ( 
            "flatte_"    + name ,
            "Flatte(%s)" % name ,
            self.mass    ,
            self.m0_980  ,
            self.m0g1    ,
            self.g2og1   ,
            self.flatte  )
        
models.append ( Flatte_pdf )                          
# =============================================================================
## @class Flatte2_pdf
#  Flatte function to describe dikaon system near threshold
#  S.M.Flatte, 
#    "Coupled-channel analysis of the \f$\pi\eta\f$ 
#    and \f$K\bar{K}\f$ systems near \f$K\bar{K}\f$ threshold  
#    Phys. Lett. B63, 224 (1976)
#  Well suitable for \f$\f_0(980)\rightarrow K^+ K^-\f$
#  @see http://www.sciencedirect.com/science/article/pii/0370269376906547
#  @see Ostap::Models::Flatte2
#  @see Ostap::Math::Flatte2
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-01-18
class Flatte2_pdf(Flatte_pdf) :
    """Flatte function:
    S.M.Flatte, ``Coupled-channel analysis of the (pi eta) and (KbarK) systems near (KbarK) threshold'' 
    Phys. Lett. B63, 224 (1976
    http://www.sciencedirect.com/science/article/pii/0370269376906547
    
    Typical case:    f0 -> K+ K- shape 
    """
    def __init__ ( self              ,
                   name              ,
                   flatte            , ## Ostap::Math::Flatte or Flatte2 
                   mn       = None   ,
                   mx       = None   ,
                   mass     = None   ,
                   m0_980   = None   ,    ## mass  of f0(980) resonance
                   m0g1     = 165000 ,    ## m0(f0(980))*gamma_1
                   g2og1    = 4.21   ) :  ## gamma2/gamma1 
        
        #
        ## initialize the base
        # 
        Flatte_pdf.__init__  ( self     , name , flatte ,
                               mn       , mx ,
                               mass     ,
                               m0_980   ,
                               m0g1     ,
                               g2og1    )

        ## delete the creatd pdf (not-needed) 
        del self.pdf
        
        ## build the actual pdf
        # 
        ## create PDF 
        self.pdf = Ostap.Models.Flatte2 ( 
            "flatte2_"    + name ,
            "Flatte2(%s)" % name ,
            self.mass    ,
            self.m0_980  ,
            self.m0g1    ,
            self.g2og1   ,
            self.flatte  ) 
        
models.append ( Flatte2_pdf )                          
# =============================================================================
## @class LASS_pdf
#  The LASS parameterization (Nucl. Phys. B296, 493 (1988))
#  describes the 0+ component of the Kpi spectrum.
#  It consists of the K*(1430) resonance together with an
#  effective range non-resonant component
#  @see Ostap::Models::LASS
#  @see Ostap::Math::LASS
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class LASS_pdf(MASS) :
    """Kappa pole:
    The LASS parameterization (Nucl. Phys. B296, 493 (1988))
    describes the 0+ component of the Kpi spectrum.
    It consists of the K*(1430) resonance together with an
    effective range non-resonant component
    """
    def __init__ ( self               ,
                   name               ,
                   mn       = None    ,
                   mx       = None    ,
                   mass     = None    ,
                   m0_1430  = None    ,     ## mass  of K*(1430)
                   g0_1430  = None    ,     ## width of K*(1430)
                   a_lass   = 1.94e-3 , 
                   r_lass   = 1.76e-3 ,
                   e_lass   = 1.0     ,    ## elasticity                    
                   mKaon    = 493.7   ,    ## kaon mass 
                   mPion    = 139.6   ) :  ## pion mass 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name     ,
                         mn      , mx       ,
                         mass    ,
                         m0_1430 , g0_1430 )
        
        self.gamma = self.sigma
        sname  = self.gamma.GetName  ()
        stitle = self.gamma.GetTitle ()
        gname  = sname .replace ( 'sigma' , 'gamma_1430' )
        gtitle = stitle.replace ( 'sigma' , 'Gamma_1440' )
        self.gamma.SetName  ( gname  ) 
        self.gamma.SetTitle ( gtitle )

        self.g0_1430 = self.gamma 
        
        self.m0_1430 = self.mean 
        sname  = self.mean.GetName  ()
        stitle = self.mean.GetTitle ()
        gname  = sname .replace ( 'mean' , 'm0_1430' )
        gtitle = stitle.replace ( 'mean' , 'm0_1440' )
        self.m0_1430.SetName  ( gname  ) 
        self.m0_1430.SetTitle ( gtitle ) 
        
        self.a_lass = makeVar ( a_lass             ,
                                'aLASS_%s'  % name ,
                                "aLASS(%s)" % name , a_lass , 
                                1.94e-3            ,
                                1.94e-3            ,
                                1.94e-3            ) 
        self.r_lass = makeVar ( r_lass             ,
                                'rLASS_%s'  % name ,
                                "rLASS(%s)" % name , r_lass , 
                                1.76e-3            ,
                                1.76e-3            ,
                                1.76e-3            ) 
        self.e_lass = makeVar ( e_lass             ,
                                'eLASS_%s'  % name ,
                                "eLASS(%s)" % name , e_lass ,
                                1.0                , 
                                1.0                ,
                                1.0                )

        self.mKaon = mKaon
        self.mPion = mPion
        
        ## create PDF 
        self.pdf = Ostap.Models.LASS ( 
            "lass_"    + name ,
            "LASS(%s)" % name ,
            self.mass    ,
            self.m0_1430 ,
            self.g0_1430 ,
            self.a_lass  ,
            self.r_lass  ,
            self.e_lass  ,
            self.mKaon   ,
            self.mPion   ) 

models.append ( LASS_pdf )                          
# =============================================================================
## @class Bugg_pdf
#  The parameterization of sigma pole by B.S.Zou and D.V.Bugg, Phys.Rev. D48 (1993) R3948.
#  @see http://dx.doi.org/10.1103/PhysRevD.48.R3948
#  @see Ostap::Models::Bugg
#  @see Ostap::Math::Bugg
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Bugg_pdf(MASS) :
    """Sigma pole
    The parameterization of sigma pole by B.S.Zou and D.V.Bugg, Phys.Rev. D48 (1993) R3948.
    http://dx.doi.org/10.1103/PhysRevD.48.R3948
    """
    def __init__ ( self              ,
                   name              ,
                   mn       = None   ,
                   mx       = None   ,
                   mass     = None   ,
                   bugg_m   = 0.9264 , 
                   bugg_g2  = 0.0024 , ## g2-parameter
                   bugg_b1  = 0.5848 , ## b1-parameter [GeV]
                   bugg_b2  = 1.6663 , ## b2-parameter [GeV^-1]
                   bugg_a   = 1.082  , ##  a-parameter [GeV^2]
                   bugg_s1  = 2.8    , ## s1-parameter [GeV^2]
                   bugg_s2  = 3.5    , ## s2-parameter
                   mPion    = 0.1396 ) :  ## pion mass 
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name     ,
                         mn      , mx       ,
                         mass    ,
                         bugg_m  , bugg_g2  ) 
        
        self.bugg_g2 = self.sigma
        self.gamam   = self.sigma
        ##
        sname  = self.gamma.GetName  ()
        stitle = self.gamma.GetTitle ()
        gname  = sname .replace ( 'sigma' , 'gamma2_Bugg' )
        gtitle = stitle.replace ( 'sigma' , 'Gamma2_Bugg' )
        self.bugg_g2.SetName  ( gname  ) 
        self.bugg_g2.SetTitle ( gtitle )
        self.gamma = self.bugg_g2 
        
        self.bugg_m = self.mean 
        sname  = self.mean.GetName  ()
        stitle = self.mean.GetTitle ()
        gname  = sname .replace ( 'mean' , 'm_Bugg' )
        gtitle = stitle.replace ( 'mean' , 'm_Bugg' )
        self.bugg_m.SetName  ( gname  ) 
        self.bugg_m.SetTitle ( gtitle ) 
        
        self.bugg_b1 = makeVar ( bugg_b1             ,
                                 'b1Bugg_%s'  % name ,
                                 "b1Bugg(%s)" % name , bugg_b1 ,
                                 0.5848 , 
                                 0 , 2  )
        
        self.bugg_b2 = makeVar ( bugg_b2             ,
                                 'b2Bugg_%s'  % name ,
                                 "b2Bugg(%s)" % name , bugg_b2 , 
                                 1.6663 ,
                                 1 , 2  ) 
        
        self.bugg_a  = makeVar ( bugg_a             ,
                                 'aBugg_%s'  % name ,
                                 "aBugg(%s)" % name , bugg_a , 
                                 1.082    ,
                                 0.5 , 5  ) 

        self.bugg_s1  = makeVar ( bugg_s1           ,
                                  's1Bugg_%s'  % name ,
                                  "s1Bugg(%s)" % name , bugg_s1 , 
                                  2.8              ,
                                  1 , 5            ) 
        
        self.bugg_s2  = makeVar ( bugg_s2           ,
                                  's2Bugg_%s'  % name ,
                                  "s2Bugg(%s)" % name , bugg_s2 , 
                                  3.5              ,
                                  1 , 5            ) 

        self.mPion = mPion
        
        ## create PDF 
        self.pdf = Ostap.Models.Bugg ( 
            "bugg_"    + name ,
            "Bugg(%s)" % name ,
            self.mass    ,
            self.bugg_m  ,
            self.bugg_g2 ,
            self.bugg_b1 ,
            self.bugg_b2 ,
            self.bugg_a  ,
            self.bugg_s1 ,
            self.bugg_s2 ,
            self.mPion   ) 
        
models.append ( Bugg_pdf )

# =============================================================================
## @class Swanson_pdf
#  S-wave cusp
#  @see LHCb-PAPER-2016-019, Appendix
#  @see E. S. Swanson, Cusps and exotic charmonia, arXiv:1504.07952
#  @see http://arxiv.org/abs/1504.07952
#  @see Ostap::Models::Swanson
#  @see Ostap::Math::Swanson
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2016-06-12
class Swanson_pdf(MASS) :
    """ S-wave cusp
    - LHCb-PAPER-2016-019, Appendix
    - E. S. Swanson, Cusps and exotic charmonia, arXiv:1504.07952
    - http://arxiv.org/abs/1504.07952
    """
    def __init__ ( self              ,
                   name              ,
                   swanson           , ## Ostap::Math::Swanson objects 
                   mn       = None   ,
                   mx       = None   ,
                   mass     = None   ,
                   beta0    = None   ) : 
        #
        ## initialize the base
        # 
        MASS.__init__  ( self    , name    ,
                         mn      , mx      ,
                         mass    ,
                         None    , None    ) 

        del self.mean 
        del self.sigma

        self.swanson = swanson 
        beta_max     = max ( swanson.mmin() , swanson.cusp() )
        self.beta0   = makeVar ( beta0 , 
                                 'b0_swanson_%s'   % name ,
                                 'b0_swanson(%s)'  % name ,
                                 beta0 , 
                                 0 , beta_max )
        ## create PDF 
        self.pdf = Ostap.Models.Swanson ( 
            "Swanson_"    + name ,
            "Swanson(%s)" % name ,
            self.mass    ,
            self.beta0   ,
            self.swanson )

models.append ( Swanson_pdf )

# =============================================================================
if '__main__' == __name__ : 
    
    from ostap.logger.line import line 
    logger.info ( __file__  + '\n' + Ostap.Line.line  ) 
    logger.info ( 80*'*'   )
    logger.info ( __doc__  )
    logger.info ( 80*'*' )
    logger.info ( ' Author  : %s' %         __author__    ) 
    logger.info ( ' Version : %s' %         __version__   ) 
    logger.info ( ' Date    : %s' %         __date__      )
    logger.info ( ' Symbols : %s' %  list ( __all__     ) )
    logger.info ( 80*'*' )
    for m in models : logger.info ( 'Model %s: %s' % ( m.__name__ ,  m.__doc__  ) ) 
    logger.info ( 80*'*' ) 
    
# =============================================================================
# The END 
# =============================================================================
