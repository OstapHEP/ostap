#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
  - double     Gauissian
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
    'DoubleGauss_pdf'        , ## double Gauss
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
                   mass             ,
                   mean      = None ,
                   sigma     = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , sigma  )
        
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
                   mass             ,
                   mean     = None  , 
                   sigma    = None  ,
                   alpha    = None  ,
                   n        = None  ) : 
                   
        #
        ## initialize the base
        #
        MASS.__init__ ( self  , name , mass  , mean , sigma    )
        
        self.__alpha = makeVar ( alpha ,
                                 'alpha_%s'        % name ,
                                 '#alpha_{CB}(%s)' % name ,  alpha  ,
                                 2.0 , 0  , 10 )        
        self.__n     = makeVar ( n   ,
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
        
    @property
    def alpha ( self ) :
        """Alpha-parameter for Crystal Ball tail"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        self.__alpha.setVal ( value )
        return self.__alpha.getVal ()
    
    @property
    def n ( self ) :
        """N-parameter for Crystal Ball tail"""
        return self.__n
    @n.setter
    def n ( self, value ) :
        self.__n.setVal ( value )
        return self.__n.getVal ()

    
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
                   mass     = None   ,
                   mean     = None   , 
                   sigma    = None   ,
                   alpha    = None   ,
                   n        = None   ) : 
                   
        
        MASS.__init__ ( self  , name  , mass , mean  , sigma )
        
        self.__alpha = makeVar ( alpha ,
                                 'alpha_%s'          % name ,
                                 '#alpha_{CBRS}(%s)' % name , alpha , 
                                 2.0 , 0  , 10      )        
        self.__n     = makeVar ( n     ,
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

    @property
    def alpha ( self ) :
        """Alpha-parameter for Crystal Ball tail"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        self.__alpha.setVal ( value )
        return self.__alpha.getVal ()
    
    @property
    def n ( self ) :
        """N-parameter for Crystal Ball tail"""
        return self.__n
    @n.setter
    def n ( self, value ) :
        self.__n.setVal ( value )
        return self.__n.getVal ()

        
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
                   mass              , 
                   mean      = None  ,
                   sigma     = None  ,
                   alphaL    = None  ,
                   alphaR    = None  ,
                   nL        = None  ,
                   nR        = None  ) : 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , sigma )
        #
        ## treat the specific parameters
        #
        self.__aL    = makeVar ( alphaL                  ,
                                 "aL_%s"          % name ,
                                 "#alpha_{L}(%s)" % name , alphaL    , 2.0 , 0 , 10 )
        self.__nL    = makeVar ( nL                      ,                     
                                 "nL_%s"          % name ,
                                 "n_{L}(%s)"      % name , nL        , 1   , 0 , 20 )
        self.__aR    = makeVar ( alphaR ,
                                 "aR_%s"          % name ,
                                 "#alpha_{R}(%s)" % name , alphaR    , 2.0 , 0 , 10 )
        self.__nR    = makeVar ( nR                      ,
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

    @property
    def aL ( self ) :
        """(left) Alpha-parameter for Crystal Ball tail"""
        return self.__aL
    @aL.setter
    def aL ( self, value ) :
        self.__aL.setVal ( value )
        return self.__aL.getVal ()

    @property
    def nL ( self ) :
        """(left) N-parameter for Crystal Ball tail"""
        return self.__nL
    @nL.setter
    def nL ( self, value ) :
        self.__nL.setVal ( value )
        return self.__nL.getVal ()

    @property
    def aR ( self ) :
        """(right) Alpha-parameter for Crystal Ball tail"""
        return self.__aR
    @aR.setter
    def aR ( self, value ) :
        self.__aR.setVal ( value )
        return self.__aR.getVal ()

    @property
    def nR ( self ) :
        """(right) N-parameter for Crystal Ball tail"""
        return self.__nR
    @nR.setter
    def nR ( self, value ) :
        self.__nR.setVal ( value )
        return self.__nR.getVal ()

        
models.append ( CB2_pdf )    
# =============================================================================
## @class Needham_pdf
#  Needham function: specific parameterisation of Crystal Ball function with
#   - \f$ n = 1 \f$  
#   - \f$ \alpha(\sigma) = a_0 + \sigma\times (a_1+\sigma \times a_2) \f$ 
#  The function is very well sutable to fit
#  \f$J/\psi \rightarrow \mu^+\mu^-\f$,
#  \f$\psi^{\prime} \rightarrow \mu^+\mu^-\f$ and
#  \f$\Upsilon \rightarrow \mu^+\mu^-\f$ signals and
#  is has been used with great success for all LCHb papers
#  on quarkonia production in dimuon final states
#  -thanks to  Matthew Needham
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
                   mass                 ,
                   mean     = 3.096     ,   ## GeV  
                   sigma    = 0.013     ,   ## GeV 
                   a0       = 1.975     ,
                   a1       = -0.0011   ,   ## GeV^-1
                   a2       = -0.00018  ) : ## GeV^-2  
        
        MASS.__init__ ( self , name  , mass , mean , sigma )
        
        #
        unit = 1000
        #
        if   3.096 in self.mass : unit = 1000 
        elif 3096  in self.mass : unit = 1
        elif 9.460 in self.mass : unit = 1000 
        elif 9460  in self.mass : unit = 1
        #
        self.__a0 = makeVar ( a0                  ,
                              "a0_%s"     % name  ,
                              "a_{0}(%s)" % name  , a0 , 
                              1.975               ,   0           , 10           )
        self.__a1 = makeVar ( a1                  ,
                              "a1_%s"     % name  ,
                              "a_{1}(%s)" % name  , a1 , 
                              -0.0011   * unit    , -10 * unit    , 10 * unit    )
        self.__a2 = makeVar ( a2                  ,
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
        
    @property
    def a0 ( self ) :
        """A0-parameter for Needham function"""
        return self.__a0
    @a0.setter
    def a0 ( self, value ) :
        self.__a0.setVal ( value )
        return self.__a0.getVal ()

    @property
    def a1 ( self ) :
        """A1-parameter for Needham function"""
        return self.__a1
    @a1.setter
    def a1 ( self, value ) :
        self.__a1.setVal ( value )
        return self.__a1.getVal ()

    @property
    def a2 ( self ) :
        """A2-parameter for Needham function"""
        return self.__a2
    @a2.setter
    def a2 ( self, value ) :
        self.__a2.setVal ( value )
        return self.__a2.getVal ()

        
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
                   mass                    ,
                   mean      = None        ,
                   sigma     = None        ,
                   alpha     = None        ,
                   n         = None        ,
                   b         = None        ) : 
                   
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , sigma ) 
        
        self.__alpha = makeVar ( alpha                     ,
                                 'alpha_%s'         % name ,
                                 '#alpha_{Apo}(%s)' % name , alpha , 
                                 2.0   , 0 , 10 )
        
        self.__n     = makeVar ( n                    ,
                                 'n_%s'        % name ,
                                 'n_{Apo}(%s)' % name , n ,
                                 2.0   , 0 , 20 )
        
        self.__b     = makeVar ( b                    ,
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

    @property
    def alpha ( self ) :
        """Alpha-parameter for Apolonious tail"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        self.__alpha.setVal ( value )
        return self.__alpha.getVal ()
    
    @property
    def n ( self ) :
        """N-parameter for Apolonios tail"""
        return self.__n
    @n.setter
    def n ( self, value ) :
        self.__n.setVal ( value )
        return self.__n.getVal ()

    @property
    def b ( self ) :
        """B-parameter for Apolonios function"""
        return self.__b
    @b.setter
    def b ( self, value ) :
        self.__b.setVal ( value )
        return self.__b.getVal ()


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
#          \frac{x-\mu}{\sigma_r} & \text{for} & x \ge \mu \\
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
                   mass               ,
                   mean      = None   ,
                   sigma     = None   ,
                   asymmetry = None   ,
                   beta      = None   ) : 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , sigma  )
        
        self.__asym = makeVar ( asymmetry                  ,
                                'asym_%s'           % name ,
                                '#kappa_{Apo2}(%s)' % name ,
                                asymmetry , -1 , 1  ) 
        
        self._lst_R   = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.__sigmaR = ROOT.RooFormulaVar (
            "sigmaR_%s"     % name   ,
            "sigma_{R}(%s)" % name   ,
            "%s*(1-%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self._lst_R   )
        
        self._lst_L   = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.__sigmaL = ROOT.RooFormulaVar (
            "sigmaL_%s"     % name   ,
            "sigma_{L}(%s)" % name   ,
            "%s*(1+%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self._lst_L   )
        
        self.__beta    = makeVar ( beta ,
                                   'beta_%s'          % name  ,
                                   '#beta_{Apo2}(%s)' % name  ,
                                   beta , 0.001  , 10000 ) 
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

    @property
    def asym ( self ) :
        """Asymmetry-parameter for Apolonious-2 function"""
        return self.__asym
    @asym.setter
    def asym ( self, value ) :
        value = float ( value ) 
        if  not -1 <= value <= 1 :
            raise AttributeError('Asymmetry parameter is out of range -1<%s<1, adjust it!' % value )
        self.__asym.setVal ( value )
        return self.__asym.getVal ()

    @property
    def beta ( self ) :
        """Beta-parameter for Apolonious-2 function"""
        return self.__beta
    @beta.setter
    def beta ( self, value ) :
        self.__beta.setVal ( value )
        return self.__beta.getVal ()

    @property
    def sigmaL ( self ) :
        """(left) sigma-parameter for Apolonious-2 function"""
        return self.__sigmaL
    
    @property
    def sigmaR ( self ) :
        """(right) sigma-parameter for Apolonious-2 function"""
        return self.__sigmaR

    


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
                   mass                  ,
                   mean      = None      ,
                   sigma     = None      ,
                   asymmetry = None      ) : 
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , sigma ) 
        
        self.__asym = makeVar ( asymmetry               ,
                                'asym_%s'        % name ,
                                '#xi_{asym}(%s)' % name ,
                                asymmetry , -1 , 1  ) 
        
        self._lst_R   = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.__sigmaR = ROOT.RooFormulaVar (
            "sigmaR_%s"     % name   ,
            "sigma_{R}(%s)" % name   ,
            "%s*(1+%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self._lst_R   )
        
        self._lst_L   = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.__sigmaL = ROOT.RooFormulaVar (
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

    @property
    def asym ( self ) :
        """Asymmetry-parameter for Bifurcated Gaussian"""
        return self.__asym
    @asym.setter
    def asym ( self, value ) :
        value = float ( value ) 
        if  not -1 <= value <= 1 :
            logger.warning('Asymmetry parameter is out of range -1<%s<1, adjust it!' % value )
            value = max  ( -1 , min ( value , 1.0 ) )  
        self.__asym.setVal ( value )
        return self.__asym.getVal ()

    @property
    def sigmaL ( self ) :
        """(left) sigma-parameter for Bifurcated Gaussian"""
        return self.__sigmaL
    
    @property
    def sigmaR ( self ) :
        """(right) sigma-parameter for Bifurcated Gaussian"""
        return self.__sigmaR



models.append ( BifurcatedGauss_pdf )    


# =============================================================================
## @class DoubleGauss_pdf
#  simple wrapper over double gaussian
#  @see RooGaussian
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-07-12
class DoubleGauss_pdf(MASS) :
    """Double Gauss :

    f(x;mu,sigma,scale,fraction) =
    fraction     * G(x;mu,sigma) +
    (1-fraction) * G(x;mu;sigma*scale) 
        
    """
    def __init__ ( self                  ,
                   name                  ,
                   mass                  ,
                   mean      = None      ,
                   sigma     = None      ,
                   fraction  = None      ,
                   scale     = None      ) :  
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , sigma )
        
        self.__scale = makeVar (
            scale ,
            'SigmaScale'        + name ,
            'SigmaScale(%s)'    % name , scale , 1 , 10 ) 
        
        ## the fraction 
        self.__fraction = makeVar (
            fraction                   , 
            'CoreFraction'      + name ,
            'CoreFraction(%s)'  % name , fraction , 0 , 1 ) 
        
        self.pdf = Ostap.Models.DoubleGauss (           
            "f2gau_"          + name ,
            "DoubleGauss(%s)" % name ,
            self.mass     ,
            self.sigma    ,
            self.fraction ,
            self.scale    ,
            self.mean    
            )

    @property
    def scale ( self ) :
        """Scale-parameter for doble Gaussian"""
        return self.__scale
    @scale.setter
    def scale ( self, value ) :
        value  = float ( value )
        if value  <= 1 : raise AttributeError ( "Scale parameter must be >1") 
        self.__scale.setVal ( abs ( value ) ) 
        return self.__scale.getVal ()

    @property
    def fraction ( self ) :
        """Fraction-parameter for Bifurcated Gaussian"""
        return self.__fraction
    @fraction.setter
    def fraction ( self, value ) :
        value = float ( value ) 
        if  not 0  <= value <= 1 : raise AttributeError ( "Fraction parameter must be 0<=f<=1") 
        self.__fraction.setVal ( value )
        return self.__fraction.getVal ()

models.append ( DoubleGauss_pdf )    
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
                   mass             ,
                   mean      = None ,
                   alpha     = None ,
                   beta      = 2    ) :  ## beta=2 is gaussian distribution 
        
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self  , name , mass , mean  , alpha ) 
        
        #
        ## rename it!
        #
        self.__alpha = self.sigma
        sname  = self.alpha.GetName  ()
        stitle = self.alpha.GetTitle ()
        gname  = sname .replace ( 'sigma' , 'alpha' )
        gtitle = stitle.replace ( 'sigma' , 'alpha' )
        self.alpha.SetName  ( gname  ) 
        self.alpha.SetTitle ( gtitle )
        
        self.__beta  = makeVar ( beta ,
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

    @property
    def alpha ( self ) :
        """alpha-parameter for Generalized Gaussian"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        value  = float ( value )
        if value  <= 0 : raise AttributeError ( "Alpha parameter must be positive!") 
        self.__alpha.setVal ( abs ( value ) ) 
        return self.__alpha.getVal ()

    @property
    def beta ( self ) :
        """beta-parameter for Generalized  Gaussian"""
        return self.__beta
    @beta.setter
    def beta ( self, value ) :
        value  = float ( value )
        if value  <= 0 : raise AttributeError ( "Beta parameter must be positive!") 
        self.__beta.setVal ( abs ( value ) ) 
        return self.__beta.getVal ()
        

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
                   mass               ,
                   mean        = None ,
                   alpha       = None ,
                   kappa       = 0    ) : ## 0 corresponds to gaussian distribution 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , alpha ) 
        
        #
        ## rename it!
        #
        self.__alpha = self.sigma        
        sname      = self.alpha.GetName  ()
        stitle     = self.alpha.GetTitle ()
        gname      = sname .replace ( 'sigma' , 'alpha' )
        gtitle     = stitle.replace ( 'sigma' , 'alpha' )
        self.alpha.SetName  ( gname  ) 
        self.alpha.SetTitle ( gtitle )
        
        self.__xi   = self.mean 
        ## sname     = self.xi.GetName  ()
        ## stitle    = self.xi.GetTitle ()
        ## gname     = sname .replace   ( 'mean' , 'xi' )
        ## gtitle    = stitle.replace   ( 'mean' , 'xi' )
        ## self.xi.SetName              ( gname  ) 
        ## self.xi.SetTitle             ( gtitle )
        
        self.__kappa = makeVar ( kappa ,
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

    @property
    def alpha ( self ) :
        """alpha-parameter for Generalized Gaussian"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        value  = float ( value )
        if value  <= 0 : raise AttributeError ( "Alpha parameter must be positive!") 
        self.__alpha.setVal ( abs ( value ) ) 
        return self.__alpha.getVal ()

    @property
    def kappa ( self ) :
        """kappa-parameter for Generalized Gaussian"""
        return self.__kappa
    @kappa.setter
    def kappa ( self, value ) :
        self.__kappa.setVal ( value ) 
        return self.__kappa.getVal ()
    
    @property
    def xi ( self ) :
        """xi-parameter (location) for Generalized Gaussian"""
        return self.__xi
    


models.append ( GenGaussV2_pdf )    
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
                   mass            , 
                   mean     = None ,
                   sigma    = None ,
                   xi       = None ,
                   rhoL     = None ,
                   rhoR     = None ) :

        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , sigma )
        #
        ## treat the specific parameters
        #
        # asymmetry 
        self.__xi    = makeVar ( xi                    ,
                                 "xi_%s"        % name ,
                                 "#xi(%s)"      % name , xi      , 0  , -1 , 1    )        
        self.__rhoL  = makeVar ( rhoL                  ,
                                 "rhol_%s"      % name ,
                                 "#rho_{L}(%s)" % name , rhoL    , 0  ,  -10 , 10 )        
        self.__rhoR  = makeVar ( rhoR                  ,
                                 "rhor_%s"      % name ,
                                 "#rho_{R}(%s)" % name , rhoR    , 0  ,  -10 , 10 )
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
            self.rhoL  ,
            self.rhoR  )

    @property
    def xi ( self ) :
        """xi-parameter (asymmetry) for Bukin function"""
        return self.__xi
    @xi.setter
    def xi ( self, value ) :
        value  = float ( value )
        if not -1 <= value  <= 1 : raise AttributeError ( "Asymmetry must be in range  -1<=xi<=1!") 
        self.__xi.setVal ( value ) 
        return self.__xi.getVal ()

    @property
    def rhoL ( self ) :
        """Rho-parameter (left tail) for Bukin function"""
        return self.__rhoL
    @rhoL.setter
    def rhoL ( self, value ) :
        value  = float ( value )
        if value < 0 : raise AttributeError ( "Rho-parameter must be non-negative!") 
        self.__rhoL.setVal ( value ) 
        return self.__rhoL.getVal ()

    @property
    def rhoR ( self ) :
        """Rho-parameter (right tail) for Bukin function"""
        return self.__rhoR
    @rhoR.setter
    def rhoR ( self, value ) :
        value  = float ( value )
        if value < 0 : raise AttributeError ( "Rho-parameter must be non-negative!") 
        self.__rhoR.setVal ( value ) 
        return self.__rhoR.getVal ()

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
                   mass             ,
                   mean      = None ,
                   sigma     = None ,
                   n         = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , sigma ) 
        
        # 
        self.__n  = makeVar ( n                    ,
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

    @property
    def n ( self ) :
        """N-parameter for Student-T function"""
        return self.__n
    @n.setter
    def n ( self, value ) :
        value  = float ( value )
        if value <= 0 : raise AttributeError ( "N-parameter must be positive!") 
        self.__n.setVal ( value ) 
        return self.__n.getVal ()

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
    for x > mu   n=n_r and dx = ( x - mu )  / sigma_r
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
                   mass             ,
                   mean      = None ,
                   sigma     = None ,
                   asymmetry = None ,
                   nL        = None , 
                   nR        = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , sigma )
        
        self.__asym = makeVar ( asymmetry                  ,
                                'asym_%s'        % name ,
                                '#xi_{asym}(%s)' % name , asymmetry , 
                                0 , -1 , 1  ) 
        
        self._lst_R   = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.__sigmaR = ROOT.RooFormulaVar (
            "sigmaR_stt_%s"     % name   ,
            "sigma_{R}(%s)" % name   ,
            "%s*(1+%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self._lst_R   )
        
        self._lst_L   = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.__sigmaL = ROOT.RooFormulaVar (
            "sigmaL_stt_%s"     % name   ,
            "sigma_{L}(%s)" % name   ,
            "%s*(1-%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self._lst_L   )
        
        # 
        self.__nL =  makeVar ( nL                     ,
                               'nL_%s'         % name ,
                               '#nL_{BST}(%s)' % name , nL , 
                               2  , 0 , 100  )
        
        self.__nR =  makeVar ( nR                    ,
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

    @property
    def nL ( self ) :
        """N-parameter (left) for Bifurcated Student-T function"""
        return self.__nL
    @nL.setter
    def nL ( self, value ) :
        value  = float ( value )
        if value <= 0 : raise AttributeError ( "N-parameter must be positive!") 
        self.__nL.setVal ( value ) 
        return self.__nL.getVal ()

    @property
    def nR ( self ) :
        """N-parameter (right) for Bifurcated Student-T function"""
        return self.__nR
    @nR.setter
    def nR ( self, value ) :
        value  = float ( value )
        if value <= 0 : raise AttributeError ( "N-parameter must be positive!") 
        self.__nR.setVal ( value ) 
        return self.__nR.getVal ()
        
    @property
    def asym ( self ) :
        """Asymmetry -parameter for Bifurcated Student-T function"""
        return self.__asym
    @asym.setter
    def asym ( self, value ) :
        value  = float ( value )
        if not -1 <= value <= 1 : raise AttributeError ( "Asymmetry must be inn range -1<=a<=1!") 
        self.__asym.setVal ( value ) 
        return self.__asym.getVal ()

    @property
    def sigmaL( self ) :
        """Sigma-parameter (left) for Bifurcated Student-T function"""
        return self.__sigmaL

    @property
    def sigmaR( self ) :
        """Sigma-parameter (right) for Bifurcated Student-T function"""
        return self.__sigmaR


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
                   mass             ,
                   mean      = None , ## mu 
                   sigma     = None ,
                   epsilon   = 0    ,
                   delta     = 1    ) :

        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , sigma )

        self.__mu      = self.mean
        self.__epsilon = makeVar ( epsilon ,
                                   'epsilon_%s'   % name ,
                                   '#epsilon(%s)' % name , epsilon ,
                                   0 , -1000 , +1000 )
        self.__delta   = makeVar ( delta ,
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

    @property
    def epsilon( self ) :
        """Epsilon-parameter for Sinh-Asinh function"""
        return self.__epsilon
    @epsilon.setter
    def epsilon ( self, value ) :
        self.__epsilon.setVal ( value ) 
        return self.__epsilon.getVal ()

    @property
    def delta ( self ) :
        """Delta-parameter for Sinh-Asinh function"""
        return self.__delta
    @delta.setter
    def delta ( self, value ) :
        value = float ( value )
        if value <= 0 : raise AttributeError ( "Delta must be positive")         
        self.__delta.setVal ( value ) 
        return self.__delta.getVal ()
        
    @property
    def mu ( self ) :
        """Mu-parameter (location) for Sinh-Asinh function"""
        return self.__mu
        

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
                   mass             ,
                   xi        = None ,   ## related to mean 
                   lam       = None ,   ## related to sigma 
                   delta     = 1    ,
                   gamma     = 0    ) :  

        #
        ## initialize the base
        #
        
        MASS.__init__  ( self , name , mass , xi , lam  ) ## mean    , sigma  )

        self.__xi = self.mean
        if None == xi : ## newly created 
            oname   = self.xi.GetName ()
            otitle  = self.xi.GetTitle()
            nname   = oname .replace ( 'mean' , 'xi' )
            ntitle  = otitle.replace ( 'mean' , 'xi' )
            self.xi.SetName  ( nname  )
            self.xi.SetTitle ( ntitle )
            
        self.__lam = self.sigma
        if None == lam : ## newly created 
            oname    = self.lambd.GetName ()
            otitle   = self.lambd.GetTitle()
            nname    = oname .replace ( 'sigma' , 'lambda' )
            ntitle   = otitle.replace ( 'sigma' , 'lambda' )
            self.lambd.SetName  ( nname  )
            self.lambd.SetTitle ( ntitle )
            self.lambd.setMax ( self.lambd.getMax() * 10 ) ## adjust it! 


        ## provide backup name 
        self.lambda_ = self.lam
        
        self.__delta   = makeVar ( delta                 ,
                                   'delta_%s'     % name ,
                                   '#delta(%s)'   % name , delta ,
                                   1 , 1.e-6 , 1000   )
        
        self.__gamma   = makeVar ( gamma               ,
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
            self.lambd     ,
            self.delta     ,
            self.gamma     ) 

    @property
    def delta ( self ) :
        """Delta-parameter for Johnson-SU function"""
        return self.__delta
    @delta.setter
    def delta ( self, value ) :
        value = float ( value )
        if value <= 0 : raise AttributeError ( "Delta must be positive")         
        self.__delta.setVal ( value ) 
        return self.__delta.getVal ()

    @property
    def gamma ( self ) :
        """Gamma-parameter for Johnson-SU function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        self.__gamma.setVal ( value ) 
        return self.__gamma.getVal ()
        
    @property
    def xi ( self ) :
        """Xi-parameter (location) for Johnson-SU function"""
        return self.__xi
    @xi.setter
    def xi ( self, value ) :
        value = float ( value )
        self.__xi.setVal ( value ) 
        return self.__xi.getVal ()

    @property
    def lambd ( self ) :
        """Lambda-parameter (location) for Johnson-SU function"""
        return self.__lam
    @lambd.setter
    def lambd ( self, value ) :
        value = float ( value )
        self.__lambd.setVal ( value ) 
        return self.__lambd.getVal ()
        

models.append ( JohnsonSU_pdf )      
# =============================================================================
## @class Atlas_pdf
#  Modified gaussian with exponential tails
#  \f[ f(x) \propto \exp( -frac{\delta x^{1+\frac{1}{1+\delta x/2}}}{2})\f],
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
    r"""Modified gaussian with exponential tails
    \f$  f(x) \propto \exp( -frac{\delta x^{1+\frac{1}{1+\delta x/2}}}{2})\f$,
    where \f$\delta x = \left| x - \mu \right|/\sigma\f$
    Function is taken from http://arxiv.org/abs/arXiv:1507.07099    
    """
    def __init__ ( self             ,
                   name             ,
                   mass             ,
                   mean      = None ,   ## related to mean 
                   sigma     = None ) : ## related to sigma 

        #
        ## initialize the base
        #
        
        MASS.__init__  ( self , name  , mass , mean , sigma  )

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
#  \f[ f(x,\mu,\sigma) \propto \frac{1}{2} \mathrm{sech} ( \frac{\pi}{2}\frac{x-\mu}{\sigma} )\f] 
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
                   mass             ,
                   mean      = None ,   ## related to mean 
                   sigma     = None ) : ## related to sigma 

        #
        ## initialize the base
        #
        
        MASS.__init__  ( self , name , mass , mean , sigma  )

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
    r""" Logistic, aka ``sech-square'' PDF
    \f$ f(x;\mu;s) = \frac{1}{4s}sech^2\left(\frac{x-\mu}{2s}\right)\f$, 
    where
    \f$  s = \sigma \frac{\sqrt{3}}{\pi}\f$
    - see https://en.wikipedia.org/wiki/Logistic_distribution
    - see Ostap::Math::Logistic
    - see Ostap::Models::Logistic
    """
    def __init__ ( self             ,
                   name             ,
                   mass             ,
                   mean      = None ,   ## related to mean 
                   sigma     = None ) : ## related to sigma 

        #
        ## initialize the base
        #
        
        MASS.__init__  ( self , name , mass , mean , sigma  )

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
                   mass             ,
                   mean      = None ,
                   sigma     = None ,
                   gamma     = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , sigma ) 

        dm = self.mass.getmax() -  self.mass.getMin()
        self.__gamma  = makeVar ( gamma               ,
                                  'gamma_%s'   % name ,   
                                  '#gamma(%s)' % name , gamma , 
                                  1.e-5*dm ,  0.3 * dm )
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

    @property
    def gamma ( self ) :
        """Gamma-parameter for Voigt function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        if value <= 0 : raise AttributeError ( "Gamma must be positive")         
        self.__gamma.setVal ( value ) 
        return self.__gamma.getVal ()


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
class PseudoVoigt_pdf(Voigt_pdf) :
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
                   mass             ,
                   mean      = None ,
                   sigma     = None ,
                   gamma     = None ) :
        #
        ## initialize the base
        # 
        Voigt_pdf.__init__  ( self , name , mass , mean , sigma , gamma ) 

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
                   mass               ,
                   mean        = None , 
                   gamma       = None ,
                   convolution = None ,
                   useFFT      = True ) :
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , gamma )
        
        self.__gamma = self.sigma
        sname  = self.gamma.GetName  ()
        stitle = self.gamma.GetTitle ()
        gname  = sname .replace ( 'sigma' , 'gamma' )
        gtitle = stitle.replace ( 'sigma' , 'Gamma' )
        self.gamma.SetName  ( gname  ) 
        self.gamma.SetTitle ( gtitle )
        
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
            
    @property
    def gamma ( self ) :
        """Gamma-parameter for Breit-Wigner function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        if value <= 0 : raise AttributeError ( "Gamma must be positive")         
        self.__gamma.setVal ( value ) 
        return self.__gamma.getVal ()


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
                   mass               ,
                   mean        = None , 
                   gamma       = None ) : 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , mean , gamma )
        
        self.__gamma = self.sigma
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
        
    @property
    def gamma ( self ) :
        """Gamma-parameter for Breit-Wigner function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        if value <= 0 : raise AttributeError ( "Gamma must be positive")         
        self.__gamma.setVal ( value ) 
        return self.__gamma.getVal ()


models.append ( BW23L_pdf )                          
# =============================================================================
## @class Flatte_pdf
#  Flatte function
#  S.M.Flatte, "Coupled-channel analysis of the \f$\pi\eta\f$ 
#    and \f$K\bar{K}\f$ systems near \f$K\bar{K}\f$ threshold  
#    Phys. Lett. B63, 224 (1976)
#  Well suitable for \f$f_0(980)\rightarrow \pi^+ \pi^-\f$
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
                   mass              ,
                   m0_980   = None   ,    ## mass  of f0(980) resonance
                   m0g1     = 165000 ,    ## m0(f0(980))*gamma_1 
                   g2og1    = 4.21   ) :  ## gamma2/gamma1 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , m0_980 , None )
        
        self.flatte = flatte 

        self.__m0_980 = self.mean 
        sname  = self.mean.GetName  ()
        stitle = self.mean.GetTitle ()
        gname  = sname .replace ( 'mean' , 'm0_980' )
        gtitle = stitle.replace ( 'mean' , 'm0_980' )
        self.m0_980.SetName  ( gname  ) 
        self.m0_980.SetTitle ( gtitle ) 
        
        self.__m0g1 = makeVar  ( m0g1                          ,
                                 'm0g1_%s'              % name ,
                                 'm_{0}*\gamma_{1}(%s)' % name , m0g1 ,
                                 165                           ,
                                 1.e-5                         ,
                                 1.e+5                         )
        
        self.__g2og1 = makeVar ( g2og1    ,
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
        
    @property
    def m0_980 ( self ) :
        """m0-parameter for Flatte-function"""
        return self.__m0_980
    @m0_980.setter
    def m0_980 ( self, value ) :
        value = float ( value )
        if value <= 0 : raise AttributeError ( "Gamma must be positive")         
        self.__m0_980.setVal ( value ) 
        return self.__m0_980.getVal ()

    @property
    def m0g1 ( self ) :
        """m0*g1-parameter for Flatte-function"""
        return self.__m0g1
    @m0_980.setter
    def m0g1 ( self, value ) :
        value = float ( value )
        self.__m0g1.setVal ( value ) 
        return self.__m0g1.getVal ()

    @property
    def g2og1 ( self ) :
        """g2/g1-parameter for Flatte-function"""
        return self.__g2og1
    @g2og1.setter
    def g2og1 ( self, value ) :
        value = float ( value )
        self.__g2og1.setVal ( value ) 
        return self.__g2og1.getVal ()


models.append ( Flatte_pdf )                          
# =============================================================================
## @class Flatte2_pdf
#  Flatte function to describe dikaon system near threshold
#  S.M.Flatte, 
#    "Coupled-channel analysis of the \f$\pi\eta\f$ 
#    and \f$K\bar{K}\f$ systems near \f$K\bar{K}\f$ threshold  
#    Phys. Lett. B63, 224 (1976)
#  Well suitable for \f$f_0(980)\rightarrow K^+ K^-\f$
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
                   mass              ,
                   m0_980   = None   ,    ## mass  of f0(980) resonance
                   m0g1     = 165000 ,    ## m0(f0(980))*gamma_1
                   g2og1    = 4.21   ) :  ## gamma2/gamma1 
        
        #
        ## initialize the base
        # 
        Flatte_pdf.__init__  ( self     , name , flatte ,
                               mass     ,
                               m0_980   ,
                               m0g1     ,
                               g2og1    )

        ## delete the created pdf (not-needed) 
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
                   mass               ,
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
        MASS.__init__  ( self , name , mass , m0_1430 , g0_1430 )
        
        self.__gamma = self.sigma
        sname  = self.gamma.GetName  ()
        stitle = self.gamma.GetTitle ()
        gname  = sname .replace ( 'sigma' , 'gamma_1430' )
        gtitle = stitle.replace ( 'sigma' , 'Gamma_1440' )
        self.gamma.SetName  ( gname  ) 
        self.gamma.SetTitle ( gtitle )

        self.g0_1430 = self.gamma 
        
        self.__m0_1430 = self.mean 
        sname  = self.mean.GetName  ()
        stitle = self.mean.GetTitle ()
        gname  = sname .replace ( 'mean' , 'm0_1430' )
        gtitle = stitle.replace ( 'mean' , 'm0_1440' )
        self.m0_1430.SetName  ( gname  ) 
        self.m0_1430.SetTitle ( gtitle ) 
        
        self.__a_lass = makeVar ( a_lass             ,
                                  'aLASS_%s'  % name ,
                                  "aLASS(%s)" % name , a_lass , 
                                  1.94e-3            ,
                                  1.94e-3            ,
                                  1.94e-3            ) 
        self.__r_lass = makeVar ( r_lass             ,
                                  'rLASS_%s'  % name ,
                                  "rLASS(%s)" % name , r_lass , 
                                  1.76e-3            ,
                                  1.76e-3            ,
                                  1.76e-3            ) 
        self.__e_lass = makeVar ( e_lass             ,
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

    @property
    def gamma ( self ) :
        """Gamma-parameter for LASS-function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        self.__gamma.setVal ( value ) 
        return self.__gamma.getVal ()

    @property
    def g0_1430 ( self ) :
        """Gamma-parameter for LASS-function"""
        return self.__gamma
    @g0_1430.setter
    def g0_1430 ( self, value ) :
        value = float ( value )
        self.__gamma.setVal ( value ) 
        return self.__gamma.getVal ()

    @property
    def m0_1430 ( self ) :
        """M0-parameter for LASS-function"""
        return self.__gamma
    @m0_1430.setter
    def m0_1430 ( self, value ) :
        value = float ( value )
        self.__m0_1430.setVal ( value ) 
        return self.__m0_1430.getVal ()

    @property
    def a ( self ) :
        """A-parameter for LASS-function"""
        return self.__a_lass
    @a.setter
    def a ( self, value ) :
        value = float ( value )
        self.__a_lass.setVal ( value ) 
        return self.__a_lass.getVal ()

    @property
    def r ( self ) :
        """R-parameter for LASS-function"""
        return self.__r_lass
    @r.setter
    def r ( self, value ) :
        value = float ( value )
        self.__r_lass.setVal ( value ) 
        return self.__r_lass.getVal ()

    @property
    def e ( self ) :
        """E-parameter for LASS-function"""
        return self.__e_lass
    @e.setter
    def e ( self, value ) :
        value = float ( value )
        self.__e_lass.setVal ( value ) 
        return self.__e_lass.getVal ()


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
                   mass              ,
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
        MASS.__init__  ( self , name , mass , bugg_m , bugg_g2 ) 
        
        self.__bugg_g2 = self.sigma
        self.__gamam   = self.sigma
        ##
        sname  = self.gamma.GetName  ()
        stitle = self.gamma.GetTitle ()
        gname  = sname .replace ( 'sigma' , 'gamma2_Bugg' )
        gtitle = stitle.replace ( 'sigma' , 'Gamma2_Bugg' )
        self.bugg_g2.SetName  ( gname  ) 
        self.bugg_g2.SetTitle ( gtitle )
        self.gamma = self.bugg_g2 
        
        self.__bugg_m = self.mean 
        sname  = self.mean.GetName  ()
        stitle = self.mean.GetTitle ()
        gname  = sname .replace ( 'mean' , 'm_Bugg' )
        gtitle = stitle.replace ( 'mean' , 'm_Bugg' )
        self.bugg_m.SetName  ( gname  ) 
        self.bugg_m.SetTitle ( gtitle ) 
        
        self.__bugg_b1 = makeVar ( bugg_b1             ,
                                 'b1Bugg_%s'  % name ,
                                 "b1Bugg(%s)" % name , bugg_b1 ,
                                 0.5848 , 
                                 0 , 2  )
        
        self.__bugg_b2 = makeVar ( bugg_b2             ,
                                   'b2Bugg_%s'  % name ,
                                   "b2Bugg(%s)" % name , bugg_b2 , 
                                   1.6663 ,
                                 1 , 2  ) 
        
        self.__bugg_a  = makeVar ( bugg_a             ,
                                   'aBugg_%s'  % name ,
                                   "aBugg(%s)" % name , bugg_a , 
                                   1.082    ,
                                   0.5 , 5  ) 
        
        self.__bugg_s1  = makeVar ( bugg_s1           ,
                                    's1Bugg_%s'  % name ,
                                    "s1Bugg(%s)" % name , bugg_s1 , 
                                    2.8              ,
                                    1 , 5            ) 
        
        self.__bugg_s2  = makeVar ( bugg_s2           ,
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

    @property
    def gamma ( self ) :
        """Gamma-parameter (g2) for Bugg function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        self.__gamma.setVal ( value ) 
        return self.__gamma.getVal ()

    @property
    def g2 ( self ) :
        """g2-parameter for Bugg function"""
        return self.__bugg_g2
    @g2.setter
    def g2 ( self, value ) :
        value = float ( value )
        self.__bugg_g2.setVal ( value ) 
        return self.__bugg_g2.getVal ()    

    @property
    def b1 ( self ) :
        """b1-parameter for Bugg function"""
        return self.__bugg_b1
    @b1.setter
    def b1 ( self, value ) :
        value = float ( value )
        self.__bugg_b1.setVal ( value ) 
        return self.__bugg_b1.getVal ()    

    @property
    def b2 ( self ) :
        """b2-parameter for Bugg function"""
        return self.__bugg_b2
    @b2.setter
    def b1 ( self, value ) :
        value = float ( value )
        self.__bugg_b2.setVal ( value ) 
        return self.__bugg_b2.getVal ()    

    @property
    def a ( self ) :
        """a-parameter for Bugg function"""
        return self.__bugg_a
    @a.setter
    def a ( self, value ) :
        value = float ( value )
        self.__bugg_a.setVal ( value ) 
        return self.__bugg_a.getVal ()    

    @property
    def s1 ( self ) :
        """s1-parameter for Bugg function"""
        return self.__bugg_s1
    @s1.setter
    def s1 ( self, value ) :
        value = float ( value )
        self.__bugg_s1.setVal ( value ) 
        return self.__bugg_s1.getVal ()    

    @property
    def s2 ( self ) :
        """s2-parameter for Bugg function"""
        return self.__bugg_s2
    @s2.setter
    def s2 ( self, value ) :
        value = float ( value )
        self.__bugg_s2.setVal ( value ) 
        return self.__bugg_s2.getVal ()    


        
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
                   mass              ,
                   beta0    = None   ) : 
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , mass , None , None    ) 

        del self.__mean
        del self.__sigma        
        mean  = None   
        sigma = None 

        self.swanson = swanson 
        beta_max     = max ( swanson.mmin() , swanson.cusp() )
        self.__beta0   = makeVar ( beta0 , 
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
        
    @property
    def beta0( self ) :
        """beta0-parameter for Swanson function"""
        return self.__beta0
    @beta0.setter
    def s2 ( self, value ) :
        value = float ( value )
        self.__beta0.setVal ( value ) 
        return self.__beta0.getVal ()    

models.append ( Swanson_pdf )

# =============================================================================
if '__main__' == __name__ : 
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )
    
# =============================================================================
# The END 
# =============================================================================
