#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/signals.py
#
#  Set of useful PDFs for various 'signal' 1D and 2D fits
#  It includes
#  - soeme empricial PDFs to describe narrow peaks: Gauss, CrystalBall, ....
#  - some PDF to describe "wide" peaks: BreitWigner,LASS, Bugg, Flatte, ...
#  - some useful PDFs to describe smooth background: phase space ;
#    expo times polynomial; phase space times polynomial, ...
#  - set of smooth non-facrorizable model for 2D fits 
#
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
# 
# =============================================================================
"""Set of useful PDFs for various 'signal' 1D and 2D fits

It includes

Empricial PDFs to describe narrow peaks

  - Gauss 
  - Crystal Ball
  - right-side Crystal Ball
  - double-side Crystal Ball
  - Needham function for J/psi, psi' and Y peaks
  - Apollonios
  - Apollonios2 (bifurcated Apollonios)
  - bifurcated Gauissian
  - double     Gauissian
  - generalized normal v1 
  - generalized normal v2
  - skew Gaussian  (temporarily removed)
  - Novosibirsk,
  - Bukin,
  - Student-T
  - bifurcated Student-T
  - SinhAsinh_pdf   
  - JohnsonSU_pdf   
  - Atlas_pdf
  - Slash_pdf
  - AsymmetricLaplace_pdf  
  - Sech_pdf   
  - Losev_pdf   
  - Logistic_pdf   
  - RaisingCosine_pdf
  - QGaussian_pdf
  - Hyperbolic_pdf
  - GenHyperbolic_pdf
  - Das_pdf
  - ExGauss_pdf
  - NormalLaplace_pdf
  - Hypatia_pdf
  - PearsonIV_pdf
  
PDF to describe 'wide' peaks

  - BreitWigner
  - BreitWigner with interference 
  - LASS
  - Bugg
  - Flatte
  - Swanson's S=wave cusp 
  - ...

Special stuff:

  - Voigt & PseudoVoigt
  - BW3L

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
    'Apollonios_pdf'         , ## Apollonios  function         
    'Apollonios2_pdf'        , ## Apollonios2 function         
    'BifurcatedGauss_pdf'    , ## bifurcated Gauss
    'DoubleGauss_pdf'        , ## double Gauss
    'GenGaussV1_pdf'         , ## generalized normal v1  
    'GenGaussV2_pdf'         , ## generalized normal v2 
    'SkewGauss_pdf'          , ## skewed gaussian (temporarily removed)
    'Novosibirsk_pdf'        , ## Novosibirsk PDF
    'Bukin_pdf'              , ## generic Bukin PDF: skewed gaussian with exponential tails
    'StudentT_pdf'           , ## Student-T function 
    'BifurcatedStudentT_pdf' , ## bifurcated Student-T function
    'PearsonIV_pdf'          , ## Pearson Type IV pdf  
    'SinhAsinh_pdf'          , ## "Sinh-arcsinh distributions". Biometrika 96 (4): 761
    'JohnsonSU_pdf'          , ## JonhsonSU-distribution 
    'Atlas_pdf'              , ## modified gaussian with exponenital tails 
    'Slash_pdf'              , ## symmetric peakk wot very heavy tails 
    'RaisingCosine_pdf'      , ## Raising  Cosine distribution
    'QGaussian_pdf'          , ## Q-gaussian distribution
    'Hyperbolic_pdf'         , ## Hyperbolic distribution
    'GenHyperbolic_pdf'      , ## Generalised Hyperbolic distribution
    'Das_pdf'                , ## Das: Gaussian with exponentrial tails 
    'Hypatia_pdf'            , ## Generalised Hyperbolic distribution
    'ExGauss_pdf'            , ## ExGauss distribution 
    'NormalLaplace_pdf'      , ## Normal Laplace distribution 
    'AsymmetricLaplace_pdf'  , ## asymmetric laplace 
    'Sech_pdf'               , ## hyperbolic secant  (inverse-cosh) 
    'Losev_pdf'              , ## asymmetric hyperbolic secant
    'Logistic_pdf'           , ## Logistic aka "sech-squared"   
    #
    ## pdfs for "wide" peaks, to be used with care - phase space corrections are large!
    # 
    'BreitWigner_pdf'        , ## (relativistic) 2-body Breit-Wigner
    'BWI_pdf'                , ## (relativistic) Breit-Wigner with interference 
    'Flatte_pdf'             , ## Flatte-function  (pipi/KK)
    'LASS_pdf'               , ## kappa-pole
    'Bugg_pdf'               , ## sigma-pole
    ##
    'Voigt_pdf'              , ## Voigt-profile
    'PseudoVoigt_pdf'        , ## PseudoVoigt-profile
    'BWMC_pdf'               , ## BWMC
    'BWPS_pdf'               , ## BWPS
    'BW3L_pdf'               , ## BW3L
    'FlattePS_pdf'           , ## Flatte-PS 
    #
    )
# =============================================================================
import ROOT, math
# =============================================================================
from   ostap.core.core          import Ostap  
from   ostap.fitting.pdfbasic   import PDF1 , all_args
from   ostap.fitting.fit1d      import PEAK , PEAKMEAN , CheckMean
from   ostap.fitting.fithelpers import Phases
import ostap.math.dalitz 
# =============================================================================
from   ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.models_signal' )
else                       : logger = getLogger ( __name__                )
# =============================================================================
models = [] 
# =============================================================================
## @class Gauss_pdf
#  simple wrapper over Gaussian-pdf
#  @see http://en.wikipedia.org/wiki/Normal_distribution
#  @see RooGaussian
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Gauss_pdf(PEAK) :
    """Trivial Gaussian function:
    http://en.wikipedia.org/wiki/Normal_distribution
    """
    def __init__ ( self               ,
                   name               ,
                   xvar               ,
                   mean        = None ,
                   sigma       = None ,
                   mean_name   = ''   , 
                   mean_title  = ''   ,
                   sigma_name  = ''   , 
                   sigma_title = ''   ) : 

        #
        ## initialize the base
        # 
        PEAK.__init__  ( self ,
                         name        = name        ,
                         xvar        = xvar        ,
                         mean        = mean        ,
                         sigma       = sigma       ,
                         mean_name   = mean_name   , 
                         mean_title  = mean_title  ,
                         sigma_name  = sigma_name  , 
                         sigma_title = sigma_title )              
        #
        ## build pdf
        # 
        self.pdf = ROOT.RooGaussian (
            self.roo_name ( 'gauss_' ) ,
            "Gauss %s" % self.name ,
            self.xvar  ,
            self.mean  ,
            self.sigma )
        ## save the configuration
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mean'  : self.mean  ,
            'sigma' : self.sigma ,
            }

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
class CrystalBall_pdf(PEAK) :
    """Crystal Ball function
    http://en.wikipedia.org/wiki/Crystal_Ball_function
    
    - T. Skwarnicki,
    'A study of the radiative cascade transitions between the Upsilon-Prime
    and Upsilon resonances', Ph.D Thesis, DESY F31-86-02(1986), Appendix E.
    http://inspirehep.net/record/230779/files/f31-86-02.pdf
    - J. E. Gaiser,
    'Charmonium Spectroscopy from Radiative Decays of the J/Psi and Psi-Prime',
    Ph.D. Thesis, SLAC-R-255 (1982), Appendix F, p 178
    http://www.slac.stanford.edu/cgi-wrap/getdoc/slac-r-255.pdf
    - M. J. Oreglia,
    'A Study of the Reactions psi prime --> gamma gamma psi',
    Ph.D. Thesis, SLAC-R-236 (1980), Appendix D.
    http://www.slac.stanford.edu/pubs/slacreports/slac-r-236.html
    
    Note:
    - Unlike the original definition parameter 'n' here is shifted by 1: n <- |n| + 1
    - Typical value of parameter alpha for 'physical' peaks is 1.5<alpha<3.0,
    - For large alpha (e.g. alpha>3), there is no sensitivity for n;
    similarly in the limit of large n, sensitivity for alpha is minimal
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mean     = None  , 
                   sigma    = None  ,
                   alpha    = None  ,
                   n        = None  ) : 
                   
        #
        ## initialize the base
        #
        PEAK.__init__ ( self , name , xvar , mean , sigma )
        
        self.__alpha = self.make_var ( alpha ,
                                       'alpha_%s'        % name ,
                                       '#alpha_{CB}(%s)' % name ,
                                       None , 2.0 , 0.05 , 5 )
        
        self.__n     = self.make_var ( n   ,
                                       'n_%s'            % name ,
                                       'n_{CB}(%s)'      % name ,
                                       None , 5.0 , 1.e-6 , 100 )
        
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.CrystalBall (
            self.roo_name ( 'cb_'  )       , 
            'Crystal Ball %s' % self.name  ,
            self.xvar  ,
            self.mean  ,
            self.sigma ,
            self.alpha ,
            self.n     )
        ## save the configuration
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mean'  : self.mean  ,
            'sigma' : self.sigma ,
            'alpha' : self.alpha ,
            'n'     : self.n     ,
            }

    @property
    def alpha ( self ) :
        """Alpha-parameter for Crystal Ball tail"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        self.set_value ( self.__alpha , value ) 
    
    @property
    def n ( self ) :
        """N-parameter for Crystal Ball tail (actually 'n+1' is used)"""
        return self.__n
    @n.setter
    def n ( self, value ) :
        self.set_value ( self.__n , value ) 

models.append ( CrystalBall_pdf )    
# =============================================================================
## @class CrystalBallRS_pdf
#  Crystal Ball function with the right side tail.
#  to be rewritten
#  @see Ostap::Models::CrystalBallRS
#  @see Ostap::Math::CrystalBallRS
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class CrystalBallRS_pdf(CrystalBall_pdf) :
    """Right-side CrystalBall    
    """
    def __init__ ( self              ,
                   name              ,
                   xvar              ,
                   mean     = None   , 
                   sigma    = None   ,
                   alpha    = None   ,
                   n        = None   ) : 
                   
        ##  the base 
        CrystalBall_pdf.__init__  ( self , name , xvar , mean , sigma , alpha , n )
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.CrystalBallRS (
            self.roo_name ( 'cbrs_' )  , 
            'Right-sided Crystal Ball %s' % self.name ,
            self.xvar  ,
            self.mean  ,
            self.sigma ,
            self.alpha ,
            self.n     )
        
        ## save the configuration
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mean'  : self.mean  ,
            'sigma' : self.sigma ,
            'alpha' : self.alpha ,
            'n'     : self.n     ,
            }
      
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
class CB2_pdf(PEAK) :
    """Double sided Crystal Ball function with both left and rigth sides
    It appears to be very powerful and is used for many LHCb papers to describe
    B-hadron mass signals, especially for B->J/psi X final states 
    
    Note:
    - Similar to CrystalBall_pdf and unlike the original definition,
    the parameters 'n' here are shifted by 1: n <- |n| + 1
    - Typical value of parameters alpha for 'physical' peaks is 1.5<alpha<3.0,
    - For large alpha (e.g. alpha>3), there is no sensitivity for n;
    similarly in the limit of large n, sensitivity for alpha is minimal
    """
    def __init__ ( self              ,
                   name              ,
                   xvar              , 
                   mean      = None  ,
                   sigma     = None  ,
                   alphaL    = None  ,
                   alphaR    = None  ,
                   nL        = None  ,
                   nR        = None  ) : 
        
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar , mean , sigma )
        #
        ## treat the specific parameters
        #
        self.__aL    = self.make_var ( alphaL                  ,
                                       "aL_%s"          % name ,
                                       "#alpha_{L}(%s)" % name ,
                                       None , 2.0 ,  0.05 , 5 )
        self.__nL    = self.make_var ( nL                      ,                     
                                       "nL_%s"          % name ,
                                       "n_{L}(%s)"      % name ,
                                       None  , 5 , 1.e-6 , 100 )
        self.__aR    = self.make_var ( alphaR ,
                                       "aR_%s"          % name ,
                                       "#alpha_{R}(%s)" % name ,
                                       None , 2.0 ,  0.05 , 5 )
        self.__nR    = self.make_var ( nR                      ,
                                       "nR_%s"          % name ,
                                       "n_{R}(%s)"      % name ,
                                       None  , 5 , 1.e-6 , 100 )
        
        self.pdf = Ostap.Models.CrystalBallDS(
            self.roo_name ( 'cb2_' ) , 
            "double-sided Crystal Ball %s" % self.name ,
            self.xvar    ,
            self.mean    ,
            self.sigma   ,
            self.aL      ,
            self.nL      ,
            self.aR      ,
            self.nR      )

        ## save the configuration
        self.config = {
            'name'   : self.name  ,
            'xvar'   : self.xvar  ,
            'mean'   : self.mean  ,
            'sigma'  : self.sigma ,
            'alphaL' : self.aL    ,
            'alphaR' : self.aR    ,
            'nL'     : self.nL    ,
            'nR'     : self.nR    ,
            }

    @property
    def aL ( self ) :
        """(left) Alpha-parameter for Crystal Ball tail"""
        return self.__aL
    @aL.setter
    def aL ( self, value ) :
        self.set_value ( self.__aL , value ) 

    @property
    def nL ( self ) :
        """(left) N-parameter for Crystal Ball tail"""
        return self.__nL
    @nL.setter
    def nL ( self, value ) : 
        self.set_value ( self.__nL , value ) 

    @property
    def aR ( self ) :
        """(right) Alpha-parameter for Crystal Ball tail"""
        return self.__aR
    @aR.setter
    def aR ( self, value ) :
        self.set_value ( self.__aR , value ) 

    @property
    def nR ( self ) :
        """(right) N-parameter for Crystal Ball tail"""
        return self.__nR
    @nR.setter
    def nR ( self, value ) :
        self.set_value ( self.__nR , value ) 
        
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
#  @see CrystalBall_pdf 
#  @see Ostap::Models::Needham 
#  @see Ostap::Math::Needham 
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Needham_pdf(PEAK) :
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
                   xvar                 ,
                   mean     = 3.096     ,   ## GeV  
                   sigma    = 0.013     ,   ## GeV 
                   a0       = 1.975     ,
                   a1       = -0.0011   ,   ## GeV^-1
                   a2       = -0.00018  ) : ## GeV^-2  
        
        PEAK.__init__ ( self , name , xvar , mean , sigma )
        
        #
        unit = 1000
        #
        if   3.096 in self.xvar : unit = 1000 
        elif 3096  in self.xvar : unit = 1
        elif 9.460 in self.xvar : unit = 1000 
        elif 9460  in self.xvar : unit = 1
        #
        self.__a0 = self.make_var ( a0                  ,
                                    "a0_%s"     % name  ,
                                    "a_{0}(%s)" % name  ,
                                    True , 1.975 , 0 , 10  )
        
        self.__a1 = self.make_var ( a1                  ,
                                    "a1_%s"     % name  ,
                                    "a_{1}(%s)" % name  ,
                                    True , -0.0011   * unit , -10 * unit , 10 * unit )
        
        self.__a2 = self.make_var ( a2                  ,
                                    "a2_%s"     % name  ,
                                    "a_{2}(%s)" % name  ,
                                    True , -0.00018  * unit**2 , -10 * unit**2 , 10 * unit**2 )
        #
        self.pdf = Ostap.Models.Needham (
            self.roo_name ( 'needham_' ) , 
            'Needham function %s' % self.name ,
            self.xvar  ,
            self.mean  ,
            self.sigma ,
            self.a0    ,
            self.a1    ,
            self.a2
            )
        
        ## save the configuration
        self.config = {
            'name'   : self.name  ,
            'xvar'   : self.xvar  ,
            'mean'   : self.mean  ,
            'sigma'  : self.sigma ,
            'a0'     : self.a0    ,
            'a1'     : self.a1    ,
            'a2'     : self.a2    ,
            }
        
    @property
    def a0 ( self ) :
        """'a0'-parameter for Needham function"""
        return self.__a0
    @a0.setter
    def a0 ( self, value ) :
        self.set_value ( self.__a0 , value ) 

    @property
    def a1 ( self ) :
        """'a1'-parameter for Needham function"""
        return self.__a1
    @a1.setter
    def a1 ( self, value ) :
        self.set_value ( self.__a1 , value ) 

    @property
    def a2 ( self ) :
        """'a2'-parameter for Needham function"""
        return self.__a2
    @a2.setter
    def a2 ( self, value ) :
        self.set_value ( self.__a2 , value ) 

            
models.append ( Needham_pdf )    
# =============================================================================
## @class Apollonios_pdf
#  simple wrapper over Apollonios PDF 
#  @see Ostap::Models::Apollonios 
#  The function is proposed by Diego Martinez Santos 
#  @see http://arxiv.org/abs/1312.5000
#  Here a bit modified version is used with redefined parameter <code>n</code>
#  to be coherent with local definitions of Crystal Ball
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Apollonios_pdf(PEAK) :
    """Apollonios function
    http://arxiv.org/abs/1312.5000
    
    The function is proposed by Diego Martinez Santos 
    https://indico.itep.ru/getFile.py/access?contribId=2&resId=1&materialId=slides&confId=262633
    
    Note :
    - similar to CrystalBall case,  parameter n is redefined:   n <- |n|+1
    to be coherent with local definitions of Crystal Ball
    - unfortuately neither sigma nor b parameters allows easy interpretation
    - typical value of parameters alpha for 'physical' peaks is 1.5<alpha<2.1,
    - for large alpha (e.g. alpha>3), there is no sensitivity for n;
    similarly in the limit of large n, sensitivity for alpha is minimal
    """
    def __init__ ( self                    ,
                   name                    ,
                   xvar                    ,
                   mean      = None        ,
                   sigma     = None        ,
                   alpha     = None        ,
                   n         = None        ,
                   b         = None        ) : 
                   
        
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar , mean , sigma ) 
        
        self.__alpha = self.make_var ( alpha                     ,
                                       'alpha_%s'         % name ,
                                       '#alpha_{Apo}(%s)' % name ,
                                       None , 2.0 , 0.01 , 5 )
        
        self.__n     = self.make_var ( n                    ,
                                       'n_%s'        % name ,
                                       'n_{Apo}(%s)' % name ,
                                       None , 5.0 , 1.e-6 , 100  )
        
        self.__b     = self.make_var ( b                    ,
                                       'b_%s'        % name ,
                                       'b_{Apo}(%s)' % name ,
                                       None , 1.e-5 , 10000 ) 
        #
        ## finally build PDF
        #
        self.pdf  = Ostap.Models.Apollonios (
            self.roo_name ( 'apo_' ) , 
            "Apollonios %s" % self.name   ,
            self.xvar   ,
            self.mean   ,
            self.sigma  ,
            self.alpha  ,
            self.n      ,
            self.b      ) 

        ## save the configuration
        self.config = {
            'name'   : self.name  ,
            'xvar'   : self.xvar  ,
            'mean'   : self.mean  ,
            'sigma'  : self.sigma ,
            'alpha'  : self.alpha ,
            'n'      : self.n     ,
            'b'      : self.b     ,
            }
        
    @property
    def alpha ( self ) :
        """'alpha'-parameter for Apollonios tail"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        self.set_value ( self.__alpha , value ) 
    
    @property
    def n ( self ) :
        """'n'-parameter for Apollonios tail"""
        return self.__n
    @n.setter
    def n ( self, value ) :
        self.set_value ( self.__n , value ) 

    @property
    def b ( self ) :
        """'b'-parameter for Apollonios function"""
        return self.__b
    @b.setter
    def b ( self, value ) :
        self.set_value ( self.__b , value ) 

models.append ( Apollonios_pdf )    
# =============================================================================
## @class Apollonios2_pdf
#  "Bifurcated Apollonios"
#  Gaussian with exponential (asymmetrical) tails
#
#  A convinient reparameterization is applied to keep reduce 
#  the correlations between "sigma"s and "beta" 
#  \f[ f(x;\mu,\sigma_l,\sigma_r,\beta) \propto 
#         \mathrm{e}^{\left|\beta\right|( \left|\beta\right| -
#         \sqrt{ \beta^2+\left(\delta x\right)^2}} \f] 
#   where 
#  \f[ \delta x  = \left\{ \begin{array}{ccc}
#          \dfrac{x-\mu}{\sigma_l} & \text{for} & x \le \mu \\
#          \dfrac{x-\mu}{\sigma_r} & \text{for} & x >   \mu \\
#          \end{array}
#          \right.\f]
#  Large betas corresponds to gaussian 
#      
#  @see Ostap::Models::Apollonios2 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-08-20
class Apollonios2_pdf(PEAK) :
    """Bifurcated Apollonios:
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
    Reparameterization reduces the correlation between 'sigma' and 'beta',
    allowing their easy interpretation.
    Large 'beta' and small 'asymmetry' corresponds to gaussian
    
    """
    def __init__ ( self               ,
                   name               ,
                   xvar               ,
                   mean      = None   ,
                   sigma     = None   ,
                   asymmetry = None   ,
                   beta      = None   ) : 
        
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar , mean , sigma  )
        
        ## asymmetry parameter 
        self.__asym = self.make_var ( asymmetry                 ,
                                      'asym_%s'          % name ,
                                      '#asym_{Apo2}(%s)' % name ,
                                      None , 0 , -1 , 1 )
        
        ## construct left and right sigmas 
        self.__sigmaL , self.__sigmaR = self.vars_from_asymmetry (
            self.sigma                                    , ## mean/average sigma 
            self.asym                                     , ## asymmetry parametet
            v1name  =  self.roo_name ( 'sigmaL' , self.name ) ,
            v2name  =  self.roo_name ( 'sigmaR' , self.name ) ,
            v1title = '#sigma_{L}: #sigma #times (1+#kappa)'        , 
            v2title = '#sigma_{R}: #sigma #times (1-#kappa)'        )

        self.__beta    = self.make_var ( beta ,
                                         'beta_%s'          % name ,
                                         '#beta_{Apo2}(%s)' % name ,
                                         None , 1 , 1.e-5 , 10000  ) 
        #
        ## finally build PDF
        #
        self.pdf  = Ostap.Models.Apollonios2 (
            self.roo_name ( 'apo2_' ) , 
            "Apollonios2 %s" % self.name ,
            self.xvar   ,
            self.mean   ,
            self.sigmaL ,
            self.sigmaR ,
            self.beta   )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            'asymmetry' : self.asym  ,
            'beta'      : self.beta  ,
            }

    @property
    def asym ( self ) :
        """'asym'- asymmetry parameter for Apollonios2 function (same as 'kappa')"""
        return self.__asym
    @asym.setter
    def asym ( self, value ) :
        self.set_value ( self.__asym , value ) 

    @property
    def kappa ( self ) :
        """'kappa'-parameter for Apollonios2 function (same as 'asym')"""
        return self.__asym
    @kappa.setter
    def kappa ( self, value ) :
        self.set_value ( self.__asym , value ) 

    @property
    def beta ( self ) :
        """'beta'-parameter for Apollonios-2 function"""
        return self.__beta
    @beta.setter
    def beta ( self, value ) :
        self.set_value ( self.__beta , value ) 

    @property
    def sigmaL ( self ) :
        """(left) sigma-parameter for Apollonios-2 function"""
        return self.__sigmaL
    
    @property
    def sigmaR ( self ) :
        """(right) sigma-parameter for Apollonios2 function"""
        return self.__sigmaR

    
models.append ( Apollonios2_pdf )    
# =============================================================================
## @class BifurcatedGauss_pdf
#  simple wrapper over bifurcated-gaussian
#  @see RooGaussian
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class BifurcatedGauss_pdf(PEAK) :
    """Bifurcated Gauss function
    
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
                   xvar                  ,
                   mean      = None      ,
                   sigma     = None      ,
                   asymmetry = None      ) : 
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar , mean , sigma ) 
        
        ## asymmetry parameter  
        self.__asym = self.make_var ( asymmetry                 ,
                                      'asym_%s'          % name ,
                                      '#asym_{asym}(%s)' % name ,
                                      None , 0 , -1 , 1  )
        
        ## constreuct left and right sigmas 
        self.__sigmaL , self.__sigmaR = self.vars_from_asymmetry (
            self.sigma                                    , ## mean/average sigma 
            self.asym                                     , ## asymmetry parametet
            v1name  =  self.roo_name ( 'sigmaL' , self.name ) ,
            v2name  =  self.roo_name ( 'sigmaR' , self.name ) ,
            v1title = '#sigma_{L}: #sigma #times (1+#kappa)'        , 
            v2title = '#sigma_{R}: #sigma #times (1-#kappa)'        )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.BifurcatedGauss (
            self.roo_name ( "bfgauss_" )  , 
            "Bifurcated Gauss %s" % self.name ,
            self.xvar   ,
            self.mean   ,
            self.sigmaL ,
            self.sigmaR )

        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            'asymmetry' : self.asym  ,
            }

    @property
    def asym ( self ) :
        """'asymmetry'-parameter for Bifurcated Gaussian"""
        return self.__asym
    @asym.setter
    def asym ( self, value ) :
        self.set_value ( self.__asym , value ) 

    @property
    def kappa ( self ) :
        """'kappa'-parameter for Bifurcated Gaussian function (same as 'asym'"""
        return self.__asym
    @kappa.setter
    def kappa ( self, value ) :
        self.set_value ( self.__asym , value ) 

    @property
    def sigmaL ( self ) :
        """(left)'sigma'-parameter for Bifurcated Gaussian"""
        return self.__sigmaL
    
    @property
    def sigmaR ( self ) :
        """(right)'sigma'-parameter for Bifurcated Gaussian"""
        return self.__sigmaR

models.append ( BifurcatedGauss_pdf )


# =============================================================================
## @class DoubleGauss_pdf
#  simple wrapper over double gaussian
#  @see RooGaussian
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2017-07-12
class DoubleGauss_pdf(PEAK) :
    """Double Gauss :

    f(x;mu,sigma,scale,fraction) =
    fraction     * G(x;mu,sigma) +
    (1-fraction) * G(x;mu;sigma*scale) 
        
    """
    def __init__ ( self                  ,
                   name                  ,
                   xvar                  ,
                   mean      = None      ,
                   sigma     = None      ,
                   fraction  = None      ,
                   scale     = None      ) :  
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar , mean , sigma )
        
        self.__scale    = self.make_var ( scale ,
                                          'SigmaScale'        + name ,
                                          'SigmaScale(%s)'    % name ,
                                          None , 1.5 , 1 , 10 ) 
        
        ## the fraction 
        self.__fraction = self.make_var ( fraction                   , 
                                          'CoreFraction'      + name ,
                                          'CoreFraction(%s)'  % name ,
                                          None , 0.75 , 0 , 1  ) 
        
        self.pdf = Ostap.Models.DoubleGauss (
            self.roo_name ( 'gauss2_' ) , 
            "double Gauss %s" % self.name ,
            self.mass     ,
            self.sigma    ,
            self.fraction ,
            self.scale    ,
            self.mean    
            )

        ## save the configuration
        self.config = {
            'name'      : self.name     ,
            'xvar'      : self.xvar     ,
            'mean'      : self.mean     ,
            'sigma'     : self.sigma    ,
            'fraction'  : self.fraction ,
            'scale'     : self.scale    ,
            }

    @property
    def scale ( self ) :
        """Scale-parameter for double Gaussian"""
        return self.__scale
    @scale.setter
    def scale ( self, value ) :
        self.set_value ( self.__scale , value ) 

    @property
    def fraction ( self ) :
        """Fraction-parameter for Bifurcated Gaussian"""
        return self.__fraction
    @fraction.setter
    def fraction ( self, value ) :
        self.set_value ( self.__fraction , value ) 

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
#
#  @see M. T. Subbotin, “On the Law of Frequency of Error”, Mat. Sb., 31:2 (1923), 296–301
#  @see http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=sm&paperid=6854&option_lang=eng
#  @see Nadarajah, Saralees (September 2005). "A generalized normal distribution".
#       Journal of Applied Statistics. 32 (7): 685–694. doi:10.1080/02664760500079464.
#  @see https://doi.org/10.1080%2F02664760500079464
class GenGaussV1_pdf(PEAK) :
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

    - see M. T. Subbotin, “On the Law of Frequency of Error”, Mat. Sb., 31:2 (1923), 296–301
    - see http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=sm&paperid=6854&option_lang=eng
    - see Nadarajah, Saralees (September 2005). "A generalized normal distribution".
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
                   mean      = None ,
                   alpha     = None ,
                   beta      = 2    ) :  ## beta=2 is gaussian distribution 
        
        
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self  , name , xvar               ,
                              mean        = mean                ,
                              sigma       = alpha               ,
                              sigma_name  = 'alpha_%s'   % name ,
                              sigma_title = '#alpha(%s)' % name )
        
        #
        ## rename it!
        #
        self.__alpha = self.sigma
        
        self.__beta  = self.make_var ( beta ,
                                       'beta_%s'        % name  ,
                                       '#beta_{v1}(%s)' % name  , beta , 
                                       2 , 1.e-4  , 1.e+6 ) 
        #
        ## finally build PDF
        #
        self.pdf = Ostap.Models.GenGaussV1 (
            self.roo_name ( 'gaussgv1_' ) , 
            "generalized Gauss-V1 %s" % self.name ,
            self.xvar   ,
            self.mean   ,
            self.alpha  ,
            self.beta   )

        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'alpha'     : self.alpha ,
            'beta'      : self.beta  ,
            }

    @property
    def alpha ( self ) :
        """'alpha'-parameter for Generalized V1 Gaussian (the same as 'sigma')"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        self.set_value ( self.__alpha , value ) 

    @property
    def beta ( self ) :
        """'beta'-parameter for Generalized  V1 Gaussian"""
        return self.__beta
    @beta.setter
    def beta ( self, value ) :
        self.set_value ( self.__beta , value ) 
        
models.append ( GenGaussV1_pdf )    
# =============================================================================
## @class GenGaussV2_pdf
#  Generalized Normal distribution v2
#  @see http://en.wikipedia.org/wiki/Generalized_normal_distribution#Version_2
#  @see Ostap::Models::GenGaussV2
#  @see Ostap::Math::GenGaussV2 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2013-12-01
class GenGaussV2_pdf(PEAK) :
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
                   xvar               ,
                   mean        = None ,
                   alpha       = None ,
                   kappa       = 0    ) : ## 0 corresponds to gaussian distribution 
        
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar ,
                         mean        = mean  ,
                         sigma       = alpha ,
                         sigma_name  = 'alpha_%s'   % name ,
                         sigma_title = '#alpha(%s)' % name )
        #
        ## rename it!
        #
        self.__alpha = self.sigma        
        self.__xi    = self.mean 
        self.__kappa = self.make_var ( kappa ,
                                       'kappa_%s'        % name  ,
                                       '#kappa_{v2}(%s)' % name  ,
                                       None , 0 , -5 , 5 )
        #
        ## finally build PDF
        #
        self.pdf = Ostap.Models.GenGaussV2 (
            self.roo_name ( 'gaussgv2_' ) , 
            "generalized Gauss-V2 %s" % self.name ,
            self.xvar   ,
            self.xi     ,
            self.alpha  ,
            self.kappa  )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'alpha'     : self.alpha ,
            'kappa'     : self.kappa ,
            }
        
    @property
    def alpha ( self ) :
        """alpha-parameter for Generalized V2 Gaussian (same as 'sigma')"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        self.set_value ( self.__alpha , value ) 

    @property
    def kappa ( self ) :
        """kappa-parameter for Generalized V2 Gaussian"""
        return self.__kappa
    @kappa.setter
    def kappa ( self, value ) :
        self.set_value ( self.__kappa , value ) 
    
    @property
    def xi ( self ) :
        """xi-parameter (location) for Generalized V2 Gaussian (same  as 'mean')"""
        return self.__xi
 
models.append ( GenGaussV2_pdf )    
# =============================================================================
## @class SkewGauss_pdf
#  Simple class that implements the skew normal distribution
#  @see http://en.wikipedia.org/wiki/Skew_normal_distribution
#  @see Ostap::Models::SkewGauss 
#  @see Ostap::Math::SkewGauss 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class SkewGauss_pdf(PEAK) :
    """Skew normal distribution
    see http://en.wikipedia.org/wiki/Skew_normal_distribution
    
    The skew normal distribution is a continuous probability distribution that
    generalises the normal distribution to allow for non-zero skewness.

    Parameters:
    - location
    - omega>0   : scale
    - alpha     : shape   (alpha=0 corresponds to gaussian distribuition)
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mean      = None ,
                   omega     = None ,
                   alpha     = 0    ) : ## alpha=0 correspond to Gaussian 
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar , 
                         mean        = mean  ,
                         sigma       = omega ,
                         sigma_name  = 'omega_%s'   % name ,
                         sigma_title = '#omega(%s)' % name )
        
        #
        ## rename it!
        #
        
        self.__omega = self.sigma        
        self.__xi    = self.mean         
        self.__alpha = self.make_var ( alpha ,
                                       'alpha_%s'   % name  ,
                                       '#alpha(%s)' % name  ,
                                       None , 0 , -1000 , 1000  ) 
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.SkewGauss (
            self.roo_name ( 'gausssk_' ) , 
            "skew Gauss %s" % self.name ,
            self.xvar   ,
            self.mean   ,
            self.omega  ,
            self.alpha  )

        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'omega'     : self.omega ,
            'alpha'     : self.alpha ,
            }

    @property
    def xi ( self ) :
        """xi-parameter (location) for Skew Gaussian (Same  as 'mean')"""
        return self.__xi
    @xi.setter
    def xi ( self, value ) :
        self.set_value ( self.__xi , value ) 

    @property
    def omega ( self ) :
        """omega-parameter (scale/width) for Skew Gaussian (same as 'sigma')"""
        return self.__omega
    @omega.setter
    def omega ( self, value ) :
        self.set_value ( self.__omega , value ) 

    @property
    def alpha( self ) :
        """alpha-parameter (shape/skew) for Skew Gaussian"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        self.set_value ( self.__alpha , value ) 

models.append ( SkewGauss_pdf )


# =============================================================================
## @class Novosibirsk_pdf
#  Novosibirsk-function for description of gaussian with tails
#  @see H.Ikeda et al., 'A detailed test of the CsI(Tl) calorimeter 
#      for BELLE with photon beams of energy between 20MeV and 5.4 GeV',
#       Nucl. Instrum. Meth. A441, (2000) 401.
#  @see DOI: 10.1016/S0168-9002(99)00992-4
#  @see https://inspirehep.net/literature/508223 
#  @see https://doi.org/10.1016/S0168-9002(99)00992-4
#
#  \f$ f(x;\mu,\sigma,\tau) = \frac{1}{\sqrt{2\pi}\sigma}
#  \mathrm{e}^{  -\frac{1}{2} \frac { \log^2 \left( 1 + \Lambda \tau \delta \right) }{\tau^2} 
#                -\frac{\tau^2}{2} } \f$
#  where
#  - \f$ \delta  = \frac{ x - \mu}{\sigma}\f$ 
#  - \f$ \Lambda = \frac{  \sinh{ \tau \sqrt{\log 4}} }{\tau\sqrt{\log 4 }}\f$ 
class Novosibirsk_pdf(PEAK) :
    """Novosibirsk-function for description of gaussian with tails
    - see H.Ikeda et al., 'A detailed test of the CsI(Tl) calorimeter 
    for BELLE with photon beams of energy between 20MeV and 5.4 GeV',
    Nucl. Instrum. Meth. A441, (2000) 401.
    - see DOI: 10.1016/S0168-9002(99)00992-4
    - see https://inspirehep.net/literature/508223 
    - see https://doi.org/10.1016/S0168-9002(99)00992-4
    """
    def __init__ ( self          ,
                   name          ,
                   xvar          , 
                   mean   = None ,
                   sigma  = None ,
                   tau    = 0    ) : ## tails/&asymmety
        
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name  , xvar , mean , sigma )
        #
        ## treat the specific parameters

        ## tail&asymmetry 
        self.__tau = self.make_var ( tau                    ,
                                     "tau_%s"        % name ,
                                     "#tau(%s)"      % name ,
                                     True , 0  , -2 , 2 )
        # 
        ## create PDF
        # 
        self.pdf = Ostap.Models.Novosibirsk  (
            self.roo_name ( 'novo_' ) , 
            "Novosibirsk %s" % self.name ,
            self.xvar  ,
            self.mean  ,
            self.sigma ,
            self.tau   )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            'tau'       : self.tau
            }

    @property
    def tau ( self ) :
        """'tau'  : tails&asymmetry for Novosibirsk function"""
        return self.__tau
    @tau.setter
    def tau ( self, value ) :
        self.set_value ( self.__tau , value ) 

models.append ( Novosibirsk_pdf )      

# =============================================================================
## @class Bukin_pdf
#  Bukin function, aka 'modified Novosibirsk function'
#  - asymmetrical gaussian-like core
#  - exponential (optionally gaussian) asymmetrical tails
#  @see http://journals.aps.org/prd/abstract/10.1103/PhysRevD.84.112007
#  @see http://arxiv.org/abs/1107.5751
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
#  @see https://doi.org/10.1007/JHEP06(2012)141     
#  @see Ostap::Math::Bukin
#  @see Analusis::Models::Bukin
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Bukin_pdf(PEAK) :
    """Bukin function, aka 'modified Novosibirsk function':
    - asymmetrical gaussian-like core
    - exponential (optionally gaussian) asymmetrical tails
    see http://journals.aps.org/prd/abstract/10.1103/PhysRevD.84.112007
    see http://arxiv.org/abs/1107.5751
    see https://doi.org/10.1007/JHEP06(2012)141     
    Here small reparameterization is applied to achieve more stable fits.
    
    It is very well suitable to describe high statistic charm meson peaks,
    as well some asymmetric distributions, e.g. log( chi^2(IP) ) 

    Parameters:
    - mean      : peak position 
    - sigma     : defines full width at half heigh 
    - xi        : asymmetry
    - rho_L>=0  : define gaussian component to left  tail
    - rho_R>=0  : define gaussian component to right tail
    Note: 
    - for rho_{L,R}=0 left/right tails are exponential
    - for large asymmetry parameter function has weird shape    
    """
    def __init__ ( self            ,
                   name            ,
                   xvar            , 
                   mean     = None ,
                   sigma    = None ,
                   xi       = None ,
                   rhoL     = None ,
                   rhoR     = None ) :

        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name  , xvar , mean , sigma )
        #
        ## treat the specific parameters
        #
        
        ## asymmetry 
        self.__xi    = self.make_var ( xi                    ,
                                       "xi_%s"        % name ,
                                       "#xi(%s)"      % name ,
                                       None , 0  , -1 , 1    )
        ## left tail
        self.__rhoL  = self.make_var ( rhoL                  ,
                                       "rhoL_%s"      % name ,
                                       "#rho_{L}(%s)" % name ,
                                       None , 0  ,  -1 , 100 )        
        ## right tail
        self.__rhoR  = self.make_var ( rhoR                  ,
                                       "rhoR_%s"      % name ,
                                       "#rho_{R}(%s)" % name ,
                                       None , 0  ,  -1 , 100 )
        # 
        ## create PDF
        # 
        self.pdf = Ostap.Models.Bukin (
            self.roo_name ( 'bukin_' ) , 
            "Bukin %s" % self.name ,
            self.xvar  ,
            self.mean  ,
            self.sigma ,
            self.xi    ,
            self.rhoL  ,
            self.rhoR  )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            'xi'        : self.xi    ,
            'rhoL'      : self.rhoL  ,
            'rhoR'      : self.rhoR  ,
            }

    @property
    def xi ( self ) :
        """'xi'-parameter (asymmetry) for Bukin function"""
        return self.__xi
    @xi.setter
    def xi ( self, value ) :
        self.set_value ( self.__xi , value ) 

    @property
    def rhoL ( self ) :
        """'rho'-parameter (left tail) for Bukin function"""
        return self.__rhoL
    @rhoL.setter
    def rhoL ( self, value ) :
        self.set_value ( self.__rhoL , value ) 

    @property
    def rhoR ( self ) :
        """'rho'-parameter (right tail) for Bukin function"""
        return self.__rhoR
    @rhoR.setter
    def rhoR ( self, value ) :
        self.set_value ( self.__rhoR , value ) 


models.append ( Bukin_pdf )      
# =============================================================================
## @class StudentT_pdf
#  Student-T distribution
#  @see http://en.wikipedia.org/wiki/Student%27s_t-distribution
#
#  \f[  f(y) = \dfrac{1}{\sqrt{\pi n}} \dfrac { \Gamma( \dfrac{n+1}{2}) } { \Gamma( \dfrac{n}{2}  ) }
#  \left( 1 + \dfrac{y^2}{n} \right)^{ -\dfrac{n+1}{2}} \f], 
#  where \f$ y = \dfrac{x - \mu}{\sigma} \f$  
# 
#  @attention Unlike the original definition parameter 'n' here is shifted by 1: n <- |n| + 1
#  @see Ostap::Models::StudentT
#  @see Ostap::Math::StudentT
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class StudentT_pdf(PEAK) :
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
                   xvar             ,
                   mean      = None ,
                   sigma     = None ,
                   n         = None ) :
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar , mean , sigma ) 
        
        # 
        self.__n  = self.make_var ( n                    ,
                                    'n_%s'        % name ,
                                    '#n_{ST}(%s)' % name ,
                                    None , 2 , 1.e-8 , 100  ) 
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.StudentT (
            self.roo_name ( 'studentt_' ) , 
            "Student's t %s" % self.name ,
            self.xvar   ,
            self.mean   ,
            self.sigma  ,
            self.n      )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            'n'         : self.n     ,
            }

    @property
    def n ( self ) :
        """'n'-parameter for Student-T function (well, actually it is 'n+1' is used)"""
        return self.__n
    @n.setter
    def n ( self, value ) :
        self.set_value ( self.__n , value ) 
    
models.append ( StudentT_pdf )

# =============================================================================
## @class BifurcatedStudentT_pdf
#  bifurcated Student-T distribution
# 
#  @see Ostap::Models::BifurcatedStudentT
#  @see Ostap::Math::BifurcatedStudentT
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class BifurcatedStudentT_pdf(PEAK) :
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
                   xvar             ,
                   mean      = None ,
                   sigma     = None ,
                   asymmetry = None ,
                   nL        = None , 
                   nR        = None ) :
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar , mean , sigma )
        
        ##  asymmetry 
        self.__asym = self.make_var ( asymmetry               ,
                                      'asym_%s'        % name ,
                                      '#xi_{asym}(%s)' % name ,
                                      None , 0 , -1 , 1  ) 
        
        ## construct left and right sigmas 
        self.__sigmaL , self.__sigmaR = self.vars_from_asymmetry (
            self.sigma                                    , ## mean/average sigma 
            self.asym                                     , ## asymmetry parametet
            v1name  =  self.roo_name ( 'sigmaL' , self.name ) ,
            v2name  =  self.roo_name ( 'sigmaR' , self.name ) ,
            v1title = '#sigma_{L}: #sigma #times (1+#kappa)'        , 
            v2title = '#sigma_{R}: #sigma #times (1-#kappa)'        )

        ## left exponent 
        self.__nL =  self.make_var ( nL                     ,
                                     'nL_%s'         % name ,
                                     '#nL_{BST}(%s)' % name ,
                                     None , 2 , 1.e-6 , 200 )
        ## right exponent 
        self.__nR =  self.make_var ( nR                    ,
                                     'nR_%s'         % name ,
                                     '#nR_{BST}(%s)' % name ,
                                     None ,  2 , 1.e-6 , 200 ) 
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.BifurcatedStudentT (
            self.roo_name ( 'sttbf_' ) , 
            "Bifurcated Student's t %s" % self.name ,
            self.xvar   ,
            self.mean   ,
            self.sigmaL ,
            self.sigmaR ,
            self.nL     ,
            self.nR     )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            'asymmetry' : self.asym  ,
            'nL'        : self.nL    ,
            'nR'        : self.nR    ,
            }

    @property     
    def nL ( self ) :
        """'n'-parameter (left) for Bifurcated Student-T function  (actually 'n+1' is used)"""
        return self.__nL
    @nL.setter
    def nL ( self, value ) :
        self.set_value ( self.__nL , value ) 

    @property
    def nR ( self ) :
        """'n'-parameter (right) for Bifurcated Student-T function (actually 'n+1' is used)"""
        return self.__nR
    @nR.setter
    def nR ( self, value ) :
        self.set_value ( self.__nR , value ) 
        
    @property
    def asym ( self ) :
        """'asymmetry'-parameter for Bifurcated Student-T function"""
        return self.__asym
    @asym.setter
    def asym ( self, value ) :
        self.set_value ( self.__asym , value ) 
        
    @property
    def kappa ( self ) :
        """'kappa'-parameter for Apollonios2 function (same as 'asym'"""
        return self.__asym
    @kappa.setter
    def kappa ( self, value ) :
        self.set_value ( self.__asym , value ) 

    @property
    def sigmaL( self ) :
        """(left)'sigma'-parameter for Bifurcated Student-T function"""
        return self.__sigmaL

    @property
    def sigmaR( self ) :
        """(right)'sigma'-parameter for Bifurcated Student-T function"""
        return self.__sigmaR
   
models.append ( BifurcatedStudentT_pdf )      


# =============================================================================
## @class PearsonIV_pdf
#  Pearson Type IV pdf
#  @see Ostap::Models::PearsonIV 
#  @see Ostap::Math::PearsonIV 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2022-07-10
class PearsonIV_pdf(PEAK) :
    """Pearson Type IV distribution:

    - see Ostap.Models.PearsonIV
    - see Ostap.Math.PearsonIV
    """
    def __init__ ( self         ,
                   name         ,
                   xvar         ,
                   mu           , 
                   varsigma     ,
                   n        = 2 , 
                   kappa    = 0 ) :
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self                   ,
                              name        = name     ,
                              xvar        = xvar     ,
                              mean        = mu       ,
                              sigma       = varsigma ,
                              mean_name   = 'mu_%s'               % name ,
                              mean_title  = '#mu_{PIV}(%s)'       % name ,
                              sigma_name  = 'varsigma_%s'         % name ,
                              sigma_title = '#varsigma_{PIV}(%s)' % name )
        
        ## location parameter 
        self.__mu       = self.mean
        ## width parameter 
        self.__varsigma = self.sigma
        
        ## n parameter 
        self.__n     = self.make_var ( n                    ,
                                       'n_%s'        % name ,
                                       'n_{PIV}(%s)' % name ,
                                       None , 1 , 1.e-6 , 200 ) 
        
        ## asymmetry parameter 
        self.__kappa = self.make_var ( kappa                     , 
                                       'kappa_%s'         % name ,
                                       '#kappa_{PIV}(%s)' % name ,
                                       None , 0 , -200 , 200 ) 
        
        ## ditto 
        self.__nu = self.__kappa

        #  make the final PDF 
        self.pdf = Ostap.Models.PearsonIV (
            self.roo_name ( 'pIV_' ) , 
            "Pearson Type IV %s" % self.name ,
            self.xvar   ,
            self.mu     ,
            self.sigma  ,
            self.n      ,
            self.kappa  )
        
        ## save the configuration
        self.config = {
            'name'      : self.name     ,
            'xvar'      : self.xvar     ,
            'mu'        : self.mu       ,
            'varsigma'  : self.varsigma ,
            'n'         : self.n        ,
            'kappa'     : self.kappa    ,
            }

    @property     
    def mu ( self ) :
        """'mu'-parameter (location) for Pearon Type IV distribution (same as 'mean')"""
        return self.__mu
    @mu.setter
    def mu ( self, value ) :
        self.set_value ( self.__mu , value ) 

    @property     
    def varsigma ( self ) :
        """'varsigma'-parameter for Pearon Type IV distribution (same as 'sigma')"""
        return self.__varsigma
    @varsigma.setter
    def varsigma ( self, value ) :
        self.set_value ( self.__varsigma , value ) 

    @property
    def n  ( self ) :
        """'n'-parameter for Pearson Type IV distribution"""
        return self.__n
    @n.setter
    def n ( self, value ) :
        self.set_value ( self.__n , value ) 

    @property
    def kappa ( self ) :
        """'kappa'-parameter (asymmetry) for Pearson Type IV distribution"""
        return self.__kappa
    @kappa.setter
    def kappa ( self, value ) :
        self.set_value ( self.__kappa , value ) 

    @property     
    def nu ( self ) :
        """'nu'-parameter (asymmetry) for Pearon Type IV distribution (same as 'kappa')"""
        return self.__nu
    @nu.setter
    def nu ( self, value ) :
        self.set_value ( self.__nu , value ) 

models.append ( PearsonIV_pdf )      

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
#    - \f$\epsilon\f$ parameter controls the skewness 
#    - \f$\delta\f$   parameter controls the kurtosis 
#   Normal distribution reappears as \f$\epsilon=0\f$ and \f$\delta=1\f$ 
#  The heavy tails correspond to \f$\delta<1\f$, 
#  light tails correpond to \f$\delta>1\f$
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2014-08-02
class SinhAsinh_pdf(PEAK) :
    """SinhAsinh-function: 
    see Jones, M. C.; Pewsey, A. (2009).
    'Sinh-arcsinh distributions'. Biometrika 96 (4): 761. 
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
                   xvar             ,
                   mean      = None , ## mu 
                   sigma     = None ,
                   epsilon   = 0    ,
                   delta     = 1    ) :

        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name  , xvar , mean , sigma )
        
        ##
        self.__mu      = self.mean
        self.__epsilon = self.make_var ( epsilon ,
                                         'epsilon_%s'   % name ,
                                         '#epsilon(%s)' % name ,
                                         None , 0 , -1000 , +1000 )
        self.__delta   = self.make_var ( delta ,
                                         'delta_%s'   % name ,
                                         '#delta(%s)' % name , 
                                         None , 1 , 1.e-6 , 1000 )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.SinhAsinh (
            self.roo_name ( 'sinasinh_' ) , 
            "Sinh-Asinh %s" % self.name ,
            self.xvar      ,
            self.mean      ,
            self.sigma     ,
            self.epsilon   ,
            self.delta     )
        
        ## save the configuration
        self.config = {
            'name'      : self.name    ,
            'xvar'      : self.xvar    ,
            'mean'      : self.mean    ,
            'sigma'     : self.sigma   ,
            'epsilon'   : self.epsilon ,
            'delta'     : self.delta   ,
            }

    @property
    def epsilon( self ) :
        """'epsilon'-parameter for Sinh-Asinh function"""
        return self.__epsilon
    @epsilon.setter
    def epsilon ( self, value ) :
        self.set_value ( self.__epsilon , value ) 

    @property
    def delta ( self ) :
        """'delta-parameter' for Sinh-Asinh function"""
        return self.__delta
    @delta.setter
    def delta ( self, value ) :
        self.set_value ( self.__delta , value ) 

    @property
    def mu ( self ) :
        """'mu'-parameter (location) for Sinh-Asinh function (same as 'mean')"""
        return self.__mu
    @mu.setter
    def mu (  self , value ) :
        self.set_value ( self.__mu , value ) 
                
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
#  \f$ z = \gamma + \delta \sinh^{-1}\dfrac{ x - \xi}{\lambda} \f$
#  follows normal distribtion with mean 0 and sigma 1.
#
#  Note:
#  Symmetric case of JonhsonSU distribution is 
#  recovered by \f$\delta\rightarrow0\f$ for 
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
class JohnsonSU_pdf(PEAK) :
    """Johnson, N. L. (1949) 
    Systems of frequency curves generated by methods of translation
    Biometrika 36: 149–176 JSTOR 2332539
    
    https://en.wikipedia.org/wiki/Johnson_SU_distribution
    
    When variable x follows Johnson-SU distribution, 
    the variable  z = gamma +  delta  sinh^{-1} ( x - xi)/ lambda 
    follows normal distribtion with mean 0 and sigma 1.
    
    Note:
    Symmetric case of JonhsonSU distribution is 
    recovered by delta -> 0  for  'sinh-asinh' distribution,
    see Jones, M. C.; Pewsey, A. (2009). 
    'Sinh-arcsinh distributions'. Biometrika 96 (4): 761. 
    doi:10.1093/biomet/asp053
    http://oro.open.ac.uk/22510
    
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   xi        = None ,   ## related to mean 
                   lambd     = None ,   ## related to sigma 
                   delta     = 1    ,
                   gamma     = 0    ) :  

        #
        ## initialize the base
        #
        
        PEAK.__init__  ( self    , name , xvar ,
                              mean        = xi                   ,
                              sigma       = lambd                ,
                              mean_name   = 'xi_%s'       % name ,
                              mean_title  = '#xi(%s)'     % name ,
                              sigma_name  = 'lambda_%s'   % name ,
                              sigma_title = '#lambda(%s)' % name )
        
        self.__xi    = self.mean
        self.__lambd = self.sigma
        self.lambd.setMax ( self.lambd.getMax() * 10 ) ## adjust it! 
        self.__delta   = self.make_var ( delta                 ,
                                         'delta_%s'     % name ,
                                         '#delta(%s)'   % name , 
                                         None , 1 , 1.e-6 , 1000 )
        self.__gamma   = self.make_var ( gamma               ,
                                         'gamma_%s'   % name ,
                                         '#gamma(%s)' % name , 
                                         None , 0 , -1000 , +1000 )
        
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.JohnsonSU (
            self.roo_name ( 'jsu_' ) , 
            "Johnson's SU %s" % self.name ,
            self.xvar      ,
            self.xi        ,
            self.lambd     ,
            self.delta     ,
            self.gamma     )

        ## save the configuration
        self.config = {
            'name'      : self.name    ,
            'xvar'      : self.xvar    ,
            'xi'        : self.xi      ,
            'lambd'     : self.lambd   ,
            'delta'     : self.delta   ,
            'gamma'     : self.gamma   ,
            }

    @property
    def delta ( self ) :
        """'delta'-parameter for Johnson-SU function"""
        return self.__delta
    @delta.setter
    def delta ( self, value ) :
        self.set_value ( self.__delta , value ) 

    @property
    def gamma ( self ) :
        """'gamma'-parameter for Johnson-SU function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
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

models.append ( JohnsonSU_pdf )      
# =============================================================================
## @class Atlas_pdf
#  Modified gaussian with exponential tails
#  \f{displaymath} f(x) \propto \exp(-\dfrac{\delta x^{1+\dfrac{1}{1+\delta x/2}}}{2})\f},
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
class Atlas_pdf(PEAK) :
    r"""Modified gaussian with exponential tails
    \f$  f(x) \propto \exp( -frac{\delta x^{1+\dfrac{1}{1+\delta x/2}}}{2})\f$,
    where \f$\delta x = \left| x - \mu \right|/\sigma\f$
    Function is taken from http://arxiv.org/abs/arXiv:1507.07099    
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mean      = None ,   ## related to mean 
                   sigma     = None ) : ## related to sigma 

        #
        ## initialize the base
        #        
        PEAK.__init__  ( self , name , xvar , mean, sigma  )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Atlas (
            self.roo_name ( 'atlas_' ) , 
            "ATLAS/ZEUS %s" % self.name ,
            self.xvar      ,
            self.mean      ,
            self.sigma     )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            }

models.append ( Atlas_pdf )      

# =============================================================================
## @class Slash_pdf
#  Symmetric function with very heavy tails.
#  @see https://en.wikipedia.org/wiki/Slash_distribution
#  The tails  are so heavy that moments does not exists
#  @see Ostap::Math::Slash
#  @see Ostap::Models::Slash
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-02-004
class Slash_pdf(PEAK) :
    """Symmetric function with very heavy tails.
    - see https://en.wikipedia.org/wiki/Slash_distribution
    The tails  are so heavy that moments does not exists
    - see Ostap::Math::Slash
    - see Ostap::Models::Slash
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mean    = None ,   ## related to mean 
                   scale   = None ) : ## related to scale

        ## initialize the base
        PEAK.__init__  ( self , name , xvar  ,
                              mean        = mean  ,
                              sigma       = scale ,
                              sigma_name  = 'scale_%s'   % name ,
                              sigma_title = '#scale(%s)' % name )
        
        self.__scale = self.sigma

        ## finally build pdf
        self.pdf = Ostap.Models.Slash (
            self.roo_name ( 'slash_' ) , 
            "Slash %s" % self.name ,
            self.xvar      ,
            self.mean      ,
            self.scale     )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'scale'     : self.scale ,
            }
    @property
    def mu ( self ) :
        """'mu'  - location parameter, the same as 'mean' or 'location'"""
        return self.mean
    @property
    def location ( self ) :
        """'location'  - location parameter, the same as 'mean' or 'mu'"""
        return self.mean
    @property
    def scale ( self ) :
        """'scale'  - scale parameter, the same as 'sigma'"""
        return self.__scale
    @scale.setter
    def scale ( self , value ) :
        self.set_value ( self.__scale , value ) 
    
models.append ( Slash_pdf )      


# =============================================================================
## @class AsymmetricLaplace_pdf
#  Asymmetric version of Laplace distribution
#  @see https://en.wikipedia.org/wiki/Asymmetric_Laplace_distribution
#  @see Ostap::Math::AsymmetricLaplace
#  @see Ostap::Models::Laplace
#  \f$  f(x) \propto \exp ( \pm \dfrac{x-\mu}{ \lambda_{L,R}} ) \f$ 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-02-004
class AsymmetricLaplace_pdf(PEAK) :
    r"""Asymmetric version of Laplace distribution:
    
    f(x) \propto \exp ( \pm \dfrac { x- \mu } { \lambda_{ L , R } } )
    
    - see https://en.wikipedia.org/wiki/Asymmetric_Laplace_distribution
    - see Ostap::Math::AsymmetricLaplace
    - see Ostap::Models::Laplace
    
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mean      = None ,   ## related to mean
                   slope     = None ,
                   asymmetry = 0    ) : ## 0 corresponds to symmetric laplace 

        ## initialize the base
        PEAK.__init__  ( self , name , xvar ,
                              mean        = mean  ,
                              sigma       = slope ,
                              sigma_name  = 'slope_%s'   % name ,
                              sigma_title = '#slope(%s)' % name )
        
        self.__asym = self.make_var ( asymmetry               ,
                                      'asym_%s'        % name ,
                                      '#asym_{AL}(%s)' % name ,
                                      None , 0 , -1 , 1  ) 
        
        self.__slope = self.sigma

        ## constreuct left and right lambdas
        self.__lambdaL , self.__lambdaR = self.vars_from_asymmetry (
            self.slope                                    , ## mean/average sigma 
            self.asym                                     , ## asymmetry parametet
            v1name  =  self.roo_name ( 'lambdaL' , self.name ) ,
            v2name  =  self.roo_name ( 'lambdaR' , self.name ) ,
            v1title = '#lambda_{L}: slope #times (1+#kappa)'        , 
            v2title = '#lambda_{R}: slope #times (1-#kappa)'        )

        ## finally build pdf
        self.pdf = Ostap.Models.AsymmetricLaplace (
            self.roo_name ( 'alaplace_' ) , 
            "Asymmetrical Laplace  %s" % self.name ,
            self.xvar      ,
            self.mean      ,
            self.lambdaL   ,
            self.lambdaR   )                   
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'slope'     : self.slope ,
            'asymmetry' : self.asym  
            }
        
    @property
    def mu ( self ) :
        """'mu'  - location parameter, the same as 'mean' or 'location'"""
        return self.mean
    @property
    def location ( self ) :
        """'location'  - location parameter, the same as 'mean' or 'mu'"""
        return self.mean

    @property
    def slope ( self ) :
        """'slope'-parameter the mean exponential slope,  the same as 'sigma'"""
        return self.__slope
    @slope.setter
    def slope ( self, value ) :
        self.set_value ( self.__slope , value ) 
    
    @property
    def asym ( self ) :
        """'asymmetry'-parameter for Asymmetric Laplace"""
        return self.__asym
    @asym.setter
    def asym ( self, value ) :
        self.set_value ( self.__asym , value ) 

    @property
    def kappa ( self ) :
        """'kappa'-parameter for Apollonios2 function (same as 'asym'"""
        return self.__asym
    @kappa.setter
    def kappa ( self, value ) :
        self.set_value ( self.__asym , value ) 

    @property
    def lambdaL ( self ) :
        """(left)'lambda'-parameter (exponential slope) for Asymmetric Laplace"""
        return self.__lambdaL
    
    @property
    def lambdaR ( self ) :
        """(right)'lambda'-parameter (exponential slope) for Asymmetric Laplace"""
        return self.__lambdaR
    

models.append ( AsymmetricLaplace_pdf )      
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
#  \f{displaymath} f(x,\mu,\sigma)
#   \propto \dfrac{1}{2} \mathrm{sech}
#   \left( \dfrac{\pi}{2}\dfrac{x-\mu}{\sigma}\right) \f}    
#   @see https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
#
#  @see Ostap::Math::Sech
#  @see Ostap::Models::Sech
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-04-025
class Sech_pdf(PEAK) :
    r"""Hyperbolic secant distribution or 'inverse-cosh' distribution
    
    The hyperbolic secant distribution shares many properties with the 
    standard normal distribution: 
    - it is symmetric with unit variance and zero mean, median and mode
    -its pdf is proportional to its characteristic function. 
     
    However, the hyperbolic secant distribution is leptokurtic; 
    that is, it has a more acute peak near its mean, and heavier tails, 
    compared with the standard normal distribution.
    
    f(x,\mu,\sigma) \propto \dfrac{1}{2} sech ( \dfrac{\pi}{2}\dfrac{x-\mu}{\sigma} )
    
    @see https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mean      = None ,   ## related to mean 
                   sigma     = None ) : ## related to sigma 

        #
        ## initialize the base
        #        
        PEAK.__init__  ( self , name , xvar , mean , sigma  )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Sech (
            self.roo_name ( 'sech_' ) , 
            "Sech %s" % self.name ,
            self.xvar      ,
            self.mean      ,
            self.sigma     ) 
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            }
    
models.append ( Sech_pdf )      

#==============================================================================
## @class Losev_pdf
#  Asymmetric variant of hyperbolic secant distribution
#   \f[ f(x;\mu,\alpha,\beta) \equiv 
#      \frac{A}{\mathrm{e}^{-\left|\alpha\right| (x-\mu)} + 
#                           \mathrm{e}^{\left|\beta\right|(x-mu)}}, \f]
#  where \f$ A = \frac{\left|\alpha\right|+\left|\beta\right|}{\pi}
#  \sin \frac{\pi\left| \beta\right| }{\left|\alpha\right|+\left|\beta\right|}\f$ 
#   - Leptokurtic distribution with exponential tails 
#   @see Losev, A., "A new lineshape for fitting x‐ray photoelectron peaks", 
#           Surf. Interface Anal., 14: 845-849. doi:10.1002/sia.740141207
#   @see  https://doi.org/10.1002/sia.740141207
#   @see  https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
#   @see Ostap::Models::Losev 
#   @see Ostap::Math::Losev 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2019-11-026
class Losev_pdf(PEAKMEAN) :
    """ Asymmetric variant of hyperbolic secant distribution
    - Leptokurtic distribution with exponential tails
    see Losev, A., 'A new lineshape for fitting x‐ray photoelectron peaks', 
    Surf. Interface Anal., 14: 845-849. doi:10.1002/sia.740141207
    see https://doi.org/10.1002/sia.740141207
    see https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
    see Ostap::Models::Losev 
    see Ostap::Math::Losev 
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mean      = None ,   ## related to mean/location  
                   alpha     = None ,   ## left tail 
                   beta      = None ) : ## rigth tail
        
        #
        ## initialize the base
        #
        
        PEAKMEAN.__init__  ( self , name , xvar , 
                             mean       = mean  ,
                             mean_name  = 'mu_%s'   % name ,
                             mean_title = '#mu(%s)' % name )


        ## left tail 
        self.__alpha = self.make_var ( alpha ,
                                       'alpha_%s'           % name ,
                                       '#alpha_{Losev}(%s)' % name , 
                                       None , 1.0 , 1.e-3 , 1000 )

        ## right tail 
        self.__beta  = self.make_var ( beta   ,
                                       'beta_%s'            % name ,
                                       '#beta_{Losev}(%s)'  % name , 
                                       None , 1.0 , 1.e-3 , 1000 )
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Losev (
            self.roo_name ( 'losev_' ) , 
            "Losev %s" % self.name ,
            self.xvar      ,
            self.mean      ,
            self.alpha     ,
            self.beta      )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'alpha'     : self.alpha ,
            'beta'      : self.beta  ,
            }
        
    @property
    def mu ( self ) :
        """'mu'- location parameter for Losev distribution (same as 'mean')
        """
        return self.mean
    @mu.setter 
    def mu ( self , value ) :
        self.mean = value 

    @property
    def alpha ( self ) :
        """`alpha'- parameter for Losev distribution (left tail)
        """
        return self.__alpha
    @alpha.setter 
    def alpha ( self , value ) :
        self.set_value ( self.__alpha , value ) 

    @property
    def beta ( self ) :
        """`beta'- parameter for Losev distribution (right tail)
        """
        return self.__beta
    @beta.setter 
    def beta ( self , value ) :
        self.set_value ( self.__beta , value ) 
    
models.append ( Losev_pdf )      

# =============================================================================
## @class Logistic_pdf
#  Logistic, aka "sech-square" PDF
#  \f$ f(x;\mu;s) = \dfrac{1}{4s}sech^2\left(\dfrac{x-\mu}{2s}\right)\f$, 
#   where
#   \f$  s = \sigma \dfrac{\sqrt{3}}{\pi}\f$
#  @see https://en.wikipedia.org/wiki/Logistic_distribution
#  @see Ostap::Math::Logistic
#  @see Ostap::Models::Logistic
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2016-06-14
class Logistic_pdf(PEAK) :
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
                   mean      = None ,   ## related to mean 
                   sigma     = None ) : ## related to sigma 

        #
        ## initialize the base
        #
        
        PEAK.__init__  ( self , name , xvar , mean , sigma  )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Logistic (
            self.roo_name ( 'logistic_' ) , 
            "Logistic %s" % self.name ,
            self.xvar      ,
            self.mean      ,
            self.sigma     ) 
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            }

    
models.append ( Logistic_pdf )      

# =============================================================================
## @class RaisingCosine_pdf 
#  "Raising cosine" distribution
#  \f{displaymath} f(x,\mu,s) = \dfrac{1}{2s} \left( 1 + \cos \pi y \right)  \f}, 
#  where \f$  y \equiv \dfrac{x-\mu}{s}\f$ 
#  @see https://en.wikipedia.org/wiki/Raised_cosine_distribution
#  @see Ostap::Models::RaisingCosine 
#  @see Ostap::Math::RaisingCosine 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-02-27
class RaisingCosine_pdf(PEAK) :
    r"""'Raising cosine' distribution
    (x,\mu,s) = \dfrac{1}{2s}   \left( 1   +\cos \pi y \right), 
    where y  \equiv = \dfrac{x-\mu}{s} 
    - see https://en.wikipedia.org/wiki/Raised_cosine_distribution
    - see Ostap::Models::RaisnCosine 
    - see Ostap::Math::RaisnCosine
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mean      = None ,   ## related to mean 
                   scale     = None ) : ## related to sigma 

        #
        ## initialize the base
        #
        
        PEAK.__init__  ( self , name , xvar ,
                         mean        = mean  ,
                         sigma       = scale ,
                         sigma_name  = 'scale_%s'   % name ,
                         sigma_title = '#scale(%s)' % name )

        
        self.__scale = self.sigma

        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.RaisingCosine (
            self.roo_name ( 'rcos_' ) , 
            "Raising Cosine %s" % self.name ,
            self.xvar      ,
            self.mean      ,
            self.scale     ) 
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'scale'     : self.scale ,
            }
        
    @property
    def scale ( self ) :
        """'scale'-parameter, the same as 'sigma'"""
        return self.__scale
    @scale.setter
    def scale ( self, value ) :
        self.set_value ( self.__scale , value ) 
    
models.append ( RaisingCosine_pdf )      

# =============================================================================
## @class QGaussian_pdf 
#  q-Gaussian distribution:
#  \f{displaymath} f(x) = \dfrac{ \sqrt{\beta}}{C_q} e_q (-\beta (x-\mu)^2)\f}, 
#  where \f$ e_q (x) = \left( 1 + (1-q)x\right)^{\dfrac{1}{1-q}}\f$ 
#  @see https://en.wikipedia.org/wiki/Q-Gaussian_distribution
#  If is equal to 
#   - scaled version of Student' t-distribution for 1<q<3
#   - Gaussian distribution for q = 1 
#   - has finite  support for q<1 
#  @see Ostap::Math::QGaussian
#  Here we use \f$ \beta = \dfrac{1}{2\sigma^2}\f$
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-02-27
class QGaussian_pdf(PEAK) :
    r"""q-Gaussian distribution:
    \f$ f(x) = \dfrac{ \sqrt{\beta}}{C_q} e_q (-\beta (x-\mu)^2)$, 
    where  \f$ e_q (x) = \left( 1 + (1-q)x\right)^{\dfrac{1}{1-q}}\f$ 
    - see https://en.wikipedia.org/wiki/Q-Gaussian_distribution
    If is equal to 
    - scaled version of Student' t-distribution for 1<q<3
    - Gaussian distribution for q = 1 
    - has finite  support for q<1 
    Note: Here we use \f$ \beta = \dfrac{1}{2\sigma^2}\f$
    - see Ostap::Math::QGaussian
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mean      = None ,   ## related to mean 
                   q         = 1    ,   ## q-value 
                   scale     = 1    ) : ## related to sigma 

        #
        ## initialize the base
        #        
        PEAK.__init__  ( self , name , xvar  , 
                              mean        = mean  ,
                              sigma       = scale,
                              sigma_name  = 'scale_%s'  % name ,
                              sigma_title = 'scale(%s)' % name )
        

        self.__scale = self.sigma
        
        ## Q 
        self.__q = self.make_var ( q               ,
                                   'q_%s'   % name ,
                                   '#q(%s)' % name ,
                                   None , 1 , -1000 , 3-1.e-6 ) 
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.QGaussian (
            self.roo_name ( 'qgauss_' ) , 
            "q-Gaussian %s" % self.name ,
            self.xvar      ,
            self.mean      ,
            self.q         ,
            self.scale     ) 
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'q'         : self.q     ,
            'scale'     : self.scale ,
            }
        
    @property
    def scale ( self ) :
        """'scale'-parameter, the same as 'sigma'"""
        return self.__scale
    @scale.setter
    def scale ( self, value ) :
        self.set_value ( self.__scale , value ) 
        
    @property
    def q ( self ) :
        """'q'-parameter"""
        return self.__q
    @q.setter
    def q ( self, value ) :
        self.set_value ( self.__q , value ) 

models.append ( QGaussian_pdf )      



# =============================================================================
## @class Hyperbolic_pdf 
#  Hyperbolic disribtion
#  @see  https://en.wikipedia.org/wiki/Hyperbolic_distribution
#  @see  Barndorff-Nielsen, Ole, 
#    "Exponentially decreasing distributions for the logarithm of particle size". 
#     Proceedings of the Royal Society of London. Series A,
#     Mathematical and Physical Sciences. 
#     The Royal Society. 353 (1674): 401–409
#     doi:10.1098/rspa.1977.0041. JSTOR 79167.
#
#  \f[  f(x;\mu, \beta, \delta, \gamma) = 
#  \frac{\gamma}{2\alpha\delta K_1(\delta \gamma)}
#  \mathrm{e}^{ -\sqrt{ \alpha^2\delta^2 + \alpha^2 (x-\mu)^2 } + \beta ( x - \mu)}
#  \f]
#  where 
#  - \f$ \alpha^2 = \beta^2\f + \gamma^2$
#  - \f$ K_1\f$ is a modified Bessel function of the second kind 
#  
# In the code we adopt parameterisation in terms of
#  - location parameter \f$\mu\f$
#  - parameter               \f$\sigma \gt  0 \f$, related to the width;
#  - dimensionless parameter \f$\kappa\f$,         related to the asymmetry;
#  - dimensionless parameter \f$\zeta   \ge 0 \f$, related to the kurtosis 
#
# The parameters are defined as:
# \f[\begin{array}{lcl}
#     \sigma^2 & \equiv & \gamma^{-2} \zeta \frac{K_2(\zeta)}{\zetaK_1(zeta)} \\
#     \kappa   & \equiv & \frac{\beta}{\sigma} \                   \
#     \zeta\equiv\delta \gamma \end{array} \f]
# - For \f$ \beta=0 (\kappa=0)\f$,  \f$\sigma^2\f$ is a variance of the distribution.
# - Large values of \f$\zeta\f$ distribtionhas small kurtosis 
# - For small \f$ \zeta \f$ distribution shows kurtosis of 3 
#
# The inverse transformation is:
# \f[ \begin{array}{lcl}
#     \beta    & = & \frac{\kappa}{\sigma}            \\
#     \delta   & = & \frac{\zeta}{\gamma}             \\
#     \gamma   & = & \frac{\sqrt{A^*(\zeta)}}{\sigma} \\
#     \alpha   & = & \sqrt { \beta^2 + \gamma^2} \end{array} \f]
# where \f$ A^{*}(\zeta) = \frac{\zeta K^*_2(\zeta)}{K^*_1(zeta)} \f$. 
# It is largely inspired by NIM A764 (2014) 150, arXiv:1312.5000, 
# but has much better properties when \f$ \zeta \rigtarrow 0 \f$ 
#  @see D. Martinez Santos and F. Dupertuis,
#          "Mass distributions marginalized over per-event errors",
#          NIM A764 (2014) 150, arXiv:1312.5000
#          DOI: 10.1016/j.nima.2014.06.081",
#
#  The final form of the distribution is 
#  \f[  f(x;\mu,\sigma,\zeta,\kappa) = 
#      \frac{ A^*(\zeta) } { 2 \sigma \sqrt{\kappa^2+A^*(\zeta)} \zeta K^*_1(\zeta) } 
#      \mathrm{e}^{\zeta - \sqrt{ (\kappa^2+A(\zeta))  \left( \frac{\zeta^2}{A(\zeta)}  +  
#      \left( \frac{x-\mu}{\sigma}\right)^2  \right) } } 
#    \f]
#  where \f$ K^*_n(x)\f$ is a scaled modified Bessel functon to th eseodn kind 
#   \f$ K^*_n(x) = \mathrm{e}^{x}K_1(x) \f$ 
#
#  In all expressions \f$ \left| \sigma \right|\f$ and 
#  \f$ \left| \zeta \right|\f$ are used instead of \f$\sigma\f$ and \f$\zeta\f$ 
# 
#  @see Ostap::Models::Hyperbolic
#  @see Ostap::Math::Hyperbolic
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-02-27
class Hyperbolic_pdf(PEAK) :
    r"""Hyperbolic distribution
    - see  https://en.wikipedia.org/wiki/Hyperbolic_distribution
    - see  Barndorff-Nielsen, Ole, 
    #    `Exponentially decreasing distributions for the logarithm of particle size'. 
    #     Proceedings of the Royal Society of London. Series A,
    #     Mathematical and Physical Sciences. 
    #     The Royal Society. 353 (1674): 401–409
    #     doi:10.1098/rspa.1977.0041. JSTOR 79167.
    - see Ostap::Math::Hyperbolic
    - see Ostap::Models::Hyperbolic
    
    Parameters are different from 'canonical'
    - mu     : related to location   (equal to mean/mode for kappa=0) 
    - sigma  : relates to width      (equal to RMS           for kappa=0)
    - zeta   : related to kurtosis   (kurtosis varies from 3 to 0 when zeta varies from 0 to infinity)
    - kappa  : related  to asymmetry 
    
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mu        = None ,   ## related to mean
                   sigma     = 1    ,   ## relatd  to width  
                   zeta      = 0    ,   ## related to kurtosis 
                   kappa     = 0    ) : ## related to asymmetry
        
        #
        ## initialize the base
        #        
        PEAK.__init__  ( self , name , xvar               , 
                              mean        = mu                 ,
                              sigma       = sigma              ,
                              mean_name   = 'mu_%s'     % name ,
                              mean_title  = '#mu(%s)'   % name )
        
        
        self.__mu    = self.mean 
        
        ## Zeta
        self.__zeta  = self.make_var ( zeta                ,
                                       'zeta_%s'    % name ,
                                       '#zeta(%s)'  % name ,
                                       None ,  1 , 0  , 100 ) 
        ## kappa  
        self.__kappa = self.make_var ( kappa               ,
                                       'kappa_%s'   % name ,
                                       '#kappa(%s)' % name ,
                                       None ,  0 , -50 , 50 ) 
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Hyperbolic (
            self.roo_name ( 'hyperbolic_' ) , 
            "Hyperbolic %s" % self.name ,
            self.xvar      ,
            self.mu        ,
            self.sigma     ,
            self.zeta      ,
            self.kappa     )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mu'        : self.mu    ,
            'sigma'     : self.sigma ,
            'zeta'      : self.zeta  ,
            'kappa'     : self.kappa }
        
    @property
    def mu ( self ) :
        """'mu' : location parameter, same as 'mean')"""
        return self.__mu
    @mu.setter
    def mu ( self , value ) :    
        self.set_value ( self.__mu , value )

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
        self.set_value ( self.__kappa , value )

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

        
models.append ( Hyperbolic_pdf )      



# =============================================================================
## @class GenHyperbolic_pdf 
#  Generalized Hyperbolic distribution
#  @see https://en.wikipedia.org/wiki/Generalised_hyperbolic_distribution 
#  
#  \f[ f(x;\lambda, \alpha,\beta,\gamma,\delta,\mu)= 
#  \frac{  (\gamma/\delta)^{\lambda} }{ \sqrt{2\pi} K_{\lambda}(\delta\gamma) }
#   \mathrm{e}^{\beta (x-\mu)}
#   \frac{  K_{\lambda -1/2} ( \alpha \sqrt{ \delta^2 + (x-\mu)^2} ) }
#   { (  \sqrt{ \delta^2 + (x-\mu)^{2} }  /\alpha)^{1/2-\lambda} } \f]
# where 
#  - $\alpha=\sqrt{\beta^2+\gamma^2}$
#
#  In the code we adopt parameterisation in terms of
#  - location parameter      \f$\mu\f$
#  - shape parameter         \f$\lambda\f$
#  - parameter               \f$\sigma \gt  0 \f$, related to the width;
#  - dimensionless parameter \f$\kappa\f$,         related to the asymmetry;
#  - dimensionless parameter \f$\zeta   \ge 0 \f$, related to the shape 
#  
# The parameters are defined as:
# \f[\begin{array}{lcl}
#     \sigma^2 & \equiv & \gamma^{-2} \zeta \frac{K_{\lambda+1}(\zeta)}{\zetaK_{\lambda}(zeta)} \\
#     \kappa   & \equiv & \frac{\beta}{\sigma} \\
#     \zeta    & \equiv & \delta \gamma \end{array} \f]
# - For \f$ \beta=0 (\kappa=0)\f$,  \f$\sigma^2\f$ is a variance of the distribution.
# - Large values of \f$\zeta\f$ distribtionhas small kurtosis 
# - For small \f$\zeta\f$ distribution shows kurtosis of 3 
#
# The inverse transformation is:
# \f[ \begin{array}{lcl}
#     \beta    & = & \frac{\kappa}{\sigma}                      \\
#     \delta   & = & \frac{\zeta}{\gamma}                       \\
#     \gamma   & = & \frac{\sqrt{A_{\lambda}^*(\zeta)}}{\sigma} \\
#     \alpha   & = & \sqrt { \beta^2 + \gamma^2} \end{array} \f]
#
# In general it has exponential tails for \f$ \lambda >0 \f$ and Gaussian core.
# For negative \f$ \lambda \f$ tails are more heavy..
#
# @see Ostap::Math::Hyperbolic
#  Useful subclasses 
#  - \f$ \lambda=1\f$ : Hyperbolic distribution  
#  - \f$ \lambda = -\frac{1}{2}\f$ : Normal Inverse Gaussian distribution
#  - \f$ \lambda=-\frac{n}{2}, \zeta\rightarrow+0\f$ : Student's t-distibution 
#  - \f$ \lambda \rightarrow \pm\infty, \kappa=0\f$ : Gaussian distribution 
#  - \f$ \zeta \rightarrow +\infty, \kappa=0\f$ : Gaussian distribution 
#
#  @see Ostap::Models::Hyperbolic
#  @see Ostap::Math::Hyperbolic
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2018-02-27
class GenHyperbolic_pdf(Hyperbolic_pdf) :
    r"""Generalised Hyperbolic distribution
    @see https://en.wikipedia.org/wiki/Generalised_hyperbolic_distribution 
    
    - see Ostap::Math::GenHyperbolic
    - see Ostap::Models::GenHyperbolic
    
    Parameters are different from 'canonical'
    - mu     : related to location   (equal to mean/mode for kappa=0) 
    - sigma  : relates to width      (equal to RMS           for kappa=0)
    - zeta   : related to kurtosis   (kurtosis varies from 3 to 0 when zeta varies from 0 to infinity)
    - kappa  : related  to asymmetry 
    - lambda : related to shape
    
    Useful subclasses 
    - \f$ \lambda = 1\f$ : Hyperbolic distributiobn  
    - \f$ \lambda = -\frac{1}{2}\f$ : Normal Inverse Gaussian distribution
    - \f$ \lambda=-\frac{n}{2}, \zeta\rightarrow+0\f$ : Stundent's t-distibution 
    - \f$ \lambda \rightarrow \pm\infty, \kappa=0\f$ : Gaussian distribution 
    - \f$ \zeta \rightarrow +\infty, \kappa=0\f$ : Gaussian distribution 

    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mu        = None ,    ## related to mean
                   sigma     = 1    ,    ## relatd  to width  
                   zeta      = 1    ,    ## related to shape 
                   kappa     = 0    ,    ## related to asymmetry
                   lambd     = 1    ) :  ## related to shape 
        
        #
        ## initialize the base
        #
        
        Hyperbolic_pdf.__init__  ( self ,
                                   name   = name  ,
                                   xvar   = xvar  , 
                                   mu     = mu    ,
                                   sigma  = sigma ,
                                   zeta   = zeta  ,
                                   kappa  = kappa ) 
        
        ## lambda 
        self.__lambda = self.make_var ( lambd               ,
                                        'lambda_%s'   % name ,
                                        '#lambda(%s)' % name ,
                                        True , -2 , -100 , 100 ) 
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.GenHyperbolic (
            self.roo_name ( 'genhyperbolic_' ) , 
            "GenHyperbolic %s" % self.name ,
            self.xvar      ,
            self.mu        ,
            self.sigma     ,
            self.zeta      ,
            self.kappa     ,
            self.lambd     )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mu'        : self.mu    ,
            'sigma'     : self.sigma ,
            'zeta'      : self.zeta  ,
            'kappa'     : self.kappa ,
            'lambd'     : self.lambd }
        
    @property
    def lambd ( self ) :
        """'lambd' : dimensionless parameter, related to shape """
        return self.__lambda
    @lambd.setter
    def lambd ( self , value ) :    
        self.set_value ( self.__lambd , value )

        
models.append ( GenHyperbolic_pdf )      

# =============================================================================
## @class Hypatia_pdf
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
class Hypatia_pdf(GenHyperbolic_pdf) :
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
    
    - see GenHyperbolic_pdf
    - see Ostap.Math.GenHyperbolic
    - see Ostap.Models.GenHyperbolic
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mu        = None ,   ## related to mean
                   sigma     =  1   ,   ## relatd  to width  
                   zeta      =  1   ,   ## related to shape 
                   kappa     =  0   ,   ## related to asymmetry
                   lambd     = -2   ,   ## related to shape 
                   sigma0    = None ,   ## width of the 'offset' Gaussian 
                   cnvpars   = {}   ) : ## convolution parameters 
        # 
        ## initialize the base
        #        
        GenHyperbolic_pdf.__init__  ( self ,
                                      name        = name  ,
                                      xvar        = xvar  , 
                                      mu          = mu    ,                                      
                                      sigma       = sigma ,
                                      zeta        = zeta  ,
                                      kappa       = kappa ,
                                      lambd       = lambd )


        self.__genhyp = self.pdf
        
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
            'mu'      : self.mu      ,
            'sigma'   : self.sigma   ,
            'zeta'    : self.zeta    ,
            'kappa'   : self.kappa   ,
            'lambd'   : self.lambd   ,
            'sigma0'  : self.sigma0  ,
            'cnvpars' : self.cnvpars 
            }
        
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
        """'sigma0' : width for the 'offset' Gaussian"""
        return self.__resolution.sigma
    @sigma0.setter
    def sigma0    ( self , value ) :
        self.__resolution.sigma = value
        
    @property
    def cnvpars ( self ) :
        """'cnvpars' : parameters for convolution"""
        return self.__cnvpars 


models.append ( Hypatia_pdf )      


# =============================================================================
## @class ExGauss_pdf
#  Exponentially modified Gaussian function, EMG
#  @see https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
#
#  It is a distibution for the varibale that is a 
#  sum (or difference for negative \f$ k\f$) 
#  of a Gaussian and exponential variables: \f$ X \sim Y + sign(k) Z \f$,  
#  where 
#  - \f$ Y \sim N(\mu,\sigma) \f$
#  - \f$ Z \sim  \frac{1}{k\sigma}\mathrm{e}^{-\frac{x}{k\sigma}} \f$ 
#  
#  For \f$ k=0\f$ one gets a Gaussian distrobution
#  - \f$ k>0\f$ corresponds to the rigth tail  
#  - \f$ kM0\f$ corresponds to the left tail  
#
# It can be considered as "single-tail" version of the Normal Laplace distribution:
#  - \f$ k = 0 \f$ corresponds to Gaussian distribution
#  - \f$ k > 0 \f$ corresponds to Normal Laplace \f$ NL(\mu,\sigma,0,k)\f$ 
#  - \f$ k < 0 \f$ corresponds to Normal Laplace \f$ NL(\mu,\sigma,\left|k\right|,0)\f$ 
#
#  @see Reed, W.J, "The Normal-Laplace Distribution and Its Relatives". 
#       In: Balakrishnan, N., Sarabia, J.M., Castillo, E. (eds) 
#       "Advances in Distribution Theory, Order Statistics, and Inference. 
#       Statistics for Industry and Technology". Birkhäuser Boston. 
#  @see https://doi.org/10.1007/0-8176-4487-3_4
#  @see Ostap::Math::NormalLaplace 
#  @see Ostap::Mdoels::NormalLaplace 
#  @see Ostap::Mdoels::ExGauss
class ExGauss_pdf(PEAK) :
    """ Exponentially modified Gaussian function, EMG
    - see https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
    
    It is a distibutiin for the varibale that is a 
    sum (or difference for negative k)  
    of a Gaussian and exponential variables.

    - k = 0 : Gaussian distribution
    - k>0  corresponds to the rigth tail  
    - k<0  corresponds to the left tail  
    
    - see Reed, W.J, "The Normal-Laplace Distribution and Its Relatives". 
    In: Balakrishnan, N., Sarabia, J.M., Castillo, E. (eds) 
    "Advances in Distribution Theory, Order Statistics, and Inference. 
    Statistics for Industry and Technology". Birkhäuser Boston.
    
    - see https://doi.org/10.1007/0-8176-4487-3_4
    - see Ostap::Math::NormalLaplace 
    - see Ostap::Mdoels::NormalLaplace 
    - see Ostap::Mdoels::ExGauss
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mu        = None ,   ## related to mean
                   varsigma  =  1   ,   ## relatd  to width
                   k         =  0   ) :
        # 
        ## initialize the base
        #        
        PEAK.__init__  ( self , name , xvar                  , 
                              mean        = mu                    ,
                              sigma       = varsigma              ,
                              mean_name   = 'mu_%s'        % name ,
                              mean_title  = '#mu(%s)'      % name ,
                              sigma_name  = 'varsigma_%s'  % name ,
                              sigma_title = 'varsigma(%s)' % name )
        
        self.__mu       = self.mean 
        self.__varsigma = self.sigma
                
        ## k 
        self.__k= self.make_var ( k              ,
                                  'k_%s'  % name ,
                                  'k(%s)' % name ,
                                  None , 0 , -200 , 200 ) 
        
        ## create PDF 
        self.pdf = Ostap.Models.ExGauss (
            self.roo_name ( 'exgauss_' ) , 
            "ExGauss %s" % self.name ,
            self.xvar      ,
            self.mu        ,
            self.varsigma  , 
            self.k         ) 
        

        self.config = {
            'name'     : self.name     ,
            'xvar'     : self.xvar     ,
            'mu'       : self.mu       ,
            'varsigma' : self.varsigma ,
            'k'        : self.k        }

    @property
    def varsigma    ( self ) :
        """'varsigma' : varsigma parameter for ExGauss function (same as 'sigma')"""
        return self.sigma
    @varsigma.setter
    def varsigma    ( self , value ) :
        self.set_value ( self.__varfsigma , value )

    @property
    def mu ( self ) :
        """'mu' : location parameter (same as 'mean')"""
        return self.__mu
    @mu.setter
    def mu ( self , value ) :    
        self.set_value ( self.__mu , value )

    @property 
    def k  ( self ) :
        """'k' :  (dimensioneless) k-parameter"""
        return self.__k
    @k.setter  
    def k ( self , value ) :
        self.set_value ( self.__k , value )
    
models.append ( ExGauss_pdf )      


# =============================================================================
## @class NormalLaplace_pdf 
#  Distribution for a sum of Gaussian and (asymmertric) Laplace variables 
#  It behaves line core Gaussian with exponential tails 
#  @see Wiliam J. Reed, "The Normal-Laplace Distribution Relatives", 
#  October, 2004
#  @see https://www.math.uvic.ca/faculty/reed/NL.draft.1.pdf
#  @see Reed, W.J, "The Normal-Laplace Distribution and Its Relatives". 
#       In: Balakrishnan, N., Sarabia, J.M., Castillo, E. (eds) 
#       "Advances in Distribution Theory, Order Statistics, and Inference. 
#       Statistics for Industry and Technology". Birkhäuser Boston. 
#  @see https://doi.org/10.1007/0-8176-4487-3_4
#
#   \f$ f(x; \mu, \sigma, k_L , k_R ) = 
#   \frac{1}{\sigma ( k_L + k_R) } 
#   \phi ( z ) \left( R ( \frac{1}{k_R} - z ) + 
#                     R ( \frac{1}{k_L} + z ) \right) 
#   \f$, where
#   - \f$ k_L,k_R \ge 0 \f$ 
#   - \f$ z = \frac{x-\mu}{\sigma} \f$ 
#   - \f$ \phi(z) \f$ is Gaussian PDF  
#   - \f$  R(x)   \f$ is Mill's ratio 
#  @see Ostap::Math::nills_normal 
#  @see Ostap::Math::nills_normal 
#  @see Ostap::Math::NormalLaplace
#
class NormalLaplace_pdf(PEAK) :
    """Distribution for a sum of Gaussian and (asymmertric) Laplace variables 
    It behaves line core Gaussian with exponential tails 
    - see Wiliam J. Reed, "The Normal-Laplace Distribution Relatives", October, 2004
    - see https://www.math.uvic.ca/faculty/reed/NL.draft.1.pdf
    - see Reed, W.J, "The Normal-Laplace Distribution and Its Relatives". 
    In: Balakrishnan, N., Sarabia, J.M., Castillo, E. (eds) 
    "Advances in Distribution Theory, Order Statistics, and Inference. 
    Statistics for Industry and Technology". Birkhäuser Boston. 
    - see https://doi.org/10.1007/0-8176-4487-3_4
    
    - see Ostap::Math::nills_normal 
    - see Ostap::Math::nills_normal 
    - see Ostap::Math::NormalLaplace
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   mu        = None ,   ## related to mean
                   varsigma  =  1   ,   ## relatd  to width
                   kL        =  0   ,
                   kR        =  0   ) :
        # 
        ## initialize the base
        #        
        PEAK.__init__  ( self , name , xvar                  , 
                              mean        = mu                    ,
                              sigma       = varsigma              ,
                              mean_name   = 'mu_%s'        % name ,
                              mean_title  = '#mu(%s)'      % name ,
                              sigma_name  = 'varsigma_%s'  % name ,
                              sigma_title = 'varsigma(%s)' % name )
        
        self.__mu       = self.mean 
        self.__varsigma = self.sigma
        
        
        ## kL 
        self.__kL = self.make_var ( kL                 ,
                                    'kL_%s'     % name ,
                                    'k_{L}(%s)' % name ,
                                    None , 0 , 0 , 200 ) 
        ## kR 
        self.__kR = self.make_var ( kR                 ,
                                    'kR_%s'     % name ,
                                    'k_{R}(%s)' % name ,
                                    None , 0 , 0 , 200 ) 
        
        ## create PDF 
        self.pdf = Ostap.Models.NormalLaplace (
            self.roo_name ( 'normlapl_' ) , 
            "NormalLaplace %s" % self.name ,
            self.xvar      ,
            self.mu        ,
            self.varsigma  , 
            self.kL        ,
            self.kR        ) 
        

        self.config = {
            'name'     : self.name     ,
            'xvar'     : self.xvar     ,
            'mu'       : self.mu       ,
            'varsigma' : self.varsigma ,
            'kL'       : self.kL       ,
            'kR'       : self.kR       }

    @property
    def varsigma    ( self ) :
        """'varsigma' : varsigma parameter for Normal Laplace function (same as 'sigma')"""
        return self.sigma
    @varsigma.setter
    def varsigma    ( self , value ) :
        self.set_value ( self.__varsigma , value )

    @property
    def mu ( self ) :
        """'mu' : location parameter (same as 'mean')"""
        return self.__mu
    @mu.setter
    def mu ( self , value ) :    
        self.set_value ( self.__mu , value )

    @property 
    def kL ( self ) :
        """'kL' :  (dimensioneless) kL-parameter"""
        return self.__kL
    @kL.setter  
    def kL ( self , value ) :
        self.set_value ( self.__kL , value )
        
    @property 
    def kR ( self ) :
        """'kR' :  (dimensioneless) kR-parameter"""
        return self.__kR
    @kR.setter  
    def kR ( self , value ) :
        self.set_value ( self.__kR , value )
        
models.append ( NormalLaplace_pdf )      

# =============================================================================
# @class Das_pdf
#  Simple gaussian function with exponential tails.
#  It corresponds to <code>ExpGaussExp</code> function, 
#  \f[ 
#   f (x ; \mu, \sigma, k_L, k_R ) = \frac{1}{\sqrt{2\pi}\sigma}
#   \left\{ \begin{array}[lcl}
#  \mathrm{e}^{  \frac{k_L^2}{2} + k_L\left(\frac{x-mu}{\sigma}\right) }
#   & \mathrm{for}  &  \left(\frac{x-\mu}{\sigma}\right) < -k_L \\   
#  \mathrm{e}^{ \frac{1}{s} \left( \frac{x-\mu}{\sigma}\right)^2}
#   & \mathrm{for}  &  -k_L < \left(\frac{x-\mu}{\sigma}\right) < k_R \\    
#  \mathrm{e}^{  \frac{k_R^2}{2} - k_R\left(\frac{x-mu}{\sigma}\right) }
#   & \mathrm{for}  &  \left(\frac{x-\mu}{\sigma}\right)> k_R   
#  \end{array} \right. \f]
#  - \f$ k_L \ge 0\f$
#  - \f$ k_R \ge 0\f$
#
#  @see Souvik Das, "A simple alternative to Crystall Ball fnuction"
#                   arXiv:1603.08591  [hep-ex]
#  @see https://arxiv.org/abs/1603.08591
#  @attention - the function is not normalized! 
#  Function was used in 
#  @see CMS collaboration, V.Khachatryan, 
#       "Search for resonant pair production of Higgs bosons decaying 
#        to two bottom quark\textendash{}antiquark pairs 
#        in proton-proton collisions at 8 TeV}",
#        Phys. Lett. B749 (2015) 560 
# @see https://arxiv.org/abs/1503.04114 
# @see https://doi.org/10.1016/j.physletb.2015.08.047 
# - Gaussian function is restored when \f$k_L,k_R \rigtharrow +\infty\f$
#
# @see Ostap::Math::Das
# @see Ostap::Models::Das
class Das_pdf(PEAK) :
    r"""Simple gaussian function with exponential tails.
    It corresponds to `ExpGaussExp` function from ref below
    
    - see Souvik Das, 'A simple alternative to Crystall Ball fnuction'
    arXiv:1603.08591  [hep-ex]
    
    - see https://arxiv.org/abs/1603.08591

    Function was used in 
    - see CMS collaboration, V.Khachatryan, 
    'Search for resonant pair production of Higgs bosons decaying 
    to two bottom quark\textendash{}antiquark pairs 
    in proton-proton collisions at 8 TeV',
    Phys. Lett. B749 (2015) 560
    
    - see https://arxiv.org/abs/1503.04114 
    - see https://doi.org/10.1016/j.physletb.2015.08.047
    
    Gaussian function is restored when \f$k_L,k_R \rigtharrow +\infty\f$
    
    - see Ostap.Math.Das
    - see Ostap.Models.Das
    """
    def __init__ ( self         ,
                   name         ,    ## the name of PDF
                   xvar         ,    ## observable
                   mu    = None ,    ## location parameter
                   sigma = None ,    ## width parameter
                   kL    = None ,    ## left tail parameter
                   kR    = None ,    ## right tail parameter
                   k     = None ,    ## tail parameer  (alternaike)
                   kappa = None ) :  ## tail asymmetry (alternative) 
        
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar ,
                              mean        = mu                ,
                              sigma       = sigma             ,
                              mean_name   = 'mu_%s'    % name ,
                              mean_title  = '#mu_(%s)' % name ) 
        
        
        if k is None and kappa is None :
            
            self.__kL    = self.make_var ( kL ,
                                           'kL_%s'    % name ,   
                                           '#k_L(%s)' % name ,
                                           None , 1 , 1.e-6 , 1000 )
            self.__kR    = self.make_var ( kR ,
                                           'kR_%s'    % name ,   
                                           '#k_R(%s)' % name ,
                                           None , 1 , 1.e-6 , 1000 )

            ## name of k-variable 
            k = Ostap.MoreRooFit.Addition2 ( self.roo_name ( 'k' , self.name ) ,
                                            '0.5*(kL+kR)' , 
                                             self.kL ,
                                             self.kR ,
                                             0.5     ,
                                             0.5     )
            
            self.__k     = self.make_var ( k , k.name , k.title , None , 1 , 1.e-6 , 1000 )
            
            kk = Ostap.MoreRooFit.Asymmetry ( self.roo_name ('kappa' , self.name ) ,
                                              '0.5*(kL-rR)/(kL+kR)' ,
                                              self.kL ,
                                              self.kR ,
                                              0.5     ) 
            
            self.__kappa = self.make_var ( kk , kk.name , kk.title , None , 0 , -1 , 1 )
 
        elif kL is None and kR is None :
            
            self.__k     = self.make_var ( k                    ,
                                           'k_%s'      % name   ,  
                                           'k(%s)'     % name   ,
                                           None , 1 , 1.e-6 , 1000 )
            self.__kappa = self.make_var ( kappa                ,
                                           'kappa_%s'   % name  ,  
                                           '#kappa(%s)' % name  ,
                                           None , 0 , -1 , 1 )

            
            kl = Ostap.MoreRooFit.Combination ( self.roo_name ( 'kL' , self.name ) ,
                                                'k*(1+kappa)' , 
                                                self.k        ,
                                                self.kappa    ,
                                                1 , 1 , 1     )
            kr = Ostap.MoreRooFit.Combination ( self.roo_name ( 'kR' , self.name ) ,
                                                'k*(1-kappa)' , 
                                                self.k        ,
                                                self.kappa    ,
                                                1 , 1 , -1    )
            
            self.__kL     = self.make_var ( kl , kl.name , kl.title , None , 1 , 1.e-6 , 1000 )
            self.__kR     = self.make_var ( kr , kr.name , kr.title , None , 1 , 1.e-6 , 1000 )

        else :
            
            raise TypeError( 'Invalid setting for k/kappa/kL/kR!' )


        ## build PDF
        self.pdf = Ostap.Models.Das (
            self.roo_name ( 'das_' ) ,
            'Das %s' % self.name     ,
            self.xvar                ,
            self.mu                  ,
            self.sigma               ,
            self.kL                  ,
            self.kR                  )

        ## save configuration
        self.config = {
            'name'  : self.name  ,
            'xvar'  : self.xvar  ,
            'mu'    : self.mu    ,
            'sigma' : self.sigma ,
            'kL'    : self.kL    if ( k  is None and kappa is None ) else None ,
            'kR'    : self.kR    if ( k  is None and kappa is None ) else None ,
            'k'     : self.k     if ( kL is None and kR    is None ) else None ,
            'kappa' : self.kappa if ( kL is None and kR    is None ) else None ,
            }

    @property
    def mu ( self ) :
        """'mu' : peak location, same as 'mean'
        """
        return self.mean
    @mu.setter 
    def mu ( self ,value ) :
        self.mean = value 
    
    @property
    def kL  ( self ) :
        """'kL' : left tail parameter
        """
        return self.__kL
    @kL.setter
    def kL ( self , value ) :
        self.setValue ( self.__kL , value )

    @property
    def kR  ( self ) :
        """'kR' : left tail parameter
        """
        return self.__kR
    @kR.setter
    def kR ( self , value ) :
        self.setValue ( self.__kR , value )

    @property
    def k   ( self ) :
        """'k' : (kL+kR)/2 parameter
        """
        return self.__k
    @k.setter
    def k ( self , value ) :
        self.setValue ( self.__k , value )

    @property
    def kappa ( self ) :
        """'k' : 0.5*(kL-kR)/(kL+kR) parameter
        """
        return self.__kappa
    @kappa.setter
    def kappa ( self , value ) :
        self.setValue ( self.__kappa , value )


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
class Voigt_pdf(PEAK) :
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
                   xvar             ,
                   m0        = None ,
                   sigma     = None ,
                   gamma     = None ) :
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar ,
                              mean        = m0                  ,
                              sigma       = sigma               ,
                              mean_name   = 'm0_%s'      % name ,
                              mean_title  = '#m_{0}(%s)' % name ) 
                         
        limits_gamma = ()
        if  self.xminmax() :
            mn , mx = self.xminmax() 
            dm = mx - mn
            limits_gamma = 1.e-5 * dm , dm
            
        self.__gamma  = self.make_var ( gamma                ,
                                        'gamma_%s'   % name  ,   
                                        '#gamma(%s)' % name  ,
                                        None , *limits_gamma )
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Voigt (
            self.roo_name ( 'voigt_' ) , 
            "Voigt %s" % self.name ,
            self.xvar   ,
            self.m0     ,
            self.gamma  ,
            self.sigma  )

        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'm0'        : self.m0    ,
            'sigma'     : self.sigma ,
            'gamma'     : self.gamma ,
            }
    
    @property
    def m0 ( self ) :
        """'m0' : m_0 parameter for Breit-Wigner function (alias for 'mean')"""
        return self.mean
    @m0.setter
    def m0 ( self , value ) :
        self.mean = value 

    @property
    def gamma ( self ) :
        """'gamma'-parameter for Voigt function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        assert 0 < value , "'gamma'-parameter must be positive"
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
#  @see https://doi.org/10.1107/S0021889800010219
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2016-06-15
class PseudoVoigt_pdf(Voigt_pdf) :
    """Voigt function:
    Convolution of non-relativistic Breit-Wigner with gaussian resolution
    
    CPU-efficient Approximation of Voight profile
    -@see T. Ida, M. Ando and H. Toraya, 
    'Extended pseudo-Voigt function for approximating the Voigt profile'
    J. Appl. Cryst. (2000). 33, 1311-1316
    - see doi:10.1107/S0021889800010219
    - see https://doi.org/10.1107/S0021889800010219
    
    Parameters
    - mean  : location 
    - gamma : gamma for Breight-Wigner pole
    - sigma : resolution parameter for gaussian 
    """
    def __init__ ( self             ,
                   name             ,
                   xvar             ,
                   m0        = None ,
                   sigma     = None ,
                   gamma     = None ) :
        #
        ## initialize the base
        # 
        Voigt_pdf.__init__  ( self , name , xvar , m0 , sigma , gamma ) 

        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.PseudoVoigt (
            self.roo_name ( 'pvoigt_' ) , 
            "Pseudo Voigt %s" % self.name ,
            self.xvar   ,
            self.mean   ,
            self.gamma  ,
            self.sigma  )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'm0'        : self.mean  ,
            'sigma'     : self.sigma ,
            'gamma'     : self.gamma ,
            }
    
models.append ( PseudoVoigt_pdf )                          
# =============================================================================
## @class BreitWigner_pdf 
#  Relativistic Breit-Wigner function using Jackson's parameterization
#  J.D.Jackson,
#  "Remarks on the Phenomenological Analysis of Resonances",
#  In Nuovo Cimento, Vol. XXXIV, N.6
#  http://www.springerlink.com/content/q773737260425652/
#  - Blatt-Weisskopf forfactors  are also possible
#  @see Ostap::Models::BreitWigner 
#  @see Ostap::Math::BreitWigner
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class BreitWigner_pdf(PEAK) :
    """Relativistic Breit-Wigner function using Jackson's parameterization
    J.D.Jackson, 'Remarks on the Phenomenological Analysis of Resonances',
    In Nuovo Cimento, Vol. XXXIV, N.6
    http://www.springerlink.com/content/q773737260425652/

    >>> bw    = Ostap.Math.BreitWigner( m_X , g_X , m_Jpsi , m_pipi , 0 )
    >>> breit = Models.BreitWigner_pdf ( 'BW'          ,
    ...                                  bw            ,
    ...                                  xvar  = mass  ,
    ...                                  m0    = m_X   ,
    ...                                  gamma = g_X   )
    
    Parameters:
    - m0          : location Breigt-Wigner function
    - gamma       : width of Breigt-Wigner function
    
    """
    def __init__ ( self               ,
                   name               ,
                   breitwigner        , ## Ostap::Math::BreitWeigner object
                   xvar               ,
                   m0          = None , 
                   gamma       = None ) :        
        #
        ## initialize the base
        #
        PEAK.__init__  ( self  , name  , xvar ,
                              mean        = m0                  ,
                              sigma       = gamma               ,
                              mean_name   = 'm0_%s'      % name ,
                              mean_title  = '#m_{0}(%s)' % name ,                         
                              sigma_name  = 'gamma_%s'   % name ,
                              sigma_title = '#Gamma(%s)' % name )
        
        bw = breitwigner
        assert isinstance ( bw , Ostap.Math.BW ), \
               'Invalid  type of the Breit-Wigner object: %s/%s' % ( bw   , type ( bw ) )
        #
        ## define the actual BW-shape using
        #      Ostap::Math::BreitWeigner object
        self.__breitwigner = breitwigner  ## Ostap::Math::BreitWeigner object

        ## create PDF 
        self.pdf = Ostap.Models.BreitWigner ( 
            self.roo_name ( 'rbw_' ) , 
            "Relativistic Breit-Wigner %s" % self.name ,
            self.xvar        ,
            self.m0          ,
            self.gamma       ,
            self.breitwigner )

        ## save the configuration
        self.config = {
            'name'        : self.name          ,
            'breitwigner' : self.breitwigner   ,
            'xvar'        : self.xvar          ,
            'm0'          : self.m0            ,
            'gamma'       : self.gamma         ,
            }

    @property
    def m0 ( self ) :
        """'m0' : m_0 parameter for Breit-Wigner function (alias for 'mean')"""
        return self.mean
    @m0.setter
    def m0 ( self , value ) :
        self.mean = value 

    @property
    def gamma ( self ) :
        """'gamma'-parameter for Breit-Wigner function (alias for 'sigma')"""
        return self.sigma 
    @gamma.setter
    def gamma ( self, value ) :
        self.sigma = value 
    
    @property
    def Gamma ( self ) :
        """'Gamma'-parameter for Breit-Wigner function (alias for 'sigma')"""
        return self.sigma 
    @Gamma.setter
    def Gamma ( self, value ) :
        self.sigma = value 
        
    @property
    def breitwigner ( self ) :
        """The Breit-Wigner function  itself"""
        return self.__breitwigner

    # =========================================================================
    ## prepare Argand plot as <code>TGraph</code>
    #  @code
    #  bw = ...
    #  argand = bw.argand ( npx = 1000 )
    #  argand.draw ( 'al')  
    #  @endcode
    #  @see  TGraph 
    def argand ( self , x_min =  None , x_max = None , npx = 1000 ) :
        """ prepare Argand plot as `TGraph`
        >>> bw = ...
        >>> argand = bw.argand ( npx = 1000 )
        >>> argand.draw ( 'al')  
        """
        bw = self.pdf.function()
        xmnmx = self.xminmax() 
        if x_min is None and xmnmx : x_min = xmnmx [ 0 ] 
        if x_max is None and xmnmx : x_max = xmnmx [ 1 ] 
        ## make Argand plot 
        return bw.argand ( xmin = x_min , xmax = x_max , npx = npx ) 
            
models.append ( BreitWigner_pdf )


# =============================================================================
## @class BWMC_pdf
#  Multi-channel version of Breit-Wigner function
#  @see Ostap::Models::BreitWignerMC
#  @see Ostap::Math::BreitWignerMC
#  @see Ostap::Math::Channel
#  @code
#  m_Kp  =  493.677 * MeV
#  m_Kz  =  497.614 * MeV
#  m_phi = 1019.462 * MeV
#  g_phi =    4.249 * MeV 
#  br_pm = 0.492
#  br_00 = 0.340 
#  br_xx = 1 - br_pm - br_00
#  ## define three  channels 
#  ch_pm = Ostap.Math.Channel ( g_phi * br_pm , m_Kp , m_Kp , 1 , 1 )
#  ch_00 = Ostap.Math.Channel ( g_phi * br_00 , m_Kz , m_Kz , 1 , 1 )
#  ch_xx = Ostap.Math.Channel ( g_phi * br_xx , 0    , 0    )
#  ## define the Breit-wigner function
#  bw    = Ostap.Math.BreitWignerMC ( m_phi , ch_pm , ch_00 , ch_xx )
#
#  mKK   = ROOT.RooRealVar ( ... ) 
#  pdf = BWMC_pdf ( 'BW' , breitwigner = bw ,
#                    xvar = mKK , mean = m_phi , gamma = g_phi ,
#                    fractions = [ br_pm , br_00 , br_xx ] ) 
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2018-11-25
class BWMC_pdf(PEAK) :
    """Multi-channel version of Relativistic Breit-Wigner function

    >>> m_Kp  =  493.677 * MeV    ## mass of K+ 
    >>> m_Kz  =  497.614 * MeV    ## mass of K0
    >>> m_phi = 1019.462 * MeV    ## mass of phi(1020)
    >>> g_phi =    4.249 * MeV    ## width of phi(1020)
    >>> br_pm = 0.492             ## Br ( phi -> K+K-) 
    >>> br_00 = 0.340             ## Br ( phi -> K0s K0L  )
    >>> br_xx = 1 - br_pm - br_00 ## Br ( phi ->  anything else )
    
    Define three decay channels:
    
    >>> ch_pm = Ostap.Math.Channel ( g_phi * br_pm , m_Kp , m_Kp , 1 , 1 )  ## K+ K- 
    >>> ch_00 = Ostap.Math.Channel ( g_phi * br_00 , m_Kz , m_Kz , 1 , 1 )  ## K0L K0L
    >>> ch_xx = Ostap.Math.Channel ( g_phi * br_xx , 0    , 0    )          ## light stuff 
    
    Define the Breit-Wigner function:
    
    >>> bw    = Ostap.Math.BreitWignerMC ( m_phi , ch_pm , ch_00 , ch_xx )

    >>> mKK   = ROOT.RooRealVar ( ... ) 
    >>> pdf = BWMC_pdf ( 'BW' , breitwigner = bw ,
    ...                  xvar = mKK , mean = m_phi , gamma = g_phi ,
    ...                  fractions = [ br_pm , br_00 , br_xx ] ) 
        
    Parameters:
    - xvar        : fitting variable/observable  
    - mean        : location of Breit-Wigner pole
    - gamma       : total width of Breit-Wigner pole
    - widths      : partial width for the channels  (mutually exclusive with gamma and fractions)
    - fractions   : branching fractions 
    
    """
    def __init__ ( self               ,
                   name               ,
                   breitwigner        , ## Ostap::Math::BreitWignerMC object
                   xvar               ,
                   m0          = None , 
                   gamma       = None ,
                   widths      = []   ,
                   fractions   = []   ) : 

        ## correct type of Breit-Wigner function?
        bw = breitwigner 
        assert isinstance ( bw , Ostap.Math.BreitWignerMC ), \
               'Invalid  type of the Breit-Wigner object: %s/%s' % ( bw   , type ( bw ) )

        ## number of channels 
        nc      = bw.nChannels    ()
        
        self.__brfrs  = ROOT.RooArgList()
        self.__widths = ROOT.RooArgList()
        
        case         = None
        self.__trash = []  ## keep the trash

        ## Valid  cases: 
        if   widths    and nc == len ( width     ) and gamma is None and not fractions : case = 1
        elif fractions and nc == len ( fractions )                   and not widths    : case = 2
        else : raise TypeError ('Gamma/widths/fraction mismatch!')

        ## partial  widths are specified:
        if 1 == case :
            
            for i in range ( nc ) : 
                gi = widths[i] 
                gg = self.make_var ( gi ,
                                     'gamma_%d_%s'     % ( i+1 , name ) ,
                                     '#Gamma_{%d}(%s)' % ( i+1 , name ) ,
                                     None )
                
                self.__trash.append ( gg )
                self.widths.add     ( gg )

            self.__trash.append ( self.widths ) 
            ## construct the total gamma
            formula = '%s ' % self.widths[0].GetName()
            for i in range ( 1 , nc ) : formula += ' + %s' % self.widths[i].GetName()  
            ## for g in self.widths[1:] :

            ## create gamma 
            gamma = Ostap.FormulaVar ( 'gamma_%s'    % name ,
                                       '#Gamma_(%s)' % name , formula , self.widths )
            self.__trash.append ( gamma )
            
        # =====================================================================
        ## initialize the base 
        # =====================================================================
        PEAK.__init__ ( self , name , xvar                ,
                        mean        =  m0                 ,
                        siga        =  gamma              ,
                        mean_name   = 'm0_%s'      % name ,
                        mean_title  = '#m_{0}(%s)' % name ,
                        sigma_name  = 'gamma_%s'   % name ,
                        sigma_title = '#Gamma(%s)' % name )
                
        ## create branching fractions 
        if 1 == case :
            
            for i in range ( nc ) :
                
                gi  = self.widths[i]
                lst = ROOT.RooArgList ( self.gamma ,  gi )
                self.__formulas_lists.append ( lst ) 
                ## br  = ROOT.RooFormulaVar ( 'brfr_%d_%s'  % ( i + 1 , name ) ,
                ##                            'Br_{%d}(%s)' % ( i + 1 , name ) ,
                ##                            '%s / %s'     % ( gi.GetName() , self.gamma.GetName() ) , lst )
                br  =  Ostap.MoreRooFit.Division ( 'brfr_%d_%s'  % ( i + 1 , name ) ,
                                                   'Br_{%d}(%s)' % ( i + 1 , name ) , self.gamma , gi )
                self.brfrs.add      ( br )
                self.__trash.append ( br )
                
        ##  branching fractions are specified 
        elif 2 == case : 
            
            for i in range ( nc ) :
                
                bi = fractions [i] 
                br = self.make_var       ( bi ,
                                           'brfr_%d_%s'  % ( i+1 , name ) ,
                                           'Br_{%d}(%s)' % ( i+1 , name ) ,None ) 
                self.brfrs.add       ( br )    
                self.__trash.append  ( br ) 
                
                ls = ROOT.RooArgList ( self.gamma ,  br )
                self.__trash.append  ( ls ) 

                ## gg  = ROOT.RooFormulaVar ( 'gamma_%d_%s'     % ( i + 1 , name ) ,
                ##                           '#Gamma_{%d}(%s)' % ( i + 1 , name ) ,
                ##                           '%s * %s'         % ( br.GetName() , self.gamma.GetName() ) , ls )
                
                gg  = Ostap.MoreRooFit.Product ( 'gamma_%d_%s'     % ( i + 1 , name ) ,
                                                 '#Gamma_{%d}(%s)' % ( i + 1 , name ) , self.gamma , br ) 
                self.widths.add      ( gg )
                self.__trash.append  ( gg ) 
            

        self.__gammas    = tuple ( [ i for i in self.widths ] )
        self.__fractions = tuple ( [ i for i in self.brfrs  ] )
                       
        #
        ## define the actual BW-shape using
        #      Ostap::Math::BreitWeignerMC object
        #
        self.__breitwigner = breitwigner  ## Ostap::Math::BreitWeignerMC object
        
        ## create PDF 
        self.pdf = Ostap.Models.BreitWignerMC ( 
            self.roo_name ( 'rbwmc_' ) , 
            "Multi-channel relativistic Breit-Wigner %s" % self.name ,
            self.xvar        ,
            self.mean        ,
            self.widths      , 
            self.breitwigner )
            
        ## save the configuration
        self.config = {
            'name'        : self.name          ,
            'breitwigner' : self.breitwigner   ,
            'xvar'        : self.xvar          ,
            'mean'        : self.mean          }
        
        if   1 == case : self.config.update (  { 'widths'    : self.gammas    ,
                                                 'fractions' : ()             } )
        elif 2 == case : self.config.update (  { 'fractions' : self.fractions ,
                                                 'gamma'     : self.gamma     ,
                                                 'widths'    : ()             } )
    @property
    def m0 ( self ) :
        """'m0' : m_0 parameter for Breit-Wigner function (alias for 'mean')"""
        return self.mean
    @m0.setter
    def m0 ( self , value ) :
        self.mean = value 

    @property
    def gamma ( self ) :
        """'gamma'-parameter for Breit-Wigner function (alias for 'sigma')"""
        return self.sigma 
    @gamma.setter
    def gamma ( self, value ) :
        self.sigma = value 
    
    @property
    def Gamma ( self ) :
        """'Gamma'-parameter for Breit-Wigner function (alias for 'sigma')"""
        return self.gamma 
    @Gamma.setter
    def Gamma ( self, value ) :
        self.sigma = value 

    @property
    def breitwigner ( self ) :
        """The Breit-Wigner function  itself"""
        return self.__breitwigner

    @property
    def widths  ( self ) :
        """'widths'  : partial widths for different decay channels"""
        return self.__widths

    @property
    def brfrs   ( self ) :
        """'brfrs'  : branching fractions for different decay channels"""
        return self.__brfrs 

    @property
    def gammas ( self ) :
        """'gammas'  : partial widths for different decay channels"""
        return self.__gammas  

    @property
    def fractions ( self ) :
        """'fractions'  : branching fractions for different decay channels"""
        return self.__fractions

    # =========================================================================
    ## prepare Argand plot as <code>TGraph</code>
    #  @code
    #  bw = ...
    #  argand = bw.argand ( npx = 1000 )
    #  argand.draw ( 'al')  
    #  @endcode
    #  @see  TGraph 
    def argand ( self , x_min =  None , x_max = None , npx = 1000 ) :
        """ prepare Argand plot as `TGraph`
        >>> bw = ...
        >>> argand = bw.argand ( npx = 1000 )
        >>> argand.draw ( 'al')  
        """
        bw = self.pdf.function()
        xmnmx = self.xminmax() 
        if x_min is None and xmnmx : x_min = xmnmx [ 0 ] 
        if x_max is None and xmnmx : x_max = xmnmx [ 1 ] 
        ## make Argand plot 
        return bw.argand ( xmin = x_min , xmax = x_max , npx = npx ) 
            

models.append ( BWMC_pdf )

# =============================================================================
## @class BWI_pdf
#  (relativistic) Breit-Wigner function + some interference
#   Breit-Wigner with some embedded interference: 
#   \f[ f(x) = \left| \upalpha b(x) + A(x)_{\mathrm{BW}} \right|^2 \f], 
#   where \f$b(x)\f$ - any smooth function and 
#   \f$ A(x)_{\mathrm{BW}} \f$ is Breit-Wigner amplitude 
#  @see Ostap.Models.BWI
class BWI_pdf (BreitWigner_pdf) :
    """ (Relativistic) Breit-Wigner function + some interference
    Breit-Wigner with some embedded interference
    - see Ostap.Models.BWI
    """
    def __init__ ( self         ,
                   name         ,
                   breitwigner  ,
                   xvar         ,
                   m0    = None ,
                   gamma = None ,
                   bkg   = -1   ,   ## background function 
                   a     = None ,   ## background scale 
                   phi   = 0    ) : ## bakcgrouns phase 
        
        ## initialize the base 
        BreitWigner_pdf.__init__ ( self , name , breitwigner , xvar , m0 , gamma )

        self.__bw = self.pdf

        if   isinstance ( bkg , PDF1            ) :             ## PDF ? 
            ## PDF?
            self.__bkg = bkg
            self.__b   = self.__bkg.pdf
        elif isinstance ( bkg , ROOT.RooAbsPdf  ) :              ## PDF? 
            ## PDF? 
            from ostap.fitting.basic import Generic1D_pdf as G1D 
            self.__bkg = G1D ( bkg , self.xvar ) 
            self.__b   = self.__bkg.pdf
        elif isinstance ( bkg , ROOT.RooRealVar ) or isinstance ( bkg , ROOT.RooConstVar ) :
            ## constant?
            self.__bkg = bkg
            self.__b   = self.make_var ( bkg              ,
                                         "bkg_%s"  % name ,
                                         "bkg(%s)" % name , None )
        elif isinstance ( bkg , ROOT.RooAbsReal ) :               ## function ?
            ## function? 
            from ostap.fitting.basic import Generic1D_pdf as G1D 
            self.__bkg = G1D ( bkg , self.xvar , special = True ) 
            self.__b   = self.__bkg.pdf
        else :                                                    ## function?
            ## function ? 
            from ostap.fitting.background import make_bkg as MKB
            self.__bkg = MKB ( bkg , 'B_4'+ self.name , self.xvar )
            self.__b   = self.__bkg.pdf 
            
        ## create background magnitude 
        self.__a  = self.make_var  ( a              ,
                                     "a_%s"    % name ,
                                     "a(%s)"   % name ,
                                     None , 0 , 1.e+5 )
        
        ## create background phase 
        self.__phi = self.make_var ( phi              ,
                                     "phi_%s"  % name ,
                                     "phi(%s)" % name ,
                                     None , -12 , 12  )

        ## finally create PDF
        self.pdf = Ostap.Models.BWI (
            self.roo_name ( 'bwi_' ) ,
            "Breit-Wigner with interference %s" % self.name  ,
            self.bw                  ,
            self.b                   ,
            self.a                   ,
            self.phi                 ) 

        ## save configuration
        self.config = {
            'name'        : self.name        ,
            'xvar'        : self.xvar        , 
            'breitwigner' : self.breitwigner , 
            'm0'          : self.m0          ,
            'gamma'       : self.gamma       ,
            'bkg'         : self.bkg         ,
            'a'           : self.a           , 
            'phi'         : self.phi         } 
                        
    @property
    def m0 ( self ) :
        """'m0' : m_0 parameter for Breit-Wigner function (alias for 'mean')"""
        return self.mean
    @m0.setter
    def m0 ( self , value ) :
        self.mean = value 


    @property
    def bw          ( self ) :
        """The Breit-Wigner PDF itself"""
        return self.__bw
    @property
    def bkg          ( self ) :
        """The background"""
        return self.__bkg
    @property
    def b           ( self ) :
        """The background"""
        return self.__b
    @property
    def a           ( self ) :
        """The background factor"""
        return self.__a
    @property
    def phi         ( self ) :
        """The background phase"""
        return self.__phi

# =============================================================================
## @class BWPS
#  Breit-Wigner function modulated with extra phase-space and polynomial factors
#  @see Ostap::Models::BWPS
#  @see Ostap::Math::BWPS
class BWPS_pdf(BreitWigner_pdf,Phases) :
    """Breit-Wigner function modulated with extra phase-space and polynomial factors
    - see Ostap.Models.BWPS
    - see Ostap.Math.BWPS
    """
    
    def __init__ ( self             ,
                   name             ,
                   breitwigner      , ## Ostap::Math::BWPS object
                   xvar             ,
                   m0        = None ,
                   gamma     = None ,
                   the_phis  = None ) :

        if   isinstance ( breitwigner , Ostap.Math.BWPS ) : pass
        elif isinstance ( breitwigner , tuple ) :
            breitwigner = Ostap.Math.BWPS  ( *breitwigner )
        else :
            raise ArgumentError("BWPS_pdf: Invalidd type of breitwigner") 
        
        ## initialize the base classes 
        BreitWigner_pdf.__init__  ( self ,
                                    name ,
                                    breitwigner = breitwigner.breit_wigner () , 
                                    xvar        = xvar   ,
                                    m0          = m0     ,
                                    gamma       = gamma  )
        
        Phases.__init__ ( self , breitwigner.npars () , the_phis ) 

        ## make "original" BW-pdf 
        self.__bw_pdf = BreitWigner_pdf ( name        = self.name + '_orig' ,
                                          breitwigner = self.breitwigner    ,
                                          xvar        = self.xvar           ,
                                          m0          = self.m0             ,
                                          gamma       = self.gamma          )
        self.__bwps = breitwigner
        
        ## finally create PDF        
        self.pdf = Ostap.Models.BWPS (
            self.roo_name ( 'bwps_' ) ,
            "Breit-Wigner with phase space %s" % self.name  ,
            self.xvar         ,
            self.m0           ,
            self.gamma        ,
            self.phi_list     ,
            self.bwps         ) 
            
        ## save configuration
        self.config = {
            'name'        : self.name  ,
            'xvar'        : self.xvar  , 
            'breitwigner' : ( self.bwps.breit_wigner () ,
                              self.bwps.phase_space  () ,
                              self.bwps.use_rho      () ,
                              self.bwps.use_N2       () ) ,                               
            'm0'          : self.mean  ,
            'gamma'       : self.gamma ,
            'the_phis'    : self.phis  }
        
    @property
    def bw_pdf  ( self ) :
        """'bw_pdf' : 'original' Breit-Wigner pdf (no additional phase space  factors)"""
        return self.__bw_pdf
    
    @property
    def bwps ( self ) :
        """The Breit-Wigner function (BWPS) itself"""
        return self.__bwps

    
models.append ( BWPS_pdf )



# =============================================================================
## @class BW3L_pdf
#  Breit-Wigner function modulated with \f$ p^{2L+1}\f$ factor
#   - it can approximate the mass distrbition from 3-body decays
# e.g.  \f$ \eta^{\prime)  \rigtharrow \left(\rho^0 
#               \rigtharrow \pi^+ \pi^-\right)\gamma \f$~decays
#    or similar  configurations  
# 
# \f[ f(x) \equiv F_{\mathrm{BW}}(x) p(x|M_0,m_3)^{2L+1} \f]
#  - \f$ p(x|M,m_3) \f$ is a momentumm of the 3rd particle, \f$P_3\f$ 
#       in the \f$ P \rightarrow \left( P_{\mathrm{BW}} \rightharrow 
#      P_1 P_2 \right) P_3 \f$ decay chain
#  - \f$ M \f$ is a (fixed) mass of "mother" particle \f$P\f$
#  - \f$ m_1\f$ is a (fixed) mass of 1st particle \f$P_1\f$
#  - \f$ m_2\f$ is a (fixed) mass of 2nd particle \f$P_2\f$
#  - \f$ m_3\f$ is a (fixed) mass of 3rd particle \f$P_3\f$
#  - \f$ x \equiv m_{23} \f$ is a mass intermediate Breit-Wigner particle \f$P_{\mathrm{BW}}\f$
#  - \f$ L \f$  is an orbital momentum between \f$ P_{\mathrm{BW}}\f$ and \f$ P_3\f$
# 
#  It is assumed that  \f$ m_1\f$  and \f$ m_2\f$ parameters 
#  are in agreement with the Breit-Wigner definition 
#  @see Ostap::Models::BW3L
#  @see Ostap::Math::BW3L
class BW3L_pdf(BreitWigner_pdf) :
    """ Breit-Wigner function modulated with  p^{2L+1} factor
    - it can approximate the mass distrbition from 3-body decays
    e.g. ( eta'  ->  ( rho0 -> pi+ pi- ) gamma ) decay or similar  configurations
    
    - see Ostap.Models.BW3L
    - see Ostap.Math.BW3L
    """
    
    def __init__ ( self             ,
                   name             ,
                   breitwigner      , ## Ostap::Math::BW3L object
                   xvar             ,
                   m0        = None ,
                   gamma     = None ) : 

        if   isinstance ( breitwigner , Ostap.Math.BW3L ) : pass
        elif isinstance ( breitwigner , tuple ) :
            breitwigner = Ostap.Math.BW3L  ( *breitwigner )
        else :
            raise ArgumentError("BW3L_pdf: Invalidd type of breitwigner") 
        
        ## initialize the base classes 
        BreitWigner_pdf.__init__  ( self ,
                                    name ,
                                    breitwigner = breitwigner.breit_wigner () , 
                                    xvar        = xvar   ,
                                    m0          = m0     ,
                                    gamma       = gamma  )
        
        ## make "original" BW-pdf 
        self.__bw_pdf = BreitWigner_pdf ( name        = self.name + '_orig' ,
                                          breitwigner = self.breitwigner    ,
                                          xvar        = self.xvar           ,
                                          m0          = self.m0             ,
                                          gamma       = self.gamma          )
        self.__bw3l = breitwigner
        
        ## finally create PDF        
        self.pdf = Ostap.Models.BW3L (
            self.roo_name ( 'bw3l_' ) ,
            "Breit-Wigner form 3-body decay  %s" % self.name  ,
            self.xvar         ,
            self.m0           ,
            self.gamma        ,
            self.bw3l         ) 
            
        ## save configuration
        self.config = {
            'name'        : self.name  ,
            'xvar'        : self.xvar  , 
            'breitwigner' : ( self.bw3l.breit_wigner () ,
                              self.bw3l.M  () ,
                              self.bw3l.m1 () ,
                              self.bw3l.m2 () ,
                              self.bw3l.m3 () ,
                              self.bw3l.L  () ) , 
            'm0'          : self.mean  ,
            'gamma'       : self.gamma }
        
    @property
    def bw_pdf  ( self ) :
        """'bw_pdf' : 'original' Breit-Wigner pdf (no additional factors)"""
        return self.__bw_pdf
    
    @property
    def bw3l ( self ) :
        """The Breit-Wigner function (BW3L) itself"""
        return self.__bw3l

    
models.append ( BW3L_pdf )


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
#  @see Ostap::Math::Flatte2
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-01-18
class Flatte_pdf(PEAKMEAN) :
    """Flatte function:
    S.M.Flatte, 'Coupled-channel analysis of the (pi eta)
    and (KbarK) systems near (KbarK) threshold' 
    Phys. Lett. B63, 224 (1976
    http://www.sciencedirect.com/science/article/pii/0370269376906547

    Typical case:    f0(980) -> pi+ pi- & K+ K- shapes 
    """
    def __init__ ( self              ,
                   name              ,
                   flatte            ,    ## Ostap::Math::Flatte/Flatte2
                   xvar              ,
                   m0       = None   ,    ## the pole 
                   m0g1     = None   ,    ## m0*g1 
                   g2og1    = None   ,    ## g2/g1 
                   g1       = None   ,    ## g1 
                   g2       = None   ,    ## g2 
                   gamma0   = None   ) :  ## gamma0 

        assert isinstance ( flatte , Ostap.Math.Flatte ), \
               'Invalid type for flatte %s' %  type ( flatte )

        ## initialize the base
        with CheckMean ( False ) :
            # for Flatte-function m0 can be outside the interesting interval 
            PEAKMEAN.__init__  ( self , name , xvar ,
                                 mean       = m0  ,
                                 mean_name  = 'm0_%s'      % name ,
                                 mean_title = '#m_{0}(%s)' % name )

        self.__my_case = 0 
        if   all_args ( self.mean , m0g1 , g2og1 ) : self.__my_case = 1 
        elif all_args ( self.mean ,   g1 , g2    ) : self.__my_case = 2
            
        assert self.case in  ( 1 , 2 ), 'Invalid combination of (m0g1,g2og1:g1,g2)arguments!'
        
        self.__flatte = flatte
            
        self.__gamma0 = self.make_var  ( gamma0                  ,
                                         'gamma0_%s'      % name ,
                                         '#Gamma_{0}(%s)' % name ,
                                         None   , 0 , gamma0 )
        if  1 == self.case : 
            
            self.__m0g1 = self.make_var  ( m0g1                     ,
                                           'm0g1_%s'          % name ,
                                           '#m_{0}#g_{1}(%s)' % name ,
                                           None  )
            
            self.__g2og1 = self.make_var ( g2og1 ,
                                           'g2og1_%s'          % name ,
                                           '#g_{2}/#g_{1}(%s)' % name ,
                                           None , 0.001 , 200  ) 

            self.__g1 = self.vars_divide   ( self.m0g1  , self.m0 , name = 'g1_%s' % name , title = "g_1(%s)" % name )
            self.__g2 = self.vars_multiply ( self.g2og1 , self.g1 , name = 'g2_%s' % name , title = "g_2(%s)" % name )
            
        elif 2 == self.case :
            
            self.__g1 =  self.make_var  ( g1                 ,
                                          'g1_%s'     % name ,
                                          'g_{1}(%s)' % name ,
                                          g1                 , None  )
            self.__g2 =  self.make_var  ( g2                 ,
                                          'g2_%s'     % name ,
                                          'g_{2}(%s)' % name ,
                                          g2                 , None )
            
            self.__m0g1  = self.vars_multiply ( self.m0 , self.g1 , name = 'm0g1_%s'  % name , title = "m_0g_1(%s)"  % name )
            self.__g2og1 = self.vars_divide   ( self.g2 , self.g1 , name = 'g2og1_%s' % name , title = "g_2/g_1(%s)" % name )
                
        ## create PDF 
        self.pdf = Ostap.Models.Flatte (
            self.roo_name ( 'flatte_' ) ,
            "Flatte %s" % self.name  ,
            self.xvar    ,
            self.m0      ,
            self.g1      ,
            self.g2      ,
            self.gamma0  ,
            self.flatte  )

        ## save the configuration
        cnf = {
            'name'        : self.name    ,
            'flatte'      : self.flatte  ,
            'xvar'        : self.xvar    ,
            'm0'          : self.m0      ,
            'gamma0'      : self.gamma0  ,
            }
        
        if   1 == self.case : cnf.update ( { 'm0g1' : self.m0g1 , 'g2og1' : self.g2og1 } )
        elif 2 == self.case : cnf.update ( { 'g1'   : self.g1   , 'g2'    : self.g2    } )

        self.config = cnf
        
    @property
    def case  (self ) :
        """'case' : How the input argument are  specified: 1: (m0g1,g2og1) vs 2: (g1,g2) """
        return self.__my_case
    
    @property
    def m0 ( self ) :
        """'m0'-parameter for Flatte-function (same as 'mean')"""
        return self.mean 
    @m0.setter
    def m0  ( self, value ) :
        self.mean = value

    @property
    def m0g1 ( self ) :
        """'m0*g1'-parameter for Flatte-function"""
        return self.__m0g1
    @m0g1.setter
    def m0g1 ( self, value ) :
        assert gasattr ( self.__m0g1 , 'setVal' ),"'m0g1'-parameter can't be set!"
        value = float ( value )
        self.__m0g1.setVal ( value ) 

    @property
    def g2og1 ( self ) :
        """'g2/g1'-parameter for Flatte-function"""
        return self.__g2og1
    @g2og1.setter
    def g2og1 ( self, value ) :
        assert hasattr ( self.__g2og1 , 'setVal'),"'g2og1'-parameter can't be set!"        
        value = float ( value )
        assert 0 < value, "'g2/g1'-parameter for Flatte-function must be positive"
        self.__g2og1.setVal ( value )

    @property
    def g1 ( self ) :
        "'g1'-parameter for Flatte-function"
        return self.__g1
    @g1.setter
    def g1 ( self , value ) :
        assert hasattr ( self.__g1 , 'setVal' ),"'g1'-parameter can't be set!"
        value = float ( value )
        self.__g1.setVal ( value ) 

    @property
    def g2 ( self ) :
        "'g2'-parameter for Flatte-function"
        return self.__g2
    @g2.setter
    def g2 ( self , value ) :
        assert hasattr ( self.__g2 , 'setVal' ),"'g2'-parameter can't be set!"
        value = float ( value )
        self.__g2.setVal ( value ) 

    @property
    def gamma0 ( self ) :
        "'gamma0'-parameter for Flatte-function"
        return self.__gamma0
    @gamma0.setter
    def gamma0 ( self , value ) :
        value = float ( value )
        self.__gamma0.setVal ( value ) 

    @property
    def flatte ( self ) :
        """The Flatte function itself"""
        return self.__flatte

    # =========================================================================
    ## prepare Argand plot as <code>TGraph</code>
    #  @code
    #  bw = ...
    #  argand = bw.argand ( npx = 1000 )
    #  argand.draw ( 'al')  
    #  @endcode
    #  @see  TGraph 
    def argand ( self , x_min =  None , x_max = None , npx = 1000 ) :
        """ prepare Argand plot as `TGraph`
        >>> bw = ...
        >>> argand = bw.argand ( npx = 1000 )
        >>> argand.draw ( 'al')  
        """
        bw    = self.pdf.function()
        xmnmx = self.xminmax() 
        if x_min is None and xmnmx : x_min = xmnmx [ 0 ] 
        if x_max is None and xmnmx : x_max = xmnmx [ 1 ] 
        ## make Argand plot 
        return bw.argand ( xmin = x_min , xmax = x_max , npx = npx ) 
            

models.append ( Flatte_pdf )                          



# ============================================================================
class FlattePS_pdf(Flatte_pdf,Phases) :
    """Flatte function:
    S.M.Flatte, 'Coupled-channel analysis of the (pi eta)
    and (KbarK) systems near (KbarK) threshold' 
    Phys. Lett. B63, 224 (1976
    http://www.sciencedirect.com/science/article/pii/0370269376906547

    Typical case:    f0(980) -> pi+ pi- & K+ K- shapes 
    """
    def __init__ ( self              ,
                   name              ,
                   flatte            ,   ## Ostap::Math::BWPS 
                   xvar              ,
                   m0       = None   ,   ## the pole 
                   m0g1     = None   ,   ## m0*g1 
                   g2og1    = None   ,   ## g2/g1 
                   g1       = None   ,   ## g1 
                   g2       = None   ,   ## g2 
                   gamma0   = None   ,   ## gamma0 
                   the_phis = None   ) : ##
        
        if   isinstance ( flatte , Ostap.Math.BWPS ) : pass
        elif isinstance ( flatte , tuple ) :
            flatte = Ostap.Math.BWPS  ( *flatte )
        else :
            raise ArgumentError("FlattePS_pdf: Invalidd type of flatte") 
        

        assert isinstance ( flatte, Ostap.Math.BWPS ),\
               'Invalid type for breitwigner %s' %  type ( flatte )
        
        ## initialize the base classes 
        Flatte_pdf.__init__  ( self ,
                               name ,
                               flatte = flatte.breit_wigner () ,
                               xvar   = xvar     , 
                               m0     = m0       ,
                               m0g1   = m0g1     ,
                               g2og1  = g2og1    ,
                               g1     = g1       ,
                               g2     = g2       ,
                               gamma0 = gamma0   )
        
        Phases.__init__ ( self , flatte.npars () , the_phis )

        ## store the "original"" Flatte PDF 
        if  1 ==  self.case : 
            self.__flatte_pdf = Flatte_pdf ( name   = self.name + "_orig" ,
                                             flatte = self.flatte ,  
                                             xvar   = self.xvar   ,
                                             m0     = self.m0     ,
                                             m0g1   = self.m0g1   ,
                                             g2og1  = self.g2og1  ,
                                             gamma0 = self.gamma0 )
        elif 2 ==  self.case : 
            self.__flatte_pdf = Flatte_pdf ( name   = self.name + "_orig" ,
                                             flatte = self.flatte ,  
                                             xvar   = self.xvar   ,
                                             m0     = self.m0     ,
                                             g1     = self.g1     ,
                                             g2     = self.g2     ,
                                             gamma0 = self.gamma0 )
            
        self.__bwps   = flatte
  
        self.__g_list = ROOT.RooArgList ( self.g1 , self.g2 , self.gamma0 )
        
        ## finally create PDF
        self.pdf = Ostap.Models.BWPS (
            self.roo_name ( 'flatteps_' ) ,
            "Flatte with phase space %s" % self.name  ,
            self.xvar         ,
            self.m0           ,
            self.g_list       ,
            self.phi_list     ,
            self.bwps         ) 
        
        ## save the configuration
        cnf = {
            'name'        : self.name    ,
            'flatte'      : ( self.bwps.breit_wigner () ,
                              self.bwps.phase_space  () ,
                              self.bwps.use_rho      () ,
                              self.bwps.use_N2       () ) ,                               
            'xvar'        : self.xvar    ,
            'm0'          : self.m0      ,
            'gamma0'      : self.gamma0  ,
            'the_phis'    : self.phis    , 
            }
        
        if   1 == self.case : cnf.update ( { 'm0g1' : self.m0g1 , 'g2og1' : self.g2og1 } )
        elif 2 == self.case : cnf.update ( { 'g1'   : self.g1   , 'g2'    : self.g2    } )

        self.config = cnf

    @property
    def flatte_pdf  ( self ) :
        """'flatte_pdf' : 'original' Flatte pdf (no additional phase space  factors)"""
        return self.__flatte_pdf
    
    @property
    def bwps ( self ) :
        """The Breit-Wigner function (BWPS) itself"""
        return self.__bwps
    
    @property
    def g_list ( self ) :
        """'g_list' list of gammas for Breit-Wigner"""
        return self.__g_list

models.append ( FlattePS_pdf )

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
class LASS_pdf(PEAK) :
    """Kappa pole:
    The LASS parameterization (Nucl. Phys. B296, 493 (1988))
    describes the 0+ component of the Kpi spectrum.
    It consists of the K*(1430) resonance together with an
    effective range non-resonant component
    """
    def __init__ ( self               ,
                   name               ,
                   xvar               ,
                   m0       = None    ,    ## mass  of K*(1430)
                   g0       = None    ,    ## width of K*(1430)
                   a        = 1.94e-3 , 
                   r        = 1.76e-3 ,
                   e        = 1.0     ,    ## elasticity                    
                   mKaon    = 493.7   ,    ## kaon mass 
                   mPion    = 139.6   ) :  ## pion mass 

        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar , 
                         mean        = m0 ,
                         sigma       = g0 ,
                         mean_name   = 'm0_%s'      % name ,
                         mean_title  = '#m_{0}(%s)' % name ,
                         sigma_name  = 'g_%s'       % name ,
                         sigma_title = 'g_0(%s)'    % name )
        
        
        self.__g0 = self.sigma
        self.__m0 = self.mean
            
        self.__a = self.make_var ( a                  ,
                                   'aLASS_%s'  % name ,
                                   "aLASS(%s)" % name ,
                                   None               , 
                                   1.94e-3            ,
                                   1.94e-3            ,
                                   1.94e-3            ) 
        self.__r = self.make_var ( r             ,
                                   'rLASS_%s'  % name ,
                                   "rLASS(%s)" % name ,
                                   None               ,
                                   1.76e-3            ,
                                   1.76e-3            ,
                                   1.76e-3            ) 
        self.__e = self.make_var ( e            ,
                                   'eLASS_%s'  % name ,
                                   "eLASS(%s)" % name ,
                                   None               , 
                                   1.0                , 
                                   0.5                ,
                                   2.0                )
        
        self.__mKaon = mKaon
        self.__mPion = mPion
        
        ## create PDF 
        self.pdf = Ostap.Models.LASS (
            self.roo_name ( 'lass_' ) ,
            "LASS/kappa %s" % self.name  ,
            self.xvar    ,
            self.m0      ,
            self.g0      ,
            self.a       ,
            self.r       ,
            self.e       ,
            self.__mKaon ,
            self.__mPion )

        ## save the configuration
        self.config = {
            'name'        : self.name    ,
            'xvar'        : self.xvar    ,
            'm0'          : self.m0      ,
            'g0'          : self.g0      ,
            'a'           : self.a       ,
            'r'           : self.r       ,
            'e'           : self.e       ,
            'mKaon'       : self.__mKaon ,
            'mPion'       : self.__mPion ,
            }
    
    @property
    def g0 ( self ) :
        """'g0'-parameter for LASS-function (same as 'sigma')"""
        return self.sigma
    @g0.setter
    def g0 ( self, value ) :
        self.sigma = value 

    @property
    def m0 ( self ) :
        """'m0'-parameter for LASS-function (same as 'mean')"""
        return self.__gamma
    @m0.setter
    def m0 ( self, value ) :
        self.mean = value 

    @property
    def a ( self ) :
        """'a'-parameter for LASS-function"""
        return self.__a_
    @a.setter
    def a ( self, value ) :
        value = float ( value )
        self.__a.setVal ( value ) 

    @property
    def r ( self ) :
        """'r'-parameter for LASS-function"""
        return self.__r
    @r.setter
    def r ( self, value ) :
        value = float ( value )
        self.__r.setVal ( value ) 

    @property
    def e ( self ) :
        """'e'-parameter for LASS-function"""
        return self.__e
    @e.setter
    def e ( self, value ) :
        value = float ( value )
        self.__e.setVal ( value ) 

models.append ( LASS_pdf )                          
# =============================================================================
## @class Bugg_pdf
#  The parameterization of sigma pole by B.S.Zou and D.V.Bugg, Phys.Rev. D48 (1993) R3948.
#  @see https://doi.org/10.1103/PhysRevD.48.R3948
#  @see Ostap::Models::Bugg
#  @see Ostap::Math::Bugg
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class Bugg_pdf(PEAK) :
    """ The parameterization of sigma pole by
    B.S.Zou and D.V.Bugg, Phys.Rev. D48 (1993) R3948.
    https://doi.org/10.1103/PhysRevD.48.R3948
    """
    def __init__ ( self           ,
                   name           ,
                   xvar           ,
                   m     = 0.9264 , 
                   g2    = 0.0024 , ## g2-parameter
                   b1    = 0.5848 , ## b1-parameter [GeV]
                   b2    = 1.6663 , ## b2-parameter [GeV^-1]
                   a     = 1.082  , ##  a-parameter [GeV^2]
                   s1    = 2.8    , ## s1-parameter [GeV^2]
                   s2    = 3.5    , ## s2-parameter
                   mPion = 0.1396 ) :  ## pion mass 
        #
        ## initialize the base
        # 
        PEAK.__init__  ( self , name , xvar ,
                         mean        = m    ,
                         sigma       = g2   ,
                         mean_name   = 'mBugg_%s'        % name ,
                         mean_title  = '#m_{Bugg}_2(%s)' % name ,
                         sigma_name  = 'gamma2_%s'       % name ,
                         sigma_title = '#gamma_2(%s)'    % name )

        
        self.__bugg_g2 = self.sigma
        self.__gamma   = self.sigma
        ##
        self.__bugg_m = self.mean

        self.__bugg_b1 = self.make_var ( b1                  ,
                                         'b1Bugg_%s'  % name ,
                                         "b1Bugg(%s)" % name ,
                                         None , 0.5848 , 0 , 2  )
        
        self.__bugg_b2 = self.make_var ( b2             ,
                                         'b2Bugg_%s'  % name ,
                                         "b2Bugg(%s)" % name ,
                                         None ,  1.6663 , 1 , 2  ) 
        
        self.__bugg_a  = self.make_var ( a             ,
                                         'aBugg_%s'  % name ,
                                         "aBugg(%s)" % name ,
                                         None , 1.082 , 0.5 , 5  ) 
        
        self.__bugg_s1  = self.make_var ( s1           ,
                                          's1Bugg_%s'  % name ,
                                          "s1Bugg(%s)" % name ,
                                          None , 2.8 , 1 , 5  ) 
        
        self.__bugg_s2  = self.make_var ( s2           ,
                                          's2Bugg_%s'  % name ,
                                          "s2Bugg(%s)" % name ,
                                          None , 3.5 , 1 , 5  ) 
        
        self.__mPion = mPion 
        ## create PDF 
        self.pdf = Ostap.Models.Bugg ( 
            self.roo_name ( 'bugg_' ) ,
            "Bugg/sigma %s" % self.name  ,
            self.xvar      ,
            self.__bugg_m  ,
            self.__bugg_g2 ,
            self.__bugg_b1 ,
            self.__bugg_b2 ,
            self.__bugg_a  ,
            self.__bugg_s1 ,
            self.__bugg_s2 ,
            self.__mPion )
        
        ## save the configuration
        self.config = {
            'name'        : self.name      ,
            'xvar'        : self.xvar      ,
            'm'           : self.__bugg_m  ,
            'g2'          : self.__bugg_g2 ,
            'b1'          : self.__bugg_b1 ,
            'b2'          : self.__bugg_b2 ,
            'a'           : self.__bugg_a  ,
            's1'          : self.__bugg_s1 ,
            's2'          : self.__bugg_s2 ,
            'mPion'       : self.__mPion   ,
            }

    @property
    def gamma ( self ) :
        """'gamma'-parameter ('g2') for Bugg function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        self.__gamma.setVal ( value ) 

    @property
    def g2 ( self ) :
        """'g2'-parameter for Bugg function"""
        return self.__bugg_g2
    @g2.setter
    def g2 ( self, value ) :
        value = float ( value )
        self.__bugg_g2.setVal ( value ) 

    @property
    def b1 ( self ) :
        """'b1'-parameter for Bugg function"""
        return self.__bugg_b1
    @b1.setter
    def b1 ( self, value ) :
        value = float ( value )
        self.__bugg_b1.setVal ( value ) 

    @property
    def b2 ( self ) :
        """'b2'-parameter for Bugg function"""
        return self.__bugg_b2
    @b2.setter
    def b1 ( self, value ) :
        value = float ( value )
        self.__bugg_b2.setVal ( value ) 

    @property
    def a ( self ) :
        """'a'-parameter for Bugg function"""
        return self.__bugg_a
    @a.setter
    def a ( self, value ) :
        value = float ( value )
        self.__bugg_a.setVal ( value ) 

    @property
    def s1 ( self ) :
        """'s1'-parameter for Bugg function"""
        return self.__bugg_s1
    @s1.setter
    def s1 ( self, value ) :
        value = float ( value )
        self.__bugg_s1.setVal ( value ) 

    @property
    def s2 ( self ) :
        """'s2'-parameter for Bugg function"""
        return self.__bugg_s2
    @s2.setter
    def s2 ( self, value ) :
        value = float ( value )
        self.__bugg_s2.setVal ( value ) 

        
models.append ( Bugg_pdf )

## # =============================================================================
## ## @class Swanson_pdf
## #  S-wave cusp
## #  @see LHCb-PAPER-2016-019, Appendix
## #  @see E. S. Swanson, Cusps and exotic charmonia, arXiv:1504.07952
## #  @see http://arxiv.org/abs/1504.07952
## #  @see Ostap::Models::Swanson
## #  @see Ostap::Math::Swanson
## #  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
## #  @date 2016-06-12
## class Swanson_pdf(PDF) :
##     """ S-wave cusp
##     - LHCb-PAPER-2016-019, Appendix
##     - E. S. Swanson, Cusps and exotic charmonia, arXiv:1504.07952
##     - http://arxiv.org/abs/1504.07952
##     """
##     def __init__ ( self              ,
##                    name              ,
##                    swanson           , ## Ostap::Math::Swanson objects 
##                    xvar              ,
##                    beta0    = None   ) : 
##         #
##         ## initialize the base
##         #
##         if   isinstance ( xvar , ROOT.TH1   ) :
##             m_title = xvar.GetTitle ()            
##             xvar    = xvar.xminmax  ()
##         elif isinstance ( xvar , ROOT.TAxis ) :
##             xvar    = xvar.GetXmin() , mass.GetXmax()
            
##         ## create the variable 
##         if isinstance ( xvar , tuple ) and 2 == len(mass) :  
##             xvar = self.make_var ( xvar       , ## var 
##                              m_name     , ## name 
##                              m_title    , ## title/comment
##                              *mass      , ## min/max 
##                              fix = None ) ## fix ? 
##         elif isinstance ( xvar , ROOT.RooAbsReal ) :
##             xvar = self.make_var ( xvar       , ## var 
##                              m_name     , ## name 
##                              m_title    , ## title/comment
##                              fix = None ) ## fix ? 
##         else :
##             raise AttributeError("Swanson: Unknown type of 'xvar' parameter %s/%s" % ( type ( xvar ) , xvar ) )

            
##         PDF.__init__  ( self , name , xvar , None , None    ) 
        
##         self.__swanson = swanson 
##         beta_max       = max ( swanson.mmin() , swanson.cusp() )
##         self.__beta0   = self.make_var ( beta0 , 
##                                    'b0_swanson_%s'   % name ,
##                                    'b0_swanson(%s)'  % name ,
##                                    beta0 , 
##                                    0 , beta_max )
##         ## create PDF 
##         self.pdf = Ostap.Models.Swanson ( 
##             "Swanson_"    + name ,
##             "Swanson(%s)" % name ,
##             self.xvar    ,
##             self.beta0   ,
##             self.swanson )
        
##         ## save the configuration
##         self.config = {
##             'name'        : self.name      ,
##             'swanson'     : self.swanson   ,
##             'xvar'        : self.xvar      ,
##             'beta0'       : self.beta0     ,
##             }
    
##     @property
##     def beta0 ( self ) :
##         """'beta0'-parameter for Swanson function"""
##         return self.__beta0
##     @beta0.setter
##     def beta0 ( self, value ) :
##         value = float ( value )
##         self.__beta0.setVal ( value ) 

##     @property
##     def swanson ( self ) :
##         """'swanson'-function itself for Swanson PDF"""
##         return self.__swanson

## models.append ( Swanson_pdf )


# =============================================================================
if '__main__' == __name__ : 
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
