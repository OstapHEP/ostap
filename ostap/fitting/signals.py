#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/signals.py
#
#  Set of useful PDFs for various ``signal'' 1D and 2D fits
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
  - skew Gaussian  (temporarily removed)
  - Bukin,
  - Student-T
  - bifurcated Student-T
  - SinhAsinh_pdf   
  - JohnsonSU_pdf   
  - Atlas_pdf
  - Slash_pdf
  - AsymmetricLaplace_pdf  
  - Sech_pdf   
  - Logistic_pdf   
  - RaisingCosine_pdf
  - QGaussian_pdf
  
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
    'SkewGauss_pdf'          , ## skewed gaussian (temporarily removed)
    'Bukin_pdf'              , ## generic Bukin PDF: skewed gaussian with exponential tails
    'StudentT_pdf'           , ## Student-T function 
    'BifurcatedStudentT_pdf' , ## bifurcated Student-T function
    'SinhAsinh_pdf'          , ## "Sinh-arcsinh distributions". Biometrika 96 (4): 761
    'JohnsonSU_pdf'          , ## JonhsonSU-distribution 
    'Atlas_pdf'              , ## modified gaussian with exponenital tails 
    'Slash_pdf'              , ## symmetric peakk wot very heavy tails 
    'RaisingCosine_pdf'      , ## Raising  Cosine distribution
    'QGaussian_pdf'          , ## Q-gaussian distribution
    'AsymmetricLaplace_pdf'  , ## asymmetric laplace 
    'Sech_pdf'               , ## hyperboilic secant  (inverse-cosh) 
    'Logistic_pdf'           , ## Logistic aka "sech-squared"   
    #
    ## pdfs for "wide" peaks, to be used with care - phase space corrections are large!
    # 
    'BreitWigner_pdf'      , ## (relativistic) 2-body Breit-Wigner
    'Flatte_pdf'           , ## Flatte-function  (pipi/KK)
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
from   ostap.fitting.basic import MASS, PDF  
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
    def __init__ ( self          ,
                   name          ,
                   xvar          ,
                   mean   = None ,
                   sigma  = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , xvar , mean , sigma  )
        
        #
        ## build pdf
        # 
        self.pdf = ROOT.RooGaussian (
            'gauss_%s'  % name ,
            "Gauss(%s)" % name ,
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
class CrystalBall_pdf(MASS) :
    """Crystal Ball function
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
                   xvar             ,
                   mean     = None  , 
                   sigma    = None  ,
                   alpha    = None  ,
                   n        = None  ) : 
                   
        #
        ## initialize the base
        #
        MASS.__init__ ( self , name , xvar , mean , sigma    )
        
        self.__alpha = self.make_var ( alpha ,
                                 'alpha_%s'        % name ,
                                 '#alpha_{CB}(%s)' % name ,  alpha  ,
                                 2.0 , 0.01  ,  5 )
        
        self.__n     = self.make_var ( n   ,
                                 'n_%s'            % name ,
                                 'n_{CB}(%s)'      % name , n       ,
                                 1.0 , 1.e-8 , 50 )
        
        #
        ## finally build PDF 
        #
        self.pdf = Ostap.Models.CrystalBall (
            'cb_%s'           % name ,
            'CrystalBall(%s)' % name ,
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
        value = float ( value )
        assert 0.1 <= value <= 5 , "``alpha'' must be between 0.1 and 5.0"
        self.__alpha.setVal ( value )
    
    @property
    def n ( self ) :
        """N-parameter for Crystal Ball tail (actually ``n+1'' is used)"""
        return self.__n
    @n.setter
    def n ( self, value ) :
        value = float ( value )
        assert 1.e-6 <= value <= 40 , "``n'' must be between 1.e-6 and 40"        
        self.__n.setVal ( value )

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
            'cbrs_%s'           % name ,
            'CrystalBallRS(%s)' % name ,
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
        MASS.__init__  ( self , name , xvar , mean , sigma )
        #
        ## treat the specific parameters
        #
        self.__aL    = self.make_var ( alphaL                  ,
                                 "aL_%s"          % name ,
                                 "#alpha_{L}(%s)" % name , alphaL    , 2.0 ,  0.01 ,  5 )
        self.__nL    = self.make_var ( nL                      ,                     
                                 "nL_%s"          % name ,
                                 "n_{L}(%s)"      % name , nL        , 1   , 1.e-8 , 50 )
        self.__aR    = self.make_var ( alphaR ,
                                 "aR_%s"          % name ,
                                 "#alpha_{R}(%s)" % name , alphaR    , 2.0 , 0.01  ,  5 )
        self.__nR    = self.make_var ( nR                      ,
                                 "nR_%s"          % name ,
                                 "n_{R}(%s)"      % name , nR        , 1   , 1.e-8 , 50 )
        
        self.pdf = Ostap.Models.CrystalBallDS(
            "cb2_"       + name ,
            "CB_{2}(%s)" % name ,
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
        value = float ( value )
        assert 0.01 < value < 5 , "alpha_L must be between 0.01 and 5" 
        self.__aL.setVal ( value )

    @property
    def nL ( self ) :
        """(left) N-parameter for Crystal Ball tail"""
        return self.__nL
    @nL.setter
    def nL ( self, value ) :
        value = float ( value )
        assert 1.e-6 < value < 40 , "N_L must be between 1.e-6 and 40" 
        self.__nL.setVal ( value )

    @property
    def aR ( self ) :
        """(right) Alpha-parameter for Crystal Ball tail"""
        return self.__aR
    @aR.setter
    def aR ( self, value ) :
        value = float ( value )
        assert 0.01 < value < 5 , "alpha_R must be between 0.01 and 5" 
        self.__aR.setVal ( value )

    @property
    def nR ( self ) :
        """(right) N-parameter for Crystal Ball tail"""
        return self.__nR
    @nR.setter
    def nR ( self, value ) :
        value = float ( value )
        assert 1.e-6 < value < 40 , "N_R must be between 1.e-6 and 40" 
        self.__nR.setVal ( value )

        
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
                   xvar                 ,
                   mean     = 3.096     ,   ## GeV  
                   sigma    = 0.013     ,   ## GeV 
                   a0       = 1.975     ,
                   a1       = -0.0011   ,   ## GeV^-1
                   a2       = -0.00018  ) : ## GeV^-2  
        
        MASS.__init__ ( self , name , xvar , mean , sigma )
        
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
                              "a_{0}(%s)" % name  , a0 , 
                              1.975               ,   0           , 10           )
        self.__a1 = self.make_var ( a1                  ,
                              "a1_%s"     % name  ,
                              "a_{1}(%s)" % name  , a1 , 
                              -0.0011   * unit    , -10 * unit    , 10 * unit    )
        self.__a2 = self.make_var ( a2                  ,
                              "a2_%s"     % name  ,
                              "a_{2}(%s)" % name  , a2 , 
                              -0.00018  * unit**2 , -10 * unit**2 , 10 * unit**2 )
        #
        self.pdf = Ostap.Models.Needham (
            'needham_%s'  % name ,
            'needham(%s)' % name ,
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
        """``a0''-parameter for Needham function"""
        return self.__a0
    @a0.setter
    def a0 ( self, value ) :
        value = float ( value ) 
        self.__a0.setVal ( value )

    @property
    def a1 ( self ) :
        """``a1''-parameter for Needham function"""
        return self.__a1
    @a1.setter
    def a1 ( self, value ) :
        value = float ( value ) 
        self.__a1.setVal ( value )

    @property
    def a2 ( self ) :
        """``a2''-parameter for Needham function"""
        return self.__a2
    @a2.setter
    def a2 ( self, value ) :
        value = float ( value ) 
        self.__a2.setVal ( value )

            
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
    - typical value of parameters alpha for ``physical'' peaks is 1.5<alpha<2.1,
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
        MASS.__init__  ( self , name , xvar , mean , sigma ) 
        
        self.__alpha = self.make_var ( alpha                     ,
                                 'alpha_%s'         % name ,
                                 '#alpha_{Apo}(%s)' % name , alpha , 
                                 2.0   , 0.01 ,   5 )
        
        self.__n     = self.make_var ( n                    ,
                                 'n_%s'        % name ,
                                 'n_{Apo}(%s)' % name , n ,
                                 2.0   , 1.e-6 , 50 )
        
        self.__b     = self.make_var ( b                    ,
                                 'b_%s'        % name ,
                                 'b_{Apo}(%s)' % name ,  b  ,
                                 1         , 0.001 , 10000 ) 
        
        #
        ## finally build PDF
        #
        self.pdf  = Ostap.Models.Apolonios (
            "apolo_"        + name ,
            "Apolonios(%s)" % name ,
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
        

    ## make a clone of this PDF    
    def clone ( self , name  = '' , xvar = None ) :
        """Make a ``clone''  of this PDF
        >>> sig1 = ...
        >>> sig2 = sig1.clone ( xvar = yvar ) 
        """
        return Apolonious_pdf ( name if name else self.name + '_copy'  ,
                                xvar if xvar else self.xvar            ,
                                self.mean  , 
                                self.sigma ,
                                self.alpha ,
                                self.n     ,
                                self.b     )
    @property
    def alpha ( self ) :
        """``alpha''-parameter for Apolonious tail"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        value = float ( value )
        assert 0.01 <=  value <= 5 , "``alpha''-parameter must be in 0.01,5 interval"
        self.__alpha.setVal ( value )        
    
    @property
    def n ( self ) :
        """``n''-parameter for Apolonios tail"""
        return self.__n
    @n.setter
    def n ( self, value ) :
        value = float ( value ) 
        assert 1.e-6 < value < 50 , "``n''-parameter must be in 1.e-6,50 interval"
        self.__n.setVal ( value )

    @property
    def b ( self ) :
        """``b''-parameter for Apolonios function"""
        return self.__b
    @b.setter
    def b ( self, value ) :
        value = float ( value )
        assert 0 < value , "``b''-parameter must be positive"
        self.__b.setVal ( value )


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
#          \dfrac{x-\mu}{\sigma_l} & \text{for} & x \le \mu \\
#          \dfrac{x-\mu}{\sigma_r} & \text{for} & x >   \mu \\
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
                   xvar               ,
                   mean      = None   ,
                   sigma     = None   ,
                   asymmetry = None   ,
                   beta      = None   ) : 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , xvar , mean , sigma  )

        
        self.__asym = self.make_var ( asymmetry                 ,
                                'asym_%s'          % name ,
                                '#asym_{Apo2}(%s)' % name , asymmetry , 0, -1 , 1  ) 
        
        self.__lst_R = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.__sigmaR = ROOT.RooFormulaVar (
            "sigmaR_%s"     % name   ,
            "sigma_{R}(%s)" % name   ,
            "%s*(1-%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self.__lst_R   )
        
        self.__lst_L = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.__sigmaL = ROOT.RooFormulaVar (
            "sigmaL_%s"     % name   ,
            "sigma_{L}(%s)" % name   ,
            "%s*(1+%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self.__lst_L   )
        
        self.__beta    = self.make_var ( beta ,
                                   'beta_%s'          % name  ,
                                   '#beta_{Apo2}(%s)' % name  ,
                                   beta , 0.01  , 1000 ) 
        #
        ## finally build PDF
        #
        self.pdf  = Ostap.Models.Apolonios2 (
            "apolo2_"        + name ,
            "Apolonios2(%s)" % name ,
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
        """``asymmetry''-parameter for Apolonious-2 function"""
        return self.__asym
    @asym.setter
    def asym ( self, value ) :
        value = float ( value ) 
        assert -1 <= value <=1 , "``asymmetry'' parameter is out of range -1,1"
        self.__asym.setVal ( value )

    @property
    def beta ( self ) :
        """``beta''-parameter for Apolonious-2 function"""
        return self.__beta
    @beta.setter
    def beta ( self, value ) :
        value =   float ( value )
        assert 0 < value , "``beta'' parameter must be positive"        
        self.__beta.setVal ( value )

    @property
    def sigmaL ( self ) :
        """(left)sigma-parameter for Apolonious-2 function"""
        return self.__sigmaL
    
    @property
    def sigmaR ( self ) :
        """(right)sigma-parameter for Apolonious-2 function"""
        return self.__sigmaR

    

models.append ( Apolonios2_pdf )    
# =============================================================================
## @class BifurcatedGauss_pdf
#  simple wrapper over bifurcated-gaussian
#  @see RooGaussian
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
class BifurcatedGauss_pdf(MASS) :
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
        MASS.__init__  ( self , name , xvar , mean , sigma ) 
        
        ## asymmetry parameter  
        self.__asym = self.make_var ( asymmetry                 ,
                                'asym_%s'          % name ,
                                '#asym_{asym}(%s)' % name ,
                                asymmetry , 0 , -1 , 1  ) 
        
        ## Right-side sigma
        self.__lst_R  = ROOT.RooArgList ( self.sigma , self.asym )
        self.__sigmaR = ROOT.RooFormulaVar (
            "sigmaR_%s"     % name   ,
            "sigma_{R}(%s)" % name   ,
            "%s*(1.0+%s)"   % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self.__lst_R   )
        
        ## Left-side sigma
        self.__lst_L  = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.__sigmaL = ROOT.RooFormulaVar (
            "sigmaL_%s"     % name   ,
            "sigma_{L}(%s)" % name   ,
            "%s*(1.0-%s)"   % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self.__lst_L   )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.BifurcatedGauss (
            "fbgau_"         + name ,
            "BifurGauss(%s)" % name ,
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
        """``asymmetry''-parameter for Bifurcated Gaussian"""
        return self.__asym
    @asym.setter
    def asym ( self, value ) :
        value = float ( value ) 
        assert -1 <= value <= 1, "``asymmetry''-parameter is out of range -1,1"
        self.__asym.setVal ( value )

    @property
    def sigmaL ( self ) :
        """(left)``sigma''-parameter for Bifurcated Gaussian"""
        return self.__sigmaL
    
    @property
    def sigmaR ( self ) :
        """(right)``sigma''-parameter for Bifurcated Gaussian"""
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
                   xvar                  ,
                   mean      = None      ,
                   sigma     = None      ,
                   fraction  = None      ,
                   scale     = None      ) :  
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , xvar , mean , sigma )
        
        self.__scale = self.make_var (
            scale ,
            'SigmaScale'        + name ,
            'SigmaScale(%s)'    % name , scale , 1 , 10 ) 
        
        ## the fraction 
        self.__fraction = self.make_var (
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
        value  = float ( value )
        assert 1 < scale, "``scale'' parameter must be >1"
        self.__scale.setVal ( value ) 

    @property
    def fraction ( self ) :
        """Fraction-parameter for Bifurcated Gaussian"""
        return self.__fraction
    @fraction.setter
    def fraction ( self, value ) :
        value = float ( value )
        assert 0 <= value <= 1 , "``fraction'' parameter must be 0<=f<=1"
        self.__fraction.setVal ( value )

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
                   xvar             ,
                   mean      = None ,
                   alpha     = None ,
                   beta      = 2    ) :  ## beta=2 is gaussian distribution 
        
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self  , name , xvar , mean  , alpha ) 
        
        #
        ## rename it!
        #
        self.__alpha = self.sigma
        if self.alpha != alpha : 
            sname  = self.alpha.GetName  ()
            stitle = self.alpha.GetTitle ()
            gname  = sname .replace ( 'sigma' , 'alpha' )
            gtitle = stitle.replace ( 'sigma' , 'alpha' )
            self.alpha.SetName  ( gname  ) 
            self.alpha.SetTitle ( gtitle )
            
        self.__beta  = self.make_var ( beta ,
                                 'beta_%s'        % name  ,
                                 '#beta_{v1}(%s)' % name  , beta , 
                                 2 , 1.e-4  , 1.e+6 ) 
        #
        ## finally build PDF
        #
        self.pdf = Ostap.Models.GenGaussV1 (
            "gengV1_"        + name ,
            "GenGaussV1(%s)" % name ,
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
        """``alpha''-parameter for Generalized V1 Gaussian (the same as ``sigma'')"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        value  = float ( value )
        assert 0 < value , "``alpha''-parameter must be positive!"
        self.__alpha.setVal ( abs ( value ) ) 

    @property
    def beta ( self ) :
        """``beta''-parameter for Generalized  V1 Gaussian"""
        return self.__beta
    @beta.setter
    def beta ( self, value ) :
        value  = float ( value )
        assert 0 < value , "``Beta''-parameter must be positive!"
        self.__beta.setVal ( abs ( value ) ) 
        
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
                   xvar               ,
                   mean        = None ,
                   alpha       = None ,
                   kappa       = 0    ) : ## 0 corresponds to gaussian distribution 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , xvar , mean , alpha ) 
        
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
        
        self.__xi    = self.mean 
        self.__kappa = self.make_var ( kappa ,
                                 'kappa_%s'        % name  ,
                                 '#kappa_{v2}(%s)' % name  , kappa , 
                                 0 , -4  , 4 ) 
        
        #
        ## finally build PDF
        #
        self.pdf = Ostap.Models.GenGaussV2 (
            "gengV2_"        + name ,
            "GenGaussV2(%s)" % name ,
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
        """alpha-parameter for Generalized V2 Gaussian (same as ``sigma'')"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        value  = float ( value )
        if value  <= 0 : raise AttributeError ( "Alpha parameter must be positive!") 
        self.__alpha.setVal ( abs ( value ) ) 

    @property
    def kappa ( self ) :
        """kappa-parameter for Generalized V2 Gaussian"""
        return self.__kappa
    @kappa.setter
    def kappa ( self, value ) :
        value  = float ( value )        
        self.__kappa.setVal ( value ) 
    
    @property
    def xi ( self ) :
        """xi-parameter (location) for Generalized V2 Gaussian (same  as ``mean'')"""
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
class SkewGauss_pdf(MASS) :
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
                   alpha     = 0    ) : ## alpha=0 correspond to gaussian 
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , xvar , mean , omega )
        
                
        #
        ## rename it!
        #
        
        self.__omega = self.sigma        
        sname      = self.sigma.GetName  ()
        stitle     = self.sigma.GetTitle ()
        gname      = sname .replace ( 'sigma' , 'omega' )
        gtitle     = stitle.replace ( 'sigma' , 'omega' )
        self.omega.SetName  ( gname  ) 
        self.omega.SetTitle ( gtitle )
        

        self.__xi   = self.mean 
        
        self.__alpha = self.make_var ( alpha ,
                                 'alpha_%s'   % name  ,
                                 '#alpha(%s)' % name  , alpha, 
                                 0 , -1000 , 1000  ) 
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.SkewGauss (
            "skewg_"         + name ,
            "SkewGauss(%s)" % name ,
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
        """xi-parameter (location) for Skew Gaussian (Same  as ``mean'')"""
        return self.__xi
    @xi.setter
    def xi ( self, value ) :
        value  = float ( value )
        self.__xi.setVal ( abs ( value ) ) 

    @property
    def omega ( self ) :
        """omega-parameter (scale/width) for Skew Gaussian (same as ``sigma'')"""
        return self.__omega
    @omega.setter
    def omega ( self, value ) :
        value  = float ( value )
        assert 0 < omega, 'Omega-parameter must be positive!'
        self.__omega.setVal ( value ) 

    @property
    def alpha( self ) :
        """alpha-parameter (shape/skew) for Skew Gaussian"""
        return self.__alpha
    @alpha.setter
    def alpha ( self, value ) :
        value  = float ( value )
        self.__alpha.setVal ( value ) 
        return self.__alpha.getVal ()
      

models.append ( SkewGauss_pdf )      
# =============================================================================
## @class Bukin_pdf
#  Bukin function, aka ``modified Novosibirsk function''
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
#  @see http://dx.doi.org/10.1007/JHEP06(2012)141     
#  @see Ostap::Math::Bukin
#  @see Analusis::Models::Bukin
#  @author Vanya BELYAEV Ivan.Belyaeve@itep.ru
#  @date 2011-07-25
class Bukin_pdf(MASS) :
    """Bukin function, aka ``modified Novosibirsk function'':
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
        MASS.__init__  ( self , name  , xvar , mean , sigma )
        #
        ## treat the specific parameters
        #
        
        ## asymmetry 
        self.__xi    = self.make_var ( xi                    ,
                                 "xi_%s"        % name ,
                                 "#xi(%s)"      % name , xi   , 0  , -1 , 1    )
        ## left tail
        self.__rhoL  = self.make_var ( rhoL                  ,
                                 "rhoL_%s"      % name ,
                                 "#rho_{L}(%s)" % name , rhoL , 0  ,  -1 , 10 )        
        ## right tail
        self.__rhoR  = self.make_var ( rhoR                  ,
                                 "rhoR_%s"      % name ,
                                 "#rho_{R}(%s)" % name , rhoR  , 0  ,  -1 , 10 )
        # 
        ## create PDF
        # 
        self.pdf = Ostap.Models.Bukin (
            "bkn_"      + name ,
            "Bukin(%s)" % name ,
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
        """``xi''-parameter (asymmetry) for Bukin function"""
        return self.__xi
    @xi.setter
    def xi ( self, value ) :
        value  = float ( value )
        assert -1 <= value <=1 , "Asymmetry must be in range  -1,1"
        self.__xi.setVal ( value ) 

    @property
    def rhoL ( self ) :
        """``rho''-parameter (left tail) for Bukin function"""
        return self.__rhoL
    @rhoL.setter
    def rhoL ( self, value ) :
        value  = float ( value )
        assert 0 <= value , "(left)``rho''-parameter must be non-negative"
        self.__rhoL.setVal ( value ) 

    @property
    def rhoR ( self ) :
        """``rho''-parameter (right tail) for Bukin function"""
        return self.__rhoR
    @rhoR.setter
    def rhoR ( self, value ) :
        value  = float ( value )
        assert 0 <= value , "(right)``rho''-parameter must be non-negative"
        self.__rhoR.setVal ( value ) 


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
                   xvar             ,
                   mean      = None ,
                   sigma     = None ,
                   n         = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , xvar , mean , sigma ) 
        
        # 
        self.__n  = self.make_var ( n                    ,
                              'n_%s'        % name ,
                              '#n_{ST}(%s)' % name , n , 
                              2 , 1.e-8 , 100  ) 
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.StudentT (
            "stT_"         + name ,
            "StudentT(%s)" % name ,
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
        """``n''-parameter for Student-T function (well, actually it is ``n+1'' is used)"""
        return self.__n
    @n.setter
    def n ( self, value ) :
        value  = float ( value )
        assert 1.e-6 < value , "``n''-parameter must be larger then 1.e-6"
        self.__n.setVal ( value ) 
    

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
                   xvar             ,
                   mean      = None ,
                   sigma     = None ,
                   asymmetry = None ,
                   nL        = None , 
                   nR        = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , xvar , mean , sigma )
        
        ##  asymmetry 
        self.__asym = self.make_var ( asymmetry               ,
                                'asym_%s'        % name ,
                                '#xi_{asym}(%s)' % name ,
                                asymmetry , 0 , -1 , 1  ) 
        
        self.__lst_R  = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.__sigmaR = ROOT.RooFormulaVar (
            "sigmaR_stt_%s"     % name   ,
            "sigma_{R}(%s)" % name   ,
            "%s*(1+%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self.__lst_R   )
        
        self.__lst_L  = ROOT.RooArgList ( self.sigma , self.asym ) 
        self.__sigmaL = ROOT.RooFormulaVar (
            "sigmaL_stt_%s"     % name   ,
            "sigma_{L}(%s)" % name   ,
            "%s*(1-%s)"     % ( self.sigma.GetName() , self.asym.GetName() ) ,
            self.__lst_L   )
        
        ## left exponent 
        self.__nL =  self.make_var ( nL                     ,
                               'nL_%s'         % name ,
                               '#nL_{BST}(%s)' % name , nL , 
                               2  , 1.e-6 , 100  )
        ## right exponent 
        self.__nR =  self.make_var ( nR                    ,
                               'nR_%s'         % name ,
                               '#nR_{BST}(%s)' % name , nR , 
                               2  , 1.e-6 , 100  ) 
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.BifurcatedStudentT (
            "bstT_"         + name ,
            "BStudentT(%s)" % name ,
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
        """``n''-parameter (left) for Bifurcated Student-T function  (actually ``n+1'' is used)"""
        return self.__nL
    @nL.setter
    def nL ( self, value ) :
        value  = float ( value )
        assert 1.e-6 < value , "(left)``n''-parameter must be lager than 1.e-6"
        self.__nL.setVal ( value ) 

    @property
    def nR ( self ) :
        """``n''-parameter (right) for Bifurcated Student-T function (actually ``n+1'' is used)"""
        return self.__nR
    @nR.setter
    def nR ( self, value ) :
        value  = float ( value )
        assert 1.e-6 < value , "(right)``n''-parameter must be lager than 1.e-6"
        self.__nR.setVal ( value ) 
        
    @property
    def asym ( self ) :
        """``asymmetry''-parameter for Bifurcated Student-T function"""
        return self.__asym
    @asym.setter
    def asym ( self, value ) :
        value  = float ( value )
        assert -1<= value <=1 , "Asymmetry must be in he  range -1,1" 
        self.__asym.setVal ( value ) 

    @property
    def sigmaL( self ) :
        """(left)``sigma''-parameter for Bifurcated Student-T function"""
        return self.__sigmaL

    @property
    def sigmaR( self ) :
        """(right)``sigma''-parameter for Bifurcated Student-T function"""
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
    """SinhAsing-function: 
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
                   xvar             ,
                   mean      = None , ## mu 
                   sigma     = None ,
                   epsilon   = 0    ,
                   delta     = 1    ) :

        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name  , xvar , mean , sigma )

        ##
        self.__mu      = self.mean
        self.__epsilon = self.make_var ( epsilon ,
                                   'epsilon_%s'   % name ,
                                   '#epsilon(%s)' % name , epsilon ,
                                   0 , -1000 , +1000 )
        self.__delta   = self.make_var ( delta ,
                                   'delta_%s'   % name ,
                                   '#delta(%s)' % name , delta ,
                                   1 , 1.e-6 , 1000   )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.SinhAsinh (
            "sinhaT_"        + name ,
            "SinhAsinh(%s)" % name ,
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
        """``epsilon''-parameter for Sinh-Asinh function"""
        return self.__epsilon
    @epsilon.setter
    def epsilon ( self, value ) :
        value = float ( vaalue ) 
        self.__epsilon.setVal ( value ) 

    @property
    def delta ( self ) :
        """``delta-parameter'' for Sinh-Asinh function"""
        return self.__delta
    @delta.setter
    def delta ( self, value ) :
        value = float ( value )
        assert 0  < value, "``delta''-parameter must be positive"         
        self.__delta.setVal ( value ) 
        
    @property
    def mu ( self ) :
        """``mu''-parameter (location) for Sinh-Asinh function (same as ``mean'')"""
        return self.__mu
    @mu.setter
    def mu (  self , value ) :
        value = float ( value )
        self.__mu.setVal ( value )
                
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
        
        MASS.__init__  ( self    , name , xvar , xi  , lambd   ) ## mean    , sigma  )

        self.__xi = self.mean
        if not self.xi is xi : ## newly created 
            oname   = self.xi.GetName ()
            otitle  = self.xi.GetTitle()
            nname   = oname .replace ( 'mean' , 'xi' )
            ntitle  = otitle.replace ( 'mean' , 'xi' )
            self.xi.SetName  ( nname  )
            self.xi.SetTitle ( ntitle )
            
        self.__lambd = self.sigma
        if not self.lambd is lambd : ## newly created 
            oname    = self.lambd.GetName ()
            otitle   = self.lambd.GetTitle()
            nname    = oname .replace ( 'sigma' , 'lambda' )
            ntitle   = otitle.replace ( 'sigma' , 'lambda' )
            self.lambd.SetName  ( nname  )
            self.lambd.SetTitle ( ntitle )
            self.lambd.setMax ( self.lambd.getMax() * 10 ) ## adjust it! 

            
        self.__delta   = self.make_var ( delta                 ,
                                   'delta_%s'     % name ,
                                   '#delta(%s)'   % name , delta ,
                                   1 , 1.e-6 , 1000   )
        
        self.__gamma   = self.make_var ( gamma               ,
                                   'gamma_%s'   % name ,
                                   '#gamma(%s)' % name , gamma ,
                                   0 , -1000 , +1000 )
        
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.JohnsonSU (
            "jSU_"          + name ,
            "JohnsonSU(%s)" % name ,
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
        """``delta''-parameter for Johnson-SU function"""
        return self.__delta
    @delta.setter
    def delta ( self, value ) :
        value = float ( value )
        assert   0 < value, "``delta''-parameter must be positive"
        self.__delta.setVal ( value ) 

    @property
    def gamma ( self ) :
        """``gamma''-parameter for Johnson-SU function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        self.__gamma.setVal ( value ) 
        
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
class Atlas_pdf(MASS) :
    """Modified gaussian with exponential tails
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
        
        MASS.__init__  ( self , name , xvar , mean, sigma  )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Atlas (
            "atlas_"    + name ,
            "ATLAS(%s)" % name ,
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
class Slash_pdf(MASS) :
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
        MASS.__init__  ( self , name , xvar , mean, scale )
        
        self.__scale = self.sigma
        if self.scale != scale : 
            sname  = self.scale.GetName  ()
            stitle = self.scale.GetTitle ()
            gname  = sname .replace ( 'sigma' , 'scale' )
            gtitle = stitle.replace ( 'sigma' , 'scale' )
            self.scale.SetName  ( gname  ) 
            self.scale.SetTitle ( gtitle )


        ## finally build pdf
        self.pdf = Ostap.Models.Slash (
            "slash_"    + name ,
            "Slash(%s)" % name ,
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
        """``mu''  - location parameter, the same as ``mean'' or ``location''"""
        return self.mean
    @property
    def location ( self ) :
        """``location''  - location parameter, the same as ``mean'' or ``mu''"""
        return self.mean
    @property
    def scale ( self ) :
        """``scale''  - scale parameter, the same as ``sigma''"""
        return self.__scale
    @scale.setter
    def scale ( self , value ) :
        value =  float ( value )
        assert 0 < scale , "``scale''-parameter must be positive"
        self.__scale.setVal ( value ) 
    
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
class AsymmetricLaplace_pdf(MASS) :
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
        MASS.__init__  ( self , name , xvar , mean , slope )
        
        self.__asym = self.make_var ( asymmetry               ,
                                'asym_%s'        % name ,
                                '#asym_{AL}(%s)' % name , asymmetry , 0 , -1 , 1  ) 

        self.__slope = self.sigma

        ## rename it if needed 
        if self.__slope != slope :
            sname  = self.slope.GetName  ()
            stitle = self.slope.GetTitle ()
            gname  = sname .replace ( 'sigma' , 'slope' )
            gtitle = stitle.replace ( 'sigma' , 'slope' )
            self.slope.SetName  ( gname  ) 
            self.slope.SetTitle ( gtitle )          
            
        ## Right-side lambda 
        self.__lst_R   = ROOT.RooArgList ( self.slope , self.asym )
        self.__lambdaR = ROOT.RooFormulaVar (
            "lambdaR_%s"     % name   ,
            "lambda_{R}(%s)" % name   ,
            "%s*(1.0+%s)"    % ( self.slope.GetName() , self.asym.GetName() ) ,
            self.__lst_R   )
        
        ## Left-side lambda 
        self.__lst_L   = ROOT.RooArgList ( self.slope , self.asym ) 
        self.__lambdaL = ROOT.RooFormulaVar (
            "lambdaL_%s"     % name   ,
            "lambda_{L}(%s)" % name   ,
            "%s*(1.0-%s)"   % ( self.slope.GetName() , self.asym.GetName() ) ,
            self.__lst_L   )
        
        ## finally build pdf
        self.pdf = Ostap.Models.AsymmetricLaplace (
            "alaplace_"    + name ,
            "ALaplace(%s)" % name ,
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
        """``mu''  - location parameter, the same as ``mean'' or ``location''"""
        return self.mean
    @property
    def location ( self ) :
        """``location''  - location parameter, the same as ``mean'' or ``mu''"""
        return self.mean

    @property
    def slope ( self ) :
        """``slope''-parameter the mean exponential slope,  the same as ``sigma''"""
        return self.__slope
    @slope.setter
    def slope ( self, value ) :
        value = float ( value ) 
        assert 0 < value , "``slope''-parameter must be positive"
        self.__slope.setVal ( value )
    
    @property
    def asym ( self ) :
        """``asymmetry''-parameter for Asymmetric Laplace"""
        return self.__asym
    @asym.setter
    def asym ( self, value ) :
        value = float ( value ) 
        assert -1 <= value <= 1, "``asymmetry''-parameter is out of range -1,1"
        self.__asym.setVal ( value )

    @property
    def lambdaL ( self ) :
        """(left)``lambda''-parameter (exponential slope) for Asymmetric Laplace"""
        return self.__lambdaL
    
    @property
    def lambdaR ( self ) :
        """(right)``lambda''-parameter (exponential slope) for Asymmetric Laplace"""
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
class Sech_pdf(MASS) :
    """Hyperbolic secant distribution or ``inverse-cosh'' distribution
    
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
        
        MASS.__init__  ( self , name , xvar , mean , sigma  )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Sech (
            "sech_"    + name ,
            "SECH(%s)" % name ,
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
class Logistic_pdf(MASS) :
    """ Logistic, aka ``sech-square'' PDF
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
        
        MASS.__init__  ( self , name , xvar , mean , sigma  )
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Logistic (
            "logistic_"    + name ,
            "Logistic(%s)" % name ,
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
class RaisingCosine_pdf(MASS) :
    r"""``Raising cosine'' distribution
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
        
        MASS.__init__  ( self , name , xvar , mean , scale  )
        
        self.__scale = self.sigma
        if self.__scale != scale : 
            sname  = self.sigma.GetName  ()
            stitle = self.sigma.GetTitle ()
            gname  = sname .replace ( 'sigma' , 'scale' )
            gtitle = stitle.replace ( 'sigma' , 'scale' )
            self.__scale.SetName  ( gname  ) 
            self.__scale.SetTitle ( gtitle )

        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.RaisingCosine (
            "rcos_"    + name ,
            "RCos(%s)" % name ,
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
        """``scale''-parameter, the same as ``sigma''"""
        return self.__scale
    @scale.setter
    def scale ( self, value ) :
        value = float ( value ) 
        assert 0 < value , "``scale''-parameter must be positive"
        self.__scale.setVal ( value ) 
    
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
class QGaussian_pdf(MASS) :
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
        
        MASS.__init__  ( self , name , xvar , mean , scale  )

        self.__scale = self.sigma
        if self.__scale != scale : 
            sname  = self.sigma.GetName  ()
            stitle = self.sigma.GetTitle ()
            gname  = sname .replace ( 'sigma' , 'scale' )
            gtitle = stitle.replace ( 'sigma' , 'scale' )
            self.__scale.SetName  ( gname  ) 
            self.__scale.SetTitle ( gtitle )

        ## Q 
        self.__q = self.make_var ( q               ,
                             'q_%s'   % name ,
                             '#q(%s)' % name , q , 1 , -1e+9 , 3-1.e-6 ) 
        
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.QGaussian (
            "qgauss_"    + name ,
            "QGauss(%s)" % name ,
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
        """``scale''-parameter, the same as ``sigma''"""
        return self.__scale
    @scale.setter
    def scale ( self, value ) :
        value = float ( value ) 
        assert 0 < value , "``scale''-parameter must be positive"
        self.__scale.setVal ( value )
        
    @property
    def q ( self ) :
        """``q''-parameter"""
        return self.__q
    @q.setter
    def q ( self, value ) :
        value = float ( value ) 
        assert 3 > value , "``q''-parameter must be <3"
        self.__q.setVal ( value ) 
    
models.append ( QGaussian_pdf )      

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
                   xvar             ,
                   mean      = None ,
                   sigma     = None ,
                   gamma     = None ) :
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , xvar , mean , sigma ) 

        limits_gamma = ()
        if  self.xminmax() :
            mn , mx = self.xminmax() 
            dm = mx - mn
            limits_gamma = 1.e-5 * dm , dm
            
        self.__gamma  = self.make_var ( gamma               ,
                                  'gamma_%s'   % name ,   
                                  '#gamma(%s)' % name , gamma ,
                                  *limits_gamma )
        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.Voigt (
            "vgt_"       + name ,
            "Voigt(%s)" % name ,
            self.xvar   ,
            self.mean   ,
            self.gamma  ,
            self.sigma  )

        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
            'sigma'     : self.sigma ,
            'gamma'     : self.gamma ,
            }
    
    @property
    def gamma ( self ) :
        """``gamma''-parameter for Voigt function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        assert 0 < value , "``gamma''-parameter must be positive"
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
                   xvar             ,
                   mean      = None ,
                   sigma     = None ,
                   gamma     = None ) :
        #
        ## initialize the base
        # 
        Voigt_pdf.__init__  ( self , name , xvar , mean , sigma , gamma ) 

        #
        ## finally build pdf
        # 
        self.pdf = Ostap.Models.PseudoVoigt (
            "pvgt_"           + name ,
            "PseudoVoigt(%s)" % name ,
            self.xvar   ,
            self.mean   ,
            self.gamma  ,
            self.sigma  )
        
        ## save the configuration
        self.config = {
            'name'      : self.name  ,
            'xvar'      : self.xvar  ,
            'mean'      : self.mean  ,
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
    ...                                  xvar  = mass  ,
    ...                                  mean  = m_X   ,
    ...                                  gamma = g_X   )
    
    Optional convolution with the resolution function is possible
    e.g. use gaussian resoltuion and fast-fourier convolution method:
    
    >>> breit = Models.BreitWigner_pdf ( 'BW'          ,
    ...                                  bw            ,
    ...                                  xvar  = mass  ,
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
                   xvar               ,
                   mean        = None , 
                   gamma       = None ,
                   convolution = None ,
                   useFFT      = True ) :
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name, xvar , mean , gamma )
        
        if gamma is not self.sigma :
            sname  = self.sigma.GetName  ()
            stitle = self.sigma.GetTitle ()
            gname  = sname .replace ( 'sigma' , 'Gamma' )
            gtitle = stitle.replace ( 'sigma' , 'Gamma' )
            self.sigma.SetName  ( gname  ) 
            self.sigma.SetTitle ( gtitle )
            if  self.xminmax () and self.sigma.minmax() : 
                mn  , mx  = self.xminmax() 
                dm = mx - mn
                smn , smx = self.sigma.minmax()                 
                self.sigma.setMin ( max (  1.e-5 * dm , smn ) )
                self.sigma.setMax ( min (  2     * dm , smx ) )
        
        #
        ## define the actual BW-shape using
        #      Ostap::Math::BreitWeigner object
        #
        self.__breitwigner = breitwigner  ## Ostap::Math::BreitWeigner object

        ## create PDF 
        self.__breit = Ostap.Models.BreitWigner ( 
            "rbw_"    + name ,
            "RBW(%s)" % name ,
            self.xvar        ,
            self.mean        ,
            self.gamma       ,
            self.breitwigner )

        
        self.__convolution = convolution  
        self.__useFFT      =  useFFT 
        if  None is convolution : self.pdf = self.__breit
        else :
            from ostap.fitting.convolution import Convolution            
            self.conv = Convolution ( 'RBW' + name  ,
                                      self.__breit  , self.xvar ,
                                      convolution   , useFFT    ) 
            self.pdf  = self.conv.pdf
            
        ## save the configuration
        self.config = {
            'name'        : self.name          ,
            'breitwigner' : self.breitwigner   ,
            'xvar'        : self.xvar          ,
            'mean'        : self.mean          ,
            'gamma'       : self.gamma         ,
            'convolution' : self.__convolution ,
            'useFFT'      : self.__useFFT      ,
            }

    @property
    def gamma ( self ) :
        """``gamma''-parameter for Breit-Wigner function (alias for ``sigma'')"""
        return self.sigma 
    @gamma.setter
    def gamma ( self, value ) :
        self.sigma = value 
    
    @property
    def Gamma ( self ) :
        """``Gamma''-parameter for Breit-Wigner function (alias for ``sigma'')"""
        return self.gamma 
    @Gamma.setter
    def Gamma ( self, value ) :
        self.sigma = value 

    @property
    def breitwigner ( self ) :
        """The Breit-Wigner function  itself"""
        return self.__breitwigner

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
                   xvar               ,
                   mean        = None , 
                   gamma       = None ) : 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , xvar , mean , gamma )
        
        self.__gamma = self.sigma
        if  self.gamma  != gamma : 
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
        self.__breitwigner = breitwigner  ## Ostap::Math::BW23L object
        
        ## create PDF 
        self.pdf = Ostap.Models.BW23L ( 
            "rbw23_"    + name ,
            "RBW23(%s)" % name ,
            self.xvar          ,
            self.mean          ,
            self.gamma         ,
            self.breitwigner   )

        ## save the configuration
        self.config = {
            'name'        : self.name          ,
            'breitwigner' : self.breitwigner   ,
            'xvar'        : self.xvar          ,
            'mean'        : self.mean          ,
            'gamma'       : self.gamma         ,
            }

    @property
    def gamma ( self ) :
        """Gamma-parameter for Breit-Wigner ``2-from-3'' function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        assert 0 < gamma , "``gamma'' must be positive"
        self.__gamma.setVal ( value ) 

    @property
    def breitwigner ( self ) :
        """The Breit-Wigner function itself"""
        return self.__breitwigner

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
#  @see Ostap::Math::Flatte2
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2014-01-18
class Flatte_pdf(MASS) :
    """Flatte function:
    S.M.Flatte, ``Coupled-channel analysis of the (pi eta)
    and (KbarK) systems near (KbarK) threshold'' 
    Phys. Lett. B63, 224 (1976
    http://www.sciencedirect.com/science/article/pii/0370269376906547

    Typical case:    f0(980) -> pi+ pi- & K+ K- shapes 
    """
    def __init__ ( self              ,
                   name              ,
                   flatte            ,    ## Ostap::Math::Flatte/Flatte2
                   xvar              ,
                   m0       = None   ,    ## the pole 
                   m0g1     = None   ,    ## m0*gamma_1 
                   g2og1    = None   ,    ## gamma2/gamma1 
                   gamma1   = None   ,    ## gamma1 
                   gamma2   = None   ,    ## gamma2 
                   gamma0   = None   ) :  ## gamma0 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , xvar , m0 , None )

        self.__flatte = flatte
            
        self.__m0 = self.mean
        if self.m0 != m0 :
            
            sname  = self.mean.GetName  ()
            stitle = self.mean.GetTitle ()
            gname  = sname .replace ( 'mean' , 'm0' )
            gtitle = stitle.replace ( 'mean' , 'm0' )
            self.m0.SetName  ( gname  ) 
            self.m0.SetTitle ( gtitle ) 
            
        self.__gamma0 = self.make_var  ( gamma0                  ,
                                         'gamma0_%s'      % name ,
                                         '#gamma_{0}(%s)' % name ,
                                         gamma0 , 0 , 0 , 5 * self.sigma.getVal() )
        
        if  gamma1 is None and gamma2 is None :
            
            vmin = 0.2 * self.mean.getMin () * self.gamma.getMin ()
            vmax = 2.0 * self.mean.getMax () * self.gamma.getMax ()
            
            self.__m0g1 = self.make_var  ( m0g1                          ,
                                           'm0g1_%s'             % name ,
                                           'm_{0}#gamma_{1}(%s)' % name ,
                                           m0g1 , m0g1 , vmin , xmax )
            
            self.__g2og1 = self.make_var ( g2og1    ,
                                           'g2og1_%s'                  % name ,
                                           '#gamma_{2}/#gamma_{1}(%s)' % name ,
                                           g2og1    ,  1  ,  0.01  , 100  ) 
            
            self.__lst1   = ROOT.RooArgList ( self.m0g1  , self.m0     ) 
            self.__gamma1 = ROOT.RooRealVar ( 
                'g1_%s'          % name ,
                '#gamma_{1}(%s)' % name ,
                '%s / %s '  % ( self.m0g1.GetName() , self.m0.GetName() ) , 
                self.__lst1  )
            self.__lst2   = ROOT.RooArgList ( self.g2og1 , self.gamma1 ) 
            self.__gamma2 = ROOT.RooRealVar ( 
                'g2_%s'          % name ,
                '#gamma_{2}(%s)' % name ,
                '%s * %s '  % ( self.g2og1.GetName() , self.gamma1.GetName() ) , 
                self.__lst2 )
            
        elif gamma1 is None : raise TypeError ( 'Flatte_pdf: gamma1 is not specified!' ) 
        elif gamma2 is None : raise TypeError ( 'Flatte_pdf: gamma2 is not specified!' ) 
        else :
            
            self.__gamma1 =  self.make_var  ( gamma1                   ,
                                              'g1_%s'          % name ,
                                              '#gamma_{1}(%s)' % name ,
                                              gamma1               ,
                                              self.gamma.getVal () ,
                                              self.gamma.getMin () ,
                                              self.gamma.getMax () )            
            self.__gamma2 =  self.make_var  ( gamma2                   ,
                                              'g2_%s'          % name ,
                                              '#gamma_{2}(%s)' % name ,
                                              gamma2   ,
                                              self.gamma.getVal () ,
                                              self.gamma.getMin () ,
                                              self.gamma.getMax () )
            
            self.__lst1  = ROOT.RooArgList ( self.m0     , self.gamma1 ) 
            self.__m0g1  = ROOT.RooFormulaVar (
                'm0g1_%s'             % name ,
                'm_{0}#gamma_{1}(%s)' % name ,
                '%s * %s ' % ( self.m0.GetName() , self.gamma1.GetName() ) ,
                self.__lst1 )
            self.__lst2  = ROOT.RooArgList ( self.gamma2 , self.gamma1 ) 
            self.__g2og1 = ROOT.RooFormulaVar ( 
                'g2og1_%s'                  % name ,
                '#gamma_{2}/#gamma_{1}(%s)' % name , g2og1 ,
                '%s / %s '  % ( self.gamma2.GetName() , self.gamma1.GetName() ) ,
                self.__lst2 )
                
        ## create PDF 
        self.pdf = Ostap.Models.Flatte ( 
            "flatte_"    + name ,
            "Flatte(%s)" % name ,
            self.xvar    ,
            self.m0      ,
            self.m0g1    ,
            self.g2og1   ,
            self.gamma0  ,
            self.flatte  )

        ## save the configuration
        self.config = {
            'name'        : self.name    ,
            'flatte'      : self.flatte  ,
            'xvar'        : self.xvar    ,
            'm0'          : self.m0      ,
            'gamma0'      : self.gamma0  ,
            }
        
        if gamma1 is None and gamma2 is None : 
            self.config.update ( { 'm0g1'   : self.m0g1   , 
                                   'g2og1'  : self.g2og1  } )
        else : 
            self.config.update ( { 'gamma1' : self.gamma1 ,
                                   'gamma2' : self.gamma2 } )
            
    @property
    def m0 ( self ) :
        """``m0''-parameter for Flatte-function (same as ``mean'')"""
        return self.__m0
    @m0.setter
    def m0  ( self, value ) :
        value = float ( value )
        self.__m0.setVal ( value ) 

    @property
    def m0g1 ( self ) :
        """``m0*g1''-parameter for Flatte-function"""
        return self.__m0g1
    @m0g1.setter
    def m0g1 ( self, value ) :
        assert not isinstance ( self.__m0g1 , ROOT.RooFormulaVar ),\
               "``m0g1''-parameter can't be set!"
        value = float ( value )
        self.__m0g1.setVal ( value ) 

    @property
    def g2og1 ( self ) :
        """``g2/g1''-parameter for Flatte-function"""
        return self.__g2og1
    @g2og1.setter
    def g2og1 ( self, value ) :
        assert not isinstance ( self.__g2og1 , ROOT.RooFormulaVar ),\
               "``g2og1''-parameter can't be set!"        
        value = float ( value )
        assert 0 < value, "``g2/g1''-parameter for Flatte-function must be positive"
        self.__g2og1.setVal ( value )

    @property
    def gamma1 ( self ) :
        "``gamma1''-parameter for Flatte-function"
        return self.__gamma1
    @gamma1.setter
    def gamma1 ( self , value ) :
        assert not isinstance ( self.__gamma1 , ROOT.RooFormulaVar ),\
               "``gamma1''-parameter can't be set!"
        value = float ( value )
        self.__gamma1.setVal ( value ) 

    @property
    def gamma2 ( self ) :
        "``gamma2''-parameter for Flatte-function"
        return self.__gamma2
    @gamma2.setter
    def gamma2 ( self , value ) :
        assert not isinstance ( self.__gamma2 , ROOT.RooFormulaVar ),\
               "``gamma2''-parameter can't be set!"
        value = float ( value )
        self.__gamma2.setVal ( value ) 

    @property
    def gamma0 ( self ) :
        "``gamma0''-parameter for Flatte-function"
        return self.__gamma0
    @gamma0.setter
    def gamma0 ( self , value ) :
        value = float ( value )
        self.__gamma0.setVal ( value ) 

    @property
    def flatte ( self ) :
        """The Flatte function itself"""
        return self.__flatte

models.append ( Flatte_pdf )                          

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
        MASS.__init__  ( self , name , xvar , m0 , g0 )
        
        self.__g0 = self.sigma
        if self.__g0 != g0  : 
            sname  = self.g0.GetName  ()
            stitle = self.g0.GetTitle ()
            gname  = sname .replace ( 'sigma' , 'g0' )
            gtitle = stitle.replace ( 'sigma' , 'g0' )
            self.g0.SetName  ( gname  ) 
            self.g0.SetTitle ( gtitle )
            
        self.__m0 = self.mean
        if self.__m0 !=  m0 : 
            sname  = self.m0.GetName  ()
            stitle = self.m0.GetTitle ()
            gname  = sname .replace ( 'mean' , 'm0' )
            gtitle = stitle.replace ( 'mean' , 'm0' )
            self.m0.SetName  ( gname  ) 
            self.m0.SetTitle ( gtitle ) 
            
        self.__a = self.make_var ( a                  ,
                                   'aLASS_%s'  % name ,
                                   "aLASS(%s)" % name , a , 
                                   1.94e-3            ,
                                   1.94e-3            ,
                                   1.94e-3            ) 
        self.__r = self.make_var ( r             ,
                                   'rLASS_%s'  % name ,
                                   "rLASS(%s)" % name , r , 
                                   1.76e-3            ,
                                   1.76e-3            ,
                                   1.76e-3            ) 
        self.__e = self.make_var ( e            ,
                                   'eLASS_%s'  % name ,
                                   "eLASS(%s)" % name , e ,
                                   1.0                , 
                                   1.0                ,
                                   1.0                )
        
        self.__mKaon = mKaon
        self.__mPion = mPion
        
        ## create PDF 
        self.pdf = Ostap.Models.LASS ( 
            "lass_"    + name ,
            "LASS(%s)" % name ,
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
        """``g0''-parameter for LASS-function (same as ``sigma'')"""
        return self.sigma
    @g0.setter
    def g0 ( self, value ) :
        self.sigma = value 

    @property
    def m0 ( self ) :
        """``m0''-parameter for LASS-function (same as ``mean'')"""
        return self.__gamma
    @m0.setter
    def m0 ( self, value ) :
        self.mean = value 

    @property
    def a ( self ) :
        """``a''-parameter for LASS-function"""
        return self.__a_
    @a.setter
    def a ( self, value ) :
        value = float ( value )
        self.__a.setVal ( value ) 

    @property
    def r ( self ) :
        """``r''-parameter for LASS-function"""
        return self.__r
    @r.setter
    def r ( self, value ) :
        value = float ( value )
        self.__r.setVal ( value ) 

    @property
    def e ( self ) :
        """``e''-parameter for LASS-function"""
        return self.__e
    @e.setter
    def e ( self, value ) :
        value = float ( value )
        self.__e.setVal ( value ) 

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
    """ The parameterization of sigma pole by
    B.S.Zou and D.V.Bugg, Phys.Rev. D48 (1993) R3948.
    http://dx.doi.org/10.1103/PhysRevD.48.R3948
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
        MASS.__init__  ( self , name , xvar , m , g2 ) 
        
        self.__bugg_g2 = self.sigma
        self.__gamma   = self.sigma
        ##
        if self.gamma != g2 :  
            sname  = self.gamma.GetName  ()
            stitle = self.gamma.GetTitle ()
            gname  = sname .replace ( 'sigma' , 'gamma2_Bugg' )
            gtitle = stitle.replace ( 'sigma' , 'Gamma2_Bugg' )
            self.bugg_g2.SetName  ( gname  ) 
            self.bugg_g2.SetTitle ( gtitle )
            self.gamma = self.bugg_g2 
            
        self.__bugg_m = self.mean
        if self.bugg_m != m :  
            sname  = self.mean.GetName  ()
            stitle = self.mean.GetTitle ()
            gname  = sname .replace ( 'mean' , 'm_Bugg' )
            gtitle = stitle.replace ( 'mean' , 'm_Bugg' )
            self.bugg_m.SetName  ( gname  ) 
            self.bugg_m.SetTitle ( gtitle ) 

        self.__bugg_b1 = self.make_var ( b1                  ,
                                         'b1Bugg_%s'  % name ,
                                         "b1Bugg(%s)" % name , b1 ,
                                         0.5848 , 
                                         0 , 2  )
        
        self.__bugg_b2 = self.make_var ( b2             ,
                                         'b2Bugg_%s'  % name ,
                                         "b2Bugg(%s)" % name , b2 , 
                                         1.6663 ,
                                         1 , 2  ) 
        
        self.__bugg_a  = self.make_var ( a             ,
                                         'aBugg_%s'  % name ,
                                         "aBugg(%s)" % name , a , 
                                         1.082    ,
                                         0.5 , 5  ) 
        
        self.__bugg_s1  = self.make_var ( s1           ,
                                          's1Bugg_%s'  % name ,
                                          "s1Bugg(%s)" % name , s1 , 
                                          2.8              ,
                                          1 , 5            ) 
        
        self.__bugg_s2  = self.make_var ( s2           ,
                                          's2Bugg_%s'  % name ,
                                          "s2Bugg(%s)" % name , s2 , 
                                          3.5              ,
                                          1 , 5            ) 
        
        self.__mPion = mPion 
        ## create PDF 
        self.pdf = Ostap.Models.Bugg ( 
            "bugg_"    + name ,
            "Bugg(%s)" % name ,
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
        """``gamma''-parameter (``g2'') for Bugg function"""
        return self.__gamma
    @gamma.setter
    def gamma ( self, value ) :
        value = float ( value )
        self.__gamma.setVal ( value ) 

    @property
    def g2 ( self ) :
        """``g2''-parameter for Bugg function"""
        return self.__bugg_g2
    @g2.setter
    def g2 ( self, value ) :
        value = float ( value )
        self.__bugg_g2.setVal ( value ) 

    @property
    def b1 ( self ) :
        """``b1''-parameter for Bugg function"""
        return self.__bugg_b1
    @b1.setter
    def b1 ( self, value ) :
        value = float ( value )
        self.__bugg_b1.setVal ( value ) 

    @property
    def b2 ( self ) :
        """``b2''-parameter for Bugg function"""
        return self.__bugg_b2
    @b2.setter
    def b1 ( self, value ) :
        value = float ( value )
        self.__bugg_b2.setVal ( value ) 

    @property
    def a ( self ) :
        """``a''-parameter for Bugg function"""
        return self.__bugg_a
    @a.setter
    def a ( self, value ) :
        value = float ( value )
        self.__bugg_a.setVal ( value ) 

    @property
    def s1 ( self ) :
        """``s1''-parameter for Bugg function"""
        return self.__bugg_s1
    @s1.setter
    def s1 ( self, value ) :
        value = float ( value )
        self.__bugg_s1.setVal ( value ) 

    @property
    def s2 ( self ) :
        """``s2''-parameter for Bugg function"""
        return self.__bugg_s2
    @s2.setter
    def s2 ( self, value ) :
        value = float ( value )
        self.__bugg_s2.setVal ( value ) 

        
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
class Swanson_pdf(PDF) :
    """ S-wave cusp
    - LHCb-PAPER-2016-019, Appendix
    - E. S. Swanson, Cusps and exotic charmonia, arXiv:1504.07952
    - http://arxiv.org/abs/1504.07952
    """
    def __init__ ( self              ,
                   name              ,
                   swanson           , ## Ostap::Math::Swanson objects 
                   xvar              ,
                   beta0    = None   ) : 
        #
        ## initialize the base
        #
        if   isinstance ( xvar , ROOT.TH1   ) :
            m_title = xvar.GetTitle ()            
            xvar    = xvar.xminmax  ()
        elif isinstance ( xvar , ROOT.TAxis ) :
            xvar    = xvar.GetXmin() , mass.GetXmax()
            
        ## create the variable 
        if isinstance ( xvar , tuple ) and 2 == len(mass) :  
            xvar = self.make_var ( xvar       , ## var 
                             m_name     , ## name 
                             m_title    , ## title/comment
                             *mass      , ## min/max 
                             fix = None ) ## fix ? 
        elif isinstance ( xvar , ROOT.RooAbsReal ) :
            xvar = self.make_var ( xvar       , ## var 
                             m_name     , ## name 
                             m_title    , ## title/comment
                             fix = None ) ## fix ? 
        else :
            raise AttributeError("Swanson: Unknown type of ``xvar'' parameter %s/%s" % ( type ( xvar ) , xvar ) )

            
        PDF.__init__  ( self , name , xvar , None , None    ) 
        
        self.__swanson = swanson 
        beta_max       = max ( swanson.mmin() , swanson.cusp() )
        self.__beta0   = self.make_var ( beta0 , 
                                   'b0_swanson_%s'   % name ,
                                   'b0_swanson(%s)'  % name ,
                                   beta0 , 
                                   0 , beta_max )
        ## create PDF 
        self.pdf = Ostap.Models.Swanson ( 
            "Swanson_"    + name ,
            "Swanson(%s)" % name ,
            self.xvar    ,
            self.beta0   ,
            self.swanson )
        
        ## save the configuration
        self.config = {
            'name'        : self.name      ,
            'swanson'     : self.swanson   ,
            'xvar'        : self.xvar      ,
            'beta0'       : self.beta0     ,
            }
    
    @property
    def beta0 ( self ) :
        """``beta0''-parameter for Swanson function"""
        return self.__beta0
    @beta0.setter
    def beta0 ( self, value ) :
        value = float ( value )
        self.__beta0.setVal ( value ) 

    @property
    def swanson ( self ) :
        """``swanson''-function itself for Swanson PDF"""
        return self.__swanson

models.append ( Swanson_pdf )

# =============================================================================
if '__main__' == __name__ : 
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )
    
# =============================================================================
# The END 
# =============================================================================
