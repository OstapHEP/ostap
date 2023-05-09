#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/models.py
#
#  Set of useful PDFs for various 1D and 2D fits
#  It includes
#  - some empricial PDFs to describe narrow peaks: Gauss, CrystalBall, ....
#  - some PDF to describe "wide" peaks: BreitWigner,LASS, Bugg, Flatter, ...
#  - some useful PDFs to describe smooth background: phase space ;
#    expo times polynomial; phase space times polynomial, ...
#  - set of smooth non-facrorizeable model for 2D fits 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2011-07-25
# 
# =============================================================================
"""Set of useful PDFs for various 1D and 2D fits

It includes:

Empricial PDFs to describe narrow peaks : 

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
  - skew Gaussian   ## temporarily disabled 
  - Bukin,
  - Student-T
  - bifurcated Student-T
  - sinh-asinh shape
  - Johnson-SU shape
  - Atlas shape
  - Slash shape
  - RasingCosine shape
  - Q-Gaussian shape
  - Hyperbolic distribution
  - Generalised Hyperbolic distribution
  - Hypatia shape 
  - Das/Higgs shape 
  - Asymmetric Laplace shape
  - Sech  shape
  - Losev  shape
  - Logistic, aka 'sech-squared' shape
  - Generalized Logistic Type IV
  
PDF to describe 'wide' peaks : 

  - BreitWigner
  - LASS
  - Bugg
  - Flatte
  - Swanson's S-wave cusp 
  - ...
 
- some useful PDFs to describe smooth background in 1D : 

  - phase space 
  - expo times polynomial
  - phase space times polynomial
  - ...
  
- set of smooth non-facrorizeable models for 2D fits 
- ...

- generic convolution PDF 
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
    'Needham_pdf'            , ## Needham function for J/psi or Y (CB function with alpha=alpha(sigma))
    'Apollonios_pdf'         , ## Apollonios function         
    'Apollonios2_pdf'        , ## Apollonios function         
    'BifurcatedGauss_pdf'    , ## bifurcated Gauss
    'DoubleGauss_pdf'        , ## double Gauss
    'GenGaussV1_pdf'         , ## generalized normal v1  
    'GenGaussV2_pdf'         , ## generalized normal v2 
    'SkewGauss_pdf'          , ## skewed gaussian
    'Bukin_pdf'              , ## generic Bukin PDF: skewed gaussian with exponential tails     
    'StudentT_pdf'           , ## Student-T function 
    'BifurcatedStudentT_pdf' , ## bifurcated Student-T function 
    'SinhAsinh_pdf'          , ## "Sinh-arcsinh distributions". Biometrika 96 (4): 761
    'JohnsonSU_pdf'          , ## Johnson-SU distributon
    'Atlas_pdf'              , ## modified gaussian wiht exponential tails 
    'Slash_pdf'              , ## Symmetie peak with vey heavy tails 
    'RaisingCosine_pdf'      , ## Raising cosine distribution
    'QGaussian_pdf'          , ## Q-gaussian distribution
    'Hyperbolic_pdf'         , ## Hyperbolic distribution
    'GenHyperbolic_pdf'      , ## Generalised Hyperbolic distribution
    'Hypatia_pdf'            , ## Hypatia shape 
    'Das_pdf'                , ## Das/Higgsshape 
    'AsymmetricLaplace_pdf'  , ## asymmetric laplace 
    'Sech_pdf'               , ## hyperbolic secant (inverse-cosh) distribution
    'Losev_pdf'              , ## asymmetric hyperbolic secant distribution
    'Logistic_pdf'           , ## Logistic aka 'sech-squared' PDF
    'GenLogisticIV_pdf'      , ## generalized Logistic Type IV with location/scale 
    #
    ## specializations:
    # 
    'D0_pdf'  , ## PDF for D0        : Bukin 
    'Dp_pdf'  , ## PDF for D+        : Bukin
    'Ds_pdf'  , ## PDF for Ds+       : Bukin 
    'Lc_pdf'  , ## PDF for Lambda_c+ : Gauss
    #
    'B0_pdf'  , ## pdf for B0        : double-sided Crystal Ball 
    'Bd_pdf'  , ## pdf for B0        : double-sided Crystal Ball 
    'Bu_pdf'  , ## pdf for B+        : double-sided Crystal Ball 
    'Bs_pdf'  , ## pdf for Bs        : double-sided Crystal Ball 
    'Bc_pdf'  , ## pdf for Bc+       : double-sided Crystal Ball 
    #
    'Manca_pdf'   , ## Manca function to fit Y->mu mu spectrum  [Y(1S),Y(2S),Y(3S)]
    'Manca2_pdf'  , ## Manca function to fit Y->mu mu spectrum  [Y(1S),Y(2S),Y(3S)]
    'MancaX_pdf'  , ## function for 2D-fit of [Y(1S),Y(2S),Y(3S)]+X
    #
    ## pdfs for "wide" peaks, to be used with care - phase space corrections are large!
    # 
    'BreitWigner_pdf'      , ## (relativistic) 2-body Breit-Wigner
    'Flatte_pdf'           , ## Flatte-function  (pipi)
    'LASS_pdf'             , ## kappa-pole
    'Bugg_pdf'             , ## sigma-pole
    ##
    'Voigt_pdf'            , ## Voigt-profile 
    'PseudoVoigt_pdf'      , ## Voigt-profile 
    'BWI_pdf'              , ## Breit-Wigner with interference 
    'BW3L_pdf'             , ## BW3L
    'BWPS_pdf'             , ## BWPS
    'BWMC_pdf'             , ## BWMC
    #
    ## "Other" distributions 
    #
    'GammaDist_pdf'         , ## Gamma-distributuon in shape/scale parameterization
    'GenGammaDist_pdf'      , ## Generalized Gamma-distribution
    'Amoroso_pdf'           , ## another view of generalized Gamma distribution
    'LogGammaDist_pdf'      , ## Gamma-distributuon in shape/scale parameterization
    'Log10GammaDist_pdf'    , ## Gamma-distributuon in shape/scale parameterization
    'LogGamma_pdf'          , ## 
    'BetaPrime_pdf'         , ## Beta-prime distribution 
    'Landau_pdf'            , ## Landau distribution 
    'Argus_pdf'             , ## Argus distribution 
    'TwoExpos_pdf'          , ## Difference of two exponents
    'SinhAsinh_pdf'         , ## "Sinh-asinh" distribution
    'Gumbel_pdf'            , ## Gumber distribution
    'Weibull_pdf'           , ## Weibull distribution
    'Argus_pdf'             , ## Argus distribution 
    #
    ## 1D-background models
    # 
    'Bkg_pdf'              , ## Background: exponential modified by positive polynom
    'PolyPos_pdf'          , ## Background: positive polynom
    'PolyEven_pdf'         , ## Background: positive even polynom
    'Monotonic_pdf'        , ## Background: positive monotonic polynom
    'Convex_pdf'           , ## Background: positive monotonic polynom with fixed sign second derivative
    'ConvexOnly_pdf'       , ## Background: positive monotonic polynom with fixed sign second derivative
    'Sigmoid_pdf'          , ## Background: sigmoid modulated by positive polynom
    'TwoExpoPoly_pdf'      , ## Background: difference of two exponents modulated by polynom
    'PSPol_pdf'            , ## phase space modulated by positive polynomial
    'PSLeftExpoPol_pdf'    , ## phase space modulated by positive polynomial
    'PS2_pdf'              , ## 2-body phase space (no parameters)
    'PSLeft_pdf'           , ## Low  edge of N-body phase space 
    'PSRight_pdf'          , ## High edge of L-body phase space from N-body decays  
    'PSNL_pdf'             , ## L-body phase space from N-body decays  
    'PS23L_pdf'            , ## 2-body phase space from 3-body decays with orbital momenta
    'PSpline_pdf'          , ## positive spline (B-spline)
    'MSpline_pdf'          , ## positive monotonic spline 
    'CSpline_pdf'          , ## positive monotonic convex or concave spline 
    'CPSpline_pdf'         , ## positive convex or concave spline
    ##
    'PSSmear_pdf'          , ## smeared PhaseSpace-based PDF 
    ##
    'Linear_pdf'           , ## positive linear polynom 
    'Parabolic_pdf'        , ## positive parabolic polynom  
    ## the native RooFit background shapes
    'RooPoly_pdf'          , ## wrapper for RooPolynomial 
    'RooCheb_pdf'          , ## wrapper for RooChebyshev 
    ##
    'CutOff_pdf'           , ## Cut-off 
    'CutOffGauss_pdf'      , ## Gaussian cut-off 
    'CutOffStudent_pdf'    , ## Student's t-like /power law cut-off 
    #
    ## 2D non-factorazable models
    #
    'PolyPos2D_pdf'   , ## A positive polynomial in 2D  
    'PSPol2D_pdf'     , ## Product of phase spaces, modulated with 2D polynomial
    'PSPol2D2_pdf'    , ## Product of phase spaces, modulated with 2D polynomial
    'ExpoPSPol2D_pdf' , ## Exponential times  phase space times positive 2D-polynomial
    'ExpoPol2D_pdf'   , ## Product of exponents times positive 2D-polynomial
    'Spline2D_pdf'    , ## generic 2D positive spline 
    #
    ## 2D non-factorazable symmetric models
    #
    'PolyPos2Dsym_pdf', ## A positive symmetric polynomial in 2D
    'PSPol2Dsym_pdf'  , ## Symmetric product of phase spaces, modulated with 2D polynomial
    'PSPol2D2sym_pdf' , ## Symmetric product of phase spaces, modulated with 2D polynomial
    'ExpoPol2Dsym_pdf', ## Symmetric version of above
    'Spline2Dsym_pdf' , ## Symmetric 2D positive spline
    ##
    'Gauss2D_pdf'     , ## 2D gauss 
    #
    ## models for Pt-spectra fitting
    #
    'Tsallis_pdf'       , ## useful model for fitting pT-spectra 
    'QGSM_pdf'          , ## useful model for fitting pT-spectra
    #
    ## Non-parametric PDFs
    #
    'RooKeys1D_pdf'     , ## 1D-wrapper for RooNDKeysPdf 
    'RooKeys2D_pdf'     , ## 2D-wrapper for RooNDKeysPdf 
    'RooKeys3D_pdf'     , ## 3D-wrapper for RooNDKeysPdf 
    #
    ## simultaneous  fit
    'SimFit'            , ##    Simultaneous fit 
    'Sim1D'             , ## 1D-simultaneous fit 
    # 
    ## helpers
    #
    'Fit1D'             , ## generic model for                1D-fit
    'Fit2D'             , ## generic model for                2D-fit
    'Fit2DSym'          , ## generic model for (symmetric)    2D-fit
    'Fit3D'             , ## generic model for                3D-fit
    'Fit3DSym'          , ## generic model for (symmetric)    3D-fit
    'Fit3DMix'          , ## generic model for (mix-symmetry) 3D-fit
    ##
    'PolyPos3D_pdf'     , ## A positive polynomial in 3D  
    'PolyPos3Dsym_pdf'  , ## A positive symmetric polynomial in 3D
    'PolyPos3DmixXY_pdf', ## A positive partly symmetric (x<-->y) polynomial in 3D 
    'PolyPos3DmixYZ_pdf', ## A positive partly symmetric (y<-->z) polynomial in 3D 
    'PolyPos3DmixXZ_pdf', ## A positive partly symmetric (x<-->z) polynomial in 3D 
    ## 
    'Generic1D_pdf'     , ## wrapper over imported RooFit (1D)-pdf  
    'Generic2D_pdf'     , ## wrapper over imported RooFit (2D)-pdf
    'Generic3D_pdf'     , ## wrapper over imported RooFit (2D)-pdf
    ##
    'Convolution_pdf'   , ## generic convolution PDF
    ##
    )
# =============================================================================
from   ostap.logger.logger          import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.fitting.models' )
else                       : logger = getLogger ( __name__               )
# =============================================================================
from ostap.fitting.pdfbasic      import *
logger.debug ("Import 1D-fit machinery            from 'fit1d'"         )
from ostap.fitting.fit1d         import *
logger.debug ("Import 2D-fit machinery            from 'fit2d'"         )
from ostap.fitting.fit2d         import *
logger.debug ("Import 2D background        models from 'models_2d'"     )
from ostap.fitting.models_2d     import *  
logger.debug ("Import 3D-fit machinery            from 'fit3d'"         )
from ostap.fitting.fit3d         import *
logger.debug ("Import 3D background        models from 'models_3d'"     )
from ostap.fitting.models_3d     import *  
logger.debug ("Import simultaneous fit            from 'simfit'"        )
from ostap.fitting.simfit        import Sim1D, SimFit 
logger.debug ("Import convolution          models from 'convolution'"   )
from ostap.fitting.convolution   import *  
logger.debug ("Import resoltuion           models from 'resolution'"    )
from ostap.fitting.resolution    import *  
logger.debug ("Import adjustment           models from 'adjust'"        )
from ostap.fitting.adjust        import *  
logger.debug ("Import modifiers                   from 'modifiers'"     )
from ostap.fitting.modifiers     import *
logger.debug ("Import functions                   from 'roofuncs'"      )
from ostap.fitting.roofuncs      import * 
logger.debug ("Import signal     (peaking) models from 'signals'"       )
from ostap.fitting.signals       import * 
logger.debug ("Import background (smooth)  models from 'background'"    )
from ostap.fitting.background    import * 
logger.debug ("Import specialized models          from 'specific'"      )
from ostap.fitting.specific      import *
logger.debug ("Import 'other'     models          from 'distributions'" )
from ostap.fitting.distributions import *

models = []
from ostap.fitting.signals       import models as _models 
models += _models
from ostap.fitting.background    import models as _models
models += _models
from ostap.fitting.distributions import models as _models
models += _models
from ostap.fitting.specific      import models as _models 
models += _models
from ostap.fitting.resolution    import models as _models 
models += _models
from ostap.fitting.models_2d     import models as _models 
models += _models
from ostap.fitting.models_3d     import models as _models 
models += _models


# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger , symbols = models )


# =============================================================================
##                                                                      The END 
# =============================================================================
