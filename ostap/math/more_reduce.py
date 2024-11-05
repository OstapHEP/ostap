#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/reduce.py
#  Module with some useful utilities for reducing some math objects 
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2011-12-01"
__all__     = ()
# =============================================================================
from    ostap.math.base        import Ostap, doubles 
from    ostap.math.reduce      import root_factory, poly_factory 
from    ostap.core.ostap_types import sequence_types 
import  ROOT, array 
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.more_reduce' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================

# =============================================================================
## BreitWigner and friends
# =============================================================================

# =============================================================================
## Reduce Ostap::Math::FormFactors::Jackson
#  @see Ostap::Math::FormFactors::Jackson
def _omffj_reduce_ ( ff ) :
    """Reduce `Ostap.Math.FormFactors.Jackson`
    - see `Ostap.Math.FormFactors.Jackson`
    """
    return root_factory , ( type ( ff ) , ff.rho() )

Ostap.Math.FormFactors.Jackson . __reduce__ = _omffj_reduce_

# =============================================================================
## Reduce Ostap::Math::FormFactors::BlattWeisskopf
#  @see Ostap::Math::FormFactors::BlattWeiskopf
def _omffbw_reduce_ ( ff ) :
    """Reduce `Ostap.Math.FormFactors.BlattWeiskopf`
    - see `Ostap.Math.FormFactors.BlattWeiskopf`
    """
    return root_factory , ( type ( ff ) , ff.L() , ff.breakup () )

Ostap.Math.FormFactors.BlattWeisskopf . __reduce__ = _omffbw_reduce_

# =============================================================================
## Reduce Ostap::Math::FormFactors::NoFormFactor
#  @see Ostap::Math::FormFactors::NoFormFactor
def _omffnf_reduce_ ( ff ) :
    """Reduce `Ostap.Math.FormFactors.NoFormFactor`
    - see `Ostap.Math.FormFactors.NoFormFactor`
    """
    return root_factory , ( type ( ff ) , )

Ostap.Math.FormFactors.NoFormFactor. __reduce__ = _omffnf_reduce_


# =============================================================================
## Reduce Ostap::Math::ChannelCW
#  @see Ostap::Math::ChannelCW
def _omccw_reduce_ ( ch ) :
    """Reduce `Ostap.Math.ChannelCW`
    - see `Ostap.Math.ChannelCW`
    """
    return root_factory , ( type ( ch ) , ch.gamma0() , ch.m1() , ch.m2() )

Ostap.Math.ChannelCW     . __reduce__ = _omccw_reduce_
Ostap.Math.ChannelQ      . __reduce__ = _omccw_reduce_
Ostap.Math.ChannelFlatte . __reduce__ = _omccw_reduce_

# =============================================================================
## Reduce Ostap::Math::ChannelFlatteBugg
#  @see Ostap::Math::ChannelFlatteBugg
def _omcfb_reduce_ ( ch ) :
    """Reduce `Ostap.Math.ChannelFlatteBugg`
    - see `Ostap.Math.ChannelFlatteBugg`
    """
    return root_factory , ( type ( ch )    ,
                            ch.gamma0   () ,
                            ch.mcharged () ,
                            ch.mneutral () ,                            
                            ch.alpha    () ,
                            ch.fc       () , 
                            ch.fn       () )
Ostap.Math.ChannelFlatteBugg . __reduce__ = _omcfb_reduce_

# =============================================================================
## Reduce Ostap::Math::Channel
#  @see Ostap::Math::Channel
def _omc_reduce_ ( ch ) :
    """Reduce `Ostap.Math.Channel`
    - see `Ostap.Math.Channel`
    """
    return root_factory , ( type ( ch ) , ch.gamma0() , ch.m1() , ch.m2() , ch.L() , ch.formfactor () )

Ostap.Math.Channel. __reduce__ = _omc_reduce_

# =============================================================================
## Reduce Ostap::Math::Channel0
#  @see Ostap::Math::Channel0
def _omc0_reduce_ ( ch ) :
    """Reduce `Ostap.Math.Channel0`
    - see `Ostap.Math.Channel0`
    """
    return root_factory , ( type ( ch ) , ch.gamma0() , ch.m1() , ch.m2() , ch.L() , ch.formfactor () , ch.qs () )

Ostap.Math.Channel0. __reduce__ = _omc0_reduce_

# =============================================================================
## reduce Ostap::Math::Channel23L
#  @see Ostap::Math::Channel23L
def _omc23l_reduce_ ( ch ) : 
    """Reduce `Ostap.Math.Channel23L`
    - see `Ostap.Math.Channel23L`
    """
    return root_factory , ( type ( ch ) , ch.channel() , ch.ps23L () ) 
   
Ostap.Math.Channel23L. __reduce__ = _omc23l_reduce_

# =============================================================================
## reduce Ostap::Math::ChanneNR3
#  @see Ostap::Math::ChannelNR3
def _omcnr3_reduce_ ( ch ) : 
    """Reduce `Ostap.Math.ChannelNR3`
    - see `Ostap.Math.ChannelNR3`
    """
    return root_factory , ( type ( ch ) , ch.m1() , ch.m2() , ch.m3() )
   
Ostap.Math.ChannelNR3. __reduce__ = _omcnr3_reduce_

# =============================================================================
## reduce Ostap::Math::ChanneGS
#  @see Ostap::Math::ChannelGS
def _omcgs_reduce_ ( ch ) : 
    """Reduce `Ostap.Math.ChannelGS`
    - see `Ostap.Math.ChannelGS`
    """
    return root_factory , ( type ( ch ) , ch.gamma0() , ch.mpi() )
   
Ostap.Math.ChannelGS. __reduce__ = _omcgs_reduce_


# =============================================================================
## Reduce Ostap::Math::BW
#  @see Ostap::Math::BW
def _ombw_reduce_ ( bw ) :
    """Reduce `Ostap.Math.BW`
    - see `Ostap.Math.BW`
    """
    return root_factory , ( type ( bw ) , bw.m0() , bw.channel() )

Ostap.Math.BW         . __reduce__ = _ombw_reduce_
Ostap.Math.BreitWigner. __reduce__ = _ombw_reduce_


# =============================================================================
## Reduce Ostap::Math::Rho0
#  @see Ostap::Math::Rho0
def _omr0_reduce_ ( bw ) :
    """Reduce `Ostap.Math.Rho0`
    - see `Ostap.Math.Rho0`
    """
    return root_factory , ( type ( bw ) , bw.m0() , bw.gamma ( 0 ) ,  bw.m1() ) 

# =============================================================================
## Reduce Ostap::Math::Kstar0
#  @see Ostap::Math::Kstar0
def _omks0_reduce_ ( bw ) :
    """Reduce `Ostap.Math.Kstar0`
    - see `Ostap.Math.Kstar0`
    """
    return root_factory , ( type ( bw ) , bw.m0() , bw.gamma ( 0 ) ,  bw.m1() , bw.m2() ) 


Ostap.Math.Rho0   . __reduce__ = _omr0_reduce_
Ostap.Math.Phi0   . __reduce__ = _omr0_reduce_
Ostap.Math.Kstar0 . __reduce__ = _omks0_reduce_

## ============================================================================
## factory for Ostap::Math::BreintWignerMC
#  @see Ostap::Math::BreintWignerMC
def _ombwmc_factory ( klass , m0 , channel , *channels )  :
    """Factory for Ostap::Math::BreintWignerMC
    - see Ostap::Math::BreintWignerMC
    """
    bwmc = klass ( m0 , channel )
    for c in channels : bwmc.addChannel ( c )
    return bwmc 

# =============================================================================
## Reduce Ostap::Math::BreitWignerMC 
#  @see Ostap::Math::BreintWignerMC 
def _ombwmc_reduce_ ( bw ) :
    """Reduce `Ostap.Math.BreitWignerMC`
    - see `Ostap.Math.BreitWignerMC`
    """
    content =  type ( bw ) , bw.m0()
    content += tuple ( bw.channel ( i ) for i in range ( bw.nChannels() ) )
    return _ombwmc_factory , content 
    
Ostap.Math.BreitWignerMC. __reduce__ = _ombwmc_reduce_

# =============================================================================
## Reduce Ostap::Math::Flatte
#  @see Ostap::Math::Flatte
def _omflt_reduce_ ( bw ) :
    """Reduce `Ostap.Math.Flatte`
    - see `Ostap.Math.Flatte`
    """
    return root_factory , ( type ( bw ) , bw.m0()    ,
                            bw.m0g1()   , bw.g2og1() ,
                            bw.mA1 ()   ,  bw.mA2()  , 
                            bw.mB1 ()   ,  bw.mB2()  , bw.g0 () ) 
                            
Ostap.Math.Flatte. __reduce__ = _omflt_reduce_


# =============================================================================
## Reduce Ostap::Math::FlatteBugg
#  @see Ostap::Math::FlatteBugg
def _omfltb_reduce_ ( bw ) :
    """Reduce `Ostap.Math.FlatteBugg`
    - see `Ostap.Math.FlatteBugg`
    """
    return root_factory , ( type ( bw ) ,
                            bw.m0      () ,
                            bw.g1      () ,
                            bw.g2og1   () ,
                            bw.alpha   () ,
                            bw.mpiplus () ,
                            bw.mpizero () ,
                            bw.mKplus  () ,
                            bw.mKzero  () ,
                            bw.g0      () )

Ostap.Math.FlatteBugg. __reduce__ = _omfltb_reduce_


# =============================================================================
## Reduce Ostap::Math::LASS
#  @see Ostap::Math::LASS
def _omlass_reduce_ ( bw ) :
    """Reduce `Ostap.Math.LASS`
    - see `Ostap.Math.LASS`
    """
    return root_factory , ( type ( bw ) , bw.m0() , bw.gamma() ,
                            bw.m1() , bw.m2() , bw.m3() ,
                            bw.a () , bw.b () , bw.e () ) 
                            
Ostap.Math.LASS. __reduce__ = _omlass_reduce_

# =============================================================================
## Reduce Ostap::Math::BWPS
#  @see Ostap::Math::BWPS
def _ombwps_reduce_ ( bw ) :
    """Reduce `Ostap.Math.BWPS`
    - see `Ostap.Math.BWPS`
    """
    return root_factory , ( type ( bw ) ,
                            bw.breit_wigner() ,
                            bw.phase_space () ,
                            bw.use_rho     () ,
                            bw.use_N2      () )
                            
Ostap.Math.BWPS. __reduce__ = _ombwps_reduce_

# =============================================================================
## Reduce Ostap::Math::BW3L
#  @see Ostap::Math::BW3L
def _ombw3l_reduce_ ( bw ) :
    """Reduce `Ostap.Math.BW3L`
    - see `Ostap.Math.BW3L`
    """
    return root_factory , ( type ( bw ) ,
                            bw.breit_wigner() ,
                            bw.M() , bw.m1() , bw.m2() , bw.m3() , bw.L() )
                            
Ostap.Math.BW3L. __reduce__ = _ombw3l_reduce_

# =============================================================================
## Reduce Ostap::Math::A2
#  @see Ostap::Math::A2
def _oma2_reduce_ ( bw ) :
    """Reduce `Ostap.Math.A2`
    - see `Ostap.Math.A2`
    """
    return root_factory , ( type ( bw ) , bw.bw() , bw.scale() )
                            
Ostap.Math.A2. __reduce__ = _oma2_reduce_

# =============================================================================
## Other peaks
# =============================================================================

# =============================================================================
## Reduce Ostap::Math::BifurcatedGauss 
#  @see Ostap::Math::BifurcatedGauss
def _ombfg_reduce_ ( peak ) :
    """Reduce `Ostap.Math.BifurcatedGauss`
    - see `Ostap.Math.BifurcatedGauss`
    """
    return root_factory , ( type ( peak ) , peak.m0 () , peak.sigmaL() , peak.sigmaR() )

Ostap.Math.BifurcatedGauss. __reduce__ = _ombfg_reduce_

# =============================================================================
## Reduce Ostap::Math::DoubleGauss 
#  @see Ostap::Math::DoubleGauss
def _om2g_reduce_ ( peak ) :
    """Reduce `Ostap.Math.DoubleGauss`
    - see `Ostap.Math.GoubleGauss`
    """
    return root_factory , ( type ( peak ) ,
                            peak.m0       () , peak.sigma() ,
                            peak.fraction () , peak.scale() ) 

Ostap.Math.DoubleGauss. __reduce__ = _om2g_reduce_

# =============================================================================
## Reduce Ostap::Math::Gauss 
#  @see Ostap::Math::Gauss
def _omg_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Gauss`
    - see `Ostap.Math.Gauss`
    """
    return root_factory , ( type ( peak ) , peak.m0 () , peak.sigma() )

Ostap.Math.Gauss. __reduce__ = _omg_reduce_

# =============================================================================
## Reduce Ostap::Math::GenGaussV1
#  @see Ostap::Math::GenGaussV1
def _omggv1_reduce_ ( peak ) :
    """Reduce `Ostap.Math.GenGaussV1`
    - see `Ostap.Math.GenGaussV1`
    """
    return root_factory , ( type ( peak ) , peak.mu() , peak.alpha() , peak.beta() )

Ostap.Math.GenGaussV1. __reduce__ = _omggv1_reduce_

# =============================================================================
## Reduce Ostap::Math::GenGaussV2
#  @see Ostap::Math::GenGaussV2
def _omggv2_reduce_ ( peak ) :
    """Reduce `Ostap.Math.GenGaussV2`
    - see `Ostap.Math.GenGaussV2`
    """
    return root_factory , ( type ( peak ) , peak.xi() , peak.alpha() , peak.kappa() )

Ostap.Math.GenGaussV2. __reduce__ = _omggv2_reduce_

# =============================================================================
## Reduce Ostap::Math::SkewGauss
#  @see Ostap::Math::SkewGauss
def _omskg_reduce_ ( peak ) :
    """Reduce `Ostap.Math.SkewGauss`
    - see `Ostap.Math.SkewGauss`
    """
    return root_factory , ( type ( peak ) , peak.xi() , peak.omega() , peak.alpha() )

Ostap.Math.SkewGauss. __reduce__ = _omskg_reduce_

# =============================================================================
## Reduce Ostap::Math::ExGauss & Ostap::Math::ExGauss2
#  @see Ostap::Math::ExGauss
#  @see Ostap::Math::ExGauss2
def _omexg_reduce_ ( peak ) :
    """Reduce `Ostap.Math.ExGauss` & `Ostap.Math.ExGauss2` 
    - see `Ostap.Math.ExGauss`
    - see `Ostap.Math.ExGauss2`
    """
    return root_factory , ( type ( peak ) , peak.mu() , peak.varsigma() , peak.k () )

Ostap.Math.ExGauss . __reduce__ = _omexg_reduce_
Ostap.Math.ExGauss2. __reduce__ = _omexg_reduce_

# =============================================================================
## Reduce Ostap::Math::Bukin2
#  @see Ostap::Math::Bukin2
def _ombk2_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Bukin2`
    - see `Ostap.Math.Bukin2`
    """
    return root_factory , ( type ( peak )     ,
                            peak.mu        () ,
                            peak.varsigmaA () ,
                            peak.varsigmaB () ,
                            peak.kA        () ,
                            peak.kB        () ,
                            peak.phi       () ) 

Ostap.Math.Bukin2 . __reduce__ = _ombk2_reduce_

# =============================================================================
## Reduce Ostap::Math::NormalLaplace
#  @see Ostap::Math::NormalLaplace
def _omnl_reduce_ ( peak ) :
    """Reduce `Ostap.Math.NormalLaplace`
    - see `Ostap.Math.NormalLaplace`
    """
    return root_factory , ( type ( peak ) , peak.mu() , peak.varsigma() , peak.kL() , peak.kR () )

Ostap.Math.NormalLaplace. __reduce__ = _omnl_reduce_

# =============================================================================
## Reduce Ostap::Math::Bukin
#  @see Ostap::Math::Nukin
def _ombu_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Bukin`
    - see `Ostap.Math.Bukin`
    """
    return root_factory , ( type ( peak ) , peak.m0() , peak.sigma() ,
                            peak.xi() , peak.rho_L() , peak.rho_R () )

Ostap.Math.Bukin. __reduce__ = _ombu_reduce_

# =============================================================================
## Reduce Ostap::Math::Novosibirsk
#  @see Ostap::Math::BNovosibirsk
def _omnovo_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Novosibirsk`
    - see `Ostap.Math.Novosibirsk`
    """
    return root_factory , ( type ( peak ) , peak.m0() , peak.sigma() , peak.tau() ) 

Ostap.Math.Novosibirsk. __reduce__ = _omnovo_reduce_

# =============================================================================
## Reduce Ostap::Math::CrystalBall
#  @see Ostap::Math::CristalBall
def _omcb_reduce_ ( peak ) :
    """Reduce `Ostap.Math.CrystalBall`
    - see `Ostap.Math.CrystalBall`
    """
    return root_factory , ( type ( peak ) , peak.m0() , peak.sigma() , peak.alpha() , peak.n() )

Ostap.Math.CrystalBall         . __reduce__ = _omcb_reduce_
Ostap.Math.CrystalBallRightSide. __reduce__ = _omcb_reduce_

# =============================================================================
## Reduce Ostap::Math::Needham
#  @see Ostap::Math::Needham
def _ommatt_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Needham`
    - see `Ostap.Math.Needham`
    """
    return root_factory , ( type ( peak ) , peak.m0 () , peak.sigma() ,
                            peak.a0 ()    , peak.a1 () , peak.a2   () ) 

Ostap.Math.Needham  . __reduce__ = _ommatt_reduce_

# =============================================================================
## Reduce Ostap::Math::CrystalBallDoubleSided
#  @see Ostap::Math::CristalBallDoubleSided 
def _omcb2_reduce_ ( peak ) :
    """Reduce `Ostap.Math.CrystalBallDoubleSided`
    - see `Ostap.Math.CrystalBallDoubleSided`
    """
    return root_factory , ( type ( peak ) , peak.m0() , peak.sigma() ,
                            peak.alpha_L() , peak.n_L() ,
                            peak.alpha_R() , peak.n_R() )
                            
Ostap.Math.CrystalBallDoubleSided. __reduce__ = _omcb2_reduce_

# =============================================================================
## Reduce Ostap::Math::Apollonios
#  @see Ostap::Math::Apollonios
def _omapo_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Apollonios`
    - see `Ostap.Math.Apollonios`
    """
    return root_factory , ( type ( peak ) , peak.m0() , peak.sigma() ,
                            peak.alpha () , peak.n () , peak.b    () )
                      
Ostap.Math.Apollonios. __reduce__ = _omapo_reduce_

# =============================================================================
## Reduce Ostap::Math::Apollonios2
#  @see Ostap::Math::Apollonios2
def _omapo2_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Apollonios2`
    - see `Ostap.Math.Apollonios2`
    """
    return root_factory , ( type ( peak ) , peak.m0()  ,
                            peak.sigmaL() , peak.sigmaR() , peak.beta () )
                      
Ostap.Math.Apollonios2. __reduce__ = _omapo2_reduce_

# =============================================================================
## Reduce Ostap::Math::StudentT
#  @see Ostap::Math::StudentT
def _omstt_reduce_ ( peak ) :
    """Reduce `Ostap.Math.StudentT`
    - see `Ostap.Math.StudentT`
    """
    return root_factory , ( type ( peak ) , peak.m0() , peak.sigma () , peak.n () )
                      
Ostap.Math.StudentT. __reduce__ = _omstt_reduce_

# =============================================================================
## Reduce Ostap::Math::BifurcatedStudentT
#  @see Ostap::Math::BifurcatedStudentT
def _ombfstt_reduce_ ( peak ) :
    """Reduce `Ostap.Math.BifurcatedStudentT`
    - see `Ostap.Math.BifurcatedStudentT`
    """
    return root_factory , ( type ( peak )  , peak.m0    () ,
                            peak.sigmaL () ,peak.sigmaR () ,
                            peak.nL     () , peak.nR    () )

Ostap.Math.BifurcatedStudentT. __reduce__ = _ombfstt_reduce_

# =============================================================================
## Reduce Ostap::Math::PearsonIV
#  @see Ostap::Math::PearsonIV
def _omp4_reduce_ ( peak ) :
    """Reduce `Ostap.Math.PearsonIV`
    - see `Ostap.Math.PearsonIV`
    """
    return root_factory , ( type ( peak )    , peak.mu () ,
                            peak.varsigma () , peak.n  () , peak.kappa () )

Ostap.Math.PearsonIV. __reduce__ = _omp4_reduce_


# =============================================================================
## reduce Ostap::Math::SkewGenT
#  @see Ostap::Math::SkewGenT 
def _omsgt_reduce_ ( peak ) :
    """Reduce `Ostap.Math.SkewGenT`
    - see `Ostap.Math.SkewGenT`
    """
    return root_factory , ( type ( peak )   , peak.mu   () ,
                            peak.sigma  ()  , peak.xi   () ,
                            peak.r      ()  , peak.zeta () )

Ostap.Math.SkewGenT. __reduce__ = _omsgt_reduce_

# =============================================================================
## reduce Ostap::Math::SkewGenErorr
#  @see Ostap::Math::SkewGenError 
def _omsge_reduce_ ( peak ) :
    """Reduce `Ostap.Math.SkewGenError`
    - see `Ostap.Math.SkewGenErorr`
    """
    return root_factory , ( type ( peak )   , peak.mu   () ,
                            peak. sigma ()  , peak.xi   () ,
                            peak.p      ()  )

Ostap.Math.SkewGenError. __reduce__ = _omsge_reduce_

# =============================================================================
## Reduce Ostap::Math::SinhAsinh
#  @see Ostap::Math::SinhAsinh
def _omshash_reduce_ ( peak ) :
    """Reduce `Ostap.Math.SinhAsinh`
    - see `Ostap.Math.SinhAsinh`
    """
    return root_factory , ( type ( peak ) , peak.location () ,
                            peak.scale () , peak.epsilon () , peak.delta () )

Ostap.Math.SinhAsinh. __reduce__ = _omshash_reduce_

# =============================================================================
## Reduce Ostap::Math::JohnsonSU
#  @see Ostap::Math::JohnsonSU
def _omjsu_reduce_ ( peak ) :
    """Reduce `Ostap.Math.JohnsonSU`
    - see `Ostap.Math.JohnsonSU`
    """
    return root_factory , ( type ( peak )    , peak.xi  () ,
                            peak.lambd () , peak.delta  () , peak.gamma () )

Ostap.Math.JohnsonSU. __reduce__ = _omjsu_reduce_

# =============================================================================
## Reduce Ostap::Math::Atlas
#  @see Ostap::Math::Atlas
def _omatlas_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Atlas`
    - see `Ostap.Math.Atlas`
    """
    return root_factory , ( type ( peak )    , peak.mean () , peak.sigma () )

Ostap.Math.Atlas   . __reduce__ = _omatlas_reduce_
Ostap.Math.Sech    . __reduce__ = _omatlas_reduce_
Ostap.Math.Logistic. __reduce__ = _omatlas_reduce_

# =============================================================================
## Reduce Ostap::Math::Losev
#  @see Ostap::Math::Losev
def _omlosev_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Losev`
    - see `Ostap.Math.Losev`
    """
    return root_factory , ( type ( peak ) , peak.mu () , peak.alpha () , peak.beta () )

Ostap.Math.Losev   . __reduce__ = _omlosev_reduce_


# =============================================================================
## Reduce Ostap::Math::GenLogisticIV 
#  @see Ostap::Math::GenLogistic4
def _omgl4_reduce_ ( peak ) :
    """Reduce `Ostap.Math.GenLogisticIV`
    - see `Ostap.Math.GenLogisticIV`
    """
    return root_factory , ( type ( peak ) ,
                            peak.mu    () ,
                            peak.sigma () ,
                            peak.alpha () ,
                            peak.beta  () )

Ostap.Math.GenLogisticIV   . __reduce__ = _omgl4_reduce_


# =============================================================================
## Reduce Ostap::Math::Slash
#  @see Ostap::Math::Slash
def _omslash_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Slash`
    - see `Ostap.Math.Slash`
    """
    return root_factory , ( type ( peak )    , peak.mean () , peak.scale () )

Ostap.Math.Slash  . __reduce__ = _omslash_reduce_

# =============================================================================
## Reduce Ostap::Math::AsymmetricLaplace
#  @see Ostap::Math::AsymmetricLaplace
def _omal_reduce_ ( peak ) :
    """Reduce `Ostap.Math.AsymmetricLaplace`
    - see `Ostap.Math.AsymmetricLaplace`
    """
    return root_factory , ( type ( peak ) , peak.mu () , peak.lambdaL ()  , peak.lambdaR() )

Ostap.Math.AsymmetricLaplace. __reduce__ = _omal_reduce_

# =============================================================================
## Reduce Ostap::Math::RaisngCosine 
#  @see Ostap::Math::RaisingCosine
def _omrcos_reduce_ ( peak ) :
    """Reduce `Ostap.Math.RaisingCosine`
    - see `Ostap.Math.RaisingCosine`
    """
    return root_factory , ( type ( peak ) , peak.mu () , peak.s () )

Ostap.Math.RaisingCosine. __reduce__ = _omrcos_reduce_

# =============================================================================
## Reduce Ostap::Math::QGaussian
#  @see Ostap::Math::QGaussian
def _omqg_reduce_ ( peak ) :
    """Reduce `Ostap.Math.QGaussian`
    - see `Ostap.Math.QGaussian`
    """
    return root_factory , ( type ( peak ) , peak.mean () , peak.scale () , peak.q () )

Ostap.Math.QGaussian. __reduce__ = _omqg_reduce_

# =============================================================================
## Reduce Ostap::Math::KGaussian
#  @see Ostap::Math::KGaussian
def _omkg_reduce_ ( peak ) :
    """Reduce `Ostap.Math.KGaussian`
    - see `Ostap.Math.KGaussian`
    """
    return root_factory , ( type ( peak ) , peak.mean () , peak.scale () , peak.kappa () )

Ostap.Math.KGaussian. __reduce__ = _omkg_reduce_

# =============================================================================
## Reduce Ostap::Math::Hyperbolic
#  @see Ostap::Math::Hyperbolic
def _omhyp_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Hyperbolic`
    - see `Ostap.Math.Hyperbolic`
    """
    return root_factory , ( type ( peak ) , peak.mu   () ,
                            peak.sigma()  , peak.zeta ()  , peak.kappa () )

Ostap.Math.Hyperbolic. __reduce__ = _omhyp_reduce_

# =============================================================================
## Reduce Ostap::Math::GenHyperbolic
#  @see Ostap::Math::GenHyperbolic
def _omghyp_reduce_ ( peak ) :
    """Reduce `Ostap.Math.GenHyperbolic`
    - see `Ostap.Math.GenHyperbolic`
    """
    return root_factory , ( type ( peak ) , peak.mu   () ,
                            peak.sigma()  , peak.zeta ()  , peak.kappa () , peak.lambd () )

Ostap.Math.GenHyperbolic. __reduce__ = _omghyp_reduce_

# =============================================================================
## Reduce Ostap::Math::Das
#  @see Ostap::Math::Das
def _omdas_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Das`
    - see `Ostap.Math.Das`
    """
    return root_factory , ( type ( peak ) , peak.mu () ,
                            peak.sigma()  , peak.kL ()  , peak.kR () )

Ostap.Math.Das. __reduce__ = _omdas_reduce_

# =============================================================================
## Reduce Ostap::Math::Hat
#  @see Ostap::Math::Hat
def _omhat_reduce_ ( peak ) :
    """Reduce `Ostap.Math.Hat`
    - see `Ostap.Math.Hat`
    """
    return root_factory , ( type ( peak ) , peak.mu   () , peak.varsigma() )

Ostap.Math.Hat. __reduce__ = _omhat_reduce_
Ostap.Math.Up. __reduce__ = _omhat_reduce_

# =============================================================================
## Reduce Ostap::Math::FupN
#  @see Ostap::Math::FupN
def _omfup_reduce_ ( peak ) :
    """Reduce `Ostap.Math.FupN`
    - see `Ostap.Math.FupN`
    """
    return root_factory , ( type ( peak ) , peak.N() , peak.mu   () , peak.varsigma() )

Ostap.Math.FupN. __reduce__ = _omfup_reduce_


# =============================================================================
## Models 
# =============================================================================

# =============================================================================
## Reduce Ostap::Math::Gumbel
#  @see Ostap::Math::Gumbel
def _omgum_reduce_ ( s ) :
    """Reduce `Ostap.Math.Gumbel`
    - see `Ostap.Math.Gumbel`
    """
    return root_factory , ( type ( s ) , s.mu () , s.beta () ) 

Ostap.Math.Gumbel  . __reduce__ = _omgum_reduce_

# =============================================================================
## Reduce Ostap::Math::GramCahrlierA
#  @see Ostap::Math::GramCharlierA
def _omgca_reduce_ ( s ) :
    """Reduce `Ostap.Math.GramCharlierA`
    - see `Ostap.Math.GramCharlierA`
    """
    return root_factory , ( type ( s ) , s.m0 () , s.sigma() , s.kappa3() , s.kappa4() ) 

Ostap.Math.GramCharlierA. __reduce__ = _omgca_reduce_

# =============================================================================
## Reduce Ostap::Math::PhaseSpacePol
#  @see Ostap::Math::PhaseSpacePol
def _ompspol_reduce_ ( s ) :
    """Reduce `Ostap.Math.PhaseSpacePol`
    - see `Ostap.Math.PhaseSpacePol`
    """
    return root_factory , ( type ( s ) , s.phasespace() , s.positive() )  

Ostap.Math.PhaseSpacePol. __reduce__ = _ompspol_reduce_

# =============================================================================
## Reduce Ostap::Math::PhaseSpaceLeftExpoPol
#  @see Ostap::Math::PhaseSpacePol
def _ompslexp_reduce_ ( s ) :
    """Reduce `Ostap.Math.PhaseSpaceLeftExpoPol`
    - see `Ostap.Math.PhaseSpaceLeftExpoPol`
    """
    return root_factory , ( type ( s ) , s.phasespace() , s.positive() . s.tau () )  

Ostap.Math.PhaseSpaceLeftExpoPol. __reduce__ = _ompslexp_reduce_

# =============================================================================
## Reduce Ostap::Math::GammaDist
#  @see Ostap::Math::GammaDist
def _omgdis_reduce_ ( s ) :
    """Reduce `Ostap.Math.GammaDist`
    - see `Ostap.Math.GammaDist`
    """
    return root_factory , ( type ( s ) , s.k() , s.theta() )  

Ostap.Math.GammaDist     . __reduce__ = _omgdis_reduce_
Ostap.Math.LogGammaDist  . __reduce__ = _omgdis_reduce_
Ostap.Math.Log10GammaDist. __reduce__ = _omgdis_reduce_

# =============================================================================
## Reduce Ostap::Math::GenGammaDist
#  @see Ostap::Math::GenGammaDist
def _omggdis_reduce_ ( s ) :
    """Reduce `Ostap.Math.GenGammaDist`
    - see `Ostap.Math.GenGammaDist`
    """
    return root_factory , ( type ( s ) , s.k() , s.theta() , s.p() , s.low() )  

Ostap.Math.GenGammaDist    . __reduce__ = _omggdis_reduce_

# =============================================================================
## Reduce Ostap::Math::LogGamma
#  @see Ostap::Math::LogGamma
def _omlgam_reduce_ ( s ) :
    """Reduce `Ostap.Math.LogGamma`
    - see `Ostap.Math.LogGamma`
    """
    return root_factory , ( type ( s ) , s.nu () , s.lambd () , s.alpha () )  

Ostap.Math.LogGamma . __reduce__ = _omlgam_reduce_

# =============================================================================
## Reduce Ostap::Math::BetaPrime
#  @see Ostap::Math::BetaPrime
def _ombprim_reduce_ ( s ) :
    """Reduce `Ostap.Math.BetaPrime`
    - see `Ostap.Math.BetaPrime`
    """
    return root_factory , ( type ( s ) , s.alpha () , s.beta() , s.scale() , s.shift () )

Ostap.Math.BetaPrime. __reduce__ = _ombprim_reduce_

# =============================================================================
## Reduce Ostap::Math::GenBetaPrime
#  @see Ostap::Math::GenBetaPrime
def _omgenbprim_reduce_ ( s ) :
    """Reduce `Ostap.Math.GenBetaPrime`
    - see `Ostap.Math.GenBetaPrime`
    """
    return root_factory , ( type ( s ) , s.alpha () , s.beta() , s.p () , s.q() , s.scale() , s.shift () )

Ostap.Math.GenBetaPrime. __reduce__ = _omgenbprim_reduce_

# =============================================================================
## Reduce Ostap::Math::Landau
#  @see Ostap::Math::Landau
def _omland_reduce_ ( s ) :
    """Reduce `Ostap.Math.Landau`
    - see `Ostap.Math.Landau`
    """
    return root_factory , ( type ( s ) , s.scale() , s.shift () )

Ostap.Math.Landau. __reduce__ = _omland_reduce_

# =============================================================================
## Reduce Ostap::Math::Weibull
#  @see Ostap::Math::Weibull
def _omwb_reduce_ ( s ) :
    """Reduce `Ostap.Math.Weibull`
    - see `Ostap.Math.Weibull`
    """
    return root_factory , ( type ( s ) , s.scale() , s.shape () , s.shift () )

Ostap.Math.Weibull. __reduce__ = _omwb_reduce_

# =============================================================================
## Reduce Ostap::Math::ExpoPosition
#  @see Ostap::Math::ExpoPositive
def _omepos_reduce_ ( s ) :
    """Reduce `Ostap.Math.ExpoPositive`
    - see `Ostap.Math.ExpoPositive`
    """
    return root_factory , ( type ( s ) , s.positive () , s.tau() )

Ostap.Math.ExpoPositive. __reduce__ = _omepos_reduce_

# =============================================================================
## Reduce Ostap::Math::Sigmoid
#  @see Ostap::Math::Sigmoid
def _omsigm_reduce_ ( s ) :
    """Reduce `Ostap.Math.Sigmoid`
    - see `Ostap.Math.Sigmoid`
    """
    return root_factory , ( type ( s ) , s.positive () , s.alpha() , s.x0() )

Ostap.Math.Sigmoid. __reduce__ = _omsigm_reduce_

# =============================================================================
## Reduce Ostap::Math::TwoExpos
#  @see Ostap::Math::TwoExpos
def _om2exp_reduce_ ( s ) :
    """Reduce `Ostap.Math.TwoExpos`
    - see `Ostap.Math.TwoExpos`
    """
    return root_factory , ( type ( s ) , s.alpha() , s.delta() , s.x0 ()  )

Ostap.Math.TwoExpos. __reduce__ = _om2exp_reduce_

# =============================================================================
## Reduce Ostap::Math::TwoExpos
#  @see Ostap::Math::TwoExpos
def _om2exp_reduce_ ( s ) :
    """Reduce `Ostap.Math.TwoExpos`
    - see `Ostap.Math.TwoExpos`
    """
    return root_factory , ( type ( s ) , s.alpha() , s.delta() , s.x0 ()  )

Ostap.Math.TwoExpos. __reduce__ = _om2exp_reduce_

# =============================================================================
## Reduce Ostap::Math::TwoExpoPositive
#  @see Ostap::Math::TwoExpoPositive
def _om2exppos_reduce_ ( s ) :
    """Reduce `Ostap.Math.TwoExpoPositive`
    - see `Ostap.Math.TwoExpoPositive`
    """
    return root_factory , ( type ( s ) , s.twoexpos() , s.positive() )

Ostap.Math.TwoExpoPositive. __reduce__ = _om2exppos_reduce_

# =============================================================================
## Reduce Ostap::Math::Rice
#  @see Ostap::Math::Rice
def _omrice_reduce_ ( s ) :
    """Reduce `Ostap.Math.Rice`
    - see `Ostap.Math.Rice`
    """
    return root_factory , ( type ( s ) , s.nu() , s.varsigma() , s.shift () )

Ostap.Math.Rice. __reduce__ = _omrice_reduce_

# =============================================================================
## Reduce Ostap::Math::GenInvGauss
#  @see Ostap::Math::GenInvGauss
def _omgig_reduce_ ( s ) :
    """Reduce `Ostap.Math.GenInvGauss`
    - see `Ostap.Math.GenInvGauss`
    """
    return root_factory , ( type ( s ) , s.theta() , s.eta() , s.p() , s.shift()  )

Ostap.Math.GenInvGauss. __reduce__ = _omgig_reduce_

# =============================================================================
## Reduce Ostap::Math::Argus
#  @see Ostap::Math::Argus
def _omargus_reduce_ ( s ) :
    """Reduce `Ostap.Math.Argus`
    - see `Ostap.Math.Argus`
    """
    return root_factory , ( type ( s ) , s.mu() , s.c() , s.chi ()  )

Ostap.Math.Argus. __reduce__ = _omargus_reduce_

# =============================================================================
## Reduce Ostap::Math::GenArgus
#  @see Ostap::Math::GenArgus
def _omgargus_reduce_ ( s ) :
    """Reduce `Ostap.Math.GenArgus`
    - see `Ostap.Math.GenArgus`
    """
    return root_factory , ( type ( s ) , s.mu() , s.c() , s.chi () , s.dp () )

Ostap.Math.GenArgus. __reduce__ = _omgargus_reduce_

# =============================================================================
## Reduce Ostap::Math::Tsallis
#  @see Ostap::Math::Tsalllis
def _omts_reduce_ ( s ) :
    """Reduce `Ostap.Math.Tsallis`
    - see `Ostap.Math.Tsallis`
    """
    return root_factory , ( type ( s ) , s.mass () , s.n () , s.T () )

Ostap.Math.Tsallis. __reduce__ = _omts_reduce_

# =============================================================================
## Reduce Ostap::Math::QGSM
#  @see Ostap::Math::QGSM
def _omqgsm_reduce_ ( s ) :
    """Reduce `Ostap.Math.QGSM`
    - see `Ostap.Math.QGSM`
    """
    return root_factory , ( type ( s ) , s.mass () , s.b () )

Ostap.Math.QGSM. __reduce__ = _omqgsm_reduce_

# =============================================================================
## Reduce Ostap::Math::Hagedorn
#  @see Ostap::Math::Hagedorn
def _omhage_reduce_ ( s ) :
    """Reduce `Ostap.Math.Hagedorn`
    - see `Ostap.Math.Hagedorn`
    """
    return root_factory , ( type ( s ) , s.mass() , s.beta() )

Ostap.Math.Hagedorn. __reduce__ = _omhage_reduce_


# =============================================================================
## reduce Ostap::Math::IrwinHall & Ostap::Math::Bates 
#  @see Ostap::Math::IrwinHall
#  @see Ostap::Math::Bates
def _omih_reduce_ ( s ) :
    """Reduce `Ostap.Math.IrwinHall` & `Ostap.Math.Bates` 
    - see `Ostap.Math.IrwinHall`
    - see `Ostap.Math.Bates` 
    """
    return root_factory , ( type ( s ) , s.n ()  )

Ostap.Math.IrwinHall. __reduce__ = _omih_reduce_
Ostap.Math.Bates    . __reduce__ = _omih_reduce_

# =============================================================================
## reduce Ostap::Math::BatesShape 
#  @see Ostap::Math::BatesShape
def _ombs_reduce_ ( s ) :
    """Reduce `Ostap.Math.BatesShape` 
    - see `Ostap.Math.BatesShape` 
    """
    return root_factory , ( type ( s ) , s.mu() , s.sigma() , s.n ()  )

Ostap.Math.BatesShape . __reduce__ = _ombs_reduce_

# =============================================================================
## reduce Ostap::Math::GenPareto, Ostap::Math::ExGenPareto and Ostap::Math::GEV
#  @see Ostap::Math::GenPareto
#  @see Ostap::Math::ExGenPareto
#  @see Ostap::Math::GEV
def _omgpd_reduce_ ( s ) :
    """Reduce `Ostap.Math.GenPareto`, `Ostap.Math.ExGenPareto` and `Ostap.Math.GeV`
    - see `Ostap.Math.GenPareto` 
    - see `Ostap.Math.ExGenPareto` 
    - see `Ostap.Math.GEV` 
    """
    return root_factory , ( type ( s ) , s.mu() , s.scale() , s.shape ()  )

Ostap.Math.GenPareto   . __reduce__ = _omgpd_reduce_
Ostap.Math.ExGenPareto . __reduce__ = _omgpd_reduce_
Ostap.Math.GEV         . __reduce__ = _omgpd_reduce_


# =============================================================================
## reduce Ostap::Math::Benini
#  @see Ostap::Math::Benini
def _omben_reduce_ ( s ) :
    """Reduce `Ostap.Math.Benini`
    - see `Ostap.Math.GenPareto` 
    """
    return root_factory , ( type ( s ) ,
                            s.pars  () ,
                            s.scale () , 
                            s.shift () )

Ostap.Math.Benini  . __reduce__ = _omben_reduce_

# =============================================================================
## reduce Ostap::Math::MPERT
#  @see Ostap::Math::MPERT
def _ommpert_reduce_ ( s ) :
    """Reduce `Ostap.Math.MPERT`
    - see `Ostap.Math.MPERT
    """
    return root_factory , ( type ( s ) ,
                            s.xmin  () ,
                            s.xmax  () ,
                            s.xi    () ,
                            s.gamma () )

Ostap.Math.MPERT . __reduce__ = _ommpert_reduce_

# =============================================================================
## Reduce Ostap::Math::HORNSdini
#  @see Ostap::Math::HORNSdini
def _omdini_reduce_ ( s ) :
    """Reduce `Ostap.Math.HORNSdini`
    - see `Ostap.Math.HORNSdini`
    """
    return root_factory , ( type ( s ) , s.a() , s.delta() , s.phi ()  )

Ostap.Math.HORNSdini. __reduce__ = _omdini_reduce_
Ostap.Math.HILLdini . __reduce__ = _omdini_reduce_

# =============================================================================
## Reduce Ostap::Math::CutOffGauss
#  @see Ostap::Math::CutOffGauss
def _omcgau_reduce_ ( s ) :
    """Reduce `Ostap.Math.CutOffGauss`
    - see `Ostap.Math.CutOffGauss`
    """
    return root_factory , ( type ( s ) , s.right () , s.x0 () , s.sigma()  )

Ostap.Math.CutOffGauss. __reduce__ = _omcgau_reduce_

# =============================================================================
## Reduce Ostap::Math::CutOffStudent
#  @see Ostap::Math::CutOffStudent
def _omcstt_reduce_ ( s ) :
    """Reduce `Ostap.Math.CutOffStudent`
    - see `Ostap.Math.CutOffStudent`
    """
    return root_factory , ( type ( s ) , s.right () , s.x0 () , s.n() , s.sigma()  )

Ostap.Math.CutOffStudent. __reduce__ = _omcstt_reduce_


# =============================================================================
# Models2 
# =============================================================================

# =============================================================================
## Reduce Ostap::Math::PS2DPol
#  @see Ostap::Math::PS2DPol
def _omps2dpol_reduce_ ( s ) :
    """Reduce `Ostap.Math.PS2DPol`
    - see `Ostap.Math.PS2DPol`
    """
    return root_factory , ( type ( s ) , s.positive() , s.psx() , s.psy () )

Ostap.Math.PS2DPol. __reduce__ = _omps2dpol_reduce_

# =============================================================================
## Reduce Ostap::Math::PS2DPolSym
#  @see Ostap::Math::PS2DPolSym
def _omps2dpols_reduce_ ( s ) :
    """Reduce `Ostap.Math.PS2DPolSym`
    - see `Ostap.Math.PS2DPolSym`
    """
    return root_factory , ( type ( s ) , s.positive() , s.psx() )

Ostap.Math.PS2DPolSym. __reduce__ = _omps2dpols_reduce_

# =============================================================================
## Reduce Ostap::Math::PS2DPol2
#  @see Ostap::Math::PS2DPol2
def _omps2dpol2_reduce_ ( s ) :
    """Reduce `Ostap.Math.PS2DPol2`
    - see `Ostap.Math.PS2DPol2`
    """
    return root_factory , ( type ( s ) , s.positive() , s.psx() , s.psy () , s.mmax() )

Ostap.Math.PS2DPol2. __reduce__ = _omps2dpol2_reduce_

# =============================================================================
## Reduce Ostap::Math::PS2DPol2Sym
#  @see Ostap::Math::PS2DPol2Sym
def _omps2dpol2s_reduce_ ( s ) :
    """Reduce `Ostap.Math.PS2DPol2Sym`
    - see `Ostap.Math.PS2DPol2Sym`
    """
    return root_factory , ( type ( s ) , s.positive() , s.psx() , s.mmax () )

Ostap.Math.PS2DPol2Sym. __reduce__ = _omps2dpol2s_reduce_

# =============================================================================
## Reduce Ostap::Math::PS2DPol3
#  @see Ostap::Math::PS2DPol3
def _omps2dpol3_reduce_ ( s ) :
    """Reduce `Ostap.Math.PS2DPol3`
    - see `Ostap.Math.PS2DPol3`
    """
    return root_factory , ( type ( s ) , s.psx() , s.psy () , s.mmax() )

Ostap.Math.PS2DPol3. __reduce__ = _omps2dpol3_reduce_

# =============================================================================
## Reduce Ostap::Math::PS2DPol3Sym
#  @see Ostap::Math::PS2DPol3Sym
def _omps2dpol3s_reduce_ ( s ) :
    """Reduce `Ostap.Math.PS2DPol3Sym`
    - see `Ostap.Math.PS2DPol3Sym`
    """
    return root_factory , ( type ( s ) , s.psx() , s.mmax () )

Ostap.Math.PS2DPol3Sym. __reduce__ = _omps2dpol3s_reduce_


# =============================================================================
## Reduce Ostap::Math::ExpoPS2DPol
#  @see Ostap::Math::ExpoPS2DPol
def _omeps2dpol_reduce_ ( s ) :
    """Reduce `Ostap.Math.ExpoPS2DPol`
    - see `Ostap.Math.ExpoPS2DPol`
    """
    return root_factory , ( type ( s ) , s.positive () , s.psy () , s.tau() )

Ostap.Math.ExpoPS2DPol. __reduce__ = _omeps2dpol_reduce_

# =============================================================================
## Reduce Ostap::Math::Expo2DPol
#  @see Ostap::Math::Expo2DPol
def _ome2dpol_reduce_ ( s ) :
    """Reduce `Ostap.Math.Expo2DPol`
    - see `Ostap.Math.Expo2DPol`
    """
    return root_factory , ( type ( s ) , s.positive () , s.tauX () , s.tauY() )

Ostap.Math.Expo2DPol. __reduce__ = _ome2dpol_reduce_

# =============================================================================
## Reduce Ostap::Math::Expo2DPolSym
#  @see Ostap::Math::Expo2DPolSym
def _ome2dpols_reduce_ ( s ) :
    """Reduce `Ostap.Math.Expo2DPolSym`
    - see `Ostap.Math.Expo2DPolSym`
    """
    return root_factory , ( type ( s ) , s.positive () , s.tau () )

Ostap.Math.Expo2DPolSym. __reduce__ = _ome2dpols_reduce_

# =============================================================================
## Reduce Ostap::Math::Gauss2D
#  @see Ostap::Math::Gauss2D
def _omg2d_reduce_ ( s ) :
    """Reduce `Ostap.Math.Gauss2D`
    - see `Ostap.Math.Gauss2D`
    """
    return root_factory , ( type ( s ) ,
                            s.muX   () , s.muY    () ,
                            s.sigmaX() , s.msigmaY() , s.theta () )

Ostap.Math.Gauss2D. __reduce__ = _omg2d_reduce_

# =============================================================================
##  Models3  
# =============================================================================

# =============================================================================
## Reduce Ostap::Math::Gauss3D
#  @see Ostap::Math::Gauss3D
def _omg3d_reduce_ ( s ) :
    """Reduce `Ostap.Math.Gauss3D`
    - see `Ostap.Math.Gauss3D`
    """
    return root_factory , ( type ( s ) ,
                            s.muX   () , s.muY    () , s.muZ    () ,                            
                            s.sigmaX() , s.msigmaY() , s.msigmaZ() ,
                            s.phi   () , s.theta  () , s.psi    () )

Ostap.Math.Gauss3D. __reduce__ = _omg3d_reduce_



# =============================================================================
## Reduce Ostap::Math::Workspace
#  @see Ostap::Math::Workspace
def _omws_reduce_ ( s ) :
    """Reduce `Ostap.Math.WorkSpace`
    - see `Ostap.Math.WorkSpace`
    """
    return root_factory , ( type ( s )        ,
                            s.size         () ,
                            s.size_cquad   () ,
                            s.size_romberg () )

Ostap.Math.WorkSpace. __reduce__ = _omws_reduce_

# =============================================================================
## Reduce Ostap::Math::Intefrator
#  @see Ostap::Math::Integrator
def _omi_reduce_ ( s ) :
    """Reduce `Ostap.Math.Integrator`
    - see `Ostap.Math.Integrator`
    """
    return root_factory , ( type ( s ) , s.ws () ) 

Ostap.Math.Integrator. __reduce__ = _omi_reduce_


# ============================================================================
## reduce 2D polynomial objects
#  @see Ostap::Math::Bernstein2D
#  @see Ostap::Math::Positive2D
#  @see Ostap::Math::LegendreSum2
def _b2d_reduce_ ( p ) :
    """reduce polynomial object
    - see  `Ostap.Math.Bernstein2D` 
    - see  `Ostap.Math.Positive2D` 
    - see  `Ostap.Math.LegendreSum2` 
    """
    return poly_factory, ( type ( p ) ,
                           array.array ( 'd' ,  p.pars() ) ,
                           p.nX   (),
                           p.nY   () ,  
                           p.xmin () ,
                           p.xmax () ,
                           p.ymin () ,
                           p.ymax () )

Ostap.Math.Bernstein2D  .__reduce__ = _b2d_reduce_
Ostap.Math.Positive2D   .__reduce__ = _b2d_reduce_
Ostap.Math.LegendreSum2 .__reduce__ = _b2d_reduce_


# ============================================================================
## Reduce KarlinShapley polynomials 
#  @see Ostap::Math::KarlinShapley 
def _kshp_reduce_ ( p ) :
    """Reduce KarlinShapley polymonials
    - see `Ostap.Math.KarlinShapley` 
    """
    pars = tuple ( [ p.A() ] + [ v for v in  p.phases1() ] + [ v for v in  p.phases1() ] ) 
    return poly_factory , ( type ( p ) ,
                            array.array ( 'd' ,  pars ) ,
                            p.xmin () ,
                            p.xmax () )

# ============================================================================
## Reduce KarlinStudden polynonial
#  @see Ostap::Math::KarlinStudden
def _kssp_reduce_ ( p ) :
    """Reduce KarlinStudden polymonials
    - see `Ostap.Math.KarlinStudden`
    """
    pars = tuple ( [ p.A() ] + [ v for v in  p.phases1() ] + [ v for v in  p.phases1() ] ) 
    return poly_factory , ( type ( p ) ,
                            array.array ( 'd' ,  pars ) ,
                            p.xmin  () ,
                            p.scale () )

Ostap.Math.KarlinShapley. __reduce__ = _kshp_reduce_
Ostap.Math.KarlinStudden. __reduce__ = _kssp_reduce_

# ============================================================================
## reduce symmetric polynomial objects 
#  @see Ostap::Math::Bernstein2DSym
#  @see Ostap::Math::Positive2DSym
#  @see Ostap::Math::Bernstein3DSym
#  @see Ostap::Math::Positive23Sym
def _b2ds_reduce_ ( p ) :
    """reduce `Ostap.Math.Bernstein2DSym` & `Ostap.Math.Positive2DSym` object
    - see  `Ostap.Math.Bernstein2DSym` 
    - see  `Ostap.Math.Positive2DSym` 
    - see  `Ostap.Math.Bernstein3DSym` 
    - see  `Ostap.Math.Positive3DSym` 
    """
    return poly_factory, ( type ( p ) ,
                           array.array ( 'd' ,  p.pars() ) ,
                           p.nX   () ,
                           p.xmin () ,
                           p.xmax () )

Ostap.Math.Bernstein2DSym.__reduce__ = _b2ds_reduce_
Ostap.Math.Positive2DSym .__reduce__ = _b2ds_reduce_
Ostap.Math.Bernstein3DSym.__reduce__ = _b2ds_reduce_
Ostap.Math.Positive3DSym .__reduce__ = _b2ds_reduce_

# ============================================================================
## reduce 3D polymonial obhjects objects
#  @see Ostap::Math::Bernstein3D
#  @see Ostap::Math::Positive3D
#  @see Ostap::Math::LegendreSum3
def _b3d_reduce_ ( p ) :
    """reduce 3D polynomial objects object
    - see  `Ostap.Math.Bernstein3D` 
    - see  `Ostap.Math.Positive3D` 
    - see  `Ostap.Math.LegendreSum3` 
    """
    return poly_factory, ( type ( p ) ,
                           array.array ( 'd' ,  p.pars() ) ,
                           p.nX   (),
                           p.nY   () ,  
                           p.nZ   () ,  
                           p.xmin () ,
                           p.xmax () ,
                           p.ymin () ,
                           p.ymax () ,
                           p.zmin () ,
                           p.zmax () )

Ostap.Math.Bernstein3D. __reduce__ = _b3d_reduce_
Ostap.Math.Positive3D   .__reduce__ = _b3d_reduce_
Ostap.Math.LegendreSum3 .__reduce__ = _b3d_reduce_


# ============================================================================
## reduce mized symmetry objects
#  @see Ostap::Math::Bernstein3DMix
#  @see Ostap::Math::Positive23Mix
def _b3dm_reduce_ ( p ) :
    """reduce `Ostap.Math.Bernstein2DMix` & `Ostap.Math.Positive2DMix` object
    - see  `Ostap.Math.Bernstein3DMix` 
    - see  `Ostap.Math.Positive3DMix` 
    """
    return poly_factory, ( type ( p ) ,
                           array.array ( 'd' ,  p.pars() ) ,
                           p.nX   () ,
                           p.nZ   () ,                           
                           p.xmin () ,
                           p.xmax () ,
                           p.zmin () ,
                           p.zmax () )

Ostap.Math.Bernstein3DMix.__reduce__ = _b3dm_reduce_
Ostap.Math.Positive3DMix .__reduce__ = _b3dm_reduce_

# ============================================================================
## reduce 4D polynomial objects
#  @see Ostap::Math::LegendreSum4
def _b4d_reduce_ ( p ) :
    """reduce 4D polynomial object
    - see  `Ostap.Math.LegendreSum4` 
    """
    return poly_factory, ( type ( p ) ,
                           array.array ( 'd' ,  p.pars() ) ,
                           p.nX   (),
                           p.nY   () ,  
                           p.nZ   () ,  
                           p.nU   () ,  
                           p.xmin () ,
                           p.xmax () ,
                           p.ymin () ,
                           p.ymax () ,
                           p.zmin () ,
                           p.zmax () ,
                           p.umin () ,
                           p.umax () )

Ostap.Math.LegendreSum4 .__reduce__ = _b4d_reduce_


# =============================================================================

_new_methods_ = [] 

_decorated_classes_  = (
    ## Formfactors 
    Ostap.Math.FormFactors.Jackson          ,
    Ostap.Math.FormFactors.BlattWeisskopf   ,
    Ostap.Math.FormFactors.NoFormFactor     ,
    ## Channels 
    Ostap.Math.ChannelCW                    , 
    Ostap.Math.ChannelQ                     ,
    Ostap.Math.ChannelFlatte                ,
    Ostap.Math.Channel                      ,
    Ostap.Math.Channel0                     ,
    Ostap.Math.Channel23L                   ,
    Ostap.Math.ChannelNR3                   ,
    Ostap.Math.ChannelGS                    ,
    ## Breit-Wigners 
    Ostap.Math.BW                           ,
    Ostap.Math.BreitWigner                  ,
    Ostap.Math.Rho0                         , 
    Ostap.Math.Phi0                         ,
    Ostap.Math.Kstar0                       ,
    Ostap.Math.BreitWignerMC                ,
    Ostap.Math.Flatte                       ,
    Ostap.Math.LASS                         ,
    Ostap.Math.BWPS                         ,
    Ostap.Math.BW3L                         ,
    Ostap.Math.A2                           ,
    ## peaks
    Ostap.Math.BifurcatedGauss              , 
    Ostap.Math.DoubleGauss                  , 
    Ostap.Math.Gauss                        , 
    Ostap.Math.GenGaussV1                   ,
    Ostap.Math.GenGaussV2                   , 
    Ostap.Math.SkewGauss                    , 
    Ostap.Math.ExGauss                      , 
    Ostap.Math.NormalLaplace                , 
    Ostap.Math.Bukin                        , 
    Ostap.Math.Novosibirsk                  , 
    Ostap.Math.CrystalBall                  , 
    Ostap.Math.CrystalBallRightSide         , 
    Ostap.Math.Needham                      , 
    Ostap.Math.CrystalBallDoubleSided       , 
    Ostap.Math.Apollonios                   , 
    Ostap.Math.Apollonios2                  , 
    Ostap.Math.StudentT                     , 
    Ostap.Math.BifurcatedStudentT           , 
    Ostap.Math.PearsonIV                    , 
    Ostap.Math.SkewGenT                     , 
    Ostap.Math.SinhAsinh                    , 
    Ostap.Math.JohnsonSU                    , 
    Ostap.Math.Atlas                        , 
    Ostap.Math.Sech                         , 
    Ostap.Math.Logistic                     , 
    Ostap.Math.Losev                        , 
    Ostap.Math.Slash                        ,     
    Ostap.Math.AsymmetricLaplace            , 
    Ostap.Math.RaisingCosine                , 
    Ostap.Math.QGaussian                    , 
    Ostap.Math.KGaussian                    , 
    Ostap.Math.Hyperbolic                   , 
    Ostap.Math.GenHyperbolic                , 
    Ostap.Math.Das                          , 
    Ostap.Math.Hat                          , 
    Ostap.Math.Up                           , 
    Ostap.Math.FupN                         ,
    ## models
    Ostap.Math.Gumbel                       , 
    Ostap.Math.GramCharlierA                , 
    Ostap.Math.PhaseSpacePol                , 
    Ostap.Math.PhaseSpaceLeftExpoPol        , 
    Ostap.Math.GammaDist                    , 
    Ostap.Math.LogGammaDist                 , 
    Ostap.Math.Log10GammaDist               , 
    Ostap.Math.GenGammaDist                 , 
    Ostap.Math.LogGamma                     , 
    Ostap.Math.BetaPrime                    , 
    Ostap.Math.Landau                       , 
    Ostap.Math.Weibull                      , 
    Ostap.Math.ExpoPositive                 , 
    Ostap.Math.Sigmoid                      , 
    Ostap.Math.TwoExpos                     , 
    Ostap.Math.TwoExpoPositive              , 
    Ostap.Math.Rice                         , 
    Ostap.Math.GenInvGauss                  , 
    Ostap.Math.Argus                        , 
    Ostap.Math.GenArgus                     , 
    Ostap.Math.Tsallis                      , 
    Ostap.Math.QGSM                         , 
    Ostap.Math.HORNSdini                    , 
    Ostap.Math.HILLdini                     , 
    Ostap.Math.CutOffGauss                  , 
    Ostap.Math.CutOffStudent                ,
    ## 2D-models
    Ostap.Math.PS2DPol                      , 
    Ostap.Math.PS2DPolSym                   , 
    Ostap.Math.PS2DPol2                     , 
    Ostap.Math.PS2DPol2Sym                  , 
    Ostap.Math.PS2DPol3                     , 
    Ostap.Math.PS2DPol3Sym                  , 
    Ostap.Math.ExpoPS2DPol                  , 
    Ostap.Math.Expo2DPol                    , 
    Ostap.Math.Expo2DPolSym                 , 
    Ostap.Math.Gauss2D                      ,
    ## 3D-models 
    Ostap.Math.Gauss3D                      ,
    ## others 
    Ostap.Math.WorkSpace                    ,
    ## 2D3D&4D polynomials
    Ostap.Math.Bernstein2D                  ,
    Ostap.Math.Bernstein2DSym               ,
    Ostap.Math.Positive2D                   ,
    Ostap.Math.Positive2DSym                ,
    Ostap.Math.Bernstein3D                  ,
    Ostap.Math.Bernstein3DSym               ,
    Ostap.Math.Bernstein3DMix               ,
    Ostap.Math.Positive3D                   ,
    Ostap.Math.Positive3DSym                ,
    Ostap.Math.Positive3DMix                ,
    Ostap.Math.LegendreSum2                 ,
    Ostap.Math.LegendreSum3                 ,
    Ostap.Math.LegendreSum4                 ,
    )

for t in _decorated_classes_ :
    _new_methods_.append ( t.__reduce__  )

_new_methods_ = tuple ( _new_methods_ ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
##                                                                      The END 
# =============================================================================

