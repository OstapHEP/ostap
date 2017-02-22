// $Id$
// ============================================================================
// Include files 
// ============================================================================
// STD & ST:
// ============================================================================
#include <limits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PDFs.h"
#include "Ostap/Iterator.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "RooArgSet.h"
#include "RooRealVar.h"
// ============================================================================
// Local 
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for namespace Ostap::Models
 *
 *
 *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
 *  @date   2011-11-30
 */
// ============================================================================
ClassImp(Ostap::Models::BreitWigner) ;
ClassImp(Ostap::Models::Rho0) ;
ClassImp(Ostap::Models::Kstar) ;
ClassImp(Ostap::Models::Phi) ;
ClassImp(Ostap::Models::BW23L) ;
ClassImp(Ostap::Models::Flatte) ;
ClassImp(Ostap::Models::Flatte2) ;
ClassImp(Ostap::Models::LASS) ;
ClassImp(Ostap::Models::LASS23L) ;
ClassImp(Ostap::Models::Bugg) ;
ClassImp(Ostap::Models::Bugg23L) ;
ClassImp(Ostap::Models::Voigt) ;
ClassImp(Ostap::Models::PseudoVoigt) ;
ClassImp(Ostap::Models::Swanson) ;
ClassImp(Ostap::Models::CrystalBall) ;
ClassImp(Ostap::Models::CrystalBallRS) ;
ClassImp(Ostap::Models::CrystalBallDS) ;
ClassImp(Ostap::Models::Needham) ;
ClassImp(Ostap::Models::Apolonios) ;
ClassImp(Ostap::Models::Apolonios2) ;
ClassImp(Ostap::Models::BifurcatedGauss) ;
ClassImp(Ostap::Models::GenGaussV1) ;
ClassImp(Ostap::Models::GenGaussV2) ;
// ClassImp(Ostap::Models::SkewGauss) ;
ClassImp(Ostap::Models::Bukin) ;
ClassImp(Ostap::Models::StudentT) ;
ClassImp(Ostap::Models::BifurcatedStudentT) ;
ClassImp(Ostap::Models::GramCharlierA) ;
ClassImp(Ostap::Models::PhaseSpace2) ;
ClassImp(Ostap::Models::PhaseSpaceLeft) ;
ClassImp(Ostap::Models::PhaseSpaceRight) ;
ClassImp(Ostap::Models::PhaseSpaceNL) ;
ClassImp(Ostap::Models::PhaseSpacePol) ;
ClassImp(Ostap::Models::PhaseSpace23L) ;
ClassImp(Ostap::Models::PolyPositive) ;
ClassImp(Ostap::Models::PolyPositiveEven) ;
ClassImp(Ostap::Models::PolyMonothonic) ;
ClassImp(Ostap::Models::PolyConvex) ;
ClassImp(Ostap::Models::PolyConvexOnly) ;
ClassImp(Ostap::Models::ExpoPositive) ;
ClassImp(Ostap::Models::PolySigmoid) ;
ClassImp(Ostap::Models::TwoExpoPositive) ;
ClassImp(Ostap::Models::GammaDist) ;
ClassImp(Ostap::Models::GenGammaDist) ;
ClassImp(Ostap::Models::Amoroso) ;
ClassImp(Ostap::Models::LogGammaDist) ;
ClassImp(Ostap::Models::Log10GammaDist) ;
ClassImp(Ostap::Models::LogGamma) ;
ClassImp(Ostap::Models::BetaPrime) ;
ClassImp(Ostap::Models::Landau) ;
ClassImp(Ostap::Models::SinhAsinh) ;
ClassImp(Ostap::Models::JohnsonSU) ;
ClassImp(Ostap::Models::Atlas) ;
ClassImp(Ostap::Models::Sech) ;
ClassImp(Ostap::Models::Logistic) ;
ClassImp(Ostap::Models::Argus) ;
ClassImp(Ostap::Models::Tsallis) ;
ClassImp(Ostap::Models::QGSM) ;
ClassImp(Ostap::Models::TwoExpos) ;
ClassImp(Ostap::Models::PositiveSpline) ;
ClassImp(Ostap::Models::MonothonicSpline) ;
ClassImp(Ostap::Models::ConvexOnlySpline) ;
ClassImp(Ostap::Models::ConvexSpline) ;
ClassImp(Ostap::Models::Poly2DPositive) ;
ClassImp(Ostap::Models::Poly2DSymPositive) ;
ClassImp(Ostap::Models::PS2DPol) ;
ClassImp(Ostap::Models::PS2DPolSym) ;
ClassImp(Ostap::Models::ExpoPS2DPol) ;
ClassImp(Ostap::Models::Expo2DPol) ;
ClassImp(Ostap::Models::Expo2DPolSym) ;
ClassImp(Ostap::Models::Spline2D) ;
ClassImp(Ostap::Models::Spline2DSym) ; 
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BreitWigner::BreitWigner 
( const char*          name  , 
  const char*          title ,
  RooAbsReal&          x     ,
  RooAbsReal&          mass  ,
  RooAbsReal&          width ,
  const double         m1    , 
  const double         m2    ,
  const unsigned short L     )
  : RooAbsPdf  (name ,title ) 
//
  , m_x     ( "x"  , "Observable" , this , x     ) 
  , m_mass  ( "m0" , "Peak"       , this , mass  ) 
  , m_width ( "g0" , "Width"      , this , width )
//
  , m_bw    ( 0 , 1 , m1 , m2 , L ) 
{
  setPars () ;
}
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BreitWigner::BreitWigner 
( const char*          name  , 
  const char*          title ,
  RooAbsReal&          x     ,
  RooAbsReal&          mass  ,
  RooAbsReal&          width ,
  const double         m1    , 
  const double         m2    ,
  const unsigned short L     , 
  const Ostap::Math::FormFactors::JacksonRho rho ) 
  : RooAbsPdf  ( name , title ) 
//
  , m_x     ( "x"  , "Observable" , this , x     ) 
  , m_mass  ( "m0" , "Peak"       , this , mass  ) 
  , m_width ( "g0" , "Width"      , this , width )
//
  , m_bw    ( 0 , 1 , m1 , m2 , L , rho ) 
{
  setPars() ;
}
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BreitWigner::BreitWigner 
( const char*                     name  , 
  const char*                     title ,
  RooAbsReal&                     x     ,
  RooAbsReal&                     mass  ,
  RooAbsReal&                     width ,
  const Ostap::Math::BreitWigner& bw    ) 
  : RooAbsPdf  ( name , title ) 
//
  , m_x     ( "x"  , "Observable" , this , x     ) 
  , m_mass  ( "m0" , "Peak"       , this , mass  ) 
  , m_width ( "g0" , "Width"      , this , width )
//
  , m_bw    ( bw ) 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::BreitWigner::BreitWigner 
( const Ostap::Models::BreitWigner& right , 
  const char*                          name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x     ( "x"  , this , right.m_x     ) 
  , m_mass  ( "m0" , this , right.m_mass  ) 
  , m_width ( "g0" , this , right.m_width )
//
  , m_bw    (               right.m_bw    ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::BreitWigner::~BreitWigner(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::BreitWigner*
Ostap::Models::BreitWigner::clone( const char* name ) const 
{ return new Ostap::Models::BreitWigner(*this,name) ; }
// ============================================================================
void Ostap::Models::BreitWigner::setPars () const 
{
  //
  m_bw.setM0     ( m_mass  ) ;
  m_bw.setGamma0 ( m_width ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::BreitWigner::evaluate() const 
{
  //
  setPars() ;
  //
  return m_bw ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::BreitWigner::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::BreitWigner::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  
  //
  return m_bw.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
// get the amplitude 
// ============================================================================
std::complex<double> Ostap::Models::BreitWigner::amplitude () const
{
  //
  m_bw.setM0   ( m_mass  ) ;
  m_bw.setGamma ( m_width ) ;
  //
  return m_bw.amplitude ( m_x ) ;
}

// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Rho0::Rho0 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          mass      ,
  RooAbsReal&          width     ,
  const double         pi_mass   )
  : Ostap::Models::BreitWigner 
    ( name    , 
      title   , 
      x       , 
      mass    ,
      width   ,
      pi_mass , 
      pi_mass , 
      1       , 
      Ostap::Math::FormFactors::Jackson_A5 )
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Rho0::Rho0 
( const Ostap::Models::Rho0& right , 
  const char*                   name  ) 
  : Ostap::Models::BreitWigner ( right , name ) 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Rho0::~Rho0(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Rho0*
Ostap::Models::Rho0::clone( const char* name ) const 
{ return new Ostap::Models::Rho0(*this,name) ; }



// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Kstar::Kstar 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          mass      ,
  RooAbsReal&          width     ,
  const double         k_mass    ,
  const double         pi_mass   ) 
  : Ostap::Models::BreitWigner 
    ( name    , 
      title   , 
      x       , 
      mass    ,
      width   ,
      k_mass  , 
      pi_mass , 
      1       , 
      Ostap::Math::FormFactors::Jackson_A2 )
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Kstar::Kstar 
( const Ostap::Models::Kstar& right , 
  const char*                    name  ) 
  : Ostap::Models::BreitWigner ( right , name ) 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Kstar::~Kstar(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Kstar*
Ostap::Models::Kstar::clone( const char* name ) const 
{ return new Ostap::Models::Kstar(*this,name) ; }


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Phi::Phi 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          mass      ,
  RooAbsReal&          width     ,
  const double         k_mass    )
  : Ostap::Models::BreitWigner
    ( name    , 
      title   , 
      x       , 
      mass    ,
      width   ,
      k_mass  , 
      k_mass  , 
      1       , 
      Ostap::Math::FormFactors::Jackson_A2 )
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Phi::Phi 
( const Ostap::Models::Phi& right , 
  const char*                    name  ) 
  : Ostap::Models::BreitWigner ( right , name ) 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Phi::~Phi(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Phi*
Ostap::Models::Phi::clone( const char* name ) const 
{ return new Ostap::Models::Phi(*this,name) ; }


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BW23L::BW23L
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          mass      ,
  RooAbsReal&          width     ,
  const double         m1        , 
  const double         m2        ,
  const unsigned short l         , 
  //
  const double         m3        , 
  const double         m         , 
  const double         L         ) 
  : RooAbsPdf ( name , title ) 
  , m_x     ( "x"     , "Observable" , this , x     ) 
  , m_mass  ( "mass"  , "BW/Peak"    , this , mass  ) 
  , m_width ( "wigth" , "BW/Width"   , this , width )
//
  , m_bw      ( 10 , 1 , m1 , m2 , m3 , m , l , L ) 
{
  setPars ();
}
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BW23L::BW23L
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          mass      ,
  RooAbsReal&          width     ,
  const double         m1        , 
  const double         m2        ,
  const unsigned short l                         ,
  const Ostap::Math::FormFactors::JacksonRho rho , 
  const double         m3        , 
  const double         m         , 
  const double         L         ) 
  : RooAbsPdf ( name , title ) 
  , m_x     ( "x"     , "Observable" , this , x     ) 
  , m_mass  ( "mass"  , "BW/Peak"    , this , mass  ) 
  , m_width ( "wigth" , "BW/Width"   , this , width )
//
  , m_bw      ( 10 , 1 , m1 , m2 , m3 , m , l , L , rho ) 
{
  setPars ();
}
// ============================================================================
// constructor from main parameters and "shape"
// ============================================================================
Ostap::Models::BW23L::BW23L
( const char*          name      , 
  const char*          title     , 
  RooAbsReal&          x         ,
  RooAbsReal&          mass      ,
  RooAbsReal&          width     ,
  const Ostap::Math::BW23L& bw   ) // shape 
  : RooAbsPdf ( name , title ) 
  , m_x     ( "x"     , "Observable" , this , x     ) 
  , m_mass  ( "mass"  , "BW/Peak"    , this , mass  ) 
  , m_width ( "wigth" , "BW/Width"   , this , width )
//
  , m_bw        ( bw ) 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::BW23L::BW23L
( const Ostap::Models::BW23L& right , 
  const char*                     name  ) 
  : RooAbsPdf ( right , name )
//
  , m_x      ( "x"     , this , right.m_x     ) 
  , m_mass   ( "mass"  , this , right.m_mass  ) 
  , m_width  ( "width" , this , right.m_width )
  , m_bw     ( right.m_bw )
{
  setPars () ;
}
// ============================================================================
Ostap::Models::BW23L::~BW23L(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::BW23L* 
Ostap::Models::BW23L::clone ( const char* name ) const 
{ return new Ostap::Models::BW23L( *this , name ) ; }
// ============================================================================
// get the amplitude 
// ============================================================================
std::complex<double> Ostap::Models::BW23L::amplitude () const
{
  setPars () ;
  return m_bw.amplitude ( m_x ) ;
}
// ============================================================================
void Ostap::Models::BW23L::setPars () const 
{
  //
  m_bw.setM0     ( m_mass  ) ;
  m_bw.setGamma0 ( m_width ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::BW23L::evaluate() const 
{
  //
  setPars () ;
  //
  return m_bw ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::BW23L::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::BW23L::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_bw.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================

  
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Flatte::Flatte 
( const char*                name      , 
  const char*                title     ,
  RooAbsReal&                x         ,
  RooAbsReal&                m0        ,
  RooAbsReal&                m0g1      ,
  RooAbsReal&                g2og1     ,
  const Ostap::Math::Flatte& flatte    ) 
  : RooAbsPdf ( name , title ) 
//
  , m_x      ( "x"     , "Observable" , this , x     ) 
  , m_m0     ( "m0"    , "Peak"       , this , m0    ) 
  , m_m0g1   ( "m0g1"  , "M0*G1"      , this , m0g1  )
  , m_g2og1  ( "g2og1" , "G2/G1"      , this , g2og1 )
    //
  , m_flatte ( flatte ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Flatte::Flatte 
( const Ostap::Models::Flatte& right , 
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x     ( "x"     , this , right.m_x     ) 
  , m_m0    ( "m0"    , this , right.m_m0    ) 
  , m_m0g1  ( "m0g1"  , this , right.m_m0g1  )
  , m_g2og1 ( "g2og1" , this , right.m_g2og1 )
//
  , m_flatte ( right.m_flatte ) 
{
  setPars() ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Flatte::~Flatte (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Flatte*
Ostap::Models::Flatte::clone( const char* name ) const 
{ return new Ostap::Models::Flatte(*this,name) ; }
// ============================================================================
void Ostap::Models::Flatte::setPars () const 
{
  //
  m_flatte.setM0     ( m_m0    ) ;
  m_flatte.setM0G1   ( m_m0g1  ) ;
  m_flatte.setG2oG1  ( m_g2og1 ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Flatte::evaluate() const 
{
  //
  setPars () ;
  //
  return m_flatte ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Flatte::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Flatte::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_flatte.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ===========================================================================
// get the amplitude 
// ===========================================================================
std::complex<double> Ostap::Models::Flatte::amplitude () const  
{
  //
  setPars () ;
  //
  return m_flatte.amplitude ( m_x ) ;
}


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Flatte2::Flatte2 
( const char*                name      , 
  const char*                title     ,
  RooAbsReal&                x         ,
  RooAbsReal&                m0        ,
  RooAbsReal&                m0g1      ,
  RooAbsReal&                g2og1     ,
  const Ostap::Math::Flatte& flatte    ) 
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"     , "Observable" , this , x     ) 
  , m_m0      ( "m0"    , "Peak"       , this , m0    ) 
  , m_m0g1    ( "m0g1"  , "M0*G1"      , this , m0g1  )
  , m_g2og1   ( "g2og1" , "G2/G1"      , this , g2og1 )
    //
  , m_flatte2 ( flatte ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Flatte2::Flatte2 
( const Ostap::Models::Flatte2& right , 
  const char*                      name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x       ( "x"     , this , right.m_x     ) 
  , m_m0      ( "m0"    , this , right.m_m0    ) 
  , m_m0g1    ( "m0g1"  , this , right.m_m0g1  )
  , m_g2og1   ( "g2og1" , this , right.m_g2og1 )
//
  , m_flatte2 ( right.m_flatte2 ) 
{
  setPars() ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Flatte2::~Flatte2 (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Flatte2*
Ostap::Models::Flatte2::clone( const char* name ) const 
{ return new Ostap::Models::Flatte2(*this,name) ; }
// ============================================================================
void Ostap::Models::Flatte2::setPars () const 
{
  //
  m_flatte2.setM0     ( m_m0    ) ;
  m_flatte2.setM0G1   ( m_m0g1  ) ;
  m_flatte2.setG2oG1  ( m_g2og1 ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Flatte2::evaluate() const 
{
  //
  setPars () ;
  //
  return m_flatte2 ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Flatte2::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Flatte2::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_flatte2.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ===========================================================================
// get the amplitude 
// ===========================================================================
std::complex<double> Ostap::Models::Flatte2::amplitude () const  
{
  //
  setPars () ;
  //
  return m_flatte2.amplitude ( m_x ) ;
}





// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::LASS::LASS
( const char*          name  , 
  const char*          title ,
  RooAbsReal&          x     ,
  RooAbsReal&          m1430 ,
  RooAbsReal&          g1430 ,
  RooAbsReal&          a     , 
  RooAbsReal&          r     , 
  RooAbsReal&          e     , 
  const double         m1    , 
  const double         m2    )
//
  : RooAbsPdf ( name , title ) 
//
  , m_x     ( "x"     , "Observable"      , this , x     ) 
  , m_m0    ( "m0"    , "K*(1430)-mass"   , this , m1430 ) 
  , m_g0    ( "g0"    , "K*(1430)-width"  , this , g1430 ) 
  , m_a     ( "a"     , "LASS-a"          , this , a     )
  , m_r     ( "r"     , "LASS-r"          , this , r     )
  , m_e     ( "e"     , "LASS-elasticity" , this , e     )
//
  , m_lass  ( m1      , m2      , 
              1430    , 300     , 
              1.94e-3 , 1.76e-1 , 1.0 ) 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::LASS::LASS 
( const Ostap::Models::LASS& right , 
  const char*                   name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x     ( "x"  , this , right.m_x  ) 
  , m_m0    ( "m0" , this , right.m_m0 ) 
  , m_g0    ( "g0" , this , right.m_g0 ) 
  , m_a     ( "a"  , this , right.m_a  ) 
  , m_r     ( "r"  , this , right.m_r  ) 
  , m_e     ( "e"  , this , right.m_e  ) 
//
  , m_lass  ( right.m_lass ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::LASS::~LASS (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::LASS*
Ostap::Models::LASS::clone( const char* name ) const 
{ return new Ostap::Models::LASS ( *this , name ) ; }
// ============================================================================
void Ostap::Models::LASS::setPars () const 
{
  //
  m_lass.setM0 ( m_m0 ) ;
  m_lass.setG0 ( m_g0 ) ;
  m_lass.setA  ( m_a  ) ;
  m_lass.setR  ( m_r  ) ;
  m_lass.setE  ( m_e  ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::LASS::evaluate() const 
{
  //
  setPars () ;
  //
  return m_lass ( m_x  ) ;
}
// ============================================================================
Int_t Ostap::Models::LASS::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::LASS::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  //
  return m_lass.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ===========================================================================
// get the complex amplitude 
// ===========================================================================
std::complex<double> Ostap::Models::LASS::amplitude() const 
{
  //
  setPars() ;
  //
  return m_lass.amplitude ( m_x  ) ;
}
// ============================================================================

// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::LASS23L::LASS23L
( const char*          name  , 
  const char*          title ,
  RooAbsReal&          x     ,
  RooAbsReal&          m1430 ,
  RooAbsReal&          g1430 ,
  RooAbsReal&          a     , 
  RooAbsReal&          r     , 
  RooAbsReal&          e     ,                
  const double         m1    , 
  const double         m2    ,
  const double         m3    , 
  const double         m     ,
  const unsigned short L     ) 
//
  : RooAbsPdf ( name , title ) 
//
  , m_x     ( "x"     , "Observable"      , this , x     ) 
  , m_m0    ( "m0"    , "K*(1430)-mass"   , this , m1430 ) 
  , m_g0    ( "g0"    , "K*(1430)-width"  , this , g1430 ) 
  , m_a     ( "a"     , "LASS-a"          , this , a     )
  , m_r     ( "r"     , "LASS-r"          , this , r     )
  , m_e     ( "e"     , "LASS-elasticity" , this , e     )
//
  , m_lass  ( m1      , m2      , m3   ,  m  , 
              1430    , 300     , 
              L       , 
              1.94e-3 , 1.76e-1 , 1.0    ) 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::LASS23L::LASS23L 
( const Ostap::Models::LASS23L& right , 
  const char*                      name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x     ( "x"  , this , right.m_x  ) 
  , m_m0    ( "m0" , this , right.m_m0 ) 
  , m_g0    ( "g0" , this , right.m_g0 ) 
  , m_a     ( "a"  , this , right.m_a  ) 
  , m_r     ( "r"  , this , right.m_r  ) 
  , m_e     ( "e"  , this , right.m_e  ) 
//
  , m_lass  ( right.m_lass ) 
{
  setPars   () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::LASS23L::~LASS23L (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::LASS23L*
Ostap::Models::LASS23L::clone( const char* name ) const 
{ return new Ostap::Models::LASS23L ( *this , name ) ; }
// ============================================================================
void Ostap::Models::LASS23L::setPars () const 
{
  //
  m_lass.setM0 ( m_m0 ) ;
  m_lass.setG0 ( m_g0 ) ;
  m_lass.setA  ( m_a  ) ;
  m_lass.setR  ( m_r  ) ;
  m_lass.setE  ( m_e  ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::LASS23L::evaluate() const 
{
  //
  setPars () ;
  //
  return m_lass ( m_x  ) ;
}
// ============================================================================
Int_t Ostap::Models::LASS23L::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::LASS23L::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_lass.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ===========================================================================
// get the complex amplitude 
// ===========================================================================
std::complex<double> Ostap::Models::LASS23L::amplitude() const 
{
  //
  setPars () ;
  //
  return m_lass.amplitude ( m_x  ) ;
}
// ============================================================================

// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Bugg::Bugg  
( const char*          name               , 
  const char*          title              ,
  RooAbsReal&          x                  ,
  RooAbsReal&          M                  ,   // sigma M 
  RooAbsReal&          g2                 ,   // sigma G2 
  RooAbsReal&          b1                 ,   // sigma B1 
  RooAbsReal&          b2                 ,   // sigma B2
  RooAbsReal&          a                  ,   // sigma a 
  RooAbsReal&          s1                 ,   // sigma s1 
  RooAbsReal&          s2                 ,   // sigma s2 
  const double         m1                 )   // mass of pi GeV 
  : RooAbsPdf ( name , title ) 
//
  , m_x     ( "x"     , "Observable"      , this , x     ) 
  , m_M     ( "M"     , "Bugg/M"          , this , M     ) 
  , m_g2    ( "g2"    , "Bugg/G2"         , this , g2    ) 
  , m_b1    ( "b1"    , "Bugg/b1"         , this , b1    )
  , m_b2    ( "b2"    , "Bugg/b2"         , this , b2    )
  , m_a     ( "a"     , "Bugg/a"          , this , a     )
  , m_s1    ( "s1"    , "Bugg/s1"         , this , s1    )
  , m_s2    ( "s2"    , "Bugg/s2"         , this , s2    )
//
  , m_bugg  ( 0.92 , 0.0024 , 0.5848 , 1.6663 , 1.082 , 2.8 , 3.5 , m1 )  
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Bugg::Bugg 
( const Ostap::Models::Bugg& right , 
  const char*                      name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x     ( "x"     , this , right.m_x     ) 
  , m_M     ( "M"     , this , right.m_M     ) 
  , m_g2    ( "g2"    , this , right.m_g2    ) 
  , m_b1    ( "b1"    , this , right.m_b1    )
  , m_b2    ( "b2"    , this , right.m_b2    )
  , m_a     ( "a"     , this , right.m_a     )
  , m_s1    ( "s1"    , this , right.m_s1    )
  , m_s2    ( "s2"    , this , right.m_s2    )
//
  , m_bugg  ( right.m_bugg ) 
{
  setPars () ;
}
// ============================================================================
Ostap::Models::Bugg::~Bugg(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Bugg*
Ostap::Models::Bugg::clone( const char* name ) const 
{ return new Ostap::Models::Bugg ( *this , name ) ; }
// ============================================================================
void Ostap::Models::Bugg::setPars () const 
{
  //
  m_bugg.setM  ( m_M  ) ;
  m_bugg.setG2 ( m_g2 ) ;
  m_bugg.setB1 ( m_b1 ) ;
  m_bugg.setB2 ( m_b2 ) ;
  m_bugg.setA  ( m_a  ) ;
  m_bugg.setS1 ( m_s1 ) ;
  m_bugg.setS2 ( m_s2 ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Bugg::evaluate() const 
{
  //
  setPars () ;
  //
  return m_bugg ( m_x  ) ;
}
// ============================================================================
Int_t Ostap::Models::Bugg::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Bugg::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  //
  return m_bugg.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ===========================================================================
// get the complex amplitude 
// ===========================================================================
std::complex<double> Ostap::Models::Bugg::amplitude() const 
{
  //
  setPars() ;
  //
  return m_bugg.amplitude ( m_x  ) ;
}
// ============================================================================

// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Bugg23L::Bugg23L  
( const char*          name               , 
  const char*          title              ,
  RooAbsReal&          x                  ,
  RooAbsReal&          M                  ,   // sigma M 
  RooAbsReal&          g2                 ,   // sigma G2 
  RooAbsReal&          b1                 ,   // sigma B1 
  RooAbsReal&          b2                 ,   // sigma B2
  RooAbsReal&          a                  ,   // sigma a 
  RooAbsReal&          s1                 ,   // sigma s1 
  RooAbsReal&          s2                 ,   // sigma s2 
  const double         m1                 ,   // mass of pi GeV 
  const double         m3                 ,   // mass of third particle 
  const double         m                  ,   // mass of mother  
  const unsigned short L                  )
  : RooAbsPdf ( name , title ) 
//
  , m_x     ( "x"     , "Observable"      , this , x     ) 
  , m_M     ( "M"     , "Bugg/M"          , this , M     ) 
  , m_g2    ( "g2"    , "Bugg/G2"         , this , g2    ) 
  , m_b1    ( "b1"    , "Bugg/b1"         , this , b1    )
  , m_b2    ( "b2"    , "Bugg/b2"         , this , b2    )
  , m_a     ( "a"     , "Bugg/a"          , this , a     )
  , m_s1    ( "s1"    , "Bugg/s1"         , this , s1    )
  , m_s2    ( "s2"    , "Bugg/s2"         , this , s2    )
//
  , m_bugg  ( 0.92 , 0.0024 , 0.5848 , 1.6663 , 1.082 , 2.8 , 3.5 , m1 , m3 , m , L )  
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Bugg23L::Bugg23L 
( const Ostap::Models::Bugg23L& right , 
  const char*                      name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x     ( "x"     , this , right.m_x     ) 
  , m_M     ( "M"     , this , right.m_M     ) 
  , m_g2    ( "g2"    , this , right.m_g2    ) 
  , m_b1    ( "b1"    , this , right.m_b1    )
  , m_b2    ( "b2"    , this , right.m_b2    )
  , m_a     ( "a"     , this , right.m_a     )
  , m_s1    ( "s1"    , this , right.m_s1    )
  , m_s2    ( "s2"    , this , right.m_s2    )
//
  , m_bugg  ( right.m_bugg ) 
{
  setPars() ;
}
// ============================================================================
Ostap::Models::Bugg23L::~Bugg23L(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Bugg23L*
Ostap::Models::Bugg23L::clone( const char* name ) const 
{ return new Ostap::Models::Bugg23L ( *this , name ) ; }
// ============================================================================
void Ostap::Models::Bugg23L::setPars () const 
{
  //
  m_bugg.setM  ( m_M  ) ;
  m_bugg.setG2 ( m_g2 ) ;
  m_bugg.setB1 ( m_b1 ) ;
  m_bugg.setB2 ( m_b2 ) ;
  m_bugg.setA  ( m_a  ) ;
  m_bugg.setS1 ( m_s1 ) ;
  m_bugg.setS2 ( m_s2 ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Bugg23L::evaluate() const 
{
  //
  setPars () ;
  //
  return m_bugg ( m_x  ) ;
}
// ============================================================================
Int_t Ostap::Models::Bugg23L::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Bugg23L::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_bugg.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ===========================================================================
// get the complex amplitude 
// ===========================================================================
std::complex<double> 
Ostap::Models::Bugg23L::amplitude() const 
{
  //
  setPars() ;
  //
  return m_bugg.amplitude ( m_x  ) ;
}
// ============================================================================

// ============================================================================
//         Voigt
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Voigt::Voigt 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          m0        , 
  RooAbsReal&          gamma     , 
  RooAbsReal&          sigma     )
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"       , "Observable" , this , x      ) 
  , m_m0      ( "m0"      , "m0"         , this , m0     ) 
  , m_gamma   ( "gamma"   , "gamma"      , this , gamma  )
  , m_sigma   ( "sigma"   , "sigma"      , this , sigma  )
//
  , m_voigt   ( 10 , 1 , 1 ) 
{
  //
  setPars () ;
  //
} 
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Voigt::Voigt
( const Ostap::Models::Voigt& right  , 
  const char*                    name   ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_m0     ( "m0"     , this , right.m_m0     ) 
  , m_gamma  ( "gamma"  , this , right.m_gamma  ) 
  , m_sigma  ( "sigma"  , this , right.m_sigma  ) 
//
  , m_voigt  ( right.m_voigt ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Voigt::~Voigt(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Voigt*
Ostap::Models::Voigt::clone( const char* name ) const 
{ return new Ostap::Models::Voigt(*this,name) ; }
// ============================================================================
void Ostap::Models::Voigt::setPars () const 
{
  //
  m_voigt.setM0     ( m_m0    ) ;
  m_voigt.setSigma  ( m_sigma  ) ;
  m_voigt.setGamma  ( m_gamma  ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Voigt::evaluate() const 
{
  //
  setPars() ;
  //
  return m_voigt    ( m_x     ) ;
}
// ============================================================================
Int_t Ostap::Models::Voigt::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Voigt::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  //
  return m_voigt.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================

// ============================================================================
//         PseudoVoigt
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::PseudoVoigt::PseudoVoigt 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          m0        , 
  RooAbsReal&          gamma     , 
  RooAbsReal&          sigma     )
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"       , "Observable" , this , x      ) 
  , m_m0      ( "m0"      , "m0"         , this , m0     ) 
  , m_gamma   ( "gamma"   , "gamma"      , this , gamma  )
  , m_sigma   ( "sigma"   , "sigma"      , this , sigma  )
//
  , m_voigt   ( 10 , 1 , 1 ) 
{
  //
  setPars () ;
  //
} 
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::PseudoVoigt::PseudoVoigt 
( const Ostap::Models::PseudoVoigt& right  , 
  const char*                    name   ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_m0     ( "m0"     , this , right.m_m0     ) 
  , m_gamma  ( "gamma"  , this , right.m_gamma  ) 
  , m_sigma  ( "sigma"  , this , right.m_sigma  ) 
//
  , m_voigt  ( right.m_voigt ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PseudoVoigt::~PseudoVoigt(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PseudoVoigt*
Ostap::Models::PseudoVoigt::clone( const char* name ) const 
{ return new Ostap::Models::PseudoVoigt(*this,name) ; }
// ============================================================================
void Ostap::Models::PseudoVoigt::setPars () const 
{
  //
  m_voigt.setM0     ( m_m0    ) ;
  m_voigt.setSigma  ( m_sigma  ) ;
  m_voigt.setGamma  ( m_gamma  ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PseudoVoigt::evaluate() const 
{
  //
  setPars() ;
  //
  return m_voigt    ( m_x     ) ;
}
// ============================================================================
Int_t Ostap::Models::PseudoVoigt::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PseudoVoigt::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  //
  return m_voigt.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================




// ============================================================================
// Swanson's S-wave cusp model
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Swanson::Swanson 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          beta0     , 
  const Ostap::Math::Swanson& sw ) 
  : RooAbsPdf ( name , title ) 
    //
  , m_x       ( "x"       , "Observable" , this , x      ) 
  , m_beta0   ( "beta0"   , "beta_0"     , this , beta0  ) 
    //
  , m_swanson ( sw ) 
{
  setPars () ;
} 
// ============================================================================
// Swanson's S-wave cusp model
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Swanson::Swanson 
( const char*          name          , 
  const char*          title         ,
  RooAbsReal&          x             , 
  RooAbsReal&          beta0         , 
  const double         m1_0          , 
  const double         m2_0          ,
  const Ostap::Math::BreitWigner& bw ) 
  : RooAbsPdf ( name , title ) 
    //
  , m_x       ( "x"       , "Observable" , this , x      ) 
  , m_beta0   ( "beta0"   , "beta_0"     , this , beta0  ) 
    //
  , m_swanson ( bw , m1_0 , m2_0 , 1.0 ) 
{
  setPars () ;
} 

// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Swanson::Swanson
( const Ostap::Models::Swanson& right  , 
  const char*                    name   ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x       ( "x"      , this , right.m_x      ) 
  , m_beta0   ( "beta0"  , this , right.m_beta0  ) 
//
  , m_swanson ( right.m_swanson ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Swanson::~Swanson(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Swanson*
Ostap::Models::Swanson::clone( const char* name ) const 
{ return new Ostap::Models::Swanson(*this,name) ; }
// ============================================================================
void Ostap::Models::Swanson::setPars () const 
{
  //
  m_swanson.setBeta0 ( m_beta0 ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Swanson::evaluate() const 
{
  //
  setPars() ;
  //
  return m_swanson ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Swanson::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Swanson::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  //
  return m_swanson.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// constructor form all parameters 
// ============================================================================
Ostap::Models::CrystalBall::CrystalBall
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          m0        ,
  RooAbsReal&          sigma     ,    
  RooAbsReal&          alpha     ,
  RooAbsReal&          n         ) 
  : RooAbsPdf ( name , title )
//
  , m_x       ( "x"       , "Observable"   , this , x      ) 
  , m_m0      ( "m0"      , "CB/mass"      , this , m0     ) 
  , m_sigma   ( "sigma"   , "CB/sigma"     , this , sigma  )
//
  , m_alpha   ( "alpha"   , "CB/alpha"     , this , alpha  ) 
  , m_n       ( "n"       , "CB/n"         , this , n      ) 
//
  , m_cb      ( 100 , 1 , 1 , 10 ) 
{
  //
  setPars () ;
  //
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::CrystalBall::CrystalBall
( const Ostap::Models::CrystalBall& right , 
  const char*                          name  ) 
  : RooAbsPdf ( right , name )
//
  , m_x       ( "x"       , this , right.m_x      ) 
  , m_m0      ( "m0"      , this , right.m_m0     ) 
  , m_sigma   ( "sigma"   , this , right.m_sigma  )
//
  , m_alpha   ( "alpha"   , this , right.m_alpha  ) 
  , m_n       ( "n"       , this , right.m_n      ) 
//
  , m_cb      ( 100 , 1 , 1 , 10 ) 
{
  setPars () ;
}
// ============================================================================
Ostap::Models::CrystalBall::~CrystalBall(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::CrystalBall*
Ostap::Models::CrystalBall::clone ( const char* name ) const 
{ return new Ostap::Models::CrystalBall (*this , name ) ; }
// ============================================================================
void Ostap::Models::CrystalBall::setPars () const 
{
  //
  m_cb.setM0      ( m_m0     ) ;
  m_cb.setSigma   ( m_sigma  ) ;
  m_cb.setAlpha   ( m_alpha  ) ;
  m_cb.setN       ( m_n      ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::CrystalBall::evaluate() const
{
  //
  setPars () ;
  //
  return m_cb ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::CrystalBall::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::CrystalBall::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars ();
  //
  return m_cb.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================





// ============================================================================
// constructor form all parameters 
// ============================================================================
Ostap::Models::CrystalBallRS::CrystalBallRS
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          m0        ,
  RooAbsReal&          sigma     ,    
  RooAbsReal&          alpha     ,
  RooAbsReal&          n         ) 
  : RooAbsPdf ( name , title )
//
  , m_x       ( "x"       , "Observable"   , this , x      ) 
  , m_m0      ( "m0"      , "CB/mass"      , this , m0     ) 
  , m_sigma   ( "sigma"   , "CB/sigma"     , this , sigma  )
//
  , m_alpha   ( "alpha"   , "CB/alpha"     , this , alpha  ) 
  , m_n       ( "n"       , "CB/n"         , this , n      ) 
//
  , m_cb      ( 100 , 1 , 1 , 10 ) 
{
  //
  setPars () ;
  //
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::CrystalBallRS::CrystalBallRS
( const Ostap::Models::CrystalBallRS& right , 
  const char*                          name  ) 
  : RooAbsPdf ( right , name )
//
  , m_x       ( "x"       , this , right.m_x      ) 
  , m_m0      ( "m0"      , this , right.m_m0     ) 
  , m_sigma   ( "sigma"   , this , right.m_sigma  )
//
  , m_alpha   ( "alpha"   , this , right.m_alpha  ) 
  , m_n       ( "n"       , this , right.m_n      ) 
//
  , m_cb      ( 100 , 1 , 1 , 10 ) 
{
  setPars () ;
}
// ============================================================================
Ostap::Models::CrystalBallRS::~CrystalBallRS(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::CrystalBallRS*
Ostap::Models::CrystalBallRS::clone ( const char* name ) const 
{ return new Ostap::Models::CrystalBallRS (*this , name ) ; }
// ============================================================================
void Ostap::Models::CrystalBallRS::setPars () const 
{
  //
  m_cb.setM0      ( m_m0     ) ;
  m_cb.setSigma   ( m_sigma  ) ;
  m_cb.setAlpha   ( m_alpha  ) ;
  m_cb.setN       ( m_n      ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::CrystalBallRS::evaluate() const
{
  //
  setPars() ;
  //
  return m_cb ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::CrystalBallRS::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::CrystalBallRS::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_cb.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================

// ============================================================================
// Double-sided CrystalBall
// ============================================================================
Ostap::Models::CrystalBallDS::CrystalBallDS
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          m0        ,
  RooAbsReal&          sigma     ,    
  RooAbsReal&          alphaL    ,    
  RooAbsReal&          nL        ,    
  RooAbsReal&          alphaR    ,    
  RooAbsReal&          nR        )
  : RooAbsPdf ( name , title )
//
  , m_x       ( "x"       , "Observable"                 , this , x      ) 
  , m_m0      ( "m0"      , "mass"                       , this , m0     ) 
  , m_sigma   ( "sigma"   , "sigma"                      , this , sigma  )
  , m_alphaL  ( "alphaL"  , "(left) alpha = 1 + |alpha|" , this , alphaL ) 
  , m_nL      ( "nL"      , "(left) n     = 1 + |n|"     , this ,     nL ) 
  , m_alphaR  ( "alphaR"  , "(left) alpha = 1 + |alpha|" , this , alphaR ) 
  , m_nR      ( "nR"      , "(left) n     = 1 + |n|"     , this ,     nR ) 
//  
  , m_cb2 ( 10 , 1 , 1 , 1 , 1  , 1 ) 
{
  //
  setPars () ;
  //
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Models::CrystalBallDS::CrystalBallDS
( const Ostap::Models::CrystalBallDS& right , 
  const char*                            name ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x       ( "x"       , this , right.m_x      ) 
  , m_m0      ( "m0"      , this , right.m_m0     ) 
  , m_sigma   ( "sigma"   , this , right.m_sigma  )
  , m_alphaL  ( "alphaL"  , this , right.m_alphaL ) 
  , m_nL      ( "nL"      , this , right.m_nL     ) 
  , m_alphaR  ( "alphaR"  , this , right.m_alphaR ) 
  , m_nR      ( "nR"      , this , right.m_nR     ) 
//  
  , m_cb2     ( right.m_cb2 ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::CrystalBallDS::~CrystalBallDS(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::CrystalBallDS*
Ostap::Models::CrystalBallDS::clone( const char* name ) const 
{ return new Ostap::Models::CrystalBallDS(*this,name) ; }
// ============================================================================
void Ostap::Models::CrystalBallDS::setPars () const 
{
  //
  m_cb2.setM0      ( m_m0     ) ;
  m_cb2.setSigma   ( m_sigma  ) ;
  m_cb2.setAlpha_L ( m_alphaL ) ;
  m_cb2.setAlpha_R ( m_alphaR ) ;
  m_cb2.setN_L     ( m_nL     ) ;
  m_cb2.setAlpha_R ( m_alphaR ) ;
  m_cb2.setN_R     ( m_nR     ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::CrystalBallDS::evaluate() const 
{
  //
  setPars () ;
  //
  return m_cb2     ( m_x      ) ;
}
// ============================================================================
Int_t Ostap::Models::CrystalBallDS::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::CrystalBallDS::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  //
  return m_cb2.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// Needham
// ============================================================================
Ostap::Models::Needham::Needham
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          m0        ,
  RooAbsReal&          sigma     ,    
  RooAbsReal&          a0        ,    
  RooAbsReal&          a1        ,    
  RooAbsReal&          a2        )
  : RooAbsPdf ( name , title )
//
  , m_x       ( "x"       , "Observable"                 , this , x      ) 
  , m_m0      ( "m0"      , "mass"                       , this , m0     ) 
  , m_sigma   ( "sigma"   , "sigma"                      , this , sigma  )
//
  , m_a0      ( "a0"      , "a0-parameter"               , this , a0     ) 
  , m_a1      ( "a1"      , "a1-parameter"               , this , a1     ) 
  , m_a2      ( "a2"      , "a2-parameter"               , this , a2     ) 
//
  , m_needham ( 100 , 1 , 1.9 , 0 , 0 ) 
{
  //
  setPars () ;
  //
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Models::Needham::Needham
( const Ostap::Models::Needham& right , 
  const char*                      name ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x       ( "x"       , this , right.m_x      ) 
  , m_m0      ( "m0"      , this , right.m_m0     ) 
  , m_sigma   ( "sigma"   , this , right.m_sigma  )
//
  , m_a0      ( "a0"      , this , right.m_a0     ) 
  , m_a1      ( "a1"      , this , right.m_a1     ) 
  , m_a2      ( "a2"      , this , right.m_a2     ) 
//
  , m_needham ( right.m_needham ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Needham::~Needham(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Needham*
Ostap::Models::Needham::clone( const char* name ) const 
{ return new Ostap::Models::Needham ( *this , name ) ; }
// ============================================================================
void Ostap::Models::Needham::setPars () const 
{
  //
  m_needham . setM0    ( m_m0     ) ;
  m_needham . setSigma ( m_sigma  ) ;
  //
  m_needham . setA0    ( m_a0     ) ;
  m_needham . setA1    ( m_a1     ) ;
  m_needham . setA2    ( m_a2     ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Needham::evaluate() const 
{
  //
  setPars () ;
  //
  return m_needham     ( m_x      ) ;
}
// ============================================================================
Int_t Ostap::Models::Needham::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Needham::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_needham.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
// get current alpha 
// ============================================================================
double Ostap::Models::Needham::alpha   () const 
{
  const double s  =         m_sigma  ;
  double       a  =         m_a0     ;
  a              += s *     m_a1     ;
  a              += s * s * m_a2     ;
  //
  return a ;
}
// ============================================================================



// ============================================================================
// Apolonios
// ============================================================================
Ostap::Models::Apolonios::Apolonios
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          m0        ,
  RooAbsReal&          sigma     ,    
  RooAbsReal&          alpha     ,    
  RooAbsReal&          n         ,    
  RooAbsReal&          b         ) 
  : RooAbsPdf ( name , title )
//
  , m_x       ( "x"       , "Observable"                 , this , x      ) 
  , m_m0      ( "m0"      , "mass"                       , this , m0     ) 
  , m_sigma   ( "sigma"   , "sigma"                      , this , sigma  )
  , m_alpha   ( "alpha"   , "alpha"                      , this , alpha  )
  , m_n       ( "n"       , "n-parameter"                , this , n      )
  , m_b       ( "b"       , "b-parameter"                , this , b      )
//
  , m_apo ( 1 , 1 , 1 , 1 , 1 ) 
{
  //
  setPars () ;
  //
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Models::Apolonios::Apolonios
( const Ostap::Models::Apolonios& right , 
  const char*                        name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x       ( "x"       , this , right.m_x      ) 
  , m_m0      ( "m0"      , this , right.m_m0     ) 
  , m_sigma   ( "sigma"   , this , right.m_sigma  )
  , m_alpha   ( "alpha"   , this , right.m_alpha  )
  , m_n       ( "n"       , this , right.m_n      )
  , m_b       ( "b"       , this , right.m_b      )
//
  , m_apo     ( right.m_apo ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Apolonios::~Apolonios (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Apolonios*
Ostap::Models::Apolonios::clone ( const char* name ) const 
{ return new Ostap::Models::Apolonios ( *this , name ) ; }
// ============================================================================
void Ostap::Models::Apolonios::setPars () const 
{
  //
  m_apo.setM0      ( m_m0     ) ;
  m_apo.setSigma   ( m_sigma  ) ;
  m_apo.setAlpha   ( m_alpha  ) ;
  m_apo.setN       ( m_n      ) ;
  m_apo.setB       ( m_b      ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Apolonios::evaluate() const 
{
  //
  setPars () ;
  //
  return m_apo ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Apolonios::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Apolonios::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_apo.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
// Apolonios2
// ============================================================================
Ostap::Models::Apolonios2::Apolonios2
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          m0        ,
  RooAbsReal&          sigmaL    ,    
  RooAbsReal&          sigmaR    ,    
  RooAbsReal&          beta      ) 
  : RooAbsPdf ( name , title )
//
  , m_x       ( "x"       , "Observable"                 , this , x      ) 
  , m_m0      ( "m0"      , "mass"                       , this , m0     ) 
  , m_sigmaL  ( "sigmaL"  , "sigmaL"                     , this , sigmaL )
  , m_sigmaR  ( "sigmaR"  , "sigmaR"                     , this , sigmaR )
  , m_beta    ( "beta"    , "beta-parameter"             , this , beta   )
//
  , m_apo2    ( 1 , 1 , 1 , 1 ) 
{
  //
  setPars () ;
  //
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Models::Apolonios2::Apolonios2
( const Ostap::Models::Apolonios2& right , 
  const char*                        name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x       ( "x"       , this , right.m_x      ) 
  , m_m0      ( "m0"      , this , right.m_m0     ) 
  , m_sigmaL  ( "sigmaL"  , this , right.m_sigmaL )
  , m_sigmaR  ( "sigmaR"  , this , right.m_sigmaR )
  , m_beta    ( "beta"    , this , right.m_beta   )
//
  , m_apo2     ( right.m_apo2 ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Apolonios2::~Apolonios2 (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Apolonios2*
Ostap::Models::Apolonios2::clone ( const char* name ) const 
{ return new Ostap::Models::Apolonios2 ( *this , name ) ; }
// ============================================================================
void Ostap::Models::Apolonios2::setPars () const 
{
  //
  m_apo2.setM0      ( m_m0     ) ;
  m_apo2.setSigmaL  ( m_sigmaL ) ;
  m_apo2.setSigmaR  ( m_sigmaR ) ;
  m_apo2.setBeta    ( m_beta   ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Apolonios2::evaluate() const 
{
  //
  setPars () ;
  //
  return m_apo2 ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Apolonios2::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Apolonios2::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_apo2.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================




// ============================================================================
// Bifurcated Gauss 
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BifurcatedGauss::BifurcatedGauss 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          peak      , 
  RooAbsReal&          sigmaL    , 
  RooAbsReal&          sigmaR    ) 
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"       , "Observable"   , this , x      ) 
  , m_peak    ( "peak"    , "peak"         , this , peak   ) 
  , m_sigmaL  ( "sigmaL"  , "sigma(left)"  , this , sigmaL )
  , m_sigmaR  ( "sigmaR"  , "sigma(right)" , this , sigmaR )
//
  , m_bg      ( 0 , 1 , 1 ) 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::BifurcatedGauss::BifurcatedGauss 
( const Ostap::Models::BifurcatedGauss& right , 
  const char*                              name   ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_peak   ( "peak"   , this , right.m_peak   ) 
  , m_sigmaL ( "sigmaL" , this , right.m_sigmaL ) 
  , m_sigmaR ( "sigmaR" , this , right.m_sigmaR ) 
//
  , m_bg     ( right.m_bg ) 
{
  setPars() ;
}
// ============================================================================
// desctructor
// ============================================================================
Ostap::Models::BifurcatedGauss::~BifurcatedGauss(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::BifurcatedGauss*
Ostap::Models::BifurcatedGauss::clone( const char* name ) const 
{ return new Ostap::Models::BifurcatedGauss ( *this , name ) ; }
// ============================================================================
void Ostap::Models::BifurcatedGauss::setPars () const 
{
  //
  m_bg . setPeak   ( m_peak   ) ;
  m_bg . setSigmaL ( m_sigmaL ) ;
  m_bg . setSigmaR ( m_sigmaR ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::BifurcatedGauss::evaluate() const 
{
  //
  setPars () ;
  //
  return m_bg    ( m_x      ) ;
}
// ============================================================================
Int_t Ostap::Models::BifurcatedGauss::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::BifurcatedGauss::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_bg.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
//         GenGaussV1 
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::GenGaussV1::GenGaussV1
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          mu        , 
  RooAbsReal&          alpha     , 
  RooAbsReal&          beta      )  
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"       , "Observable" , this , x      ) 
  , m_mu      ( "mu"      , "mu"         , this , mu     ) 
  , m_alpha   ( "alpha"   , "alpha"      , this , alpha  ) 
  , m_beta    ( "beta"    , "beta"       , this , beta   ) 
//
  , m_ggv1    ( 0 , 1 , 2 ) 
{
  //
  setPars () ;
  //
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::GenGaussV1::GenGaussV1 
( const Ostap::Models::GenGaussV1& right , 
  const char*                         name   ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_mu     ( "mu"     , this , right.m_mu     ) 
  , m_alpha  ( "alpha"  , this , right.m_alpha  ) 
  , m_beta   ( "beta"   , this , right.m_beta   ) 
//
  , m_ggv1   ( right.m_ggv1 ) 
{
  setPars () ;
}
// ============================================================================
// desctructor
// ============================================================================
Ostap::Models::GenGaussV1::~GenGaussV1(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::GenGaussV1*
Ostap::Models::GenGaussV1::clone( const char* name ) const 
{ return new Ostap::Models::GenGaussV1 ( *this , name ) ; }
// ============================================================================
void Ostap::Models::GenGaussV1::setPars () const 
{
  //
  m_ggv1.setMu     ( m_mu    ) ;
  m_ggv1.setAlpha  ( m_alpha ) ;
  m_ggv1.setBeta   ( m_beta  ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::GenGaussV1::evaluate() const 
{
  //
  setPars () ;
  //
  return m_ggv1    ( m_x      ) ;
}
// ============================================================================
Int_t Ostap::Models::GenGaussV1::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::GenGaussV1::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  //
  return m_ggv1.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
//         GenGaussV2
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::GenGaussV2::GenGaussV2
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          xi        , 
  RooAbsReal&          alpha     , 
  RooAbsReal&          kappa     )  
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"       , "Observable" , this , x      ) 
  , m_xi      ( "xi"      , "xi"         , this , xi     ) 
  , m_alpha   ( "alpha"   , "alpha"      , this , alpha  ) 
  , m_kappa   ( "kappa"   , "kappa"      , this , kappa  ) 
//
  , m_ggv2    ( 0 , 1 , 0 ) 
{
  //
  setPars () ;
  //
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::GenGaussV2::GenGaussV2 
( const Ostap::Models::GenGaussV2& right , 
  const char*                         name   ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_xi     ( "xi"     , this , right.m_xi     ) 
  , m_alpha  ( "alpha"  , this , right.m_alpha  ) 
  , m_kappa  ( "kappa"  , this , right.m_kappa  ) 
//
  , m_ggv2   ( right.m_ggv2 ) 
{
  setPars () ;
}
// ============================================================================
// desctructor
// ============================================================================
Ostap::Models::GenGaussV2::~GenGaussV2(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::GenGaussV2*
Ostap::Models::GenGaussV2::clone( const char* name ) const 
{ return new Ostap::Models::GenGaussV2 ( *this , name ) ; }
// ============================================================================
void Ostap::Models::GenGaussV2::setPars () const 
{
  //
  m_ggv2.setXi     ( m_xi    ) ;
  m_ggv2.setAlpha  ( m_alpha ) ;
  m_ggv2.setKappa  ( m_kappa ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::GenGaussV2::evaluate() const 
{
  //
  setPars () ;
  //
  return m_ggv2    ( m_x      ) ;
}
// ============================================================================
Int_t Ostap::Models::GenGaussV2::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::GenGaussV2::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_ggv2.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================




// // ============================================================================
// //         SkewGauss
// // ============================================================================
// // constructor from all parameters 
// // ============================================================================
// Ostap::Models::SkewGauss::SkewGauss
// ( const char*          name      , 
//   const char*          title     ,
//   RooAbsReal&          x         , 
//   RooAbsReal&          xi        , 
//   RooAbsReal&          omega     , 
//   RooAbsReal&          alpha     )  
//   : RooAbsPdf ( name , title ) 
// //
//   , m_x       ( "x"       , "Observable" , this , x      ) 
//   , m_xi      ( "xi"      , "xi"         , this , xi     ) 
//   , m_omega   ( "omega"   , "omega"      , this , omega  ) 
//   , m_alpha   ( "alpha"   , "alpha"      , this , alpha  ) 
// //
//   , m_sg    ( 0 , 1 , 0 ) 
// {
//   //
//   setPars () ;
//   //
// }
// // ============================================================================
// // "copy" constructor 
// // ============================================================================
// Ostap::Models::SkewGauss::SkewGauss
// ( const Ostap::Models::SkewGauss& right , 
//   const char*                        name   ) 
//   : RooAbsPdf ( right , name ) 
// //
//   , m_x      ( "x"      , this , right.m_x      ) 
//   , m_xi     ( "xi"     , this , right.m_xi     ) 
//   , m_omega  ( "omega"  , this , right.m_omega  ) 
//   , m_alpha  ( "alpha"  , this , right.m_alpha  ) 
// //
//   , m_sg     ( right.m_sg ) 
// {
//   setPars () ;
// }
// // ============================================================================
// // desctructor
// // ============================================================================
// Ostap::Models::SkewGauss::~SkewGauss (){}
// // ============================================================================
// // clone 
// // ============================================================================
// Ostap::Models::SkewGauss*
// Ostap::Models::SkewGauss::clone( const char* name ) const 
// { return new Ostap::Models::SkewGauss ( *this , name ) ; }
// // ============================================================================
// void Ostap::Models::SkewGauss::setPars () const 
// {
//   //
//   m_sg . setXi     ( m_xi    ) ;
//   m_sg . setOmega  ( m_omega ) ;
//   m_sg . setAlpha  ( m_alpha ) ;
//   //
// }
// // ============================================================================
// // the actual evaluation of function 
// // ============================================================================
// Double_t Ostap::Models::SkewGauss::evaluate() const 
// {
//   //
//   setPars () ;
//   //
//   return m_sg      ( m_x      ) ;
// }
// // ============================================================================
// Int_t Ostap::Models::SkewGauss::getAnalyticalIntegral
// ( RooArgSet&     allVars      , 
//   RooArgSet&     analVars     ,
//   const char* /* rangename */ ) const 
// {
//   if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
//   return 0 ;
// }
// // ============================================================================
// Double_t Ostap::Models::SkewGauss::analyticalIntegral 
// ( Int_t       code      , 
//   const char* rangeName ) const 
// {
//   assert ( code == 1 ) ;
//   if ( 1 != code ) {}
//   //
//   setPars () ;
//   //
//   return m_sg.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
// }
// // ============================================================================

// ============================================================================
//         Bukin
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Bukin::Bukin 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          peak      , 
  RooAbsReal&          sigma     , 
  RooAbsReal&          xi        ,
  RooAbsReal&          rhoL      ,
  RooAbsReal&          rhoR      )
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"       , "Observable" , this , x      ) 
  , m_peak    ( "peak"    , "peak"       , this , peak   ) 
  , m_sigma   ( "sigma"   , "sigma"      , this , sigma  )
  , m_xi      ( "xi"      , "xi"         , this , xi     )
  , m_rhoL    ( "rhoL"    , "rhoL"       , this , rhoL   )
  , m_rhoR    ( "rhoR"    , "rhoR"       , this , rhoR   )
//
  , m_bukin   ( 10 , 1 , 0 , 0 , 0 ) 
{
  //
  setPars () ;
  //
} 
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Bukin::Bukin
( const Ostap::Models::Bukin& right  , 
  const char*                    name   ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_peak   ( "peak"   , this , right.m_peak   ) 
  , m_sigma  ( "sigma"  , this , right.m_sigma  ) 
  , m_xi     ( "xi"     , this , right.m_xi     ) 
  , m_rhoL   ( "rhoL"   , this , right.m_rhoL   ) 
  , m_rhoR   ( "rhoR"   , this , right.m_rhoR   ) 
//
  , m_bukin  ( right.m_bukin ) 
{
  setPars() ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Bukin::~Bukin(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Bukin*
Ostap::Models::Bukin::clone( const char* name ) const 
{ return new Ostap::Models::Bukin(*this,name) ; }
// ============================================================================
void Ostap::Models::Bukin::setPars () const 
{
  //
  m_bukin.setPeak   ( m_peak  ) ;
  m_bukin.setSigma  ( m_sigma ) ;
  m_bukin.setXi     ( m_xi    ) ;
  m_bukin.setRho_L  ( m_rhoL  ) ;
  m_bukin.setRho_R  ( m_rhoR  ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Bukin::evaluate() const 
{
  //
  setPars () ;
  //
  return m_bukin    ( m_x     ) ;
}
// ============================================================================
Int_t Ostap::Models::Bukin::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Bukin::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  //
  return m_bukin.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================

  
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::StudentT::StudentT 
( const char*          name  , 
  const char*          title ,
  RooAbsReal&          x     ,
  RooAbsReal&          mu    ,
  RooAbsReal&          sigma ,
  RooAbsReal&          n     )
  : RooAbsPdf  (name ,title ) 
//
  , m_x     ( "x"     , "Observable" , this , x     ) 
  , m_mu    ( "mu"    , "Peak"       , this , mu    ) 
  , m_sigma ( "sigma" , "Width"      , this , sigma )
  , m_n     ( "n"     , "N"          , this , n     )
//
  , m_stt   ( 0 , 1 , 1 ) 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::StudentT::StudentT
( const Ostap::Models::StudentT& right , 
  const char*                       name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x     ( "x"     , this , right.m_x     ) 
  , m_mu    ( "mu"    , this , right.m_mu    ) 
  , m_sigma ( "sigma" , this , right.m_sigma )
  , m_n     ( "n"     , this , right.m_n     )
//
  , m_stt   (               right.m_stt    ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::StudentT::~StudentT(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::StudentT*
Ostap::Models::StudentT::clone( const char* name ) const 
{ return new Ostap::Models::StudentT(*this,name) ; }
// ============================================================================
void Ostap::Models::StudentT::setPars () const 
{
  //
  m_stt.setM     ( m_mu    ) ;
  m_stt.setSigma ( m_sigma ) ;
  m_stt.setN     ( m_n     ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::StudentT::evaluate() const 
{
  //
  setPars () ;
  //
  return m_stt   ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::StudentT::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::StudentT::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_stt.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BifurcatedStudentT::BifurcatedStudentT 
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          mu     ,
  RooAbsReal&          sigmaL ,
  RooAbsReal&          sigmaR ,
  RooAbsReal&          nL     ,
  RooAbsReal&          nR     )
  : RooAbsPdf  (name ,title ) 
//
  , m_x      ( "x"     , "Observable" , this , x      ) 
  , m_mu     ( "mu"    , "Peak"       , this , mu     ) 
  , m_sigmaL ( "sigmaL" , "Width(L)"  , this , sigmaL )
  , m_sigmaR ( "sigmaR" , "Width(R)"  , this , sigmaR )
  , m_nL     ( "nL"     , "N(L)"      , this , nL     )
  , m_nR     ( "nR"     , "N(R)"      , this , nR     )
//
  , m_stt   ( 0 , 1 , 1 , 2 , 2 ) 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::BifurcatedStudentT::BifurcatedStudentT 
( const Ostap::Models::BifurcatedStudentT& right , 
  const char*                       name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_mu     ( "mu"     , this , right.m_mu     ) 
  , m_sigmaL ( "sigmaL" , this , right.m_sigmaL )
  , m_sigmaR ( "sigmaR" , this , right.m_sigmaR )
  , m_nL     ( "nL"     , this , right.m_nL     )
  , m_nR     ( "nR"     , this , right.m_nR     )
    //
  , m_stt   (               right.m_stt    ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::BifurcatedStudentT::~BifurcatedStudentT (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::BifurcatedStudentT*
Ostap::Models::BifurcatedStudentT::clone( const char* name ) const 
{ return new Ostap::Models::BifurcatedStudentT(*this,name) ; }
// ============================================================================
void Ostap::Models::BifurcatedStudentT::setPars () const 
{
  //
  m_stt.setM      ( m_mu     ) ;
  m_stt.setSigmaL ( m_sigmaL ) ;
  m_stt.setSigmaR ( m_sigmaR ) ;
  m_stt.setNL     ( m_nL     ) ;
  m_stt.setNR     ( m_nR     ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::BifurcatedStudentT::evaluate() const 
{
  //
  setPars () ;
  //
  return m_stt   ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::BifurcatedStudentT::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::BifurcatedStudentT::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_stt.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
//         Gram-Charlier type A 
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::GramCharlierA::GramCharlierA 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          m0        , 
  RooAbsReal&          sigma     , 
  RooAbsReal&          kappa3    ,
  RooAbsReal&          kappa4    )
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"       , "Observable" , this , x      ) 
  , m_m0      ( "m0"      , "m0"         , this , m0     ) 
  , m_sigma   ( "sigma"   , "sigma"      , this , sigma  )
  , m_kappa3  ( "kappa3"  , "kappa3"     , this , kappa3 )
  , m_kappa4  ( "kappa4"  , "kappa4"     , this , kappa4 )
//
  , m_gca         ( 10 , 1 , 0 , 0 ) 
{
  //
  setPars () ;
  //
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::GramCharlierA::GramCharlierA
( const Ostap::Models::GramCharlierA& right , 
  const char*                            name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_m0     ( "m0"     , this , right.m_m0     ) 
  , m_sigma  ( "sigma"  , this , right.m_sigma  ) 
  , m_kappa3 ( "kappa3" , this , right.m_kappa3 ) 
  , m_kappa4 ( "kappa4" , this , right.m_kappa4 ) 
//
  , m_gca    ( right.m_gca ) 
{
  setPars() ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::GramCharlierA::~GramCharlierA(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::GramCharlierA*
Ostap::Models::GramCharlierA::clone( const char* name ) const 
{ return new Ostap::Models::GramCharlierA(*this,name) ; }
// ============================================================================
void Ostap::Models::GramCharlierA::setPars () const 
{
  //
  m_gca.setM0     ( m_m0     ) ;
  m_gca.setSigma  ( m_sigma  ) ;
  m_gca.setKappa3 ( m_kappa3 ) ;
  m_gca.setKappa4 ( m_kappa4 ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::GramCharlierA::evaluate() const 
{
  //
  setPars() ;
  //
  return m_gca ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::GramCharlierA::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::GramCharlierA::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  //
  return m_gca.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// Two-body phase space 
// ============================================================================
Ostap::Models::PhaseSpace2::PhaseSpace2
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  const double         m1        , 
  const double         m2        ) 
  : RooAbsPdf ( name , title ) 
//
  , m_x   ( "x" , "Observable" , this , x ) 
//  
  , m_ps2 ( m1 , m2 ) 
{}
// ============================================================================
// "copy constructor"
// ============================================================================
Ostap::Models::PhaseSpace2::PhaseSpace2
( const Ostap::Models::PhaseSpace2& right , const char* name )  
  : RooAbsPdf ( right , name )
//
  , m_x       ( "x"       , this , right.m_x      ) 
//
  , m_ps2     ( right.m_ps2 ) 
//
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PhaseSpace2::~PhaseSpace2(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PhaseSpace2*
Ostap::Models::PhaseSpace2::clone( const char* name ) const 
{ return new Ostap::Models::PhaseSpace2(*this,name) ; }
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PhaseSpace2::evaluate() const 
{ return m_ps2 ( m_x ) ; }
// ============================================================================
Int_t Ostap::Models::PhaseSpace2::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PhaseSpace2::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  return m_ps2.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
// Left-edge of N-body phase space 
// ============================================================================
Ostap::Models::PhaseSpaceLeft::PhaseSpaceLeft
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          threshold ,
  const unsigned short N         ) 
  : RooAbsPdf ( name , title ) 
//
  , m_x         ( "x"  , "Observable" , this , x         ) 
  , m_threshold ( "th" , "Threshold"  , this , threshold  ) 
//  
  , m_left ( 10 , N ) 
{
  setPars () ;  
}
// ============================================================================
// "copy constructor"
// ============================================================================
Ostap::Models::PhaseSpaceLeft::PhaseSpaceLeft
( const Ostap::Models::PhaseSpaceLeft& right , const char* name )  
  : RooAbsPdf ( right , name )
//
  , m_x         ( "x"  , this , right.m_x         ) 
  , m_threshold ( "tr" , this , right.m_threshold ) 
//
  , m_left     ( right.m_left ) 
//
{
  setPars () ;  
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PhaseSpaceLeft::~PhaseSpaceLeft(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PhaseSpaceLeft*
Ostap::Models::PhaseSpaceLeft::clone( const char* name ) const 
{ return new Ostap::Models::PhaseSpaceLeft(*this,name) ; }
// ============================================================================
void Ostap::Models::PhaseSpaceLeft::setPars () const 
{
  //
  m_left.setThreshold ( m_threshold ) ;  
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PhaseSpaceLeft::evaluate() const 
{
  //
  setPars () ;
  //
  return m_left ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::PhaseSpaceLeft::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PhaseSpaceLeft::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_left.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
// Right-edge of L-body phase space in N-body decays  
// ============================================================================
Ostap::Models::PhaseSpaceRight::PhaseSpaceRight
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          threshold ,
  const unsigned short L         ,
  const unsigned short N         ) 
  : RooAbsPdf ( name , title ) 
//
  , m_x         ( "x"  , "Observable" , this , x         ) 
  , m_threshold ( "th" , "Threshold"  , this , threshold  ) 
//  
  , m_right ( 10 , L , N ) 
{
  setPars () ;
}
// ============================================================================
// "copy constructor"
// ============================================================================
Ostap::Models::PhaseSpaceRight::PhaseSpaceRight
( const Ostap::Models::PhaseSpaceRight& right , const char* name )  
  : RooAbsPdf ( right , name )
//
  , m_x         ( "x"  , this , right.m_x         ) 
  , m_threshold ( "tr" , this , right.m_threshold ) 
//
  , m_right     ( right.m_right ) 
//
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PhaseSpaceRight::~PhaseSpaceRight(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PhaseSpaceRight*
Ostap::Models::PhaseSpaceRight::clone( const char* name ) const 
{ return new Ostap::Models::PhaseSpaceRight(*this,name) ; }
// ============================================================================
void Ostap::Models::PhaseSpaceRight::setPars () const 
{
  //
  m_right.setThreshold ( m_threshold ) ;  
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PhaseSpaceRight::evaluate() const 
{
  //
  setPars () ;
  //
  return m_right ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::PhaseSpaceRight::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PhaseSpaceRight::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_right.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::PhaseSpaceNL::PhaseSpaceNL
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          low       ,
  RooAbsReal&          high      ,    
  const unsigned short N         , 
  const unsigned short L         ) 
  : RooAbsPdf ( name  , title )
//
  , m_x     ( "x"     , "Observable" , this , x     ) 
  , m_low   ( "low"   , "m(low)"     , this , low   ) 
  , m_high  ( "high"  , "m(high)"    , this , high  ) 
//
  , m_ps    ( 1 , 2 , N , L ) 
{
  setPars () ;
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Models::PhaseSpaceNL::PhaseSpaceNL
( const Ostap::Models::PhaseSpaceNL& right , const char* name ) 
  : RooAbsPdf ( right  , name )
//
  , m_x     ( "x"     , this , right.m_x     ) 
  , m_low   ( "low"   , this , right.m_low   ) 
  , m_high  ( "low"   , this , right.m_high  ) 
//
  , m_ps    (                  right.m_ps    ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PhaseSpaceNL::~PhaseSpaceNL(){}
// ===========================================================================
// clone 
// ============================================================================
Ostap::Models::PhaseSpaceNL*
Ostap::Models::PhaseSpaceNL::clone( const char* name ) const 
{ return new Ostap::Models::PhaseSpaceNL(*this,name) ; }
// ============================================================================
void Ostap::Models::PhaseSpaceNL::setPars () const 
{
  //
  m_ps.setThresholds ( m_low , m_high ) ;  
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PhaseSpaceNL::evaluate() const 
{
  //
  setPars () ;
  //
  return m_ps ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::PhaseSpaceNL::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PhaseSpaceNL::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_ps.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================

// ============================================================================
// Two-body phase space from 3-body decays 
// ============================================================================
Ostap::Models::PhaseSpace23L::PhaseSpace23L
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  const double         m1        , 
  const double         m2        ,
  const double         m3        ,
  const double         m         ,
  const unsigned short L         , 
  const unsigned short l         ) 
  : RooAbsPdf ( name , title ) 
  , m_x       ( "x" , "Observable" , this , x ) 
  , m_ps23L   ( m1 , m2 , m3 , m , L , l ) 
{}
// ============================================================================
// "copy constructor"
// ============================================================================
Ostap::Models::PhaseSpace23L::PhaseSpace23L
( const Ostap::Models::PhaseSpace23L& right , const char* name )  
  : RooAbsPdf ( right , name )
  , m_x       ( "x"       , this , right.m_x      ) 
  , m_ps23L   ( right.m_ps23L ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PhaseSpace23L::~PhaseSpace23L(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PhaseSpace23L*
Ostap::Models::PhaseSpace23L::clone( const char* name ) const 
{ return new Ostap::Models::PhaseSpace23L(*this,name) ; }
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PhaseSpace23L::evaluate() const 
{ return m_ps23L ( m_x ) ; }
// ============================================================================
Int_t Ostap::Models::PhaseSpace23L::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PhaseSpace23L::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  return m_ps23L.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================






// ============================================================================
//  PhaseSpace x poly
// ============================================================================
Ostap::Models::PhaseSpacePol::PhaseSpacePol 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  const double         low       , 
  const double         high      ,
  const unsigned short N         , 
  const unsigned short L         , 
  RooAbsReal&          phi1      ) 
  : RooAbsPdf ( name , title ) 
//
  , m_x        ( "x"       , "Observable"   , this , x   ) 
  , m_phis     ( "phi"     , "Coefficients" , this )
    //
  , m_ps       ( low , high , L , N , 1 ) 
{
  m_phis.add ( phi1 ) ;
  //
  const RooRealVar* v = dynamic_cast<RooRealVar*>(&x) ;
  if ( 0 != v ) 
  {
    const double xmin = std::max ( std::min ( low , high ) , v -> getMin () ) ;
    const double xmax = std::min ( std::max ( low , high ) , v -> getMax () ) ;
    m_ps = Ostap::Math::PhaseSpacePol ( m_ps.phasespace() , 1 , xmin , xmax ) ;  
  }
  //
  setPars() ;
}
// ============================================================================
//  PhaseSpace x poly
// ============================================================================
Ostap::Models::PhaseSpacePol::PhaseSpacePol 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  const double         low       , 
  const double         high      ,
  const unsigned short N         , 
  const unsigned short L         , 
  RooAbsReal&          phi1      ,
  RooAbsReal&          phi2      ) 
  : RooAbsPdf ( name , title ) 
//
  , m_x        ( "x"       , "Observable"   , this , x   ) 
  , m_phis     ( "phi"     , "Coefficients" , this )
//
  , m_ps       (     low , high , L , N , 2 ) 
{
  m_phis.add ( phi1 ) ;
  m_phis.add ( phi2 ) ;
  //
  const RooRealVar* v = dynamic_cast<RooRealVar*>(&x) ;
  if ( 0 != v ) 
  {
    const double xmin = std::max ( std::min ( low , high ) , v -> getMin () ) ;
    const double xmax = std::min ( std::max ( low , high ) , v -> getMax () ) ;
    m_ps = Ostap::Math::PhaseSpacePol ( m_ps.phasespace() , 2 , xmin , xmax ) ;  
  }
  //
  setPars() ;
}
// ============================================================================
//  PhaseSpace x poly
// ============================================================================
Ostap::Models::PhaseSpacePol::PhaseSpacePol 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  const double         low       , 
  const double         high      ,
  const unsigned short N         , 
  const unsigned short L         , 
  RooAbsReal&          phi1      ,
  RooAbsReal&          phi2      , 
  RooAbsReal&          phi3      ) 
  : RooAbsPdf ( name , title ) 
//
  , m_x        ( "x"       , "Observable"   , this , x   ) 
  , m_phis     ( "phi"     , "Coefficients" , this )
//
  , m_ps       (     low , high , L , N , 3 ) 
{
  m_phis.add ( phi1 ) ;
  m_phis.add ( phi2 ) ;
  m_phis.add ( phi3 ) ;
  //
  const RooRealVar* v = dynamic_cast<RooRealVar*>(&x) ;
  if ( 0 != v ) 
  {
    const double xmin = std::max ( std::min ( low , high ) , v -> getMin () ) ;
    const double xmax = std::min ( std::max ( low , high ) , v -> getMax () ) ;
    m_ps = Ostap::Math::PhaseSpacePol ( m_ps.phasespace() , 3 , xmin , xmax ) ;  
  }
  //
  setPars() ;
}
// ============================================================================
//  PhaseSpace x poly
// ============================================================================
Ostap::Models::PhaseSpacePol::PhaseSpacePol 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  const double         low       , 
  const double         high      ,
  const unsigned short N         , 
  const unsigned short L         , 
  RooArgList&          phis      )
  : RooAbsPdf ( name , title ) 
//
  , m_x        ( "x"       , "Observable"   , this , x    ) 
  , m_phis     ( "phi"     , "Coefficients" , this ) 
//
  , m_ps       (                  low , high , L , N , phis.getSize() ) 
{
  //
  RooAbsArg*   coef = 0 ;
  unsigned     num  = 0 ;
  Ostap::Utils::Iterator   tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_ps.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  //
  const RooRealVar* v = dynamic_cast<RooRealVar*>(&x) ;
  if ( 0 != v ) 
  {
    const double xmin = std::max ( std::min ( low , high ) , v -> getMin () ) ;
    const double xmax = std::min ( std::max ( low , high ) , v -> getMax () ) ;
    m_ps = Ostap::Math::PhaseSpacePol ( m_ps.phasespace() , phis.getSize() , xmin , xmax ) ;  
  }
  //
  setPars() ;
}
// ============================================================================
//  PhaseSpace x poly
// ============================================================================
Ostap::Models::PhaseSpacePol::PhaseSpacePol 
( const char*                      name      , 
  const char*                      title     ,
  RooAbsReal&                      x         ,
  const Ostap::Math::PhaseSpaceNL& ps        , 
  RooAbsReal&                      phi1      )
  : RooAbsPdf ( name , title ) 
//
  , m_x        ( "x"       , "Observable"   , this , x   ) 
  , m_phis     ( "phi"     , "Coefficients" , this )
//
  , m_ps       ( ps , 1 ) 
{
  m_phis.add ( phi1 ) ;
  //
  const RooRealVar* v = dynamic_cast<RooRealVar*>(&x) ;
  if ( 0 != v ) 
  {
    const double xmin = std::max ( ps .  lowEdge () , v -> getMin () ) ;
    const double xmax = std::min ( ps . highEdge () , v -> getMax () ) ;
    m_ps = Ostap::Math::PhaseSpacePol ( m_ps.phasespace() , 1 , xmin , xmax ) ;  
  }
  //
  setPars() ;
}
// ============================================================================
//  PhaseSpace x poly
// ============================================================================
Ostap::Models::PhaseSpacePol::PhaseSpacePol 
( const char*                      name      , 
  const char*                      title     ,
  RooAbsReal&                      x         ,
  const Ostap::Math::PhaseSpaceNL& ps        , 
  RooAbsReal&                      phi1      ,
  RooAbsReal&                      phi2      )
  : RooAbsPdf ( name , title ) 
//
  , m_x        ( "x"       , "Observable"   , this , x   ) 
  , m_phis     ( "phi"     , "Coefficients" , this )
//
  , m_ps       ( ps , 2 ) 
{
  m_phis.add ( phi1 ) ;
  m_phis.add ( phi2 ) ;
  //
  const RooRealVar* v = dynamic_cast<RooRealVar*>(&x) ;
  if ( 0 != v ) 
  {
    const double xmin = std::max ( ps .  lowEdge () , v -> getMin () ) ;
    const double xmax = std::min ( ps . highEdge () , v -> getMax () ) ;
    m_ps = Ostap::Math::PhaseSpacePol ( m_ps.phasespace() , 2 , xmin , xmax ) ;  
  }
  //
  setPars() ;
}
// ============================================================================
//  PhaseSpace x poly
// ============================================================================
Ostap::Models::PhaseSpacePol::PhaseSpacePol 
( const char*                      name      , 
  const char*                      title     ,
  RooAbsReal&                      x         ,
  const Ostap::Math::PhaseSpaceNL& ps        , 
  RooAbsReal&                      phi1      ,
  RooAbsReal&                      phi2      ,
  RooAbsReal&                      phi3      )
  : RooAbsPdf ( name , title ) 
//
  , m_x        ( "x"       , "Observable"   , this , x   ) 
  , m_phis     ( "phi"     , "Coefficients" , this )
//
  , m_ps       ( ps , 3 ) 
{
  m_phis.add ( phi1 ) ;
  m_phis.add ( phi2 ) ;
  m_phis.add ( phi3 ) ;
  //
  const RooRealVar* v = dynamic_cast<RooRealVar*>(&x) ;
  if ( 0 != v ) 
  {
    const double xmin = std::max ( ps .  lowEdge () , v -> getMin () ) ;
    const double xmax = std::min ( ps . highEdge () , v -> getMax () ) ;
    m_ps = Ostap::Math::PhaseSpacePol ( m_ps.phasespace() , 3 , xmin , xmax ) ;  
  }
  //
  setPars() ;
}
// ============================================================================
//  PhaseSpace x poly
// ============================================================================
Ostap::Models::PhaseSpacePol::PhaseSpacePol 
( const char*                      name      , 
  const char*                      title     ,
  RooAbsReal&                      x         ,
  const Ostap::Math::PhaseSpaceNL& ps        , 
  RooArgList&                      phis      )
  : RooAbsPdf ( name , title ) 
//
  , m_x        ( "x"       , "Observable"   , this , x    ) 
  , m_phis     ( "phi"     , "Coefficients" , this )
    //
  , m_ps       ( ps , phis.getSize() ) 
{
  //
  RooAbsArg*   coef = 0 ;
  unsigned     num  = 0 ;
  Ostap::Utils::Iterator   tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_ps.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  //
  const RooRealVar* v = dynamic_cast<RooRealVar*>(&x) ;
  if ( 0 != v ) 
  {
    const double xmin = std::max ( ps .  lowEdge () , v -> getMin () ) ;
    const double xmax = std::min ( ps . highEdge () , v -> getMax () ) ;
    m_ps = Ostap::Math::PhaseSpacePol ( m_ps.phasespace() , phis.getSize()  , xmin , xmax ) ;  
  }
  //
  setPars() ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PhaseSpacePol::~PhaseSpacePol () {}
// ======================================================================
// "copy" constructor 
// ======================================================================
Ostap::Models::PhaseSpacePol::PhaseSpacePol 
( const Ostap::Models::PhaseSpacePol& right , 
  const char*                            name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x    ) 
  , m_phis     ( "phis"   , this , right.m_phis ) 
//
  , m_ps       ( right.m_ps       ) 
{
  setPars() ;
}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PhaseSpacePol* 
Ostap::Models::PhaseSpacePol::clone ( const char* name ) const 
{ return new Ostap::Models::PhaseSpacePol ( *this , name ) ; }
// ============================================================================
void Ostap::Models::PhaseSpacePol::setPars () const 
{
  //
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_ps.setPar ( k  , phi ) ;
    //
    ++k ;
  }
  //
}
// ============================================================================
// evaluate the function
// ============================================================================
Double_t Ostap::Models::PhaseSpacePol::evaluate () const 
{
  //
  setPars () ;
  //
  return m_ps ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::PhaseSpacePol::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PhaseSpacePol::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_ps.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
// generic polinomial
// ============================================================================
Ostap::Models::PolyPositive::PolyPositive 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  const RooArgList&    phis      , 
  const double         xmin      , 
  const double         xmax      ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x ) 
  , m_phis     ( "phi"     , "Coefficients" , this     )
    //
  , m_positive ( phis.getSize() , xmin , xmax ) 
{
  //
  Ostap::Utils::Iterator tmp ( phis ) ;
  RooAbsArg* coef = 0 ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PolyPositive::PolyPositive
( const Ostap::Models::PolyPositive&  right ,      
  const char*                            name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_positive ( right.m_positive ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PolyPositive::~PolyPositive() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PolyPositive*
Ostap::Models::PolyPositive::clone( const char* name ) const 
{ return new Ostap::Models::PolyPositive(*this,name) ; }
// ============================================================================
void Ostap::Models::PolyPositive::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_positive.setPar ( k  , phi ) ;
    //
    ++k ;
  }
}
//
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PolyPositive::evaluate() const 
{
  //
  setPars () ;
  //
  return m_positive ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::PolyPositive::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PolyPositive::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_positive.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
// generic even polinomial 
// ============================================================================
Ostap::Models::PolyPositiveEven::PolyPositiveEven 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  const RooArgList&    phis      , 
  const double         xmin      , 
  const double         xmax      ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x ) 
  , m_phis     ( "phi"     , "Coefficients" , this     )
    //
  , m_even     ( phis.getSize() , xmin , xmax ) 
{
  //
  Ostap::Utils::Iterator tmp ( phis ) ;
  RooAbsArg* coef = 0 ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PolyPositiveEven::PolyPositiveEven
( const Ostap::Models::PolyPositiveEven&  right ,      
  const char*                                name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_even     ( right.m_even ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PolyPositiveEven::~PolyPositiveEven() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PolyPositiveEven*
Ostap::Models::PolyPositiveEven::clone( const char* name ) const 
{ return new Ostap::Models::PolyPositiveEven(*this,name) ; }
// ============================================================================
void Ostap::Models::PolyPositiveEven::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_even.setPar ( k  , phi ) ;
    //
    ++k ;
  }
}
//
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PolyPositiveEven::evaluate() const 
{
  //
  setPars () ;
  //
  return m_even ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::PolyPositiveEven::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PolyPositiveEven::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_even.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// monothonic polinomial
// ============================================================================
Ostap::Models::PolyMonothonic::PolyMonothonic
( const char*          name       , 
  const char*          title      ,
  RooAbsReal&          x          ,
  const RooArgList&    phis       , 
  const double         xmin       , 
  const double         xmax       , 
  const bool           increasing ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x ) 
  , m_phis     ( "phi"     , "Coefficients" , this     )
//
  , m_monothonic ( phis.getSize() , xmin , xmax , increasing ) 
{
  //
  Ostap::Utils::Iterator tmp ( phis ) ;
  RooAbsArg* coef = 0 ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PolyMonothonic::PolyMonothonic
( const Ostap::Models::PolyMonothonic&  right ,      
  const char*                              name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x          ( "x"      , this , right.m_x     ) 
  , m_phis       ( "phis"   , this , right.m_phis  ) 
    //
  , m_monothonic ( right.m_monothonic ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PolyMonothonic::~PolyMonothonic() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PolyMonothonic*
Ostap::Models::PolyMonothonic::clone( const char* name ) const 
{ return new Ostap::Models::PolyMonothonic(*this,name) ; }
// ============================================================================
void Ostap::Models::PolyMonothonic::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_monothonic.setPar ( k  , phi ) ;
    //
    ++k ;
  }
}
//
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PolyMonothonic::evaluate() const 
{
  //
  setPars () ;
  //
  return m_monothonic ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::PolyMonothonic::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PolyMonothonic::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_monothonic.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================




// ============================================================================
// convex polinomial
// ============================================================================
Ostap::Models::PolyConvex::PolyConvex
( const char*          name       , 
  const char*          title      ,
  RooAbsReal&          x          ,
  const RooArgList&    phis       , 
  const double         xmin       , 
  const double         xmax       , 
  const bool           increasing ,
  const bool           convex     ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x ) 
  , m_phis     ( "phi"     , "Coefficients" , this     )
    //
  , m_convex ( phis.getSize() , xmin , xmax , increasing , convex ) 
{
  //
  RooAbsArg* coef = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PolyConvex::PolyConvex
( const Ostap::Models::PolyConvex&  right ,      
  const char*                          name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x          ( "x"      , this , right.m_x     ) 
  , m_phis       ( "phis"   , this , right.m_phis  ) 
    //
  , m_convex ( right.m_convex ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PolyConvex::~PolyConvex () { }
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PolyConvex*
Ostap::Models::PolyConvex::clone( const char* name ) const 
{ return new Ostap::Models::PolyConvex(*this,name) ; }
// ============================================================================
void Ostap::Models::PolyConvex::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_convex.setPar ( k  , phi ) ;
    //
    ++k ;
  }
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PolyConvex::evaluate() const 
{
  //
  setPars () ;
  //
  return m_convex ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::PolyConvex::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PolyConvex::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_convex.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================




// ============================================================================
// convex polinomial
// ============================================================================
Ostap::Models::PolyConvexOnly::PolyConvexOnly
( const char*          name       , 
  const char*          title      ,
  RooAbsReal&          x          ,
  const RooArgList&    phis       , 
  const double         xmin       , 
  const double         xmax       , 
  const bool           convex     ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x ) 
  , m_phis     ( "phi"     , "Coefficients" , this     )
//
  , m_convex ( phis.getSize() , xmin , xmax , convex ) 
{
  //
  RooAbsArg* coef = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PolyConvexOnly::PolyConvexOnly
( const Ostap::Models::PolyConvexOnly&  right ,      
  const char*                              name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x          ( "x"      , this , right.m_x     ) 
  , m_phis       ( "phis"   , this , right.m_phis  ) 
    //
  , m_convex ( right.m_convex ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PolyConvexOnly::~PolyConvexOnly () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PolyConvexOnly*
Ostap::Models::PolyConvexOnly::clone( const char* name ) const 
{ return new Ostap::Models::PolyConvexOnly(*this,name) ; }
// ============================================================================
void Ostap::Models::PolyConvexOnly::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  Ostap::Utils::Iterator it ( m_phis ) ;
  unsigned short k = 0 ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_convex.setPar ( k  , phi ) ;
    //
    ++k ;
  }
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PolyConvexOnly::evaluate() const 
{
  //
  setPars () ;
  //
  return m_convex ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::PolyConvexOnly::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PolyConvexOnly::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_convex.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// sigmoid polinomial
// ============================================================================
Ostap::Models::PolySigmoid::PolySigmoid
( const char*          name       , 
  const char*          title      ,
  RooAbsReal&          x          ,
  const RooArgList&    phis       , 
  const double         xmin       , 
  const double         xmax       , 
  RooAbsReal&          alpha      ,
  RooAbsReal&          x0         )
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x      ) 
  , m_phis     ( "phi"     , "Coefficients" , this          )
  , m_alpha    ( "alpha"   , "Alpha"        , this , alpha  )
  , m_x0       ( "x0"      , "X0"           , this , x0     )
    //
  , m_sigmoid  ( phis.getSize() , xmin , xmax , m_alpha , m_x0 ) 
{
  //
  RooAbsArg* coef = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PolySigmoid::PolySigmoid
( const Ostap::Models::PolySigmoid&  right ,      
  const char*                           name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x          ( "x"      , this , right.m_x     ) 
  , m_phis       ( "phis"   , this , right.m_phis  ) 
  , m_alpha      ( "alpha"  , this , right.m_alpha ) 
  , m_x0         ( "x0"     , this , right.m_x0    ) 
    //
  , m_sigmoid ( right.m_sigmoid ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PolySigmoid::~PolySigmoid(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PolySigmoid*
Ostap::Models::PolySigmoid::clone( const char* name ) const 
{ return new Ostap::Models::PolySigmoid(*this,name) ; }
// ============================================================================
void Ostap::Models::PolySigmoid::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_sigmoid.setPar ( k  , phi ) ;
    //
    ++k ;
  }
  //
  m_sigmoid.setAlpha ( m_alpha ) ;
  m_sigmoid.setX0    ( m_x0    ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PolySigmoid::evaluate() const 
{
  //
  setPars () ;
  //
  return m_sigmoid ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::PolySigmoid::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PolySigmoid::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_sigmoid.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// positive spline 
// ============================================================================
/*  constructor with the spline 
 *  @param name  the name 
 *  @param title the  title
 *  @param x     the  variable 
 *  @param spine the spline  
 *  @param phis  vector of parameters 
 */
// ============================================================================
Ostap::Models::PositiveSpline::PositiveSpline 
( const char*                        name, 
  const char*                        title     ,
  RooAbsReal&                        x         ,
  const Ostap::Math::PositiveSpline& spline    ,   // the spline 
  RooArgList&                        phis      )   // parameters
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x ) 
  , m_phis     ( "phi"     , "Coefficients" , this     )
    //
  , m_spline   ( spline ) 
{
  //
  RooAbsArg* coef = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PositiveSpline::PositiveSpline
( const Ostap::Models::PositiveSpline&  right ,      
  const char*                            name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_spline ( right.m_spline ) 
{
  setPars () ;
}

// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PositiveSpline::~PositiveSpline() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PositiveSpline*
Ostap::Models::PositiveSpline::clone( const char* name ) const 
{ return new Ostap::Models::PositiveSpline(*this,name) ; }
// ============================================================================
void Ostap::Models::PositiveSpline::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_spline.setPar ( k  , phi ) ;
    //
    ++k ;
  }
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PositiveSpline::evaluate() const 
{
  //
  setPars () ;
  //
  return m_spline ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::PositiveSpline::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PositiveSpline::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_spline.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================





// ============================================================================
// monothonic spline 
// ============================================================================
/* constructor with the spline 
 *  @param name  the name 
 *  @param title the  title
 *  @param x     the  variable 
 *  @param spine the spline  
 *  @param phis  vector of parameters 
 */
// ============================================================================
Ostap::Models::MonothonicSpline::MonothonicSpline 
( const char*                          name, 
  const char*                          title     ,
  RooAbsReal&                          x         ,
  const Ostap::Math::MonothonicSpline& spline    ,   // the spline 
  RooArgList&                          phis      )   // parameters
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x ) 
  , m_phis     ( "phi"     , "Coefficients" , this     )
    //
  , m_spline   ( spline ) 
{
  //
  RooAbsArg* coef = 0 ;
  Ostap::Utils::Iterator  tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::MonothonicSpline::MonothonicSpline 
( const Ostap::Models::MonothonicSpline&  right ,      
  const char*                                name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_spline ( right.m_spline ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::MonothonicSpline::~MonothonicSpline() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::MonothonicSpline*
Ostap::Models::MonothonicSpline::clone( const char* name ) const 
{ return new Ostap::Models::MonothonicSpline(*this,name) ; }
// ============================================================================
void Ostap::Models::MonothonicSpline::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator  it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_spline.setPar ( k  , phi ) ;
    //
    ++k ;
  }
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::MonothonicSpline::evaluate() const 
{
  //
  setPars () ;
  //
  return m_spline ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::MonothonicSpline::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::MonothonicSpline::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_spline.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================

// ============================================================================
// convex spline 
// ============================================================================
/* constructor with the spline 
 *  @param name  the name 
 *  @param title the  title
 *  @param x     the  variable 
 *  @param spine the spline  
 *  @param phis  vector of parameters 
 */
// ============================================================================
Ostap::Models::ConvexSpline::ConvexSpline 
( const char*                          name, 
  const char*                          title     ,
  RooAbsReal&                          x         ,
  const Ostap::Math::ConvexSpline& spline    ,   // the spline 
  RooArgList&                          phis      )   // parameters
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x ) 
  , m_phis     ( "phi"     , "Coefficients" , this     )
    //
  , m_spline   ( spline ) 
{
  //
  RooAbsArg* coef = 0 ;
  Ostap::Utils::Iterator tmp  ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::ConvexSpline::ConvexSpline 
( const Ostap::Models::ConvexSpline&  right ,      
  const char*                                name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_spline ( right.m_spline ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::ConvexSpline::~ConvexSpline() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::ConvexSpline*
Ostap::Models::ConvexSpline::clone( const char* name ) const 
{ return new Ostap::Models::ConvexSpline(*this,name) ; }
// ============================================================================
void Ostap::Models::ConvexSpline::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it  ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_spline.setPar ( k  , phi ) ;
    //
    ++k ;
  }
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::ConvexSpline::evaluate() const 
{
  //
  setPars () ;
  //
  return m_spline ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::ConvexSpline::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::ConvexSpline::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_spline.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================

// ============================================================================
// convex spline 
// ============================================================================
/* constructor with the spline 
 *  @param name  the name 
 *  @param title the  title
 *  @param x     the  variable 
 *  @param spine the spline  
 *  @param phis  vector of parameters 
 */
// ============================================================================
Ostap::Models::ConvexOnlySpline::ConvexOnlySpline 
( const char*                          name, 
  const char*                          title     ,
  RooAbsReal&                          x         ,
  const Ostap::Math::ConvexOnlySpline& spline    ,   // the spline 
  RooArgList&                          phis      )   // parameters
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x ) 
  , m_phis     ( "phi"     , "Coefficients" , this     )
    //
  , m_spline   ( spline ) 
{
  //
  RooAbsArg* coef = 0 ;
  Ostap::Utils::Iterator tmp  ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::ConvexOnlySpline::ConvexOnlySpline 
( const Ostap::Models::ConvexOnlySpline&  right ,      
  const char*                                name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_spline ( right.m_spline ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::ConvexOnlySpline::~ConvexOnlySpline() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::ConvexOnlySpline*
Ostap::Models::ConvexOnlySpline::clone( const char* name ) const 
{ return new Ostap::Models::ConvexOnlySpline(*this,name) ; }
// ============================================================================
void Ostap::Models::ConvexOnlySpline::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it  ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_spline.setPar ( k  , phi ) ;
    //
    ++k ;
  }
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::ConvexOnlySpline::evaluate() const 
{
  //
  setPars () ;
  //
  return m_spline ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::ConvexOnlySpline::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::ConvexOnlySpline::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_spline.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
// generic polinomial times exponent 
// ============================================================================
Ostap::Models::ExpoPositive::ExpoPositive
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          tau       , 
  const RooArgList&    phis      , 
  const double         xmin      , 
  const double         xmax      ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x   ) 
  , m_tau      ( "tau"     , "Exponential"  , this , tau )
  , m_phis     ( "phi"     , "Coefficients" , this )
//
  , m_positive ( phis.getSize() , 0 , xmin , xmax ) 
{
  //
  RooAbsArg* coef = 0 ;
  Ostap::Utils::Iterator  tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  //
  m_positive.setTau ( m_tau ) ;
  setPars () ;
  //
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::ExpoPositive::ExpoPositive
( const Ostap::Models::ExpoPositive&  right ,      
  const char*                            name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x        ( "x"      , this , right.m_x    ) 
  , m_tau      ( "tau"    , this , right.m_tau  )
  , m_phis     ( "phis"   , this , right.m_phis ) 
//
  , m_positive ( right.m_positive ) 
{
  setPars () ;
  //
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::ExpoPositive::~ExpoPositive() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::ExpoPositive*
Ostap::Models::ExpoPositive::clone( const char* name ) const 
{ return new Ostap::Models::ExpoPositive(*this,name) ; }
// ============================================================================
void Ostap::Models::ExpoPositive::setPars () const 
{
  //
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator  it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_positive.setPar ( k  , phi ) ;
    //
    ++k ;
  }
  //
  m_positive.setTau ( m_tau ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::ExpoPositive::evaluate() const 
{
  //
  setPars () ;
  //
  return m_positive ( m_x   ) ;
}
// ============================================================================
Int_t Ostap::Models::ExpoPositive::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::ExpoPositive::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_positive.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// generic polinomial times exponent 
// ============================================================================
Ostap::Models::TwoExpoPositive::TwoExpoPositive
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          alpha     , 
  RooAbsReal&          delta     , 
  RooAbsReal&          x0        , 
  const RooArgList&    phis      , 
  const double         xmin      , 
  const double         xmax      ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x     ) 
  , m_alpha    ( "alpha"   , "slope 1"      , this , alpha )
  , m_delta    ( "delta"   , "delta slope"  , this , delta )
  , m_x0       ( "x0"      , "threshold"    , this , x0    )
  , m_phis     ( "phi"     , "Coefficients" , this )
    //
  , m_2expopos ( phis.getSize() , 1 , 2 , 1 , xmin , xmax ) 
{
  //
  RooAbsArg* coef = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
  }
  //
  setPars () ;
  //
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::TwoExpoPositive::TwoExpoPositive
( const Ostap::Models::TwoExpoPositive&  right ,      
  const char*                               name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_alpha    ( "alpha"  , this , right.m_alpha )
  , m_delta    ( "delta"  , this , right.m_delta )
  , m_x0       ( "x0"     , this , right.m_x0    )
  , m_phis     ( "phis"   , this , right.m_phis  ) 
//
  , m_2expopos ( right.m_2expopos ) 
{
  setPars () ;
  //
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::TwoExpoPositive::~TwoExpoPositive() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::TwoExpoPositive*
Ostap::Models::TwoExpoPositive::clone( const char* name ) const 
{ return new Ostap::Models::TwoExpoPositive(*this,name) ; }
// ============================================================================
void Ostap::Models::TwoExpoPositive::setPars () const 
{
  //
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phi   = r->getVal ( nset ) ;
    //
    m_2expopos.setPar ( k  , phi ) ;
    //
    ++k ;
  }
  //
  m_2expopos.setAlpha ( m_alpha ) ;
  m_2expopos.setDelta ( m_delta ) ;
  m_2expopos.setX0    ( m_x0    ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::TwoExpoPositive::evaluate() const 
{
  //
  setPars () ;
  //
  return m_2expopos( m_x   ) ;
}
// ============================================================================
Int_t Ostap::Models::TwoExpoPositive::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::TwoExpoPositive::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_2expopos.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::GammaDist::GammaDist
( const char*          name  , 
  const char*          title ,
  RooAbsReal&          x     ,
  RooAbsReal&          k     ,
  RooAbsReal&          theta )
  : RooAbsPdf  (name ,title ) 
//
  , m_x       ( "x"     , "Observable" , this , x     ) 
  , m_k       ( "k"     , "Shape"      , this , k     ) 
  , m_theta   ( "theta" , "Scale"      , this , theta )
//
  , m_gamma   ( 1 , 1 ) 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::GammaDist::GammaDist
( const Ostap::Models::GammaDist& right ,
  const char*                        name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x     ( "x"     , this , right.m_x     ) 
  , m_k     ( "k"     , this , right.m_k     ) 
  , m_theta ( "theta" , this , right.m_theta )
//
  , m_gamma (                  right.m_gamma ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::GammaDist::~GammaDist () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::GammaDist*
Ostap::Models::GammaDist::clone( const char* name ) const 
{ return new Ostap::Models::GammaDist ( *this , name) ; }
// ============================================================================
void Ostap::Models::GammaDist::setPars () const 
{
  //
  m_gamma.setK     ( m_k     ) ;
  m_gamma.setTheta ( m_theta ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::GammaDist::evaluate() const 
{
  //
  setPars () ;
  //
  return m_gamma   ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::GammaDist::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::GammaDist::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  // 
  return m_gamma.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::GenGammaDist::GenGammaDist
( const char*          name  , 
  const char*          title ,
  RooAbsReal&          x     ,
  RooAbsReal&          k     ,
  RooAbsReal&          theta ,
  RooAbsReal&          p     ,
  RooAbsReal&          low   )
  : RooAbsPdf  (name ,title ) 
//
  , m_x       ( "x"     , "Observable" , this , x     ) 
  , m_k       ( "k"     , "Shape"      , this , k     ) 
  , m_theta   ( "theta" , "Scale"      , this , theta )
  , m_p       ( "p"     , "P"          , this , p     )
  , m_low     ( "low"   , "Low"        , this , low   )
//
  , m_ggamma   ( 2 , 1 , 1 , 0 ) 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::GenGammaDist::GenGammaDist
( const Ostap::Models::GenGammaDist& right ,
  const char*                           name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"     , this , right.m_x      ) 
  , m_k      ( "k"     , this , right.m_k      ) 
  , m_theta  ( "theta" , this , right.m_theta  )
  , m_p      ( "p"     , this , right.m_p      )
  , m_low    ( "low"   , this , right.m_low    )
//
  , m_ggamma (                  right.m_ggamma ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::GenGammaDist::~GenGammaDist () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::GenGammaDist*
Ostap::Models::GenGammaDist::clone( const char* name ) const 
{ return new Ostap::Models::GenGammaDist ( *this , name) ; }
// ============================================================================
void Ostap::Models::GenGammaDist::setPars () const 
{
  //
  m_ggamma.setK     ( m_k     ) ;
  m_ggamma.setTheta ( m_theta ) ;
  m_ggamma.setP     ( m_p     ) ;
  m_ggamma.setLow   ( m_low   ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::GenGammaDist::evaluate() const 
{
  //
  setPars () ;
  //
  return m_ggamma   ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::GenGammaDist::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::GenGammaDist::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_ggamma.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
//    Amoroso
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Amoroso::Amoroso
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          theta     , 
  RooAbsReal&          alpha     , 
  RooAbsReal&          beta      ,
  RooAbsReal&          a         ) 
  : RooAbsPdf ( name ,title ) 
//
  , m_x       ( "x"     , "Observable" , this , x     ) 
  , m_theta   ( "theta" , "theta"      , this , theta )
  , m_alpha   ( "alpha" , "alpha"      , this , alpha )
  , m_beta    ( "beta"  , "beta"       , this , beta  )
  , m_a       ( "a"     , "a"          , this , a     )
//
  , m_amoroso   ( 1 , 1 , 1 , 0 ) 
{
  setPars ();
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Amoroso::Amoroso
( const Ostap::Models::Amoroso& right , 
  const char*                      name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x       ( "x"     , this , right.m_x     ) 
  , m_theta   ( "theta" , this , right.m_theta )
  , m_alpha   ( "alpha" , this , right.m_alpha )
  , m_beta    ( "beta"  , this , right.m_beta  )
  , m_a       ( "a"     , this , right.m_a     )
//
  , m_amoroso ( right.m_amoroso ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Amoroso::~Amoroso () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Amoroso*
Ostap::Models::Amoroso::clone( const char* name ) const 
{ return new Ostap::Models::Amoroso ( *this , name) ; }
// ============================================================================
void Ostap::Models::Amoroso::setPars () const 
{
  //
  m_amoroso.setTheta     ( m_theta ) ;
  m_amoroso.setAlpha     ( m_alpha ) ;
  m_amoroso.setBeta      ( m_beta  ) ;
  m_amoroso.setA         ( m_a     ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Amoroso::evaluate() const 
{
  //
  setPars () ;
  //
  return m_amoroso   ( m_x     ) ;
}
// ============================================================================
Int_t Ostap::Models::Amoroso::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Amoroso::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_amoroso.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================





// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::LogGammaDist::LogGammaDist
( const char*          name  , 
  const char*          title ,
  RooAbsReal&          x     ,
  RooAbsReal&          k     ,
  RooAbsReal&          theta )
  : RooAbsPdf  (name ,title ) 
//
  , m_x       ( "x"     , "Observable" , this , x     ) 
  , m_k       ( "k"     , "Shape"      , this , k     ) 
  , m_theta   ( "theta" , "Scale"      , this , theta )
//
  , m_gamma   ( 1 , 1 ) 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::LogGammaDist::LogGammaDist
( const Ostap::Models::LogGammaDist& right ,
  const char*                           name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x     ( "x"     , this , right.m_x     ) 
  , m_k     ( "k"     , this , right.m_k     ) 
  , m_theta ( "theta" , this , right.m_theta )
//
  , m_gamma (                  right.m_gamma ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::LogGammaDist::~LogGammaDist () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::LogGammaDist*
Ostap::Models::LogGammaDist::clone( const char* name ) const 
{ return new Ostap::Models::LogGammaDist ( *this , name) ; }
// ============================================================================
void Ostap::Models::LogGammaDist::setPars () const 
{
  //
  m_gamma.setK     ( m_k     ) ;
  m_gamma.setTheta ( m_theta ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::LogGammaDist::evaluate() const 
{
  //
  setPars () ;
  //
  return m_gamma   ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::LogGammaDist::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::LogGammaDist::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_gamma.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Log10GammaDist::Log10GammaDist
( const char*          name  , 
  const char*          title ,
  RooAbsReal&          x     ,
  RooAbsReal&          k     ,
  RooAbsReal&          theta )
  : RooAbsPdf  (name ,title ) 
//
  , m_x       ( "x"     , "Observable" , this , x     ) 
  , m_k       ( "k"     , "Shape"      , this , k     ) 
  , m_theta   ( "theta" , "Scale"      , this , theta )
//
  , m_gamma   ( 1 , 1 ) 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Log10GammaDist::Log10GammaDist
( const Ostap::Models::Log10GammaDist& right ,
  const char*                           name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x     ( "x"     , this , right.m_x     ) 
  , m_k     ( "k"     , this , right.m_k     ) 
  , m_theta ( "theta" , this , right.m_theta )
//
  , m_gamma (                  right.m_gamma ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Log10GammaDist::~Log10GammaDist () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Log10GammaDist*
Ostap::Models::Log10GammaDist::clone( const char* name ) const 
{ return new Ostap::Models::Log10GammaDist ( *this , name) ; }
// ============================================================================
void Ostap::Models::Log10GammaDist::setPars () const 
{
  //
  m_gamma.setK     ( m_k     ) ;
  m_gamma.setTheta ( m_theta ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Log10GammaDist::evaluate() const 
{
  //
  setPars () ;
  //
  return m_gamma   ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Log10GammaDist::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Log10GammaDist::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_gamma.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::LogGamma::LogGamma
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          nu     ,
  RooAbsReal&          lambda ,
  RooAbsReal&          alpha  )
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"      , "Observable" , this , x      ) 
  , m_nu      ( "nu"     , "nu"         , this , nu     ) 
  , m_lambda  ( "lambda" , "lambda"     , this , lambda ) 
  , m_alpha   ( "alpha"  , "alpha"      , this , alpha  )
//
  , m_lgamma   ( 0 , 1 , 1 ) 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::LogGamma::LogGamma
( const Ostap::Models::LogGamma& right ,
  const char*                       name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_nu     ( "nu"     , this , right.m_nu     ) 
  , m_lambda ( "lambda" , this , right.m_lambda )
  , m_alpha  ( "alpha"  , this , right.m_alpha  )
//
  , m_lgamma (                   right.m_lgamma ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::LogGamma::~LogGamma () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::LogGamma*
Ostap::Models::LogGamma::clone( const char* name ) const 
{ return new Ostap::Models::LogGamma ( *this , name) ; }
// ============================================================================
void Ostap::Models::LogGamma::setPars () const 
{
  //
  m_lgamma.setNu     ( m_nu     ) ;
  m_lgamma.setLambda ( m_lambda ) ;
  m_lgamma.setAlpha  ( m_alpha  ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::LogGamma::evaluate() const 
{
  //
  setPars () ;
  //
  return m_lgamma    ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::LogGamma::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::LogGamma::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  //
  return m_lgamma.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BetaPrime::BetaPrime
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          alpha  ,
  RooAbsReal&          beta   ,
  RooAbsReal&          scale  ,
  RooAbsReal&          shift  )
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"      , "Observable" , this , x      ) 
  , m_alpha   ( "alpha"  , "alpha"      , this , alpha  ) 
  , m_beta    ( "beta"   , "beta"       , this , beta   ) 
  , m_scale   ( "scale"  , "scale"      , this , scale  ) 
  , m_shift   ( "shift"  , "shift"      , this , shift  ) 
    //
  , m_betap   ( 3 , 3 , 1 , 0 ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::BetaPrime::BetaPrime
( const Ostap::Models::BetaPrime& right ,
  const char*                       name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_alpha  ( "alpha"  , this , right.m_alpha  )
  , m_beta   ( "beta"   , this , right.m_beta   )
  , m_scale  ( "scale"  , this , right.m_scale  )
  , m_shift  ( "shift"  , this , right.m_shift  )
//
  , m_betap  (                   right.m_betap  ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::BetaPrime::~BetaPrime () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::BetaPrime*
Ostap::Models::BetaPrime::clone( const char* name ) const 
{ return new Ostap::Models::BetaPrime ( *this , name) ; }
// ============================================================================
void Ostap::Models::BetaPrime::setPars () const 
{
  //
  m_betap.setAlpha  ( m_alpha  ) ;
  m_betap.setBeta   ( m_beta   ) ;
  m_betap.setScale  ( m_scale  ) ;
  m_betap.setShift  ( m_shift  ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::BetaPrime::evaluate() const 
{
  //
  setPars () ;
  //
  return m_betap    ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::BetaPrime::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::BetaPrime::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_betap.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::SinhAsinh::SinhAsinh
( const char*          name    , 
  const char*          title   ,
  RooAbsReal&          x       ,
  RooAbsReal&          mu      ,
  RooAbsReal&          sigma   ,
  RooAbsReal&          epsilon ,
  RooAbsReal&          delta   )
  : RooAbsPdf ( name , title ) 
    //
  , m_x       ( "x"       , "Observable"    , this , x       ) 
  , m_mu      ( "mu"      , "mu/location"   , this , mu      ) 
  , m_sigma   ( "sigma"   , "sigma/scale"   , this , sigma   ) 
  , m_epsilon ( "epsilon" , "epsilon/skew"  , this , epsilon ) 
  , m_delta   ( "delta"   , "delta/tail"    , this , delta   ) 
    //
  , m_sinhasinh  ( 1 , 1 , 0 , 1 )  
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::SinhAsinh::SinhAsinh
( const Ostap::Models::SinhAsinh& right , 
  const char*                         name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x       ( "x"       , this , right.m_x       ) 
  , m_mu      ( "mu"      , this , right.m_mu      )
  , m_sigma   ( "sigma"   , this , right.m_sigma   )
  , m_epsilon ( "epsilon" , this , right.m_epsilon )
  , m_delta   ( "delta"   , this , right.m_delta   )
    //
  , m_sinhasinh ( right.m_sinhasinh ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::SinhAsinh::~SinhAsinh(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::SinhAsinh*
Ostap::Models::SinhAsinh::clone( const char* name ) const 
{ return new Ostap::Models::SinhAsinh ( *this , name) ; }
// ============================================================================
void Ostap::Models::SinhAsinh::setPars () const 
{
  //
  m_sinhasinh.setMu       ( m_mu      ) ;
  m_sinhasinh.setSigma    ( m_sigma   ) ;
  m_sinhasinh.setEpsilon  ( m_epsilon ) ;
  m_sinhasinh.setDelta    ( m_delta   ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::SinhAsinh::evaluate() const 
{
  //
  setPars () ;
  //
  return m_sinhasinh ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::SinhAsinh::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::SinhAsinh::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_sinhasinh.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::JohnsonSU::JohnsonSU
( const char*          name    , 
  const char*          title   ,
  RooAbsReal&          x       ,
  RooAbsReal&          xi      ,
  RooAbsReal&          lam     ,
  RooAbsReal&          delta   ,
  RooAbsReal&          gamma   )
  : RooAbsPdf ( name , title ) 
    //
  , m_x       ( "x"       , "Observable"    , this , x       ) 
  , m_xi      ( "xi"      , "mu/location"   , this , xi      ) 
  , m_lambda  ( "lambda"  , "lambda/scale"  , this , lam     ) 
  , m_delta   ( "delta"   , "delta/shape"   , this , delta   ) 
  , m_gamma   ( "gamma"   , "gamma/shape"   , this , gamma   ) 
    //
  , m_johnsonSU  ( 0 , 1 , 1 , 0 )  
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::JohnsonSU::JohnsonSU
( const Ostap::Models::JohnsonSU&  right , 
  const char*                         name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x       ( "x"       , this , right.m_x       ) 
  , m_xi      ( "xi"      , this , right.m_xi      )
  , m_lambda  ( "sigma"   , this , right.m_lambda  )
  , m_delta   ( "delta"   , this , right.m_delta   )
  , m_gamma   ( "gamma"   , this , right.m_gamma   )
    //
  , m_johnsonSU ( right.m_johnsonSU ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::JohnsonSU::~JohnsonSU(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::JohnsonSU*
Ostap::Models::JohnsonSU::clone( const char* name ) const 
{ return new Ostap::Models::JohnsonSU ( *this , name) ; }
// ============================================================================
void Ostap::Models::JohnsonSU::setPars () const 
{
  //
  m_johnsonSU.setXi       ( m_xi      ) ;
  m_johnsonSU.setLambda   ( m_lambda  ) ;
  m_johnsonSU.setDelta    ( m_delta   ) ;
  m_johnsonSU.setGamma    ( m_gamma   ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::JohnsonSU::evaluate() const 
{
  //
  setPars () ;
  //
  return m_johnsonSU ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::JohnsonSU::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::JohnsonSU::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_johnsonSU.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Landau::Landau
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          scale  ,
  RooAbsReal&          shift  )
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"      , "Observable" , this , x      ) 
  , m_scale   ( "scale"  , "scale"      , this , scale  ) 
  , m_shift   ( "shift"  , "shift"      , this , shift  ) 
    //
  , m_landau  ( 1 , 0 ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Landau::Landau
( const Ostap::Models::Landau& right ,
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_scale  ( "scale"  , this , right.m_scale  )
  , m_shift  ( "shift"  , this , right.m_shift  )
//
  , m_landau (                   right.m_landau ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Landau::~Landau () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Landau*
Ostap::Models::Landau::clone( const char* name ) const 
{ return new Ostap::Models::Landau ( *this , name) ; }
// ============================================================================
void Ostap::Models::Landau::setPars () const 
{
  //
  m_landau.setScale  ( m_scale  ) ;
  m_landau.setShift  ( m_shift  ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Landau::evaluate() const 
{
  //
  setPars () ;
  //
  return m_landau ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Landau::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Landau::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_landau.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Atlas::Atlas
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          mu     ,
  RooAbsReal&          sigma  )
  : RooAbsPdf ( name , title ) 
    //
  , m_x      ( "x"      , "Observable" , this , x      ) 
  , m_mu     ( "mu"     , "location"   , this , mu     ) 
  , m_sigma  ( "sigma"  , "sigma"      , this , sigma  ) 
    //
  , m_atlas  ( 0 , 1 ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Atlas::Atlas
( const Ostap::Models::Atlas&  right ,
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_mu     ( "mu"     , this , right.m_mu     )
  , m_sigma  ( "sigma"  , this , right.m_sigma  )
//
  , m_atlas  (                   right.m_atlas  ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Atlas::~Atlas () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Atlas*
Ostap::Models::Atlas::clone( const char* name ) const 
{ return new Ostap::Models::Atlas ( *this , name) ; }
// ============================================================================
void Ostap::Models::Atlas::setPars () const 
{
  //
  m_atlas.setMean  ( m_mu    ) ;
  m_atlas.setSigma ( m_sigma ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Atlas::evaluate() const 
{
  //
  setPars () ;
  //
  return m_atlas ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Atlas::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Atlas::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_atlas.integral ( m_x.min( rangeName ) , m_x.max( rangeName ) ) ;
}
// ============================================================================



// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Sech::Sech
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          mu     ,
  RooAbsReal&          sigma  )
  : RooAbsPdf ( name , title ) 
    //
  , m_x      ( "x"      , "Observable" , this , x      ) 
  , m_mu     ( "mu"     , "location"   , this , mu     ) 
  , m_sigma  ( "sigma"  , "sigma"      , this , sigma  ) 
    //
  , m_sech  ( 0 , 1 ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Sech::Sech
( const Ostap::Models::Sech&  right ,
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_mu     ( "mu"     , this , right.m_mu     )
  , m_sigma  ( "sigma"  , this , right.m_sigma  )
    //
  , m_sech  (                   right.m_sech    ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Sech::~Sech () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Sech*
Ostap::Models::Sech::clone( const char* name ) const 
{ return new Ostap::Models::Sech( *this , name) ; }
// ============================================================================
void Ostap::Models::Sech::setPars () const 
{
  //
  m_sech.setMean  ( m_mu    ) ;
  m_sech.setSigma ( m_sigma ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Sech::evaluate() const 
{
  //
  setPars () ;
  //
  return m_sech ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Sech::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Sech::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_sech.integral ( m_x.min( rangeName ) , m_x.max( rangeName ) ) ;
}
// ============================================================================



// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Logistic::Logistic
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          mu     ,
  RooAbsReal&          sigma  )
  : RooAbsPdf ( name , title ) 
    //
  , m_x       ( "x"      , "Observable" , this , x      ) 
  , m_mu      ( "mu"     , "location"   , this , mu     ) 
  , m_sigma   ( "sigma"  , "sigma"      , this , sigma  ) 
    //
  , m_logistic ( 0 , 1 ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Logistic::Logistic
( const Ostap::Models::Logistic&  right ,
  const char*                        name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_mu     ( "mu"     , this , right.m_mu     )
  , m_sigma  ( "sigma"  , this , right.m_sigma  )
    //
  , m_logistic  ( right.m_logistic ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Logistic::~Logistic () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Logistic*
Ostap::Models::Logistic::clone( const char* name ) const 
{ return new Ostap::Models::Logistic( *this , name) ; }
// ============================================================================
void Ostap::Models::Logistic::setPars () const 
{
  //
  m_logistic.setMean  ( m_mu    ) ;
  m_logistic.setSigma ( m_sigma ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Logistic::evaluate() const 
{
  //
  setPars () ;
  //
  return m_logistic ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Logistic::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Logistic::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_logistic.integral ( m_x.min( rangeName ) , m_x.max( rangeName ) ) ;
}
// ============================================================================



// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Argus::Argus
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          shape  ,
  RooAbsReal&          high   ,
  RooAbsReal&          low    )
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"      , "Observable" , this , x      ) 
  , m_shape   ( "shape"  , "shape"      , this , shape  ) 
  , m_high    ( "high"   , "high"       , this , high   ) 
  , m_low     ( "low"    , "low"        , this , low    ) 
    //
  , m_argus  ( 1 , 1 , 0 ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Argus::Argus
( const Ostap::Models::Argus&  right ,
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_shape  ( "shape"  , this , right.m_shape  )
  , m_high   ( "high"   , this , right.m_high   )
  , m_low    ( "low"    , this , right.m_low    )
//
  , m_argus  (                   right.m_argus  ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Argus::~Argus () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Argus*
Ostap::Models::Argus::clone( const char* name ) const 
{ return new Ostap::Models::Argus ( *this , name) ; }
// ============================================================================
void Ostap::Models::Argus::setPars () const 
{
  //
  m_argus.setShape  ( m_shape  ) ;
  m_argus.setLow    ( m_low    ) ;
  m_argus.setHigh   ( m_high   ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Argus::evaluate() const 
{
  //
  setPars () ;
  //
  return m_argus ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Argus::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Argus::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_argus.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Tsallis::Tsallis
( const char*           name      , 
  const char*           title     ,
  RooAbsReal&           x         , 
  RooAbsReal&           n         ,   // parameter N 
  RooAbsReal&           T         ,   // parameter T
  RooAbsReal&           mass      )   // particle mass (fixed)
  : RooAbsPdf ( name , title ) 
    //
  , m_x       ( "x"      , "Observable"  , this , x      ) 
  , m_n       ( "n"      , "shape"       , this , n      ) 
  , m_T       ( "T"      , "temperature" , this , T      ) 
  , m_mass    ( "m"      , "mass"        , this , mass   ) 
    //
  , m_tsallis  ( 0 , 10 , 1 ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Tsallis::Tsallis
( const Ostap::Models::Tsallis& right ,
  const char*                      name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x       ( "x"  , this , right.m_x       ) 
  , m_n       ( "n"  , this , right.m_n       )
  , m_T       ( "T"  , this , right.m_T       )
  , m_mass    ( "m"  , this , right.m_mass    )
    //
  , m_tsallis (               right.m_tsallis ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Tsallis::~Tsallis () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Tsallis*
Ostap::Models::Tsallis::clone( const char* name ) const 
{ return new Ostap::Models::Tsallis ( *this , name) ; }
// ============================================================================
void Ostap::Models::Tsallis::setPars () const 
{
  //
  m_tsallis.setMass  ( m_mass ) ;
  m_tsallis.setN     ( m_n    ) ;
  m_tsallis.setT     ( m_T    ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Tsallis::evaluate() const 
{
  //
  setPars () ;
  //
  return m_tsallis ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Tsallis::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Tsallis::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_tsallis.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::QGSM::QGSM
( const char*           name      , 
  const char*           title     ,
  RooAbsReal&           x         , 
  RooAbsReal&           b         ,   // parameter b 
  RooAbsReal&           mass      )   // particle mass (fixed)
  : RooAbsPdf ( name , title ) 
    //
  , m_x       ( "x"      , "Observable"  , this , x      ) 
  , m_b       ( "b"      , "slope"       , this , b      ) 
  , m_mass    ( "m"      , "mass"        , this , mass   ) 
    //
  , m_qgsm  ( 0 , 1 ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::QGSM::QGSM
( const Ostap::Models::QGSM&    right ,
  const char*                      name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x    ( "x"  , this , right.m_x    ) 
  , m_b    ( "b"  , this , right.m_b    )
  , m_mass ( "m"  , this , right.m_mass )
    //
  , m_qgsm (               right.m_qgsm ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::QGSM::~QGSM () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::QGSM*
Ostap::Models::QGSM::clone( const char* name ) const 
{ return new Ostap::Models::QGSM ( *this , name) ; }
// ============================================================================
void Ostap::Models::QGSM::setPars () const 
{
  //
  m_qgsm.setMass  ( m_mass ) ;
  m_qgsm.setB     ( m_b    ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::QGSM::evaluate() const 
{
  //
  setPars () ;
  //
  return m_qgsm ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::QGSM::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::QGSM::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_qgsm.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================




// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::TwoExpos::TwoExpos 
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          alpha  ,
  RooAbsReal&          delta  ,
  RooAbsReal&          x0     )
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"      , "Observable" , this , x     ) 
  , m_alpha   ( "alpha"  , "alpha"      , this , alpha ) 
  , m_delta   ( "delta"  , "delta"      , this , delta ) 
  , m_x0      ( "x0"     , "x0"         , this , x0    ) 
    //
  , m_2expos  ( 1 , 1 , 0 )
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::TwoExpos::TwoExpos 
( const Ostap::Models::TwoExpos&  right ,
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x      ( "x"      , this , right.m_x     ) 
  , m_alpha  ( "alpha"  , this , right.m_alpha ) 
  , m_delta  ( "delta"  , this , right.m_delta ) 
  , m_x0     ( "x0"      , this , right.m_x0    ) 
    //
  , m_2expos (                   right.m_2expos ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::TwoExpos::~TwoExpos (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::TwoExpos*
Ostap::Models::TwoExpos::clone( const char* name ) const 
{ return new Ostap::Models::TwoExpos ( *this , name) ; }
// ============================================================================
void Ostap::Models::TwoExpos::setPars () const 
{
  //
  m_2expos.setAlpha  ( m_alpha  ) ;
  m_2expos.setDelta  ( m_delta  ) ;
  m_2expos.setX0     ( m_x0     ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::TwoExpos::evaluate() const 
{
  //
  setPars () ;
  //
  return m_2expos( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::TwoExpos::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::TwoExpos::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  //
  return m_2expos.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================



// ============================================================================
// generic polinomial
// ============================================================================
Ostap::Models::Poly2DPositive::Poly2DPositive 
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const unsigned short nX        , 
  const unsigned short nY        ,
  RooArgList&          phis      ) 
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
//
  , m_positive ( nX , nY , x.getMin() , x.getMax() , y.getMin() , y.getMax() ) 
{
  //
  RooAbsArg*   coef = 0 ;
  unsigned int num  = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_positive.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
    ++num ;  
  }
  //
  if ( num != m_positive.npars() ) 
  { Ostap::throwException ( "Invalid size of parameters vector", 
                            "Ostap::Poly2DPositive"            ) ; }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Poly2DPositive::Poly2DPositive
( const Ostap::Models::Poly2DPositive&  right ,      
  const char*                              name  ) 
  : RooAbsPdf  ( right , name ) 
//
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_y        ( "y"      , this , right.m_y     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
//
  , m_positive ( right.m_positive ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Poly2DPositive::~Poly2DPositive(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Poly2DPositive*
Ostap::Models::Poly2DPositive::clone( const char* name ) const 
{ return new Ostap::Models::Poly2DPositive(*this,name) ; }
// ============================================================================
void Ostap::Models::Poly2DPositive::setPars () const 
{
  //
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phiv   = r->getVal ( nset ) ;
    //
    m_positive.setPar ( k  , phiv ) ;
    //
    ++k ;
  }
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Poly2DPositive::evaluate() const 
{
  //
  setPars () ;
  //
  return m_positive ( m_x , m_y ) ; 
}
// ============================================================================
Int_t Ostap::Models::Poly2DPositive::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if      ( matchArgs ( allVars , analVars , m_x , m_y ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars       , m_y ) ) { return 3 ; }
  //
  return 0 ;
}
//_____________________________________________________________________________
Double_t Ostap::Models::Poly2DPositive::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( 1 == code || 2 ==code || 3 == code ) ;
  //
  setPars () ;
  //
  return 
    1 == code ? m_positive.integral   (        m_x.min(rangeName) , m_x.max(rangeName) ,
                                               m_y.min(rangeName) , m_y.max(rangeName) ) : 
    2 == code ? m_positive.integrateX ( m_y  , m_x.min(rangeName) , m_x.max(rangeName) ) : 
    3 == code ? m_positive.integrateY ( m_x  , m_y.min(rangeName) , m_y.max(rangeName) ) : 0.0 ;  
}
// ============================================================================
// generic polinomial
// ============================================================================
Ostap::Models::Poly2DSymPositive::Poly2DSymPositive 
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const unsigned short n         ,
  RooArgList&          phis      ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
//
  , m_positive ( n , x.getMin() , x.getMax() ) 
{
  //
  if (  x.getMin() != y.getMin() || x.getMax() != y.getMax() )
  { Ostap::throwException ( "Non-equal x/y-edges "     , 
                            "Ostap::Poly2DSymPositive" ) ; }
  RooAbsArg* coef = 0 ;
  unsigned num = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_positive.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
    ++num ;  
  }
  //
  if ( num != m_positive.npars() ) 
  { Ostap::throwException ( "Invalid size of parameters vector" , 
                            "Ostap::Poly2DSymPositive"          ) ; }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Poly2DSymPositive::Poly2DSymPositive
( const Ostap::Models::Poly2DSymPositive&  right ,      
  const char*                                 name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_y        ( "y"      , this , right.m_y     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
//
  , m_positive ( right.m_positive ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Poly2DSymPositive::~Poly2DSymPositive() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Poly2DSymPositive*
Ostap::Models::Poly2DSymPositive::clone( const char* name ) const 
{ return new Ostap::Models::Poly2DSymPositive(*this,name) ; }
// ============================================================================
void Ostap::Models::Poly2DSymPositive::setPars () const 
{
  //
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phiv   = r->getVal ( nset ) ;
    //
    m_positive.setPar ( k  , phiv ) ;
    //
    ++k ;
  }
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Poly2DSymPositive::evaluate() const 
{
  //
  setPars () ;
  //
  return m_positive ( m_x , m_y ) ; 
}
// ============================================================================
Int_t Ostap::Models::Poly2DSymPositive::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if      ( matchArgs ( allVars , analVars , m_x , m_y ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars       , m_y ) ) { return 3 ; }
  //
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Poly2DSymPositive::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( 1 == code || 2 == code || 3 == code ) ;
  //
  setPars () ;
  //
  return 
    1 == code ? m_positive.integral   (        m_x.min(rangeName) , m_x.max(rangeName) , 
                                               m_y.min(rangeName) , m_y.max(rangeName) ) : 
    2 == code ? m_positive.integrateX ( m_y  , m_x.min(rangeName) , m_x.max(rangeName) ) : 
    3 == code ? m_positive.integrateY ( m_x  , m_y.min(rangeName) , m_y.max(rangeName) ) : 0.0 ;  
}

// ============================================================================
//  PS(x)*PS(y)*Polynom 
// ============================================================================
Ostap::Models::PS2DPol::PS2DPol
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const Ostap::Math::PhaseSpaceNL& psx , 
  const Ostap::Math::PhaseSpaceNL& psy ,
  const unsigned short nX        ,
  const unsigned short nY        ,
  RooArgList&          phis      ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_function ( psx, psy , nX , nY , 
                 x.getMin() , x.getMax() ,
                 y.getMin() , y.getMax() )
{
  //  
  RooAbsArg* coef = 0 ;
  unsigned num = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_function.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
    ++num ;  
  }
  //
  if ( num != m_function.npars() ) 
  { Ostap::throwException ( "Invalid size of parameters vector", 
                            "Ostap::PS22DPol"                  ) ; }
  //
}
// ============================================================================
// generic polinomial
// ============================================================================
Ostap::Models::PS2DPol::PS2DPol
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const Ostap::Math::PS2DPol& ps , 
  RooArgList&          phis      ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_function ( ps ) 
{
  //
  RooAbsArg* coef = 0 ;
  unsigned num = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_function.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
    ++num ;  
  }
  //
  if ( num != m_function.npars() ) 
  { Ostap::throwException ( "Invalid size of parameters vector", 
                            "Ostap::PS22DPol"                  ) ; }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PS2DPol::PS2DPol
( const Ostap::Models::PS2DPol&  right ,      
  const char*                       name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_y        ( "y"      , this , right.m_y     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_function ( right.m_function )
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PS2DPol::~PS2DPol() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PS2DPol*
Ostap::Models::PS2DPol::clone( const char* name ) const 
{ return new Ostap::Models::PS2DPol ( *this , name ) ; }
// ============================================================================
void Ostap::Models::PS2DPol::setPars () const 
{
  //
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phiv   = r->getVal ( nset ) ;
    //
    m_function.setPar ( k  , phiv ) ;
    //
    ++k ;
  }
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PS2DPol::evaluate() const 
{
  //
  setPars () ;
  //
  return m_function ( m_x , m_y ) ;
}
// ============================================================================
Int_t Ostap::Models::PS2DPol::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if      ( matchArgs ( allVars , analVars , m_x , m_y ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars       , m_y ) ) { return 3 ; }
  //
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PS2DPol::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( 1 == code || 2 == code || 3 == code ) ;
  //
  setPars () ;
  //
  return 
    1 == code ? m_function.integral   (        m_x.min(rangeName) , m_x.max(rangeName) , 
                                               m_y.min(rangeName) , m_y.max(rangeName) ) : 
    2 == code ? m_function.integrateX ( m_y  , m_x.min(rangeName) , m_x.max(rangeName) ) : 
    3 == code ? m_function.integrateY ( m_x  , m_y.min(rangeName) , m_y.max(rangeName) ) : 0.0 ;  
}
// ============================================================================
//  PS(x)*PS(y)*SymPolynom 
// ============================================================================
Ostap::Models::PS2DPolSym::PS2DPolSym
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const Ostap::Math::PhaseSpaceNL& ps ,
  const unsigned short n         ,
  RooArgList&          phis      ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_function ( ps , n , x.getMin() , x.getMax() ) 
{
  //
  RooAbsArg* coef = 0 ;
  unsigned   num  = 0 ;  
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_function.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
    ++num ;  
  }
  //
  if ( num != m_function.npars() ) 
  { Ostap::throwException ( "Invalid size of parameters vector", 
                            "Ostap::PS22DPol"                  ) ; }
  //
}
// ============================================================================
// generic polinomial
// ============================================================================
Ostap::Models::PS2DPolSym::PS2DPolSym
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const Ostap::Math::PS2DPolSym& ps , 
  RooArgList&          phis      ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_function ( ps ) 
{
  //
  RooAbsArg* coef = 0 ;
  unsigned num = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_function.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
    ++num ;  
  }
  //
  if ( num != m_function.npars() ) 
  { Ostap::throwException ( "Invalid size of parameters vector" , 
                            "Ostap::PS22DPol"                   ) ; }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PS2DPolSym::PS2DPolSym
( const Ostap::Models::PS2DPolSym& right ,      
  const char*                         name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_y        ( "y"      , this , right.m_y     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_function ( right.m_function )
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PS2DPolSym::~PS2DPolSym() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PS2DPolSym*
Ostap::Models::PS2DPolSym::clone( const char* name ) const 
{ return new Ostap::Models::PS2DPolSym ( *this , name ) ; }
// ============================================================================
void Ostap::Models::PS2DPolSym::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phiv   = r->getVal ( nset ) ;
    //
    m_function.setPar ( k  , phiv ) ;
    //
    ++k ;
  }
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PS2DPolSym::evaluate() const 
{
  //
  setPars () ;
  //
  return m_function ( m_x , m_y ) ;
}
// ============================================================================
Int_t Ostap::Models::PS2DPolSym::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if      ( matchArgs ( allVars , analVars , m_x , m_y ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars       , m_y ) ) { return 3 ; }
  //
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PS2DPolSym::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( 1 == code || 2 == code || 3 == code ) ;
  //
  setPars () ;
  //
  return 
    1 == code ? m_function.integral   (        m_x.min(rangeName) , m_x.max(rangeName) , 
                                               m_y.min(rangeName) , m_y.max(rangeName) ) : 
    2 == code ? m_function.integrateX ( m_y  , m_x.min(rangeName) , m_x.max(rangeName) ) : 
    3 == code ? m_function.integrateY ( m_x  , m_y.min(rangeName) , m_y.max(rangeName) ) : 0.0 ;  
}


// ============================================================================
//  exp(x)*PS(y)*Polynom 
// ============================================================================
Ostap::Models::ExpoPS2DPol::ExpoPS2DPol
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  RooAbsReal&          tau       ,
  const Ostap::Math::PhaseSpaceNL& psy ,
  const unsigned short nX        ,
  const unsigned short nY        ,
  RooArgList&          phis      ) 
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"       , "Observable-X"      , this , x   ) 
  , m_y        ( "y"       , "Observable-Y"      , this , y   ) 
  , m_tau      ( "tau"     , "Exponential slope" , this , tau ) 
  , m_phis     ( "phis"    , "Coefficients"      , this       )
    //
  , m_function ( psy        ,
                 x.getMin() , x.getMax() ,
                 nX         , nY         , 
                 y.getMin() , y.getMax() )
{
  //
  RooAbsArg* coef = 0 ;
  unsigned   num = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_function.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
    ++num ;  
  }
  //
  if ( num != m_function.npars() ) 
  { Ostap::throwException ( "Invalid size of parameters vector" , 
                            "Ostap::ExpoPS22DPol"               ) ; }
}
// ============================================================================
// generic polinomial
// ============================================================================
Ostap::Models::ExpoPS2DPol::ExpoPS2DPol
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  RooAbsReal&          tau       ,
  const Ostap::Math::ExpoPS2DPol& ps , 
  RooArgList&          phis      ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X"      , this , x   ) 
  , m_y        ( "y"       , "Observable-Y"      , this , y   ) 
  , m_tau      ( "tau"     , "Exponential slope" , this , tau ) 
  , m_phis     ( "phis"    , "Coefficients"      , this       )
    //
  , m_function ( ps ) 
{
  //  
  RooAbsArg* coef = 0 ;
  unsigned num = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_function.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
    ++num ;  
  }
  //
  if ( num != m_function.npars() ) 
  { Ostap::throwException ( "Invalid size of parameters vector" ,  
                            "Ostap::ExpoPS22DPol"               ) ; }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::ExpoPS2DPol::ExpoPS2DPol
( const Ostap::Models::ExpoPS2DPol&  right ,      
  const char*                       name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_y        ( "y"      , this , right.m_y     ) 
  , m_tau      ( "tau"    , this , right.m_tau   ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_function ( right.m_function )
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::ExpoPS2DPol::~ExpoPS2DPol() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::ExpoPS2DPol*
Ostap::Models::ExpoPS2DPol::clone( const char* name ) const 
{ return new Ostap::Models::ExpoPS2DPol ( *this , name ) ; }
// ============================================================================
void Ostap::Models::ExpoPS2DPol::setPars () const 
{
  //
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  unsigned short   k     = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phiv   = r->getVal ( nset ) ;
    //
    m_function.setPar ( k  , phiv ) ;
    //
    ++k ;
  }
  //
  m_function.setTau ( m_tau ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::ExpoPS2DPol::evaluate() const 
{
  //
  setPars () ;
  //
  return m_function ( m_x , m_y ) ;
}
// ============================================================================
Int_t Ostap::Models::ExpoPS2DPol::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if      ( matchArgs ( allVars , analVars , m_x , m_y ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars       , m_y ) ) { return 3 ; }
  //
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::ExpoPS2DPol::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( 1 == code || 2 == code || 3 == code ) ;
  //
  setPars () ;
  //
  return 
    1 == code ? m_function.integral   (        m_x.min(rangeName) , m_x.max(rangeName) , 
                                               m_y.min(rangeName) , m_y.max(rangeName) ) : 
    2 == code ? m_function.integrateX ( m_y  , m_x.min(rangeName) , m_x.max(rangeName) ) : 
    3 == code ? m_function.integrateY ( m_x  , m_y.min(rangeName) , m_y.max(rangeName) ) : 0.0 ;  
}



// ============================================================================
//  exp(x)*exp(y)*Polynom 
// ============================================================================
Ostap::Models::Expo2DPol::Expo2DPol
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  RooAbsReal&          taux      ,
  RooAbsReal&          tauy      ,
  const unsigned short nX        ,
  const unsigned short nY        ,
  RooArgList&          phis      ) 
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"       , "Observable-X"      , this , x    ) 
  , m_y        ( "y"       , "Observable-Y"      , this , y    ) 
  , m_taux     ( "taux"    , "Exponential slope" , this , taux ) 
  , m_tauy     ( "tauy"    , "Exponential slope" , this , tauy ) 
  , m_phis     ( "phis"    , "Coefficients"      , this        )
    //
  , m_function ( x.getMin() , x.getMax() ,
                 y.getMin() , y.getMax() ,
                 nX         , nY         )
{
  //
  RooAbsArg* coef = 0 ;
  unsigned num = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_function.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
    ++num ;  
  }
  //
  if ( num != m_function.npars() ) 
  { Ostap::throwException ( "Invalid size of parameters vector" , 
                            "Ostap::ExpoPS22DPol"               ) ; }
  //
  m_function.setTauX ( m_taux ) ;
  m_function.setTauY ( m_tauy ) ;
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Expo2DPol::Expo2DPol
( const Ostap::Models::Expo2DPol&  right ,      
  const char*                       name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_y        ( "y"      , this , right.m_y     ) 
  , m_taux     ( "taux"   , this , right.m_taux  ) 
  , m_tauy     ( "tauy"   , this , right.m_tauy  ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_function ( right.m_function )
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Expo2DPol::~Expo2DPol() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Expo2DPol*
Ostap::Models::Expo2DPol::clone( const char* name ) const 
{ return new Ostap::Models::Expo2DPol ( *this , name ) ; }
// ============================================================================
void Ostap::Models::Expo2DPol::setPars () const 
{
  //
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phiv   = r->getVal ( nset ) ;
    //
    m_function.setPar ( k  , phiv ) ;
    //
    ++k ;
  }
  //
  m_function.setTauX ( m_taux ) ;
  m_function.setTauY ( m_tauy ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Expo2DPol::evaluate() const 
{
  //
  setPars () ;
  //
  return m_function ( m_x , m_y ) ;
}
// ============================================================================
Int_t Ostap::Models::Expo2DPol::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if      ( matchArgs ( allVars , analVars , m_x , m_y ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars       , m_y ) ) { return 3 ; }
  //
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Expo2DPol::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( 1 == code || 2 == code || 3 == code ) ;
  //
  setPars () ;
  //
  return 
    1 == code ? m_function.integral   (        m_x.min(rangeName) , m_x.max(rangeName) , 
                                               m_y.min(rangeName) , m_y.max(rangeName) ) : 
    2 == code ? m_function.integrateX ( m_y  , m_x.min(rangeName) , m_x.max(rangeName) ) : 
    3 == code ? m_function.integrateY ( m_x  , m_y.min(rangeName) , m_y.max(rangeName) ) : 0.0 ;  
}

// ============================================================================
//  exp(x)*exp(y)*SymPolynom 
// ============================================================================
Ostap::Models::Expo2DPolSym::Expo2DPolSym
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  RooAbsReal&          tau       ,
  const unsigned short n         ,
  RooArgList&          phis      ) 
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"       , "Observable-X"      , this , x    ) 
  , m_y        ( "y"       , "Observable-Y"      , this , y    ) 
  , m_tau      ( "tau"     , "Exponential slope" , this , tau  ) 
  , m_phis     ( "phis"    , "Coefficients"      , this        )
    //
  , m_function ( x.getMin() , x.getMax() , n )
{
  //
  RooAbsArg* coef = 0 ;
  unsigned num = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_function.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
    ++num ;  
  }
  //
  if ( num != m_function.npars() ) 
  { Ostap::throwException ( "Invalid size of parameters vector" , 
                            "Ostap::ExpoPS22DPol"               ) ; }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Expo2DPolSym::Expo2DPolSym
( const Ostap::Models::Expo2DPolSym& right ,      
  const char*                           name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_y        ( "y"      , this , right.m_y     ) 
  , m_tau      ( "tau"    , this , right.m_tau   ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_function ( right.m_function )
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Expo2DPolSym::~Expo2DPolSym() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Expo2DPolSym*
Ostap::Models::Expo2DPolSym::clone( const char* name ) const 
{ return new Ostap::Models::Expo2DPolSym ( *this , name ) ; }
// ============================================================================
void Ostap::Models::Expo2DPolSym::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phiv   = r->getVal ( nset ) ;
    //
    m_function.setPar ( k  , phiv ) ;
    //
    ++k ;
  }
  //
  m_function.setTau ( m_tau ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Expo2DPolSym::evaluate() const 
{
  //
  setPars () ;
  //
  return m_function ( m_x , m_y ) ;
}
// ============================================================================
Int_t Ostap::Models::Expo2DPolSym::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if      ( matchArgs ( allVars , analVars , m_x , m_y ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars       , m_y ) ) { return 3 ; }
  //
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Expo2DPolSym::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( 1 == code || 2 == code || 3 == code ) ;
  //
  setPars () ;
  //
  return 
    1 == code ? m_function.integral   (        m_x.min(rangeName) , m_x.max(rangeName) , 
                                               m_y.min(rangeName) , m_y.max(rangeName) ) : 
    2 == code ? m_function.integrateX ( m_y  , m_x.min(rangeName) , m_x.max(rangeName) ) : 
    3 == code ? m_function.integrateY ( m_x  , m_y.min(rangeName) , m_y.max(rangeName) ) : 0.0 ;  
}
// ============================================================================







// ============================================================================
// generic 2D-spline
// ============================================================================
Ostap::Models::Spline2D::Spline2D
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const Ostap::Math::Spline2D& spline, 
  RooArgList&          phis      ) 
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_spline   ( spline ) 
{
  //
  RooAbsArg*   coef = 0 ;
  unsigned int num  = 0 ;
  Ostap::Utils::Iterator   tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_spline.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
    ++num ;  
  }
  //
  if ( num != m_spline.npars() ) 
  { Ostap::throwException ( "Invalid size of parameters vector" , 
                            "Ostap::Spline2D"                   ) ; }  
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Spline2D::Spline2D
( const Ostap::Models::Spline2D&  right ,      
  const char*                              name  ) 
  : RooAbsPdf  ( right , name ) 
//
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_y        ( "y"      , this , right.m_y     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
//
  , m_spline   ( right.m_spline ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Spline2D::~Spline2D(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Spline2D*
Ostap::Models::Spline2D::clone( const char* name ) const 
{ return new Ostap::Models::Spline2D(*this,name) ; }
// ============================================================================
void Ostap::Models::Spline2D::setPars () const 
{
  //
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator   it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phiv   = r->getVal ( nset ) ;
    //
    m_spline.setPar ( k  , phiv ) ;
    //
    ++k ;
  }
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Spline2D::evaluate() const 
{
  //
  setPars () ;
  //
  return m_spline ( m_x , m_y ) ; 
}
// ============================================================================
Int_t Ostap::Models::Spline2D::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if      ( matchArgs ( allVars , analVars , m_x , m_y ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars       , m_y ) ) { return 3 ; }
  //
  return 0 ;
}
//_____________________________________________________________________________
Double_t Ostap::Models::Spline2D::analyticalIntegral 
( Int_t       code       , 
  const char* rangeName  ) const 
{
  //
  assert ( 1 == code || 2 == code || 3 == code ) ;
  //
  setPars () ;
  //
  return 
    1 == code ? m_spline.integral   (        m_x.min(rangeName) , m_x.max(rangeName) , 
                                             m_y.min(rangeName) , m_y.max(rangeName) ) : 
    2 == code ? m_spline.integrateX ( m_y  , m_x.min(rangeName) , m_x.max(rangeName) ) : 
    3 == code ? m_spline.integrateY ( m_x  , m_y.min(rangeName) , m_y.max(rangeName) ) : 0.0 ;  
}

// ============================================================================
// generic 2D symmetric spline
// ============================================================================
Ostap::Models::Spline2DSym::Spline2DSym
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const Ostap::Math::Spline2DSym& spline, 
  RooArgList&          phis      ) 
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_spline   ( spline ) 
{
  //
  RooAbsArg*   coef = 0 ;
  unsigned int num  = 0 ;
  Ostap::Utils::Iterator tmp ( phis ) ;
  while ( ( coef = (RooAbsArg*) tmp.next() ) && num < m_spline.npars() )
  {
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( coef ) ;
    if ( 0 == r ) { continue ; }
    m_phis.add ( *coef ) ;
    ++num ;  
  }
  //
  if ( num != m_spline.npars() ) 
  { Ostap::throwException ( "Invalid size of parameters vector" , 
                            "Ostap::Spline2DSym"                ) ; }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Spline2DSym::Spline2DSym
( const Ostap::Models::Spline2DSym&  right ,      
  const char*                              name  ) 
  : RooAbsPdf  ( right , name ) 
//
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_y        ( "y"      , this , right.m_y     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
//
  , m_spline   ( right.m_spline ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Spline2DSym::~Spline2DSym(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Spline2DSym*
Ostap::Models::Spline2DSym::clone( const char* name ) const 
{ return new Ostap::Models::Spline2DSym(*this,name) ; }
// ============================================================================
void Ostap::Models::Spline2DSym::setPars () const 
{
  //
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_phis ) ;
  while ( ( phi = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( phi ) ;
    if ( 0 == r ) { continue ; }
    //
    const double phiv   = r->getVal ( nset ) ;
    //
    m_spline.setPar ( k  , phiv ) ;
    //
    ++k ;
  }
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Spline2DSym::evaluate() const 
{
  //
  setPars () ;
  //
  return m_spline ( m_x , m_y ) ; 
}
// ============================================================================
Int_t Ostap::Models::Spline2DSym::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if      ( matchArgs ( allVars , analVars , m_x , m_y ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars       , m_y ) ) { return 3 ; }
  //
  return 0 ;
}
//_____________________________________________________________________________
Double_t Ostap::Models::Spline2DSym::analyticalIntegral 
( Int_t       code       , 
  const char* rangeName  ) const 
{
  //
  assert ( 1 == code || 2 == code || 3 == code ) ;
  //
  setPars () ;
  //
  return 
    1 == code ? m_spline.integral   (        m_x.min(rangeName) , m_x.max(rangeName) , 
                                             m_y.min(rangeName) , m_y.max(rangeName) ) : 
    2 == code ? m_spline.integrateX ( m_y  , m_x.min(rangeName) , m_x.max(rangeName) ) : 
    3 == code ? m_spline.integrateY ( m_x  , m_y.min(rangeName) , m_y.max(rangeName) ) : 0.0 ;  
}

// ============================================================================
// The END 
// ============================================================================
