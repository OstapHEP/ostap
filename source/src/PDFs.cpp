// ============================================================================
// Include files 
// ============================================================================
// STD & ST:
// ============================================================================
#include <limits>
// ============================================================================
// Local
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
#include "local_roofit.h"
// ============================================================================
/** @file 
 *  Implementation file for namespace Ostap::Models
 *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
 *  @date   2011-11-30
 */
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
#include "BatchHelpers.h"
// ============================================================================
typedef BatchHelpers::BracketAdapter<double> BA ;
// ============================================================================
namespace 
{
  // ==========================================================================
  template<class TX, class FUN>
  void compute_X ( RooSpan<double> output , FUN& fun , TX x ) 
  {
    const int n = output.size();
    for ( int i = 0 ; i < n ; ++i ) { output [ i ] = fun ( x [ i ] ) ; }
  }
  // ==========================================================================
}
#endif
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
  , m_bw ( std::make_unique<Ostap::Math::BreitWigner>( 0 , 1 , m1 , m2 , L ) )
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
  , m_bw ( std::make_unique<Ostap::Math::BreitWigner> ( 0 , 1 , m1 , m2 , L , rho ) ) 
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
  , m_bw ( bw.clone() ) 
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
  , m_bw    (  right.m_bw->clone()        ) 
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
Ostap::Models::BreitWigner::clone ( const char* name ) const 
{ return new Ostap::Models::BreitWigner(*this,name) ; }
// ============================================================================
void Ostap::Models::BreitWigner::setPars () const 
{
  //
  m_bw->setM0     ( m_mass  ) ;
  m_bw->setGamma0 ( m_width ) ;
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
  return  ( *m_bw ) ( m_x ) ;
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
  return m_bw->integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
// get the amplitude 
// ============================================================================
std::complex<double> Ostap::Models::BreitWigner::amplitude () const
{
  //
  m_bw -> setM0    ( m_mass  ) ;
  m_bw -> setGamma ( m_width ) ;
  //
  return m_bw->amplitude ( m_x ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::BreitWigner::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , *m_bw , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================






// ============================================================================
// multi-channel BreiWigner
// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Models::BreitWignerMC::BreitWignerMC
( const char*                                 name      , 
  const char*                                 title     ,
  RooAbsReal&                                 x         ,
  RooAbsReal&                                 mass      , 
  RooArgList&                                 widths    ,
  const Ostap::Math::BreitWignerMC& bw ) 
  : RooAbsPdf ( name , title ) 
    //
  , m_x        ( "x"       , "Observable"   , this , x    ) 
  , m_mass     ( "m0"      , "PolePosition" , this , mass ) 
  , m_widths   ( "widths"  , "Widths"       , this ) 
    //
  , m_bw ( bw.clone() ) 
{
  //
  ::copy_real   ( widths , m_widths , "Invalid width parameter!" ,
                  "Ostap::Models::BreitWignerMC" ) ;
  //
  Ostap::Assert ( ::size ( m_widths ) == m_bw->nChannels() , 
                  "Widths/#channels mismatch" , 
                  "Ostap::Models::BreitWignerMC" ) ;
  //
  setPars() ;
}
// ============================================================================
// "copy" constructor
// ============================================================================
Ostap::Models::BreitWignerMC::BreitWignerMC
( const Ostap::Models::BreitWignerMC& right, const char* name  )
  : RooAbsPdf ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x      ) 
  , m_mass     ( "m0"     , this , right.m_mass   ) 
  , m_widths   ( "widths" , this , right.m_widths ) 
    //
  , m_bw       ( right.m_bw->clone() ) 
{
  setPars() ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::BreitWignerMC::~BreitWignerMC(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::BreitWignerMC*
Ostap::Models::BreitWignerMC::clone ( const char* name ) const 
{ return new Ostap::Models::BreitWignerMC( *this , name ) ; }
// ============================================================================
void Ostap::Models::BreitWignerMC::setPars () const 
{
  // set the mass 
  m_bw->setM0 ( m_mass ) ;
  // 
  // set the widths
  RooAbsArg*       g   = 0 ;
  const RooArgSet* nset  = m_widths.nset() ;
  //
  unsigned short k = 0 ;
  Ostap::Utils::Iterator it ( m_widths ) ;
  while ( ( g = (RooAbsArg*) it.next() ) )
  {
    const RooAbsReal* r = dynamic_cast<RooAbsReal*> ( g ) ;
    if ( 0 == r ) { continue ; }
    //
    const double width = r->getVal ( nset ) ;
    //
    m_bw->setGamma ( k  , width ) ;
    //
    ++k ;
  }
  //
}
// ============================================================================
// get the amplitude
// ============================================================================
std::complex<double>  Ostap::Models::BreitWignerMC::amplitude () const  
{
  setPars() ;
  return m_bw->amplitude( m_x ) ;
}
// ============================================================================
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Models::BreitWignerMC::evaluate() const 
{
  setPars() ;
  return (*m_bw)(m_x) ;
}   
// ============================================================================
Int_t Ostap::Models::BreitWignerMC::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::BreitWignerMC::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  return m_bw->integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================  
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::BreitWignerMC::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , *m_bw , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================

// ============================================================================
/// Breit-wigner with interference 
// ============================================================================
// constructor from Breit-Wigner and backround 
// ============================================================================
Ostap::Models::BWI::BWI 
( const char*                       name  , 
  const Ostap::Models::BreitWigner& bw    ,
  RooAbsReal&                       b     , 
  RooAbsReal&                       ab    , 
  RooAbsReal&                       phib  ) 
  : BreitWigner ( bw , name ) 
  , m_b     ( "b"     , "backgground"       , this , b     ) 
  , m_ab    ( "ab"    , "backgronud factor" , this , ab    ) 
  , m_phib  ( "phib"  , "background phase"  , this , phib  )
{}
// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Models::BWI::BWI
( const char*          name      ,
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          mass      ,
  RooAbsReal&          width     ,
  const double         m1        ,
  const double         m2        ,
  const unsigned short L         , 
  RooAbsReal&          b         ,
  RooAbsReal&          ab        , 
  RooAbsReal&          phib      ) 
  : BreitWigner  ( name , title , x , mass , width , m1 , m2 , L ) 
  , m_b     ( "b"     , "backgground"       , this , b     ) 
  , m_ab    ( "ab"    , "backgronud factor" , this , ab    ) 
  , m_phib  ( "phib"  , "background phase"  , this , phib  )
{}
// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Models::BWI::BWI
( const char*          name      ,
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          mass      ,
  RooAbsReal&          width     ,
  const double         m1        ,
  const double         m2        ,
  const unsigned short L                         ,
  const Ostap::Math::FormFactors::JacksonRho rho , 
  RooAbsReal&          b         ,
  RooAbsReal&          ab        , 
  RooAbsReal&          phib      ) 
  : BreitWigner  ( name , title , x , mass , width , m1 , m2 , L , rho ) 
  , m_b     ( "b"     , "backgground"       , this , b     ) 
  , m_ab    ( "ab"    , "backgronud factor" , this , ab    ) 
  , m_phib  ( "phib"  , "background phase"  , this , phib  )
{}
// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Models::BWI::BWI
( const char*          name          ,
  const char*          title         ,
  RooAbsReal&          x             ,
  RooAbsReal&          mass          ,
  RooAbsReal&          width         ,
  const Ostap::Math::BreitWigner& bw , 
  RooAbsReal&          b             ,
  RooAbsReal&          ab            , 
  RooAbsReal&          phib          ) 
  : BreitWigner  ( name , title , x , mass , width , bw    ) 
  , m_b     ( "b"     , "backgground"       , this , b     ) 
  , m_ab    ( "ab"    , "backgronud factor" , this , ab    ) 
  , m_phib  ( "phib"  , "background phase"  , this , phib  )
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::BWI::BWI
( const Ostap::Models::BWI& right , 
  const char*               name  ) 
  : BreitWigner ( right , name )
    //
  , m_b      ( "b"    , this , right.m_b    ) 
  , m_ab     ( "ab"   , this , right.m_ab   ) 
  , m_phib   ( "phib" , this , right.m_phib )
{}
// ============================================================================
Ostap::Models::BWI::~BWI(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::BWI* 
Ostap::Models::BWI::clone ( const char* name ) const 
{ return new Ostap::Models::BWI ( *this , name ) ; }
// ============================================================================
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Models::BWI::evaluate () const 
{
  //
  setPars () ;
  //
  const double b    = m_b    ;
  const double ab   = m_ab   ;
  const double phib = m_phib ;
  //
  const std::complex<double> ib  = ab * std::exp ( std::complex<double>(0,1) * phib ) ;
  //
  const double x = m_x ;
  const std::complex<double> g   = function().gamma ( x ) / function().gamma0 () ;
  //
  const std::complex<double> amp = g.real() * amplitude() + b * ib ;
  //
  return std::norm ( amp ) ;
}
// ============================================================================
Int_t Ostap::Models::BWI::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const { return 0 ; }
// ============================================================================



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
  return m_bw.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::BW23L::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_bw , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================



  
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
  RooAbsReal&                g0        ,
  const Ostap::Math::Flatte& flatte    ) 
  : RooAbsPdf ( name , title ) 
//
  , m_x      ( "x"     , "Observable"    , this , x     ) 
  , m_m0     ( "m0"    , "Peak"          , this , m0    ) 
  , m_m0g1   ( "m0g1"  , "M0*Gamma1"     , this , m0g1  )
  , m_g2og1  ( "g2og1" , "Gamma2/Gamma1" , this , g2og1 )
  , m_g0     ( "g0"    , "Gamma0"        , this , g0    )
    //
  , m_flatte ( flatte.clone () ) 
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
  , m_g0    ( "g0"    , this , right.m_g0    )
    //
  , m_flatte ( right.m_flatte->clone() ) 
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
  m_flatte->setM0     ( m_m0    ) ;
  m_flatte->setM0G1   ( m_m0g1  ) ;
  m_flatte->setG2oG1  ( m_g2og1 ) ;
  m_flatte->setG0     ( m_g0    ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Flatte::evaluate() const 
{
  setPars () ;
  return (*m_flatte) ( m_x ) ;
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
  return m_flatte->integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ===========================================================================
// get the amplitude 
// ===========================================================================
std::complex<double> Ostap::Models::Flatte::amplitude () const  
{
  setPars () ;
  return m_flatte->amplitude ( m_x ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Flatte::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , *m_flatte , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================



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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::LASS::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_lass , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================



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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::LASS23L::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_lass , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================


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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Bugg::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_bugg , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================


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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Bugg23L::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_bugg , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================

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
  return m_voigt.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Voigt::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_voigt , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================



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
  return m_voigt.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PseudoVoigt::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_voigt , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================




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
  return m_swanson.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Swanson::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_swanson , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================


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
  return m_cb.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::CrystalBall::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_cb , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================





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
  return m_cb.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::CrystalBallRS::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_cb , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_cb2.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::CrystalBallDS::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x . empty ()      ) { return {} ; }
  // 
  auto m0     = m_m0    . getValBatch ( begin , batchSize ) ;
  if ( !m0    . empty ()  ) { return {} ; }
  //
  auto sigma  = m_sigma . getValBatch ( begin , batchSize ) ;
  if ( !sigma . empty ()  ) { return {} ; }
  //
  auto alphaL = m_alphaL. getValBatch ( begin , batchSize ) ;
  if ( !alphaL. empty ()  ) { return {} ; }
  //
  auto alphaR = m_alphaR. getValBatch ( begin , batchSize ) ;
  if ( !alphaR. empty ()  ) { return {} ; }
  //
  auto nL     = m_nL    . getValBatch ( begin , batchSize ) ;
  if ( !nL    . empty ()  ) { return {} ; }
  //
  auto nR     = m_nR    . getValBatch ( begin , batchSize ) ;
  if ( !nR    . empty ()  ) { return {} ; }
  //
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_cb2 , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================










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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Needham::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto m0     = m_m0    . getValBatch ( begin , batchSize ) ;
  if ( !m0    . empty ()  ) { return {} ; }
  //
  auto sigma  = m_sigma . getValBatch ( begin , batchSize ) ;
  if ( !sigma . empty ()  ) { return {} ; }
  //
  auto  a0    = m_a0    . getValBatch ( begin , batchSize ) ;
  if ( !a0    . empty ()  ) { return {} ; }
  auto  a1    = m_a1    . getValBatch ( begin , batchSize ) ;
  if ( !a0    . empty ()  ) { return {} ; }
  auto  a2    = m_a2    . getValBatch ( begin , batchSize ) ;
  if ( !a2    . empty ()  ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_needham , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================












// ============================================================================
// Apollonios
// ============================================================================
Ostap::Models::Apollonios::Apollonios
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
Ostap::Models::Apollonios::Apollonios
( const Ostap::Models::Apollonios& right , 
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
Ostap::Models::Apollonios::~Apollonios (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Apollonios*
Ostap::Models::Apollonios::clone ( const char* name ) const 
{ return new Ostap::Models::Apollonios ( *this , name ) ; }
// ============================================================================
void Ostap::Models::Apollonios::setPars () const 
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
Double_t Ostap::Models::Apollonios::evaluate() const 
{
  //
  setPars () ;
  //
  return m_apo ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Apollonios::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Apollonios::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_apo.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Apollonios::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x     = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x     . empty ()      ) { return {} ; }
  // 
  auto  m0    = m_m0    . getValBatch ( begin , batchSize ) ;
  if ( !m0    . empty ()  ) { return {} ; }
  //
  auto  sigma = m_sigma . getValBatch ( begin , batchSize ) ;
  if ( !sigma . empty ()  ) { return {} ; }
  //
  auto  alpha = m_alpha . getValBatch ( begin , batchSize ) ;
  if ( !alpha . empty ()  ) { return {} ; }
  //
  auto  n     = m_n     . getValBatch ( begin , batchSize ) ;
  if ( !n . empty ()  ) { return {} ; }
  //
  auto  b     = m_b     . getValBatch ( begin , batchSize ) ;
  if ( !b     . empty ()  ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_apo , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================



// ============================================================================
// Apollonios2
// ============================================================================
Ostap::Models::Apollonios2::Apollonios2
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
Ostap::Models::Apollonios2::Apollonios2
( const Ostap::Models::Apollonios2& right , 
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
Ostap::Models::Apollonios2::~Apollonios2 (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Apollonios2*
Ostap::Models::Apollonios2::clone ( const char* name ) const 
{ return new Ostap::Models::Apollonios2 ( *this , name ) ; }
// ============================================================================
void Ostap::Models::Apollonios2::setPars () const 
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
Double_t Ostap::Models::Apollonios2::evaluate() const 
{
  //
  setPars () ;
  //
  return m_apo2 ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Apollonios2::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Apollonios2::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_apo2.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Apollonios2::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x     . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()  ) { return {} ; }
  // 
  auto  m0     = m_m0    . getValBatch ( begin , batchSize ) ;
  if ( !m0     . empty ()  ) { return {} ; }
  //
  auto  sigmaL = m_sigmaL . getValBatch ( begin , batchSize ) ;
  if ( !sigmaL . empty ()  ) { return {} ; }
  //
  auto  sigmaR = m_sigmaR . getValBatch ( begin , batchSize ) ;
  if ( !sigmaR . empty ()  ) { return {} ; }
  //
  auto  beta   = m_beta   . getValBatch ( begin , batchSize ) ;
  if ( !beta   . empty ()  ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_apo2 , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================




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
  return m_bg.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::BifurcatedGauss::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x       . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()  ) { return {} ; }
  // 
  auto  peak   = m_peak    . getValBatch ( begin , batchSize ) ;
  if ( !peak   . empty ()  ) { return {} ; }
  //
  auto  sigmaL = m_sigmaL . getValBatch ( begin , batchSize ) ;
  if ( !sigmaL . empty ()  ) { return {} ; }
  //
  auto  sigmaR = m_sigmaR . getValBatch ( begin , batchSize ) ;
  if ( !sigmaR . empty ()  ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_bg , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================


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
  return m_ggv1.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::GenGaussV1::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x       . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()  ) { return {} ; }
  // 
  auto  mu     = m_mu      . getValBatch ( begin , batchSize ) ;
  if ( !mu     . empty ()  ) { return {} ; }
  //
  auto  alpha  = m_alpha   . getValBatch ( begin , batchSize ) ;
  if ( !alpha  . empty ()  ) { return {} ; }
  //
  auto  beta   = m_beta    . getValBatch ( begin , batchSize ) ;
  if ( !beta   . empty ()  ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_ggv1 , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================



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
  return m_ggv2.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::GenGaussV2::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x       . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()  ) { return {} ; }
  // 
  auto  xi     = m_xi      . getValBatch ( begin , batchSize ) ;
  if ( !xi     . empty ()  ) { return {} ; }
  //
  auto  alpha  = m_alpha   . getValBatch ( begin , batchSize ) ;
  if ( !alpha  . empty ()  ) { return {} ; }
  //
  auto  kappa  = m_kappa   . getValBatch ( begin , batchSize ) ;
  if ( !kappa  . empty ()  ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_ggv2 , x ) ;
  //
  return output ;
}
// ======================================================================
#endif
// ======================================================================




// ============================================================================
//         SkewGauss
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::SkewGauss::SkewGauss
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          xi        , 
  RooAbsReal&          omega     , 
  RooAbsReal&          alpha     )  
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "x"       , "Observable" , this , x      ) 
  , m_xi      ( "xi"      , "xi"         , this , xi     ) 
  , m_omega   ( "omega"   , "omega"      , this , omega  ) 
  , m_alpha   ( "alpha"   , "alpha"      , this , alpha  ) 
//
  , m_sg    ( 0 , 1 , 0 ) 
{
  //
  setPars () ;
  //
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::SkewGauss::SkewGauss
( const Ostap::Models::SkewGauss& right , 
  const char*                        name   ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_xi     ( "xi"     , this , right.m_xi     ) 
  , m_omega  ( "omega"  , this , right.m_omega  ) 
  , m_alpha  ( "alpha"  , this , right.m_alpha  ) 
//
  , m_sg     ( right.m_sg ) 
{
  setPars () ;
}
// ============================================================================
// desctructor
// ============================================================================
Ostap::Models::SkewGauss::~SkewGauss (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::SkewGauss*
Ostap::Models::SkewGauss::clone( const char* name ) const 
{ return new Ostap::Models::SkewGauss ( *this , name ) ; }
// ============================================================================
void Ostap::Models::SkewGauss::setPars () const 
{
  //
  m_sg . setXi     ( m_xi    ) ;
  m_sg . setOmega  ( m_omega ) ;
  m_sg . setAlpha  ( m_alpha ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::SkewGauss::evaluate() const 
{
  //
  setPars () ;
  //
  return m_sg      ( m_x      ) ;
}
// ============================================================================
Int_t Ostap::Models::SkewGauss::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::SkewGauss::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_sg.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::SkewGauss::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x       . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()  ) { return {} ; }
  // 
  auto  xi     = m_xi      . getValBatch ( begin , batchSize ) ;
  if ( !xi     . empty ()  ) { return {} ; }
  //
  auto  omega  = m_omega   . getValBatch ( begin , batchSize ) ;
  if ( !omega  . empty ()  ) { return {} ; }
  //
  auto  alpha  = m_alpha   . getValBatch ( begin , batchSize ) ;
  if ( !alpha  . empty ()  ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_sg , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================


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
  return m_bukin.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Bukin::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x       . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()  ) { return {} ; }
  // 
  auto  peak   = m_peak    . getValBatch ( begin , batchSize ) ;
  if ( !peak   . empty ()  ) { return {} ; }
  //
  auto  sigma  = m_sigma   . getValBatch ( begin , batchSize ) ;
  if ( !sigma  . empty ()  ) { return {} ; }
  //
  auto  xi     = m_xi      . getValBatch ( begin , batchSize ) ;
  if ( !xi     . empty ()  ) { return {} ; }
  //
  auto  rhoL   = m_rhoL    . getValBatch ( begin , batchSize ) ;
  if ( !rhoL   . empty ()  ) { return {} ; }
  //
  auto  rhoR   = m_rhoR    . getValBatch ( begin , batchSize ) ;
  if ( !rhoR   . empty ()  ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_bukin , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_stt.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::StudentT::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x       . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()  ) { return {} ; }
  // 
  auto  mu     = m_mu      . getValBatch ( begin , batchSize ) ;
  if ( !mu     . empty ()  ) { return {} ; }
  //
  auto  sigma  = m_sigma   . getValBatch ( begin , batchSize ) ;
  if ( !sigma  . empty ()  ) { return {} ; }
  //
  auto  n      = m_n       . getValBatch ( begin , batchSize ) ;
  if ( !n      . empty ()  ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_stt , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_stt.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::BifurcatedStudentT::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x       . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()  ) { return {} ; }
  // 
  auto  mu     = m_mu      . getValBatch ( begin , batchSize ) ;
  if ( !mu     . empty ()  ) { return {} ; }
  //
  auto  sigmaL = m_sigmaL  . getValBatch ( begin , batchSize ) ;
  if ( !sigmaL . empty ()  ) { return {} ; }
  //
  auto  sigmaR = m_sigmaR  . getValBatch ( begin , batchSize ) ;
  if ( !sigmaR . empty ()  ) { return {} ; }
  //
  auto  nL     = m_nL      . getValBatch ( begin , batchSize ) ;
  if ( !nL     . empty ()  ) { return {} ; }
  //
  auto  nR     = m_nR      . getValBatch ( begin , batchSize ) ;
  if ( !nR     . empty ()  ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_stt , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_gca.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::GramCharlierA::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x       . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()  ) { return {} ; }
  // 
  auto  m0     = m_m0      . getValBatch ( begin , batchSize ) ;
  if ( !m0     . empty ()  ) { return {} ; }
  //
  auto  sigma  = m_sigma   . getValBatch ( begin , batchSize ) ;
  if ( !sigma  . empty ()  ) { return {} ; }
  //
  auto  kappa3 = m_kappa3  . getValBatch ( begin , batchSize ) ;
  if ( !kappa3 . empty ()  ) { return {} ; }
  //
  auto  kappa4 = m_kappa4  . getValBatch ( begin , batchSize ) ;
  if ( !kappa4 . empty ()  ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_gca , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PhaseSpace2::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x       . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()  ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  compute_X ( output , m_ps2 , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_left.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PhaseSpaceLeft::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  // 
  auto  t      = m_threshold . getValBatch ( begin , batchSize ) ;
  if ( !t      . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_left , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_right.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PhaseSpaceRight::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  // 
  auto  t      = m_threshold . getValBatch ( begin , batchSize ) ;
  if ( !t      . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_right , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_ps.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PhaseSpaceNL::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  // 
  auto  l      = m_low       . getValBatch ( begin , batchSize ) ;
  if ( !l      . empty ()    ) { return {} ; }
  // 
  auto  h      = m_high      . getValBatch ( begin , batchSize ) ;
  if ( !h      . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_ps , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PhaseSpace23L::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  compute_X ( output , m_ps23L , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  , m_ps       ( low , high , L , N , phis.getSize() ) 
{
  //
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , "Ostap::Models::PhaseSpacePol" );
  Ostap::Assert ( ::size ( m_phis ) == m_ps.npars() , 
                  "#phis/#npars mismatch!"            , "Ostap::Models::PhaseSpacePol" );
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
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , "Ostap::Models::PhaseSpacePol" );
  Ostap::Assert ( ::size ( m_phis ) == m_ps.npars() , 
                  "#phis/#npars mismatch!"            , "Ostap::Models::PhaseSpacePol" );
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
  ::set_pars ( m_phis ,  m_ps ) ;
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
  return m_ps.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PhaseSpacePol::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_ps , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================



// ============================================================================
//  PhaseSpaceLeft x expo x pol 
// ============================================================================
Ostap::Models::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol 
( const char*                        name      , 
  const char*                        title     ,
  RooRealVar&                        x         ,
  const Ostap::Math::PhaseSpaceLeft& ps        , 
  RooAbsReal&                        tau       ,
  RooAbsReal&                        scale     ,
  RooArgList&                        phis      )
  : RooAbsPdf ( name , title ) 
    //
  , m_x        ( "x"       , "Observable"   , this , x     ) 
  , m_tau      ( "tau"     , "Exponent"     , this , tau   )
  , m_scale    ( "scale"   , "Scale-factor" , this , scale )
  , m_phis     ( "phi"     , "Coefficients" , this )
    //
  , m_ps       ( ps , phis.getSize() , 0.0 , x.getMin() , x.getMax() ) 
{
  //
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , 
                  "Ostap::Models::PhaseSpaceLeftExpoPol" ) ;
  Ostap::Assert ( ::size ( m_phis ) == m_ps.npars() , 
                  "#phis/#npars mismatch!"            , 
                  "Ostap::Models::PhaseSpaceLeftExpoPol" ) ;
  //
  setPars() ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PhaseSpaceLeftExpoPol::~PhaseSpaceLeftExpoPol () {}
// ======================================================================
// "copy" constructor 
// ======================================================================
Ostap::Models::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol 
( const Ostap::Models::PhaseSpaceLeftExpoPol& right , 
  const char*                                 name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_tau      ( "tau"    , this , right.m_tau   )
  , m_scale    ( "scale"  , this , right.m_scale )
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_ps       ( right.m_ps       ) 
{
  setPars() ;
}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PhaseSpaceLeftExpoPol* 
Ostap::Models::PhaseSpaceLeftExpoPol::clone ( const char* name ) const 
{ return new Ostap::Models::PhaseSpaceLeftExpoPol ( *this , name ) ; }
// ============================================================================
void Ostap::Models::PhaseSpaceLeftExpoPol::setPars () const 
{
  //
  m_ps.setTau    ( m_tau   ) ;
  m_ps.setScale  ( m_scale ) ;
  //
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  ::set_pars ( m_phis , m_ps ) ;
  //
}
// ============================================================================
// evaluate the function
// ============================================================================
Double_t Ostap::Models::PhaseSpaceLeftExpoPol::evaluate () const 
{
  setPars () ;
  return m_ps ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::PhaseSpaceLeftExpoPol::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PhaseSpaceLeftExpoPol::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_ps.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PhaseSpaceLeftExpoPol::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_ps , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , 
                  "Ostap::Models::PolyPositive" ) ;
  Ostap::Assert ( ::size ( m_phis ) == m_positive.npars()  , 
                  "#phis/#npars mismatch!"             ,
                  "Ostap::Models::PolyPositive" ) ;
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
  ::set_pars ( m_phis , m_positive ) ;
  //
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
  return m_positive.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PolyPositive::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_positive , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , 
                  "Ostap::Models::PolyPositiveEven" ) ;
  Ostap::Assert ( ::size ( m_phis ) == m_even.npars() , 
                  "#phis/#npars mismatch!"            ,
                  "Ostap::Models::PolyPositiveEven" ) ;
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
  ::set_pars ( m_phis , m_even ) ;
  //
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
  return m_even.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PolyPositiveEven::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_even , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================


// ============================================================================
// monotonic polinomial
// ============================================================================
Ostap::Models::PolyMonotonic::PolyMonotonic
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
  , m_monotonic ( phis.getSize() , xmin , xmax , increasing ) 
{
  //
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , 
                  "Ostap::Models::PolyMonotonic" ) ;
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PolyMonotonic::PolyMonotonic
( const Ostap::Models::PolyMonotonic&  right ,      
  const char*                              name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x          ( "x"      , this , right.m_x     ) 
  , m_phis       ( "phis"   , this , right.m_phis  ) 
    //
  , m_monotonic ( right.m_monotonic ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::PolyMonotonic::~PolyMonotonic() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PolyMonotonic*
Ostap::Models::PolyMonotonic::clone( const char* name ) const 
{ return new Ostap::Models::PolyMonotonic(*this,name) ; }
// ============================================================================
void Ostap::Models::PolyMonotonic::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  ::set_pars ( m_phis , m_monotonic ) ;
  //
}
//
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PolyMonotonic::evaluate() const 
{
  //
  setPars () ;
  //
  return m_monotonic ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::PolyMonotonic::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PolyMonotonic::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_monotonic.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PolyMonotonic::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_monotonic , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , "Ostap::Models::PolyConvex" ) ;
  //
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
  ::set_pars ( m_phis , m_convex ) ;
  //
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
  return m_convex.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PolyConvex::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_convex , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , "Ostap::Models::PolyConvexOnly" ) ;
  //
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
  ::set_pars ( m_phis , m_convex ) ;
  //
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
  return m_convex.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PolyConvexOnly::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_convex , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , "Ostap::Models::PolySigmoid" ) ;
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
  ::set_pars ( m_phis , m_sigmoid ) ;
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
  return m_sigmoid.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PolySigmoid::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto  alpha  = m_alpha     . getValBatch ( begin , batchSize ) ;
  if ( !alpha  . empty ()    ) { return {} ; }
  // 
  auto  x0     = m_x0        . getValBatch ( begin , batchSize ) ;
  if ( !x0     . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_sigmoid , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , "Ostap::Models::PositiveSpline" ) ;
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
  ::set_pars ( m_phis , m_spline ) ;
  //
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
  return m_spline.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::PositiveSpline::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_spline , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================





// ============================================================================
// monotonic spline 
// ============================================================================
/** constructor with the spline 
 *  @param name  the name 
 *  @param title the  title
 *  @param x     the  variable 
 *  @param spine the spline  
 *  @param phis  vector of parameters 
 */
// ============================================================================
Ostap::Models::MonotonicSpline::MonotonicSpline 
( const char*                          name, 
  const char*                          title     ,
  RooAbsReal&                          x         ,
  const Ostap::Math::MonotonicSpline& spline    ,   // the spline 
  RooArgList&                          phis      )   // parameters
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable"   , this , x ) 
  , m_phis     ( "phi"     , "Coefficients" , this     )
    //
  , m_spline   ( spline ) 
{
  //
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , "Ostap::Models::MonotonicSpline" ) ;
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::MonotonicSpline::MonotonicSpline 
( const Ostap::Models::MonotonicSpline&  right ,      
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
Ostap::Models::MonotonicSpline::~MonotonicSpline() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::MonotonicSpline*
Ostap::Models::MonotonicSpline::clone( const char* name ) const 
{ return new Ostap::Models::MonotonicSpline(*this,name) ; }
// ============================================================================
void Ostap::Models::MonotonicSpline::setPars () const 
{
  RooAbsArg*       phi   = 0 ;
  const RooArgSet* nset  = m_phis.nset() ;
  //
  std::vector<double> sin2phi ;
  //
  ::set_pars ( m_phis , m_spline ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::MonotonicSpline::evaluate() const 
{
  //
  setPars () ;
  //
  return m_spline ( m_x ) ; 
}
// ============================================================================
Int_t Ostap::Models::MonotonicSpline::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::MonotonicSpline::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_spline.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::MonotonicSpline::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_spline , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================


// ============================================================================
// convex spline 
// ============================================================================
/** constructor with the spline 
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
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , "Ostap::Models::ConvexSpline" ) ;
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
  ::set_pars ( m_phis , m_spline ) ;
  //
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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::ConvexSpline::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_spline , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================

// ============================================================================
// convex spline 
// ============================================================================
/** constructor with the spline 
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
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , "Ostap::Models::ConvexOnlySpline" ) ;
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
  ::set_pars ( m_phis , m_spline ) ;
  //
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
  return m_spline.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::ConvexOnlySpline::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_spline , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , "Ostap::Models::ExpoPositive" ) ;
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
  ::set_pars ( m_phis , m_positive ) ;
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
  return m_positive.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::ExpoPositive::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto  t      = m_tau         . getValBatch ( begin , batchSize ) ;
  if ( !t      . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_positive , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  ::copy_real   ( phis , m_phis , "Invalid parameter!" , "Ostap::Models::TwoExpoPositive" ) ;
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
  ::set_pars ( m_phis , m_2expopos ) ;
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
  return m_2expopos.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::TwoExpoPositive::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto  alpha  = m_alpha     . getValBatch ( begin , batchSize ) ;
  if ( !alpha  . empty ()    ) { return {} ; }
  // 
  auto  delta  = m_delta     . getValBatch ( begin , batchSize ) ;
  if ( !delta  . empty ()    ) { return {} ; }
  // 
  auto  x0     = m_x0         . getValBatch ( begin , batchSize ) ;
  if ( !x0     . empty ()    ) { return {} ; }
  // 
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_2expopos , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_gamma.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
namespace 
{
  template<class TX , 
           class TK ,
           class TT, class FUN>
  void compute_GD ( RooSpan<double> output , FUN& fun , TX x , TK k , TT theta ) 
  {
    const int n = output.size();
    for ( int i = 0 ; i < n ; ++i ) 
    {
      fun.setK        ( k     [i] ) ;
      fun.setTheta    ( theta [i] ) ;
      output[i] = fun ( x     [i] ) ; 
    }
  }
}
// ============================================================================
RooSpan<double> 
Ostap::Models::GammaDist::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  auto  k      = m_k         . getValBatch ( begin , batchSize ) ;
  auto  theta  = m_theta     . getValBatch ( begin , batchSize ) ;
  //
  const bool bx = x    .empty()  ;
  const bool bk = k    .empty()  ;
  const bool bt = theta.empty()  ;
  //
  if ( bx &&  bk &&  bt ) { return {}; } //
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  if      (  bx && !bk && !bt ) 
  { compute_GD ( output , m_gamma ,         x   , BA ( m_k ) , BA ( m_theta ) ) ; }
  else if ( !bx &&  bk && !bt ) 
  { compute_GD ( output , m_gamma , BA  ( m_x ) ,        k   , BA ( m_theta ) ) ; }
  else if ( !bx && !bk &&  bt ) 
  { compute_GD ( output , m_gamma , BA  ( m_x ) , BA ( m_k ) ,        theta   ) ; }
  else if ( !bx &&  bk &&  bt ) 
  { compute_GD ( output , m_gamma , BA  ( m_x ) ,        k   ,        theta   ) ; }
  else if (  bx && !bk &&  bt ) 
  { compute_GD ( output , m_gamma ,         x   , BA ( m_k ) ,        theta   ) ; }
  else if (  bx &&  bk && !bt ) 
  { compute_GD ( output , m_gamma , x           ,        k   , BA ( m_theta ) ) ; }
  else if (  bx &&  bk &&  bt ) 
  { compute_GD ( output , m_gamma , x           ,        k   ,        theta   ) ; }
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_ggamma.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::GenGammaDist::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_ggamma , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_amoroso.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Amoroso::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_amoroso , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_gamma.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::LogGammaDist::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  auto  k      = m_k         . getValBatch ( begin , batchSize ) ;
  auto  theta  = m_theta     . getValBatch ( begin , batchSize ) ;
  //
  const bool bx = x    .empty()  ;
  const bool bk = k    .empty()  ;
  const bool bt = theta.empty()  ;
  //
  if ( bx &&  bk &&  bt ) { return {}; } //
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  if      (  bx && !bk && !bt ) 
  { compute_GD ( output , m_gamma ,         x   , BA ( m_k ) , BA ( m_theta ) ) ; }
  else if ( !bx &&  bk && !bt ) 
  { compute_GD ( output , m_gamma , BA  ( m_x ) ,        k   , BA ( m_theta ) ) ; }
  else if ( !bx && !bk &&  bt ) 
  { compute_GD ( output , m_gamma , BA  ( m_x ) , BA ( m_k ) ,        theta   ) ; }
  else if ( !bx &&  bk &&  bt ) 
  { compute_GD ( output , m_gamma , BA  ( m_x ) ,        k   ,        theta   ) ; }
  else if (  bx && !bk &&  bt ) 
  { compute_GD ( output , m_gamma ,         x   , BA ( m_k ) ,        theta   ) ; }
  else if (  bx &&  bk && !bt ) 
  { compute_GD ( output , m_gamma , x           ,        k   , BA ( m_theta ) ) ; }
  else if (  bx &&  bk &&  bt ) 
  { compute_GD ( output , m_gamma , x           ,        k   ,        theta   ) ; }
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_gamma.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Log10GammaDist::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  auto  k      = m_k         . getValBatch ( begin , batchSize ) ;
  auto  theta  = m_theta     . getValBatch ( begin , batchSize ) ;
  //
  const bool bx = x    .empty()  ;
  const bool bk = k    .empty()  ;
  const bool bt = theta.empty()  ;
  //
  if ( bx &&  bk &&  bt ) { return {}; } //
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  if      (  bx && !bk && !bt ) 
  { compute_GD ( output , m_gamma ,         x   , BA ( m_k ) , BA ( m_theta ) ) ; }
  else if ( !bx &&  bk && !bt ) 
  { compute_GD ( output , m_gamma , BA  ( m_x ) ,        k   , BA ( m_theta ) ) ; }
  else if ( !bx && !bk &&  bt ) 
  { compute_GD ( output , m_gamma , BA  ( m_x ) , BA ( m_k ) ,        theta   ) ; }
  else if ( !bx &&  bk &&  bt ) 
  { compute_GD ( output , m_gamma , BA  ( m_x ) ,        k   ,        theta   ) ; }
  else if (  bx && !bk &&  bt ) 
  { compute_GD ( output , m_gamma ,         x   , BA ( m_k ) ,        theta   ) ; }
  else if (  bx &&  bk && !bt ) 
  { compute_GD ( output , m_gamma , x           ,        k   , BA ( m_theta ) ) ; }
  else if (  bx &&  bk &&  bt ) 
  { compute_GD ( output , m_gamma , x           ,        k   ,        theta   ) ; }
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_lgamma.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::LogGamma::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_lgamma , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_betap.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::BetaPrime::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_betap , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_sinhasinh.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::SinhAsinh::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_sinhasinh , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_johnsonSU.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::JohnsonSU::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_johnsonSU , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_landau.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Landau::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_landau , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_atlas.integral ( m_x.min( rangeName ) , m_x.max( rangeName ) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Atlas::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_atlas , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_sech.integral ( m_x.min ( rangeName ) , m_x.max ( rangeName ) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Sech::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_sech , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================





// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Losev::Losev
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          mu     ,
  RooAbsReal&          alpha  ,
  RooAbsReal&          beta   )
  : RooAbsPdf ( name , title ) 
    //
  , m_x      ( "x"      , "Observable"  , this , x      ) 
  , m_mu     ( "mu"     , "location"    , this , mu     ) 
  , m_alpha  ( "alpha"  , "left-slope"  , this , alpha  ) 
  , m_beta   ( "beta"   , "right-slope" , this , beta   ) 
    //
  , m_losev  ( 0 , 1 , 1 ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Losev::Losev
( const Ostap::Models::Losev&  right ,
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_mu     ( "mu"     , this , right.m_mu     )
  , m_alpha  ( "alpha"  , this , right.m_alpha  )
  , m_beta   ( "beta"   , this , right.m_beta   )
    //
  , m_losev  (                   right.m_losev  ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Losev::~Losev () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Losev*
Ostap::Models::Losev::clone( const char* name ) const 
{ return new Ostap::Models::Losev( *this , name) ; }
// ============================================================================
void Ostap::Models::Losev::setPars () const 
{
  //
  m_losev.setMu    ( m_mu    ) ;
  m_losev.setAlpha ( m_alpha ) ;
  m_losev.setBeta  ( m_beta  ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Losev::evaluate() const 
{
  //
  setPars () ;
  //
  return m_losev ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Losev::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Losev::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_losev.integral ( m_x.min ( rangeName ) , m_x.max ( rangeName ) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Losev::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_losev , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_logistic.integral ( m_x.min( rangeName ) , m_x.max( rangeName ) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Logistic::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_logistic , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_argus.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Argus::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_argus , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Slash::Slash
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          mu     ,
  RooAbsReal&          scale  ) 
  : RooAbsPdf ( name , title ) 
    //
  , m_x       ( "x"      , "Observable" , this , x      ) 
  , m_mu      ( "mu"     , "location"   , this , mu     ) 
  , m_scale   ( "scale"  , "scale"      , this , scale  ) 
    //
  , m_slash   ( 0 , 1 ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Slash::Slash
( const Ostap::Models::Slash&  right ,
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "x"      , this , right.m_x      ) 
  , m_mu     ( "mu"     , this , right.m_mu     )
  , m_scale  ( "scale"  , this , right.m_scale  )
//
  , m_slash  (                   right.m_slash  ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Slash::~Slash () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Slash*
Ostap::Models::Slash::clone( const char* name ) const 
{ return new Ostap::Models::Slash( *this , name) ; }
// // ============================================================================
void Ostap::Models::Slash::setPars () const 
{  
  m_slash.setMu    ( m_mu    ) ;
  m_slash.setScale ( m_scale ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Slash::evaluate() const 
{
  setPars () ;
  return m_slash ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Slash::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Slash::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_slash.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Slash::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_slash , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::AsymmetricLaplace::AsymmetricLaplace
( const char*          name    , 
  const char*          title   ,
  RooAbsReal&          x       ,
  RooAbsReal&          mu      ,
  RooAbsReal&          lambdaL ,
  RooAbsReal&          lambdaR )
  : RooAbsPdf ( name , title ) 
    //
  , m_x       ( "x"       , "Observable"                  , this , x       ) 
  , m_mu      ( "mu"      , "location"                    , this , mu      ) 
  , m_lambdaL ( "lambdaL" , "``left'' exponential slope"  , this , lambdaL ) 
  , m_lambdaR ( "lambdaR" , "``right'' exponential slope" , this , lambdaR ) 
    //
  , m_laplace ( 0 , 1 , 1 ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::AsymmetricLaplace::AsymmetricLaplace
( const Ostap::Models::AsymmetricLaplace&  right ,
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x        ( "x"       , this , right.m_x       ) 
  , m_mu       ( "mu"      , this , right.m_mu      )
  , m_lambdaL  ( "lambdaL" , this , right.m_lambdaL )
  , m_lambdaR  ( "lambdaR" , this , right.m_lambdaR )
    //
  , m_laplace  ( right.m_laplace ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::AsymmetricLaplace::~AsymmetricLaplace(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::AsymmetricLaplace*
Ostap::Models::AsymmetricLaplace::clone( const char* name ) const 
{ return new Ostap::Models::AsymmetricLaplace( *this , name) ; }
// ============================================================================
void Ostap::Models::AsymmetricLaplace::setPars () const 
{
  m_laplace.setMu      ( m_mu      ) ;
  m_laplace.setLambdaL ( m_lambdaL ) ;
  m_laplace.setLambdaR ( m_lambdaR ) ;  
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::AsymmetricLaplace::evaluate() const 
{
  setPars () ;
  return m_laplace( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::AsymmetricLaplace::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::AsymmetricLaplace::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_laplace.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::AsymmetricLaplace::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_laplace , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_tsallis.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Tsallis::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_tsallis , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_qgsm.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::QGSM::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_qgsm , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
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
  return m_2expos.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::TwoExpos::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_2expos , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::DoubleGauss::DoubleGauss
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          sigma     ,  // narrow sigma 
  RooAbsReal&          fraction  ,  // fraction of narrow sigma 
  RooAbsReal&          scale     ,  // wide/narrow sigma ratio    
  RooAbsReal&          mean      )  // mean, presumably  fixed at 0
  //
  : RooAbsPdf( name , title ) 
  , m_x        ( "x"         , "Observable"          , this , x        ) 
  , m_sigma    ( "sigma"     , "Narrow sigma"        , this , sigma    ) 
  , m_fraction ( "fraction"  , "Fraction"            , this , fraction ) 
  , m_scale    ( "scale"     , "Scale"               , this , scale    ) 
  , m_mean     ( "mean"      , "Mean"                , this , mean     ) 
  , m_2gauss  () 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::DoubleGauss::DoubleGauss
( const Ostap::Models::DoubleGauss& right ,  
  const char*                          name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x        ( "x"        , this , right.m_x        ) 
  , m_sigma    ( "sigma"    , this , right.m_sigma    )
  , m_fraction ( "fraction" , this , right.m_fraction )
  , m_scale    ( "scale"    , this , right.m_scale    )
  , m_mean     ( "mean"     , this , right.m_mean    )
  , m_2gauss   ( right.m_2gauss ) 
{
  setPars() ;
}
// ============================================================================
// clone it!
// ============================================================================
Ostap::Models::DoubleGauss*
Ostap::Models::DoubleGauss::clone( const char* name ) const 
{ return new Ostap::Models::DoubleGauss ( *this , name ) ; }
// ============================================================================
void Ostap::Models::DoubleGauss::setPars () const 
{
  m_2gauss.setPeak      ( m_mean     ) ;
  m_2gauss.setSigma     ( m_sigma    ) ;
  m_2gauss.setScale     ( m_scale    ) ;
  m_2gauss.setFraction  ( m_fraction ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::DoubleGauss::evaluate() const 
{
  //
  setPars() ;
  return m_2gauss ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::DoubleGauss::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if (matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::DoubleGauss::analyticalIntegral
( Int_t       code      , 
  const char* rangeName ) const
{
  assert( code == 1 ) ;
  if ( 1 != code ){}
  //
  const double  xmax = m_x.max ( rangeName )  ;
  const double  xmin = m_x.min ( rangeName )  ;
  //
  setPars() ;
  return m_2gauss.integral ( xmin , xmax ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::DoubleGauss::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_2gauss , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================

// ============================================================================
/*  constructor from all parameters
 *  @param  z    the  variable 
 *  @param  eta  the shape eta-parameter 
 *  @param  b    the scale b-parameter 
 *  @param  xmin the bias  parameter
 */
// ============================================================================
Ostap::Models::Gumbel::Gumbel
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          mu        , 
  RooAbsReal&          beta      )
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"     , "Observable"           , this , x    ) 
  , m_mu       ( "mu"    , "Shift parameter/mode" , this , mu   ) 
  , m_beta     ( "beta"  , "Scale parameter"      , this , beta ) 
  , m_gumbel   () 
{
  setPars() ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Gumbel::Gumbel
( const Ostap::Models::Gumbel& right ,
  const char*                     name  ) 
  : RooAbsPdf  ( right , name ) 
  , m_x        ( "x"    , this , right.m_x    ) 
  , m_mu       ( "mu"   , this , right.m_mu   ) 
  , m_beta     ( "beta" , this , right.m_beta ) 
  , m_gumbel   ( right.m_gumbel ) 
{
  setPars() ;
}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Gumbel*
Ostap::Models::Gumbel::clone( const char* name ) const 
{ return new Ostap::Models::Gumbel(*this,name) ; }
// ============================================================================
void Ostap::Models::Gumbel::setPars () const 
{
  m_gumbel.setMu   ( m_mu   ) ;
  m_gumbel.setBeta ( m_beta ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Gumbel::evaluate() const 
{
  //
  setPars ();
  return m_gumbel ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Gumbel::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Gumbel::analyticalIntegral
( Int_t       code      , 
  const char* rangeName ) const
{
  assert ( code == 1 ) ;
  if ( 1 != code ){}
  //
  const double xmin =  m_x.min ( rangeName ) ;
  const double xmax =  m_x.max ( rangeName ) ;
  //
  setPars() ;
  return m_gumbel.integral ( xmin , xmax ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Gumbel::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_gumbel , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================


// ============================================================================
/* constructor from all parameters
 *  @param  x    the  variable 
 *  @param  scale the scale parameter
 *  @param  shape the shape parameter 
 *  @param  shift the shift parameter 
 */
// ============================================================================
Ostap::Models::Weibull::Weibull 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , // observable 
  RooAbsReal&          scale     , // scale/lambda 
  RooAbsReal&          shape     , // shape/k 
  RooAbsReal&          shift     ) // shift/x0 
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"      , "Observable"             , this , x     ) 
  , m_scale    ( "scale"  , "Scale parameter/lambda" , this , scale ) 
  , m_shape    ( "shape"  , "Shape parameter/k"      , this , shape ) 
  , m_shift    ( "shift"  , "Shift parameter/x0"     , this , shift ) 
  , m_weibull  ()  
{
  setPars () ;  
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Weibull::Weibull 
( const Ostap::Models::Weibull& right ,
  const char*                     name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"     , this , right.m_x     ) 
  , m_scale    ( "scale" , this , right.m_scale ) 
  , m_shape    ( "shape" , this , right.m_shape ) 
  , m_shift    ( "shift" , this , right.m_shift ) 
  , m_weibull  ( right.m_weibull ) 
{
  setPars () ;  
}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Weibull*
Ostap::Models::Weibull::clone( const char* name ) const 
{ return new Ostap::Models::Weibull(*this,name) ; }
// ============================================================================
void Ostap::Models::Weibull::setPars () const 
{
  m_weibull.setScale ( m_scale ) ;
  m_weibull.setShape ( m_shape ) ;
  m_weibull.setShift ( m_shift ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Weibull::evaluate() const 
{
  setPars() ;
  return m_weibull ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Weibull::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Weibull::analyticalIntegral
( Int_t       code      , 
  const char* rangeName ) const
{
  assert ( code == 1 ) ;
  if ( 1 != code ){}
  //
  const double xmin =  m_x.min ( rangeName ) ;
  const double xmax =  m_x.max ( rangeName ) ;
  //
  setPars() ;
  return m_weibull.integral ( xmin , xmax ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Weibull::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_weibull , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================



// ============================================================================
/*  constructor from all parameters
 *  @param  x      the variable 
 *  @param  mean   the mean/mode/median/location 
 *  @param  scale  the scale parameter 
 */
// ============================================================================
Ostap::Models::RaisingCosine::RaisingCosine 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , // observable 
  RooAbsReal&          mean      , // mean
  RooAbsReal&          scale     ) // scale
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"      , "Observable"               , this , x    ) 
  , m_mean     ( "mean"   , "Mean/location parameter"  , this , mean ) 
  , m_scale    ( "scale"  , "Scale parameter"          , this , scale ) 
  , m_rcos  ()  
{
  setPars () ;  
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::RaisingCosine::RaisingCosine 
( const Ostap::Models::RaisingCosine& right ,
  const char*                     name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"     , this , right.m_x     ) 
  , m_mean     ( "mean"  , this , right.m_mean  ) 
  , m_scale    ( "scale" , this , right.m_scale ) 
  , m_rcos     ( right.m_rcos ) 
{
  setPars () ;  
}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::RaisingCosine*
Ostap::Models::RaisingCosine::clone( const char* name ) const 
{ return new Ostap::Models::RaisingCosine(*this,name) ; }
// ============================================================================
void Ostap::Models::RaisingCosine::setPars () const 
{
  m_rcos.setMean  ( m_mean  ) ;
  m_rcos.setScale ( m_scale ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::RaisingCosine::evaluate() const 
{
  setPars() ;
  return m_rcos ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::RaisingCosine::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::RaisingCosine::analyticalIntegral
( Int_t       code      , 
  const char* rangeName ) const
{
  assert ( code == 1 ) ;
  if ( 1 != code ){}
  //
  const double xmin =  m_x.min ( rangeName ) ;
  const double xmax =  m_x.max ( rangeName ) ;
  //
  setPars() ;
  return m_rcos.integral ( xmin , xmax ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::RaisingCosine::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_rcos , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================




// ============================================================================
/*  constructor from all parameters
 *  @param  x      the variable 
 *  @param  mean   the mean/mode/median/location 
 *  @param  q      the q-value 
 *  @param  scale  the scale parameter 
 */
// ============================================================================
Ostap::Models::QGaussian::QGaussian
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , // observable 
  RooAbsReal&          mean      , // mean
  RooAbsReal&          q         , // q
  RooAbsReal&          scale     ) // scale
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"      , "Observable"               , this , x     ) 
  , m_mean     ( "mean"   , "Mean/location parameter"  , this , mean  ) 
  , m_q        ( "q"      , "Q-parameter"              , this , q     ) 
  , m_scale    ( "scale"  , "Scale parameter"          , this , scale ) 
  , m_qgauss   ()  
{
  setPars () ;  
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::QGaussian::QGaussian
( const Ostap::Models::QGaussian& right ,
  const char*                     name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"     , this , right.m_x     ) 
  , m_mean     ( "mean"  , this , right.m_mean  ) 
  , m_q        ( "q"     , this , right.m_q     ) 
  , m_scale    ( "scale" , this , right.m_scale ) 
  , m_qgauss   ( right.m_qgauss ) 
{
  setPars () ;  
}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::QGaussian*
Ostap::Models::QGaussian::clone( const char* name ) const 
{ return new Ostap::Models::QGaussian(*this,name) ; }
// ============================================================================
void Ostap::Models::QGaussian::setPars () const 
{
  m_qgauss.setMean  ( m_mean  ) ;
  m_qgauss.setQ     ( m_q     ) ;
  m_qgauss.setScale ( m_scale ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::QGaussian::evaluate() const 
{
  setPars() ;
  return m_qgauss ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::QGaussian::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::QGaussian::analyticalIntegral
( Int_t       code      , 
  const char* rangeName ) const
{
  assert ( code == 1 ) ;
  if ( 1 != code ){}
  //
  const double xmin =  m_x.min ( rangeName ) ;
  const double xmax =  m_x.max ( rangeName ) ;
  //
  setPars() ;
  return m_qgauss.integral ( xmin , xmax ) ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::QGaussian::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto  x      = m_x         . getValBatch ( begin , batchSize ) ;
  if (  x      . empty ()    ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars () ;
  compute_X ( output , m_qgauss , x ) ;
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================


// ============================================================================
// Flat in 1D
// ============================================================================
Ostap::Models::Uniform::Uniform
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         )
  : RooAbsPdf  ( name , title ) 
  , m_dim      ( 1 ) 
  , m_x        ( "x" , "x-observable" , this , x )
{}
// ============================================================================
// Flat in 2D
// ============================================================================
Ostap::Models::Uniform::Uniform
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          y         )
  : RooAbsPdf  ( name , title ) 
  , m_dim      ( 2 ) 
  , m_x        ( "x" , "x-observable" , this , x )
  , m_y        ( "y" , "y-observable" , this , y )
{}
// ============================================================================
// Flat in 3D
// ============================================================================
Ostap::Models::Uniform::Uniform
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          y         ,
  RooAbsReal&          z         )
  : RooAbsPdf  ( name , title ) 
  , m_dim      ( 3 ) 
  , m_x        ( "x" , "x-observable" , this , x )
  , m_y        ( "y" , "y-observable" , this , y )
  , m_z        ( "z" , "z-observable" , this , z )
{}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Uniform::Uniform
( const Ostap::Models::Uniform& right ,
  const char*                   name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_dim      ( right.m_dim ) 
  , m_x        ( "!x" , this , right.m_x ) 
  , m_y        ( "!y" , this , right.m_y ) 
  , m_z        ( "!z" , this , right.m_z ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Uniform::~Uniform(){} 
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Uniform*
Ostap::Models::Uniform::clone ( const char* name ) const 
{ return new Ostap::Models::Uniform ( *this , name ) ; }
// ============================================================================
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Models::Uniform::evaluate() const { return 1 ; }
// ============================================================================
Int_t Ostap::Models::Uniform::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  //
  if      ( 3 == m_dim && matchArgs ( allVars , analVars , m_x , m_y , m_z ) ) { return 1 ; }
  else if ( 3 == m_dim && matchArgs ( allVars , analVars , m_x , m_z       ) ) { return 2 ; }
  else if ( 3 == m_dim && matchArgs ( allVars , analVars , m_y , m_z       ) ) { return 3 ; }
  else if ( 2 <= m_dim && matchArgs ( allVars , analVars , m_x , m_y       ) ) { return 4 ; }
  else if ( 3 == m_dim && matchArgs ( allVars , analVars , m_z             ) ) { return 5 ; }
  else if ( 2 <= m_dim && matchArgs ( allVars , analVars , m_y             ) ) { return 6 ; }
  else if ( 1 <= m_dim && matchArgs ( allVars , analVars , m_x             ) ) { return 7 ; }
  //
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Uniform::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  // 3D-itegral 
  if      ( 3 == m_dim && 1 == code ) 
  {
    return 
      ( m_x.max ( rangeName ) - m_x.min ( rangeName ) ) *
      ( m_y.max ( rangeName ) - m_y.min ( rangeName ) ) *
      ( m_z.max ( rangeName ) - m_z.min ( rangeName ) ) ;
  }
  // 2D-itegral: x,z
  else if ( 3 == m_dim && 2 == code ) 
  {
    return 
      ( m_x.max ( rangeName ) - m_x.min ( rangeName ) ) *
      ( m_z.max ( rangeName ) - m_z.min ( rangeName ) ) ;    
  }
  // 2D-itegral: y,z
  else if ( 3 == m_dim && 3 == code ) 
  {
    return 
      ( m_y.max ( rangeName ) - m_y.min ( rangeName ) ) *
      ( m_z.max ( rangeName ) - m_z.min ( rangeName ) ) ;    
  }
  // 2D-itegral: x,y
  else if ( 2 <= m_dim && 4 == code ) 
  {
    return 
      ( m_x.max ( rangeName ) - m_x.min ( rangeName ) ) *
      ( m_y.max ( rangeName ) - m_y.min ( rangeName ) ) ;    
  }
  // 1D-itegral: z 
  else if ( 3 == m_dim && 5 == code ) 
  {
    return 
      ( m_z.max ( rangeName ) - m_z.min ( rangeName ) ) ;    
  }
  // 1D-itegral: y 
  else if ( 2 <= m_dim && 6 == code ) 
  {
    return 
      ( m_y.max ( rangeName ) - m_y.min ( rangeName ) ) ;    
  }
  // 1D-itegral: x 
  else if ( 1 <= m_dim && 7 == code ) 
  {
    return 
      ( m_x.max ( rangeName ) - m_x.min ( rangeName ) ) ;    
  }
  //
  assert ( 1 > 2 ) ;
  //
  return 0 ;
}
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Uniform::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  //
  auto  x = m_x . getValBatch ( begin , batchSize ) ;
  if ( ( 1 == m_dim ) && x.empty () )                    { return {} ; }
  //
  if  ( 2 <= m_dim ) 
  {
    //
    auto  y = m_y . getValBatch ( begin , batchSize ) ;
    if ( ( 2 == m_dim ) && x .empty () && y.empty() )    { return {} ; }
    //
    if ( 3 == m_dim ) 
    {
      auto  z = m_z . getValBatch ( begin , batchSize ) ;
      if ( x.empty () && y.empty() && z.empty() )        { return {} ; }
    }
  }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  for ( auto& o : output ) { o = 1.0 ; }  
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================


// ============================================================================
ClassImp(Ostap::Models::Uniform            ) 
ClassImp(Ostap::Models::BreitWigner        ) 
ClassImp(Ostap::Models::BreitWignerMC      ) 
ClassImp(Ostap::Models::BWI  ) 
ClassImp(Ostap::Models::BW23L              ) 
ClassImp(Ostap::Models::Flatte             ) 
ClassImp(Ostap::Models::LASS               ) 
ClassImp(Ostap::Models::LASS23L            ) 
ClassImp(Ostap::Models::Bugg               ) 
ClassImp(Ostap::Models::Bugg23L            ) 
ClassImp(Ostap::Models::Voigt              ) 
ClassImp(Ostap::Models::PseudoVoigt        ) 
ClassImp(Ostap::Models::Swanson            ) 
ClassImp(Ostap::Models::CrystalBall        ) 
ClassImp(Ostap::Models::CrystalBallRS      ) 
ClassImp(Ostap::Models::CrystalBallDS      ) 
ClassImp(Ostap::Models::Needham            ) 
ClassImp(Ostap::Models::Apollonios         ) 
ClassImp(Ostap::Models::Apollonios2        ) 
ClassImp(Ostap::Models::BifurcatedGauss    ) 
ClassImp(Ostap::Models::GenGaussV1         ) 
ClassImp(Ostap::Models::GenGaussV2         ) 
ClassImp(Ostap::Models::SkewGauss          ) 
ClassImp(Ostap::Models::Bukin              ) 
ClassImp(Ostap::Models::StudentT           ) 
ClassImp(Ostap::Models::BifurcatedStudentT ) 
ClassImp(Ostap::Models::GramCharlierA      ) 
ClassImp(Ostap::Models::PhaseSpace2        ) 
ClassImp(Ostap::Models::PhaseSpaceLeft     ) 
ClassImp(Ostap::Models::PhaseSpaceRight    ) 
ClassImp(Ostap::Models::PhaseSpaceNL       ) 
ClassImp(Ostap::Models::PhaseSpacePol      ) 
ClassImp(Ostap::Models::PhaseSpaceLeftExpoPol ) 
ClassImp(Ostap::Models::PhaseSpace23L      ) 
ClassImp(Ostap::Models::PolyPositive       ) 
ClassImp(Ostap::Models::PolyPositiveEven   ) 
ClassImp(Ostap::Models::PolyMonotonic      ) 
ClassImp(Ostap::Models::PolyConvex         ) 
ClassImp(Ostap::Models::PolyConvexOnly     ) 
ClassImp(Ostap::Models::ExpoPositive       ) 
ClassImp(Ostap::Models::PolySigmoid        )
ClassImp(Ostap::Models::TwoExpoPositive    ) 
ClassImp(Ostap::Models::GammaDist          ) 
ClassImp(Ostap::Models::GenGammaDist       ) 
ClassImp(Ostap::Models::Amoroso            ) 
ClassImp(Ostap::Models::LogGammaDist       ) 
ClassImp(Ostap::Models::Log10GammaDist     ) 
ClassImp(Ostap::Models::LogGamma           )
ClassImp(Ostap::Models::BetaPrime          ) 
ClassImp(Ostap::Models::Landau             ) 
ClassImp(Ostap::Models::SinhAsinh          ) 
ClassImp(Ostap::Models::JohnsonSU          ) 
ClassImp(Ostap::Models::Atlas              ) 
ClassImp(Ostap::Models::Sech               ) 
ClassImp(Ostap::Models::Losev              ) 
ClassImp(Ostap::Models::Logistic           ) 
ClassImp(Ostap::Models::Argus              ) 
ClassImp(Ostap::Models::Slash              ) 
ClassImp(Ostap::Models::AsymmetricLaplace  ) 
ClassImp(Ostap::Models::Tsallis            ) 
ClassImp(Ostap::Models::QGSM               ) 
ClassImp(Ostap::Models::TwoExpos           ) 
ClassImp(Ostap::Models::DoubleGauss        ) 
ClassImp(Ostap::Models::Gumbel             )
ClassImp(Ostap::Models::Weibull            )
ClassImp(Ostap::Models::RaisingCosine      )
ClassImp(Ostap::Models::QGaussian          )
ClassImp(Ostap::Models::PositiveSpline     ) 
ClassImp(Ostap::Models::MonotonicSpline    ) 
ClassImp(Ostap::Models::ConvexOnlySpline   )
ClassImp(Ostap::Models::ConvexSpline       )
// ============================================================================
//                                                                      The END 
// ============================================================================
