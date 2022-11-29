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
#include "RooConstVar.h"
#include "RooGlobalFunc.h"
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
//  Shape1D
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Models::Shape1D::Shape1D
( const Ostap::Models::Shape1D& right ,
  const char*                   name  )
  : RooAbsPdf ( right , name ) 
  , m_x        ( "x"  , this , right.m_x ) 
  , m_function ( right.m_function ) 
{}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Shape1D*
Ostap::Models::Shape1D::clone ( const char* name ) const 
{ return new Ostap::Models::Shape1D(*this,name) ; }
// ============================================================================
//  Shape2D
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Models::Shape2D::Shape2D
( const Ostap::Models::Shape2D& right ,
  const char*                   name  )
  : RooAbsPdf ( right , name ) 
  , m_x        ( "x"  , this , right.m_x ) 
  , m_y        ( "y"  , this , right.m_y ) 
  , m_function ( right.m_function ) 
{}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Shape2D*
Ostap::Models::Shape2D::clone ( const char* name ) const 
{ return new Ostap::Models::Shape2D(*this,name) ; }
// ============================================================================
//  Shape3D
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Models::Shape3D::Shape3D
( const Ostap::Models::Shape3D& right ,
  const char*                   name  )
  : RooAbsPdf ( right , name ) 
  , m_x        ( "x"  , this , right.m_x ) 
  , m_y        ( "y"  , this , right.m_y ) 
  , m_z        ( "z"  , this , right.m_z ) 
  , m_function ( right.m_function ) 
{}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Shape3D*
Ostap::Models::Shape3D::clone ( const char* name ) const 
{ return new Ostap::Models::Shape3D(*this,name) ; }
// ============================================================================





// ============================================================================
//  Histo1D
// ============================================================================
// constructor
// ============================================================================
Ostap::Models::Histo1D::Histo1D
( const char*                 name  , 
  const char*                 title , 
  RooAbsReal&                 x     ,
  const Ostap::Math::Histo1D& histo ) 
  : RooAbsPdf  (  name ,  title ) 
  , m_x        ( "x"   , "x-variable" , this , x ) 
  , m_histo    ( histo ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Models::Histo1D::Histo1D
( const Ostap::Models::Histo1D& right ,
  const char*                   name  )
  : RooAbsPdf ( right , name  ) 
  , m_x       ( "x"  , this , right.m_x ) 
  , m_histo   ( right.m_histo ) 
{}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Histo1D*
Ostap::Models::Histo1D::clone ( const char* name ) const 
{ return new Ostap::Models::Histo1D(*this,name) ; }
// ============================================================================
//  Histo2D
// ============================================================================
// constructor
// ============================================================================
Ostap::Models::Histo2D::Histo2D
( const char*                 name  , 
  const char*                 title , 
  RooAbsReal&                 x     ,
  RooAbsReal&                 y     ,
  const Ostap::Math::Histo2D& histo ) 
  : RooAbsPdf  (  name ,  title ) 
  , m_x        ( "x"   , "x-variable" , this , x ) 
  , m_y        ( "y"   , "y-variable" , this , y ) 
  , m_histo    ( histo ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Models::Histo2D::Histo2D
( const Ostap::Models::Histo2D& right ,
  const char*                   name  )
  : RooAbsPdf ( right , name  ) 
  , m_x       ( "x"  , this , right.m_x ) 
  , m_y       ( "y"  , this , right.m_y ) 
  , m_histo   ( right.m_histo ) 
{}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Histo2D*
Ostap::Models::Histo2D::clone ( const char* name ) const 
{ return new Ostap::Models::Histo2D(*this,name) ; }
// ============================================================================
//  Histo3D
// ============================================================================
// constructor
// ============================================================================
Ostap::Models::Histo3D::Histo3D
( const char*                 name  , 
  const char*                 title , 
  RooAbsReal&                 x     ,
  RooAbsReal&                 y     ,
  RooAbsReal&                 z     ,
  const Ostap::Math::Histo3D& histo ) 
  : RooAbsPdf  (  name ,  title ) 
  , m_x        ( "x"   , "x-variable" , this , x ) 
  , m_y        ( "y"   , "y-variable" , this , y ) 
  , m_z        ( "z"   , "z-variable" , this , z ) 
  , m_histo    ( histo ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Models::Histo3D::Histo3D
( const Ostap::Models::Histo3D& right ,
  const char*                   name  )
  : RooAbsPdf ( right , name  ) 
  , m_x       ( "x"  , this , right.m_x ) 
  , m_y       ( "y"  , this , right.m_y ) 
  , m_z       ( "z"  , this , right.m_z ) 
  , m_histo   ( right.m_histo ) 
{}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Histo3D*
Ostap::Models::Histo3D::clone ( const char* name ) const 
{ return new Ostap::Models::Histo3D(*this,name) ; }
// ============================================================================



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
  , m_x      ( "x"  , "Observable" , this , x     ) 
  , m_mass   ( "m0" , "Peak"       , this , mass  ) 
  , m_widths ( "g"  , "Widths"     , this )
    //
  , m_bw ( std::make_unique<Ostap::Math::BreitWigner>( 0 , 1 , m1 , m2 , L ) )
{
  m_widths.add (  width ) ;
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
  , m_x      ( "x"  , "Observable" , this , x     ) 
  , m_mass   ( "m0" , "Peak"       , this , mass  ) 
  , m_widths ( "g"  , "Widths"     , this )
    //
  , m_bw ( std::make_unique<Ostap::Math::BreitWigner> ( 0 , 1 , m1 , m2 , L , rho ) ) 
{
  m_widths.add (  width ) ;
}
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BreitWigner::BreitWigner 
( const char*            name  , 
  const char*            title ,
  RooAbsReal&            x     ,
  RooAbsReal&            mass  ,
  RooAbsReal&            width ,
  const Ostap::Math::BW& bw    ) 
  : RooAbsPdf  ( name , title ) 
    //
  , m_x      ( "x"  , "Observable" , this , x     ) 
  , m_mass   ( "m0" , "Peak"       , this , mass  ) 
  , m_widths ( "g"  , "Widths"     , this )
    //
  , m_bw ( bw.clone() ) 
{
  m_widths.add (  width ) ;
}
// ============================================================================
// (protected) constructor from all parameters
// ============================================================================
Ostap::Models::BreitWigner::BreitWigner
( const char* name          , 
  const char* title         ,
  RooAbsReal& x             ,
  RooAbsReal& mass          , 
  RooArgList& widths        ,
  const Ostap::Math::BW& bw ) 
  : RooAbsPdf ( name , title ) 
    //
  , m_x        ( "x"  , "Observable" , this , x    ) 
  , m_mass     ( "m0" , "Peak"       , this , mass ) 
  , m_widths   ( "g"  , "Widths"     , this ) 
    //
  , m_bw ( bw.clone() ) 
{
  //
  ::copy_real   ( widths , m_widths , "Invalid width parameter!" ,
                  "Ostap::Models::BreitWigner" ) ;
  //
  Ostap::Assert ( ::size ( m_widths ) == m_bw->nChannels() , 
                  "Widths/#channels mismatch" , 
                  "Ostap::Models::BreitWigner" ) ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::BreitWigner::BreitWigner 
( const Ostap::Models::BreitWigner& right , 
  const char*                       name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x      ( "x"  , this , right.m_x      ) 
  , m_mass   ( "m0" , this , right.m_mass   ) 
  , m_widths ( "g"  , this , right.m_widths )
    //
  , m_bw     (  right.m_bw->clone()        ) 
{
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
  m_bw -> setM0    (                 m_mass     ) ;
  m_bw -> setGamma ( ::get_par ( 0 , m_widths ) ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::BreitWigner::evaluate() const 
{ setPars() ; return  ( *m_bw ) ( m_x ) ; }
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
  return m_bw->integral ( m_x.min ( rangeName ) , m_x.max ( rangeName ) ) ;
}
// ============================================================================
// get the amplitude 
// ============================================================================
std::complex<double> Ostap::Models::BreitWigner::amplitude () const
{ setPars () ; return m_bw->amplitude ( m_x ) ; }
// ============================================================================






// ============================================================================
// multi-channel BreiWigner
// ============================================================================
// constructor from all parameters
// ============================================================================
Ostap::Models::BreitWignerMC::BreitWignerMC
( const char*                                 name   , 
  const char*                                 title  ,
  RooAbsReal&                                 x      ,
  RooAbsReal&                                 mass   , 
  RooArgList&                                 widths ,
  const Ostap::Math::BreitWignerMC&           bwmc   ) 
  : BreitWigner ( name , title , x , mass , widths , bwmc ) 
{
  Ostap::Assert ( ::size ( widths )  == bwmc.nChannels() , 
                  "Invalid number of width-parameters"   , 
                  "Ostap::Models::BreitWignerMC"         ) ;
}
// ============================================================================
// "copy" constructor
// ============================================================================
Ostap::Models::BreitWignerMC::BreitWignerMC
( const Ostap::Models::BreitWignerMC& right, const char* name  )
  : BreitWigner ( right , name )
{}
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
  m_bw -> setM0 ( m_mass ) ;
  //
  Ostap::Math::BreitWignerMC* bwmc = (Ostap::Math::BreitWignerMC*) m_bw.get() ;
  //
  const unsigned short nc = bwmc -> nChannels() ;
  for ( unsigned short i = 0 ; i < nc ; ++i ) 
  { bwmc -> setGamma ( i , ::get_par ( i , m_widths ) ) ; }
  //
}
// ============================================================================
// access to underlying function
// ============================================================================
const Ostap::Math::BreitWignerMC&  
Ostap::Models::BreitWignerMC::breit_wigner_MC () const 
{ 
  setPars() ;
  const Ostap::Math::BreitWignerMC* bwmc = 
    (const Ostap::Math::BreitWignerMC*) m_bw.get() ;
  return *bwmc ; 
}
// ============================================================================


// ============================================================================
/// Breit-wigner with interference 
// ============================================================================
// constructor from Breit-Wigner and background 
// ============================================================================
Ostap::Models::BWI::BWI 
( const char*                       name      , 
  const char*                       title     , 
  const Ostap::Models::BreitWigner& bw        ,  
  RooAbsReal&                       magnitude , 
  RooAbsReal&                       phase     , 
  RooAbsReal&                       scale1    , 
  RooAbsReal&                       scale2    ) 
  : BreitWigner ( bw , name ) 
  , m_magnitude ( "!magnitude" , "background magnitude" , this , magnitude ) 
  , m_phase     ( "!phib"      , "background phase"     , this , phase     )
  , m_scale1    ( "!scale1"    , "background scale1"    , this , scale1    ) 
  , m_scale2    ( "!scale1"    , "background scale2"    , this , scale2    ) 
  // clone the original 
  , m_original ( bw.clone ( nullptr ) )
{
  SetTitle ( title ) ;  
}
// ============================================================================
// constructor from Breit-Wigner and background 
// ============================================================================
Ostap::Models::BWI::BWI 
( const char*                       name      , 
  const char*                       title     , 
  const Ostap::Models::BreitWigner& bw        ,
  RooAbsReal&                       magnitude , 
  RooAbsReal&                       phase     , 
  RooAbsReal&                       scale1    ) 
  : BWI ( name , title , bw , magnitude , phase , 
          scale1 , RooFit::RooConst ( 1 ) )
{}
// ============================================================================
// constructor from the Breit-Wigner and background 
// ============================================================================
Ostap::Models::BWI::BWI 
( const char*                       name      , 
  const char*                       title     , 
  const Ostap::Models::BreitWigner& bw        ,
  RooAbsReal&                       magnitude , 
  RooAbsReal&                       phase     ) 
  : BWI ( name , title , bw , magnitude , phase , 
          RooFit::RooConst ( 1 ) , RooFit::RooConst ( 1 ) ) 
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::BWI::BWI
( const Ostap::Models::BWI& right , 
  const char*               name  ) 
  : BreitWigner ( right , name )
    //
  , m_magnitude ( "!magnitude" , this , right.m_magnitude ) 
  , m_phase     ( "!phase"     , this , right.m_phase     )
  , m_scale1    ( "!scale1"    , this , right.m_scale1    ) 
  , m_scale2    ( "!scale2"    , this , right.m_scale2    ) 
    // clone the original 
  , m_original ( right.m_original->clone( nullptr ) )
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
// get the amplitude
// ============================================================================
#include <iostream> 
std::complex<double> Ostap::Models::BWI::amplitude () const
{
  // setPars () ;
  //
  const double x         = m_x ;
  //
  const double magnitude = m_magnitude ; // background magnitude 
  const double phase     = m_phase     ; // background phase 
  const double scale1    = m_scale1    ; // background magnitude scale1 
  const double scale2    = m_scale2    ; // background phase     scale2 
  //
  const double b         = magnitude * scale1 ;
  const double phi       = phase     * scale2 ;
  //
  const std::complex<double> ib  = 
    std::complex<double> ( std::cos ( phi ) , std::sin ( phi ) ) ;
  //
  const std::complex<double> amp = b * ib + 
    Ostap::Models::BreitWigner::amplitude() ;
  //
  return amp ;
}
// ============================================================================
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Models::BWI::evaluate () const 
{
  const std::complex<double> amp = this->amplitude() ;
  return m_bw->breit_wigner ( m_x , amp ) ;
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
Ostap::Models::Flatte::Flatte 
( const char*                name   , 
  const char*                title  ,
  RooAbsReal&                x      ,
  RooAbsReal&                m0     ,
  RooAbsReal&                g1     ,
  RooAbsReal&                g2     ,
  RooAbsReal&                g0     ,
  const Ostap::Math::Flatte& flatte ) 
  : BreitWigner ( name ,  title , x , m0 , g1 , flatte ) 
{
  m_widths.add ( g2 ) ;
  m_widths.add ( g0 ) ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Flatte::Flatte 
( const Ostap::Models::Flatte& right , 
  const char*                     name  ) 
  : BreitWigner ( right , name ) 
{}
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
  Ostap::Math::Flatte* flatte = (Ostap::Math::Flatte*) m_bw.get() ;
  //
  flatte -> setM0   ( m_mass      ) ;
  flatte -> setG1   ( ::get_par ( 0 , m_widths ) ) ;
  flatte -> setG2   ( ::get_par ( 1 , m_widths ) ) ;
  flatte -> setGam0 ( ::get_par ( 2 , m_widths ) ) ;
}
// ============================================================================
// access to underlying function
// ============================================================================
const Ostap::Math::Flatte& 
Ostap::Models::Flatte::flatte    () const 
{
  setPars () ;
  const Ostap::Math::Flatte* flatte = (const Ostap::Math::Flatte*) m_bw.get() ;
  return *flatte ;
}
// ============================================================================

// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::FlatteBugg::FlatteBugg 
( const char*                name   , 
  const char*                title  ,
  RooAbsReal&                x      ,
  RooAbsReal&                m0     ,
  RooAbsReal&                g1     ,
  RooAbsReal&                g2     ,
  RooAbsReal&                g0     ,
  const Ostap::Math::FlatteBugg& flatte ) 
  : BreitWigner ( name ,  title , x , m0 , g1 , flatte ) 
{
  m_widths.add ( g2 ) ;
  m_widths.add ( g0 ) ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::FlatteBugg::FlatteBugg 
( const Ostap::Models::FlatteBugg& right , 
  const char*                     name  ) 
  : BreitWigner ( right , name ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::FlatteBugg::~FlatteBugg (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::FlatteBugg*
Ostap::Models::FlatteBugg::clone( const char* name ) const 
{ return new Ostap::Models::FlatteBugg(*this,name) ; }
// ============================================================================
void Ostap::Models::FlatteBugg::setPars () const 
{
  //
  Ostap::Math::FlatteBugg* flatte = (Ostap::Math::FlatteBugg*) m_bw.get() ;
  //
  flatte -> setM0   ( m_mass      ) ;
  flatte -> setG1   ( ::get_par ( 0 , m_widths ) ) ;
  flatte -> setG2   ( ::get_par ( 1 , m_widths ) ) ;
  flatte -> setGam0 ( ::get_par ( 2 , m_widths ) ) ;
}
// ============================================================================
// access to underlying function
// ============================================================================
const Ostap::Math::FlatteBugg& 
Ostap::Models::FlatteBugg::flatte_bugg    () const 
{
  setPars () ;
  const Ostap::Math::FlatteBugg* flatte = 
    (const Ostap::Math::FlatteBugg*) m_bw.get() ;
  return *flatte ;
}
// ============================================================================












// ============================================================================
// constructor from all parameters
// ===========================================================================
Ostap::Models::LASS::LASS 
( const char*                name   ,
  const char*                title  ,
  RooAbsReal&                x      ,
  RooAbsReal&                m0     ,
  RooAbsReal&                g0     ,
  RooAbsReal&                a      ,
  RooAbsReal&                b      ,
  RooAbsReal&                e      , 
  const Ostap::Math::LASS&   lass   ) 
  : BreitWigner ( name ,  title , x , m0 , g0 , lass ) 
  , m_a      ( "a"     , "a-parameter" , this , a  ) 
  , m_b      ( "b"     , "b-parameter" , this , b  ) 
  , m_e      ( "e"     , "elasticity"  , this , e  ) 
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::LASS::LASS 
( const Ostap::Models::LASS& right , 
  const char*                     name  ) 
  : BreitWigner ( right  , name ) 
    //
  , m_a  ( "a"  , this , right.m_a  ) 
  , m_b  ( "b"  , this , right.m_b  ) 
  , m_e  ( "e"  , this , right.m_e  ) 
{}
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
  Ostap::Math::LASS* lass = (Ostap::Math::LASS*) m_bw.get() ;
  //
  lass -> setM0    ( m_mass       ) ;
  lass -> setGamma ( ::get_par ( 0 , m_widths ) ) ;
  lass -> setA     ( m_a          ) ;
  lass -> setB     ( m_b          ) ;
  lass -> setE     ( m_e          ) ;
}
// ============================================================================
// access to underlying function
// ============================================================================
const Ostap::Math::LASS& 
Ostap::Models::LASS::lass () const 
{
  setPars () ;
  const Ostap::Math::LASS* lass = (const Ostap::Math::LASS*) m_bw.get() ;
  return *lass ;
}
// ============================================================================



// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BWPS::BWPS 
( const char*              name  , 
  const char*              title ,
  RooAbsReal&              x     ,
  RooAbsReal&              m0    ,
  RooAbsReal&              gamma ,
  RooArgList&              phis  , 
  const Ostap::Math::BWPS& bwps  ) 
  : RooAbsPdf  (name ,title ) 
    //
  , m_x      ( "x"    , "Observable" , this , x     ) 
  , m_m0     ( "m0"   , "Peak"       , this , m0    ) 
  , m_gamma  ( "g"    , "Width"      , this )
  , m_phis   ( "phis" , "Phis"       , this )
    //
  , m_bwps   ( bwps ) 
{
  m_gamma.add ( gamma ) ;
  ::copy_real   ( phis , m_phis , "Invalid phis parameter!" ,
                  "Ostap::Models::BWPS" ) ;
  //
  Ostap::Assert ( ::size ( m_phis  ) == m_bwps.npars() , 
                  "#phis mismatch"      , 
                  "Ostap::Models::BWPS" ) ;
  Ostap::Assert ( ::size ( m_gamma ) == m_bwps.nChannels() , 
                  "#channels mismatch"      , 
                  "Ostap::Models::BWPS" ) ;
}

// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BWPS::BWPS 
( const char*              name  , 
  const char*              title ,
  RooAbsReal&              x     ,
  RooAbsReal&              m0    ,
  RooArgList&              gamma ,
  RooArgList&              phis  , 
  const Ostap::Math::BWPS& bwps  ) 
  : RooAbsPdf  (name ,title ) 
    //
  , m_x      ( "x"    , "Observable" , this , x  ) 
  , m_m0     ( "m0"   , "Peak"       , this , m0 ) 
  , m_gamma  ( "g"    , "Width"      , this )
  , m_phis   ( "phis" , "Phis"       , this )
    //
  , m_bwps   ( bwps ) 
{
  ::copy_real ( gamma , m_gamma , "Invalid gamma parameter!" , "Ostap::Models::BWPS" ) ;
  ::copy_real ( phis  , m_phis  , "Invalid phis  parameter!" , "Ostap::Models::BWPS" ) ;
  //
  Ostap::Assert ( ::size ( m_phis  ) == m_bwps.npars      () , 
                  "#phis mismatch"      , 
                  "Ostap::Models::BWPS" ) ;
  Ostap::Assert ( ::size ( m_gamma ) == m_bwps.nChannels  () , 
                  "#channels mismatch"      , 
                  "Ostap::Models::BWPS" ) ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::BWPS::BWPS
( const Ostap::Models::BWPS& right , 
  const char*                name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x      ( "x"    , this , right.m_x      ) 
  , m_m0     ( "m0"   , this , right.m_m0     ) 
  , m_gamma  ( "g"    , this , right.m_gamma  )
  , m_phis   ( "phis" , this , right.m_phis   )
    //
  , m_bwps   ( right.m_bwps )   
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::BWPS::~BWPS(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::BWPS*
Ostap::Models::BWPS::clone ( const char* name ) const 
{ return new Ostap::Models::BWPS(*this,name) ; }
// ============================================================================
void Ostap::Models::BWPS::setPars () const 
{
  //
  m_bwps.setM0    ( m_m0  ) ;
  //
  const unsigned short np = m_bwps.npars ();
  for ( unsigned short i = 0 ; i < np ; ++i ) 
  { m_bwps.setPar   ( i , ::get_par ( i , m_phis  ) ) ; }
  //
  const unsigned short nc = m_bwps.nChannels ();
  for ( unsigned short i = 0 ; i < nc ; ++i ) 
  { m_bwps.setGamma ( i , ::get_par ( i , m_gamma ) ) ; }
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::BWPS::evaluate() const 
{ setPars() ; return  m_bwps  ( m_x ) ; }
// ============================================================================
Int_t Ostap::Models::BWPS::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::BWPS::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  return m_bwps.integral ( m_x.min ( rangeName ) , m_x.max ( rangeName ) ) ;
}
// ============================================================================
// get the amplitude 
// ============================================================================
std::complex<double> 
Ostap::Models::BWPS::amplitude () const
{ setPars () ; return m_bwps.amplitude ( m_x ) ; }
// ============================================================================




// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BW3L::BW3L 
( const char*              name  , 
  const char*              title ,
  RooAbsReal&              x     ,
  RooAbsReal&              m0    ,
  RooAbsReal&              gamma ,
  const Ostap::Math::BW3L& bw3l  ) 
  : RooAbsPdf  (name ,title ) 
    //
  , m_x      ( "x"    , "Observable" , this , x     ) 
  , m_m0     ( "m0"   , "Peak"       , this , m0    ) 
  , m_gamma  ( "g"    , "Width"      , this )
    //
  , m_bw3l   ( bw3l ) 
{
  m_gamma.add ( gamma ) ;
  //
  Ostap::Assert ( ::size ( m_gamma ) == m_bw3l.nChannels() , 
                  "#channels mismatch"  , 
                  "Ostap::Models::BW3L" ) ;
}
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::BW3L::BW3L
( const char*              name  , 
  const char*              title ,
  RooAbsReal&              x     ,
  RooAbsReal&              m0    ,
  RooArgList&              gamma ,
  const Ostap::Math::BW3L& bw3l  ) 
  : RooAbsPdf  (name ,title ) 
    //
  , m_x      ( "x"    , "Observable" , this , x  ) 
  , m_m0     ( "m0"   , "Peak"       , this , m0 ) 
  , m_gamma  ( "g"    , "Width"      , this )
    //
  , m_bw3l   ( bw3l ) 
{
  ::copy_real ( gamma , m_gamma , "Invalid gamma parameter!" , "Ostap::Models::BWPS" ) ;
  Ostap::Assert ( ::size ( m_gamma ) == m_bw3l.nChannels  () , 
                  "#channels mismatch"      , 
                  "Ostap::Models::BW3L" ) ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::BW3L::BW3L
( const Ostap::Models::BW3L& right , 
  const char*                name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x      ( "x"    , this , right.m_x      ) 
  , m_m0     ( "m0"   , this , right.m_m0     ) 
  , m_gamma  ( "g"    , this , right.m_gamma  )
    //
  , m_bw3l   ( right.m_bw3l )   
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::BW3L::~BW3L(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::BW3L*
Ostap::Models::BW3L::clone ( const char* name ) const 
{ return new Ostap::Models::BW3L(*this,name) ; }
// ============================================================================
void Ostap::Models::BW3L::setPars () const 
{
  //
  m_bw3l.setM0    ( m_m0  ) ;
  //
  const unsigned short nc = m_bw3l.nChannels ();
  for ( unsigned short i = 0 ; i < nc ; ++i ) 
  { m_bw3l.setGamma ( i , ::get_par ( i , m_gamma ) ) ; }
  //
}
// ============================================================================
// get the amplitude 
// ============================================================================
std::complex<double> 
Ostap::Models::BW3L::amplitude () const
{ setPars () ; return m_bw3l.amplitude ( m_x ) ; }
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::BW3L::evaluate() const 
{ setPars() ; return  m_bw3l  ( m_x ) ; }
// ============================================================================
Int_t Ostap::Models::BW3L::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::BW3L::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars() ;
  return m_bw3l.integral ( m_x.min ( rangeName ) , m_x.max ( rangeName ) ) ;
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
  return m_voigt.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}




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


// ============================================================================
//         ExGauss
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::ExGauss::ExGauss
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          mu        , 
  RooAbsReal&          varsigma  , 
  RooAbsReal&          k         )  
  : RooAbsPdf ( name , title ) 
//
  , m_x        ( "!x"         , "Observable" , this , x        ) 
  , m_mu       ( "!mu"        , "mu"         , this , mu       ) 
  , m_varsigma ( "!varsigma"  , "varsigma"   , this , varsigma ) 
  , m_k        ( "!k"         , "k"          , this , k        ) 
  //
  , m_eg       () 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::ExGauss::ExGauss
( const Ostap::Models::ExGauss& right , 
  const char*                        name   ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x         ( "!x"        , this , right.m_x        ) 
  , m_mu        ( "!mu"       , this , right.m_mu       ) 
  , m_varsigma  ( "!varsigma" , this , right.m_varsigma ) 
  , m_k         ( "!k"        , this , right.m_k        ) 
    //
  , m_eg     ( right.m_eg ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::ExGauss::~ExGauss (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::ExGauss*
Ostap::Models::ExGauss::clone( const char* name ) const 
{ return new Ostap::Models::ExGauss ( *this , name ) ; }
// ============================================================================
void Ostap::Models::ExGauss::setPars () const 
{
  m_eg . setMu       ( m_mu       ) ;
  m_eg . setVarsigma ( m_varsigma ) ;
  m_eg . setK        ( m_k        ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::ExGauss::evaluate() const 
{
  setPars () ;
  return m_eg ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::ExGauss::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::ExGauss::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_eg.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================




// ============================================================================
//        NormalLaplace 
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::NormalLaplace::NormalLaplace
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          mu        , 
  RooAbsReal&          varsigma  , 
  RooAbsReal&          kL        ,  
  RooAbsReal&          kR        )  
  : RooAbsPdf ( name , title ) 
//
  , m_x        ( "!x"         , "Observable" , this , x        ) 
  , m_mu       ( "!mu"        , "mu"         , this , mu       ) 
  , m_varsigma ( "!varsigma"  , "varsigma"   , this , varsigma ) 
  , m_kL       ( "!kL"        , "k-left"     , this , kL       ) 
  , m_kR       ( "!kR"        , "k-right"    , this , kR       ) 
  //
  , m_nl       () 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::NormalLaplace::NormalLaplace
( const Ostap::Models::NormalLaplace& right , 
  const char*                        name   ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x         ( "!x"        , this , right.m_x        ) 
  , m_mu        ( "!mu"       , this , right.m_mu       ) 
  , m_varsigma  ( "!varsigma" , this , right.m_varsigma ) 
  , m_kL        ( "!kL"       , this , right.m_kL       ) 
  , m_kR        ( "!kR"       , this , right.m_kR       ) 
    //
  , m_nl        ( right.m_nl ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::NormalLaplace::~NormalLaplace(){}

// ============================================================================
// clone 
// ============================================================================
Ostap::Models::NormalLaplace*
Ostap::Models::NormalLaplace::clone( const char* name ) const 
{ return new Ostap::Models::NormalLaplace( *this , name ) ; }
// ============================================================================
void Ostap::Models::NormalLaplace::setPars () const 
{
  m_nl . setMu       ( m_mu       ) ;
  m_nl . setVarsigma ( m_varsigma ) ;
  m_nl . setKL       ( m_kL       ) ;
  m_nl . setKR       ( m_kR       ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::NormalLaplace::evaluate() const 
{
  setPars () ;
  return m_nl ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::NormalLaplace::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::NormalLaplace::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_nl.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
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




// ============================================================================
//         Novosibirsk
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Novosibirsk::Novosibirsk
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , 
  RooAbsReal&          peak      , 
  RooAbsReal&          sigma     , 
  RooAbsReal&          tau       ) 
  : RooAbsPdf ( name , title ) 
//
  , m_x       ( "!x"       , "Observable" , this , x      ) 
  , m_peak    ( "!peak"    , "peak"       , this , peak   ) 
  , m_sigma   ( "!sigma"   , "sigma"      , this , sigma  )
  , m_tau     ( "!tau"     , "tau"        , this , tau    )
//
  , m_novosibirsk ( 1 , 1 , 0 ) 
{
  //
  setPars () ;
  //
} 
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Novosibirsk::Novosibirsk
( const Ostap::Models::Novosibirsk& right  , 
  const char*                       name   ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "!x"      , this , right.m_x      ) 
  , m_peak   ( "!peak"   , this , right.m_peak   ) 
  , m_sigma  ( "!sigma"  , this , right.m_sigma  ) 
  , m_tau    ( "!tau"    , this , right.m_tau     ) 
//
  , m_novosibirsk ( right.m_novosibirsk ) 
{
  setPars() ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Novosibirsk::~Novosibirsk(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Novosibirsk*
Ostap::Models::Novosibirsk::clone( const char* name ) const 
{ return new Ostap::Models::Novosibirsk(*this,name) ; }
// ============================================================================
void Ostap::Models::Novosibirsk::setPars () const 
{
  m_novosibirsk.setPeak   ( m_peak  ) ;
  m_novosibirsk.setSigma  ( m_sigma ) ;
  m_novosibirsk.setTau    ( m_tau   ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Novosibirsk::evaluate() const 
{
  setPars () ;
  return m_novosibirsk   ( m_x     ) ;
}
// ============================================================================
Int_t Ostap::Models::Novosibirsk::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Novosibirsk::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  setPars() ;
  return m_novosibirsk.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
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
  return m_stt.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================





// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::PearsonIV::PearsonIV
( const char*          name     , 
  const char*          title    ,
  RooAbsReal&          x        ,
  RooAbsReal&          mu       ,
  RooAbsReal&          varsigma ,
  RooAbsReal&          n        ,
  RooAbsReal&          kappa    )
  : RooAbsPdf  (name , title ) 
//
  , m_x         ( "!x"        , "Observable" , this , x        ) 
  , m_mu        ( "!mu"       , "Peak"       , this , mu       ) 
  , m_varsigma  ( "!varsigma" , "Width"      , this , varsigma )
  , m_n         ( "!n"        , "N"          , this , n        )
  , m_kappa     ( "!kappa"    , "kappa"      , this , kappa    )
//
  , m_p4    () 
{
  setPars () ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::PearsonIV::PearsonIV
( const Ostap::Models::PearsonIV& right , 
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x         ( "!x"        , this , right.m_x        ) 
  , m_mu        ( "!mu"       , this , right.m_mu       ) 
  , m_varsigma  ( "!varsigma" , this , right.m_varsigma )
  , m_n         ( "!n"        , this , right.m_n        )
  , m_kappa     ( "!kappa"    , this , right.m_kappa    )
    //
  , m_p4 ( right.m_p4) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::PearsonIV::~PearsonIV (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PearsonIV*
Ostap::Models::PearsonIV::clone( const char* name ) const 
{ return new Ostap::Models::PearsonIV(*this,name) ; }
// ============================================================================
void Ostap::Models::PearsonIV::setPars () const 
{
  //
  m_p4.setMu       ( m_mu       ) ;
  m_p4.setVarsigma ( m_varsigma ) ;
  m_p4.setN        ( m_n        ) ;
  m_p4.setKappa    ( m_kappa    ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PearsonIV::evaluate() const 
{
  //
  setPars () ;
  //
  return m_p4 ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::PearsonIV::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::PearsonIV::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_p4.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
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
  RooAbsReal&          scale     ,
  const Ostap::Math::PhaseSpaceLeft& left ) 
  : RooAbsPdf ( name , title ) 
    //
  , m_x         ( "x"     , "Observable" , this , x          ) 
  , m_threshold ( "th"    , "Threshold"  , this , threshold  ) 
  , m_scale     ( "scale" , "Scale"      , this , scale      ) 
    //  
  , m_left ( left ) 
{
  setPars () ;  
}
// ============================================================================
// Left-edge of N-body phase space 
// ============================================================================
Ostap::Models::PhaseSpaceLeft::PhaseSpaceLeft
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          threshold , 
  const Ostap::Math::PhaseSpaceLeft& left ) 
 : PhaseSpaceLeft ( name , title , x , threshold , RooFit::RooConst ( 1.0 ) , left )
{}
// ============================================================================
// Left-edge of N-body phase space 
// ============================================================================
Ostap::Models::PhaseSpaceLeft::PhaseSpaceLeft
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          threshold ,
  RooAbsReal&          scale     ,
  const Ostap::Math::PhaseSpace2& left ) 
  : PhaseSpaceLeft ( name , title , x , threshold , scale , Ostap::Math::PhaseSpaceLeft ( left ) )
{}
// ============================================================================
// Left-edge of N-body phase space 
// ============================================================================
Ostap::Models::PhaseSpaceLeft::PhaseSpaceLeft
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          threshold ,
  const Ostap::Math::PhaseSpace2& left ) 
  : PhaseSpaceLeft ( name , title , x , threshold , RooFit::RooConst ( 1.0 )  , Ostap::Math::PhaseSpaceLeft ( left  ) )
{}
// ============================================================================
// Left-edge of N-body phase space 
// ============================================================================
Ostap::Models::PhaseSpaceLeft::PhaseSpaceLeft
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          threshold ,
  RooAbsReal&          scale     ,
  const Ostap::Math::PhaseSpace3& left ) 
  : PhaseSpaceLeft ( name , title , x , threshold , scale , Ostap::Math::PhaseSpaceLeft ( left ) )
{}
// ============================================================================
// Left-edge of N-body phase space 
// ============================================================================
Ostap::Models::PhaseSpaceLeft::PhaseSpaceLeft
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          threshold ,
  const Ostap::Math::PhaseSpace3& left ) 
  : PhaseSpaceLeft ( name , title , x , threshold , RooFit::RooConst ( 1.0 )  , Ostap::Math::PhaseSpaceLeft ( left  ) )
{}
// ============================================================================
// Left-edge of N-body phase space 
// ============================================================================
Ostap::Models::PhaseSpaceLeft::PhaseSpaceLeft
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          threshold ,
  RooAbsReal&          scale     ,
  const Ostap::Math::PhaseSpace3s& left ) 
  : PhaseSpaceLeft ( name , title , x , threshold , scale , Ostap::Math::PhaseSpaceLeft ( left ) )
{}
// ============================================================================
// Left-edge of N-body phase space 
// ============================================================================
Ostap::Models::PhaseSpaceLeft::PhaseSpaceLeft
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          threshold ,
  const Ostap::Math::PhaseSpace3s& left ) 
  : PhaseSpaceLeft ( name , title , x , threshold , RooFit::RooConst ( 1.0 )  , Ostap::Math::PhaseSpaceLeft ( left  ) )
{}
// ============================================================================
// Left-edge of N-body phase space 
// ============================================================================
Ostap::Models::PhaseSpaceLeft::PhaseSpaceLeft
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          threshold ,
  RooAbsReal&          scale     ,
  const unsigned short N         )
  : PhaseSpaceLeft ( name , title , x , threshold , scale , Ostap::Math::PhaseSpaceLeft ( 10 , N ) )
{}
// ============================================================================
// Left-edge of N-body phase space 
// ============================================================================
Ostap::Models::PhaseSpaceLeft::PhaseSpaceLeft
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,
  RooAbsReal&          threshold ,
  const unsigned short N         )
  : PhaseSpaceLeft ( name , title , x , threshold , RooFit::RooConst ( 1.0 )  , Ostap::Math::PhaseSpaceLeft ( 10 , N ) )
{}
// ============================================================================
// "copy constructor"
// ============================================================================
Ostap::Models::PhaseSpaceLeft::PhaseSpaceLeft
( const Ostap::Models::PhaseSpaceLeft& right , const char* name )  
  : RooAbsPdf ( right , name )
    //
  , m_x         ( "x"     , this , right.m_x         ) 
  , m_threshold ( "tr"    , this , right.m_threshold ) 
  , m_scale     ( "scale" , this , right.m_scale     ) 
    //
  , m_left      ( right.m_left ) 
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
  m_left.setScale     ( m_scale     ) ;  
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


// ============================================================================
// Two-body phase space from 3-body decays 
// ============================================================================
Ostap::Models::PhaseSpace23L::PhaseSpace23L
( const char*                       name      , 
  const char*                       title     ,
  RooAbsReal&                       x         ,
  const Ostap::Kinematics::Dalitz&  dalitz    ,  
  const unsigned short              L         , 
  const unsigned short              l         ) 
  : RooAbsPdf ( name , title ) 
  , m_x       ( "x" , "Observable" , this , x ) 
  , m_ps23L   ( dalitz , L , l ) 
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
//  PhaseSpaceLeft x expo x pol 
// ============================================================================
Ostap::Models::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol 
( const char*                        name      , 
  const char*                        title     ,
  RooRealVar&                        x         ,
  const Ostap::Math::PhaseSpace2&    ps        , 
  RooAbsReal&                        tau       ,
  RooAbsReal&                        scale     ,
  RooArgList&                        phis      )
  : PhaseSpaceLeftExpoPol  ( name , title , x , 
                             Ostap::Math::PhaseSpaceLeft ( ps ) , 
                             tau  , scale , phis )
{}
// ============================================================================
//  PhaseSpaceLeft x expo x pol 
// ============================================================================
Ostap::Models::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol 
( const char*                        name      , 
  const char*                        title     ,
  RooRealVar&                        x         ,
  const Ostap::Math::PhaseSpace3&    ps        , 
  RooAbsReal&                        tau       ,
  RooAbsReal&                        scale     ,
  RooArgList&                        phis      )
  : PhaseSpaceLeftExpoPol  ( name , title , x , 
                             Ostap::Math::PhaseSpaceLeft ( ps ) , 
                             tau  , scale , phis )
{}
// ============================================================================
//  PhaseSpaceLeft x expo x pol 
// ============================================================================
Ostap::Models::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol 
( const char*                        name      , 
  const char*                        title     ,
  RooRealVar&                        x         ,
  const Ostap::Math::PhaseSpace3s&   ps        , 
  RooAbsReal&                        tau       ,
  RooAbsReal&                        scale     ,
  RooArgList&                        phis      )
  : PhaseSpaceLeftExpoPol  ( name , title , x , 
                             Ostap::Math::PhaseSpaceLeft ( ps ) , 
                             tau  , scale , phis )
{}
// ============================================================================
//  PhaseSpaceLeft x expo x pol 
// ============================================================================
Ostap::Models::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol 
( const char*                        name      , 
  const char*                        title     ,
  RooRealVar&                        x         ,
  const Ostap::Math::PhaseSpaceNL&   ps        , 
  RooAbsReal&                        tau       ,
  RooAbsReal&                        scale     ,
  RooArgList&                        phis      )
  : PhaseSpaceLeftExpoPol  ( name , title , x , 
                             Ostap::Math::PhaseSpaceLeft ( ps ) , 
                             tau  , scale , phis )
{}
// ============================================================================
//  PhaseSpaceLeft x expo x pol 
// ============================================================================
Ostap::Models::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol 
( const char*                        name      , 
  const char*                        title     ,
  RooRealVar&                        x         ,
  const unsigned short               ps        , 
  RooAbsReal&                        tau       ,
  RooAbsReal&                        scale     ,
  RooArgList&                        phis      )
  : PhaseSpaceLeftExpoPol  ( name , title , x , 
                             Ostap::Math::PhaseSpaceLeft ( 10.0 , ps ) , 
                             tau  , scale , phis )
{}
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
  ::set_pars ( m_phis , m_positive ) ;
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
{ ::set_pars ( m_phis , m_even ) ; }
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
{ ::set_pars ( m_phis , m_monotonic ) ; }
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
{ ::set_pars ( m_phis , m_convex ) ; }
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
{ ::set_pars ( m_phis , m_convex ) ; }
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
{ ::set_pars ( m_phis , m_spline ) ; }
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
{  ::set_pars ( m_phis , m_spline ) ; }
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
{ ::set_pars ( m_phis , m_spline ) ; }
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
{ ::set_pars ( m_phis , m_spline ) ; }
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



// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Argus::Argus
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          mu     ,
  RooAbsReal&          c      ,
  RooAbsReal&          chi    )
  : RooAbsPdf ( name , title ) 
    //
  , m_x       ( "!x"   , "Observable" , this , x      ) 
  , m_mu      ( "!mu"  , "mu"         , this , mu     ) 
  , m_c       ( "!c"   , "c"          , this , c      ) 
  , m_chi     ( "!chi" , "chi"        , this , chi    ) 
    //
  , m_argus  ( 1 , 1 , 0 ) 
{
  setPars() ;
}
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Argus::Argus
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          c      ,
  RooAbsReal&          chi    )
  : Argus ( name , title , x , c , c , chi )
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Argus::Argus
( const Ostap::Models::Argus&  right ,
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "!x"      , this , right.m_x      ) 
  , m_mu     ( "!mu"     , this , right.m_mu     )
  , m_c      ( "!c"      , this , right.m_c      )
  , m_chi    ( "!chi"    , this , right.m_chi    )
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
  m_argus.setMu  ( m_mu     ) ;
  m_argus.setC   ( m_c      ) ;
  m_argus.setChi ( m_chi    ) ;
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



// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::GenArgus::GenArgus
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          mu     ,
  RooAbsReal&          c      ,
  RooAbsReal&          chi    ,
  RooAbsReal&          dp     )
  : RooAbsPdf ( name , title ) 
    //
  , m_x       ( "!x"   , "Observable" , this , x      ) 
  , m_mu      ( "!mu"  , "mu"         , this , mu     ) 
  , m_c       ( "!c"   , "c"          , this , c      ) 
  , m_chi     ( "!chi" , "chi"        , this , chi    ) 
  , m_dp      ( "!do"  , "dp"         , this , dp    ) 
    //
  , m_argus  ( 1 , 1 , 0 , 1.5 ) 
{
  setPars() ;
}
// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::GenArgus::GenArgus
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,
  RooAbsReal&          c      ,
  RooAbsReal&          chi    ,
  RooAbsReal&          dp     )
  : GenArgus ( name , title , x , c , c , chi , dp )
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::GenArgus::GenArgus
( const Ostap::Models::GenArgus&  right ,
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
//
  , m_x      ( "!x"      , this , right.m_x      ) 
  , m_mu     ( "!mu"     , this , right.m_mu     )
  , m_c      ( "!c"      , this , right.m_c      )
  , m_chi    ( "!chi"    , this , right.m_chi    )
  , m_dp     ( "!dp"     , this , right.m_dp     )
//
  , m_argus  (                   right.m_argus  ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::GenArgus::~GenArgus () {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::GenArgus*
Ostap::Models::GenArgus::clone( const char* name ) const 
{ return new Ostap::Models::GenArgus ( *this , name) ; }
// ============================================================================
void Ostap::Models::GenArgus::setPars () const 
{
  //
  m_argus.setMu  ( m_mu     ) ;
  m_argus.setC   ( m_c      ) ;
  m_argus.setChi ( m_chi    ) ;
  m_argus.setDp  ( m_dp     ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::GenArgus::evaluate() const 
{
  //
  setPars () ;
  //
  return m_argus ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::GenArgus::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::GenArgus::analyticalIntegral 
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



// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Hat::Hat
( const char*          name     , 
  const char*          title    ,
  RooAbsReal&          x        ,
  RooAbsReal&          mu       ,
  RooAbsReal&          varsigma ) 
  : RooAbsPdf ( name , title ) 
    //
  , m_x        ( "!x"        , "observable" , this , x        ) 
  , m_mu       ( "!mu"       , "location"   , this , mu       ) 
  , m_varsigma ( "!varsigma" , "scale"      , this , varsigma ) 
    //
  , m_hat ( 0 , 1  ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Hat::Hat
( const Ostap::Models::Hat&  right ,
  const char*                     name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x         ( "!x"        , this , right.m_x        ) 
  , m_mu        ( "!mu"       , this , right.m_mu       )
  , m_varsigma  ( "!varsigma" , this , right.m_varsigma )
    //
  , m_hat  ( right.m_hat ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Hat::~Hat(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Hat*
Ostap::Models::Hat::clone( const char* name ) const 
{ return new Ostap::Models::Hat( *this , name) ; }
// ============================================================================
void Ostap::Models::Hat::setPars () const 
{
  m_hat.setMu       ( m_mu       ) ;
  m_hat.setVarsigma ( m_varsigma ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Hat::evaluate() const 
{
  setPars () ;
  return m_hat ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Hat::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Hat::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_hat.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::Up::Up
( const char*          name     , 
  const char*          title    ,
  RooAbsReal&          x        ,
  RooAbsReal&          mu       ,
  RooAbsReal&          varsigma ) 
  : RooAbsPdf ( name , title ) 
    //
  , m_x        ( "!x"        , "observable" , this , x        ) 
  , m_mu       ( "!mu"       , "location"   , this , mu       ) 
  , m_varsigma ( "!varsigma" , "scale"      , this , varsigma ) 
    //
  , m_up ( 0 , 1  ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::Up::Up
( const Ostap::Models::Up& right ,
  const char*              name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x         ( "!x"        , this , right.m_x        ) 
  , m_mu        ( "!mu"       , this , right.m_mu       )
  , m_varsigma  ( "!varsigma" , this , right.m_varsigma )
    //
  , m_up  ( right.m_up ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::Up::~Up(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Up*
Ostap::Models::Up::clone( const char* name ) const 
{ return new Ostap::Models::Up( *this , name) ; }
// ============================================================================
void Ostap::Models::Up::setPars () const 
{
  m_up.setMu       ( m_mu       ) ;
  m_up.setVarsigma ( m_varsigma ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Up::evaluate() const 
{
  setPars () ;
  return m_up ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Up::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Up::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_up.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Models::FupN::FupN
( const char*          name     , 
  const char*          title    ,
  RooAbsReal&          x        ,
  const unsigned short N        ,
  RooAbsReal&          mu       ,
  RooAbsReal&          varsigma ) 
  : RooAbsPdf ( name , title ) 
    //
  , m_x        ( "!x"        , "observable" , this , x        ) 
  , m_mu       ( "!mu"       , "location"   , this , mu       ) 
  , m_varsigma ( "!varsigma" , "scale"      , this , varsigma ) 
    //
  , m_fupN ( N , 0 , 1  ) 
{
  setPars() ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::Models::FupN::FupN
( const Ostap::Models::FupN& right ,
  const char*                name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_x         ( "!x"        , this , right.m_x        ) 
  , m_mu        ( "!mu"       , this , right.m_mu       )
  , m_varsigma  ( "!varsigma" , this , right.m_varsigma )
    //
  , m_fupN ( right.m_fupN ) 
{
  setPars () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::FupN::~FupN(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::FupN*
Ostap::Models::FupN::clone( const char* name ) const 
{ return new Ostap::Models::FupN( *this , name) ; }
// ============================================================================
void Ostap::Models::FupN::setPars () const 
{
  m_fupN.setMu       ( m_mu       ) ;
  m_fupN.setVarsigma ( m_varsigma ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::FupN::evaluate() const 
{
  setPars () ;
  return m_fupN ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::FupN::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::FupN::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( code == 1 ) ;
  if ( 1 != code ) {}
  //
  setPars () ;
  return m_fupN.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
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
  return m_2expos.integral ( m_x.min(rangeName) , m_x.max(rangeName) ) ;
}
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
  RooAbsReal&          scale     , // scale
  RooAbsReal&          q         ) // q
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"      , "Observable"               , this , x     ) 
  , m_mean     ( "mean"   , "Mean/location parameter"  , this , mean  ) 
  , m_scale    ( "scale"  , "Scale parameter"          , this , scale ) 
  , m_q        ( "q"      , "Q-parameter"              , this , q     ) 
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
  , m_x        ( "!x"     , this , right.m_x     ) 
  , m_mean     ( "!mean"  , this , right.m_mean  ) 
  , m_scale    ( "!scale" , this , right.m_scale ) 
  , m_q        ( "!q"     , this , right.m_q     ) 
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


// ============================================================================
/*  constructor from all parameters
 *  @param  x      the variable 
 *  @param  mean   the mean/mode/median/location 
 *  @param  kappa  the kappa-value 
 *  @param  scale  the scale parameter 
 */
// ============================================================================
Ostap::Models::KGaussian::KGaussian
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         , // observable 
  RooAbsReal&          mean      , // mean
  RooAbsReal&          scale     , // scale
  RooAbsReal&          kappa     ) // kappa
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "!x"      , "Observable"               , this , x     ) 
  , m_mean     ( "!mean"   , "Mean/location parameter"  , this , mean  ) 
  , m_scale    ( "!scale"  , "Scale parameter"          , this , scale ) 
  , m_kappa    ( "!kappa"  , "kappa-parameter"          , this , kappa ) 
  , m_kgauss   ()  
{
  setPars () ;  
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::KGaussian::KGaussian
( const Ostap::Models::KGaussian& right ,
  const char*                     name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "!x"     , this , right.m_x     ) 
  , m_mean     ( "!mean"  , this , right.m_mean  ) 
  , m_scale    ( "!scale" , this , right.m_scale ) 
  , m_kappa    ( "!kappa" , this , right.m_kappa ) 
  , m_kgauss   ( right.m_kgauss ) 
{
  setPars () ;  
}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::KGaussian*
Ostap::Models::KGaussian::clone( const char* name ) const 
{ return new Ostap::Models::KGaussian(*this,name) ; }
// ============================================================================
void Ostap::Models::KGaussian::setPars () const 
{
  m_kgauss.setMean  ( m_mean  ) ;
  m_kgauss.setKappa ( m_kappa ) ;
  m_kgauss.setScale ( m_scale ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::KGaussian::evaluate() const 
{
  setPars() ;
  return m_kgauss ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::KGaussian::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::KGaussian::analyticalIntegral
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
  return m_kgauss.integral ( xmin , xmax ) ;
}
// ============================================================================



// ============================================================================
/*  constructor from all parameters
 *  @param name  name of PDF
 *  @param title name of PDF
 *  @param x     observable 
 *  @param mu    related to location 
 *  @param beta  related to asymmetry
 *  @param gamma related to width 
 *  @param delta related to width 
 */
// ============================================================================
Ostap::Models::Hyperbolic::Hyperbolic 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,   // observable 
  RooAbsReal&          mu        ,   // location
  RooAbsReal&          sigma     ,   // respondible for asymmetry
  RooAbsReal&          zeta      ,   // 
  RooAbsReal&          kappa     )   // related to width 
  : RooAbsPdf    ( name , title ) 
  , m_x          ( "x"      , "Observable"            , this , x     ) 
  , m_mu         ( "mu"     , "Location parameter"    , this , mu    ) 
  , m_sigma      ( "sigma"  , "Sigma    parameter"    , this , sigma ) 
  , m_zeta       ( "zeta"   , "Zeta     parameter"    , this , zeta  ) 
  , m_kappa      ( "kappa"  , "Kappa    parameter"    , this , kappa ) 
  , m_hyperbolic ()  
{
  setPars () ;  
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Hyperbolic::Hyperbolic 
( const Ostap::Models::Hyperbolic& right ,
  const char*                      name  ) 
  : RooAbsPdf    ( right , name ) 
    //
  , m_x          ( "x"     , this , right.m_x     ) 
  , m_mu         ( "mu"    , this , right.m_mu    )  
  , m_sigma      ( "sigma" , this , right.m_sigma ) 
  , m_zeta       ( "zeta"  , this , right.m_zeta  ) 
  , m_kappa      ( "kappa" , this , right.m_kappa ) 
  , m_hyperbolic ( right.m_hyperbolic ) 
{
  setPars () ;  
}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Hyperbolic*
Ostap::Models::Hyperbolic::clone ( const char* name ) const 
{ return new Ostap::Models::Hyperbolic ( *this , name ) ; }
// ============================================================================
void Ostap::Models::Hyperbolic::setPars () const 
{
  m_hyperbolic.setMu    ( m_mu    ) ;
  m_hyperbolic.setSigma ( m_sigma ) ;
  m_hyperbolic.setZeta  ( m_zeta  ) ;
  m_hyperbolic.setKappa ( m_kappa ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Hyperbolic::evaluate() const 
{
  setPars() ;
  return m_hyperbolic ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Hyperbolic::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Hyperbolic::analyticalIntegral
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
  return m_hyperbolic.integral ( xmin , xmax ) ;
}
// ============================================================================




// ============================================================================
/*  constructor from all parameters
 *  @param name  name of PDF
 *  @param title name of PDF
 *  @param x      observable 
 *  @param mu     related to location 
 *  @param sigma  related to width
 *  @param zeta   related to shape 
 *  @param kappa  related to asymmetry 
 *  @param lambda related to shape 
 */
// ============================================================================
Ostap::Models::GenHyperbolic::GenHyperbolic 
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         ,   // observable 
  RooAbsReal&          mu        ,   // location
  RooAbsReal&          sigma     ,   // related to width 
  RooAbsReal&          zeta      ,   // related to shape 
  RooAbsReal&          kappa     ,   // related to asymmetry 
  RooAbsReal&          lambda    )   // related to shapoe 
  : RooAbsPdf    ( name , title ) 
  , m_x          ( "x"      , "Observable"            , this , x      ) 
  , m_mu         ( "mu"     , "Location parameter"    , this , mu     ) 
  , m_sigma      ( "sigma"  , "Sigma    parameter"    , this , sigma  ) 
  , m_zeta       ( "zeta"   , "Zeta     parameter"    , this , zeta   ) 
  , m_kappa      ( "kappa"  , "Kappa    parameter"    , this , kappa  ) 
  , m_lambda     ( "lambda" , "Lambda   parameter"    , this , lambda ) 
  , m_hyperbolic ()  
{
  setPars () ;  
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::GenHyperbolic::GenHyperbolic 
( const Ostap::Models::GenHyperbolic& right ,
  const char*                      name  ) 
  : RooAbsPdf    ( right , name ) 
    //
  , m_x          ( "x"      , this , right.m_x      ) 
  , m_mu         ( "mu"     , this , right.m_mu     )  
  , m_sigma      ( "sigma"  , this , right.m_sigma  ) 
  , m_zeta       ( "zeta"   , this , right.m_zeta   ) 
  , m_kappa      ( "kappa"  , this , right.m_kappa  ) 
  , m_lambda     ( "lambda" , this , right.m_lambda ) 
  , m_hyperbolic ( right.m_hyperbolic ) 
{
  setPars () ;  
}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::GenHyperbolic*
Ostap::Models::GenHyperbolic::clone ( const char* name ) const 
{ return new Ostap::Models::GenHyperbolic ( *this , name ) ; }
// ============================================================================
void Ostap::Models::GenHyperbolic::setPars () const 
{ 
  m_hyperbolic.setMu     ( m_mu     ) ;
  m_hyperbolic.setSigma  ( m_sigma  ) ;
  m_hyperbolic.setZeta   ( m_zeta   ) ;
  m_hyperbolic.setKappa  ( m_kappa  ) ;
  m_hyperbolic.setLambda ( m_lambda ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::GenHyperbolic::evaluate() const 
{
  setPars() ;
  return m_hyperbolic ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::GenHyperbolic::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::GenHyperbolic::analyticalIntegral
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
  return m_hyperbolic.integral ( xmin , xmax ) ;
}
// ============================================================================



// ============================================================================
/** constructor from all parameters
 *  @param name  name of PDF
 *  @param title name of PDF
 *  @param x      observable 
 *  @param mu     location 
 *  @param sigma  width for Gaussian core 
 *  @param kL     left tail 
 *  @param kR     right tail 
 */
// ============================================================================
Ostap::Models::Das::Das 
( const char*          name   , 
  const char*          title  ,
  RooAbsReal&          x      ,   // observable 
  RooAbsReal&          mu     ,   // location
  RooAbsReal&          sigma  ,   // width 
  RooAbsReal&          kL     ,   // left tail 
  RooAbsReal&          kR     )   // right tail 
  : RooAbsPdf    ( name , title ) 
  , m_x          ( "x"      , "Observable"            , this , x      ) 
  , m_mu         ( "mu"     , "Location parameter"    , this , mu     ) 
  , m_sigma      ( "sigma"  , "Sigma    parameter"    , this , sigma  ) 
  , m_kL         ( "kL"     , "Left tail"             , this , kL     ) 
  , m_kR         ( "kR"     , "Right tail"            , this , kR     ) 
  , m_das ()  
{
  setPars () ;  
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Das::Das
( const Ostap::Models::Das& right ,
  const char*                      name  ) 
  : RooAbsPdf    ( right , name ) 
    //
  , m_x          ( "x"      , this , right.m_x      ) 
  , m_mu         ( "mu"     , this , right.m_mu     )  
  , m_sigma      ( "sigma"  , this , right.m_sigma  ) 
  , m_kL         ( "kL"     , this , right.m_kL     ) 
  , m_kR         ( "kR"     , this , right.m_kR     ) 
  , m_das        ( right.m_das ) 
{
  setPars () ;  
}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Das*
Ostap::Models::Das::clone ( const char* name ) const 
{ return new Ostap::Models::Das( *this , name ) ; }
// ============================================================================
void Ostap::Models::Das::setPars () const 
{ 
  m_das.setMu     ( m_mu     ) ;
  m_das.setSigma  ( m_sigma  ) ;
  m_das.setKL     ( m_kL     ) ;
  m_das.setKR     ( m_kR     ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Das::evaluate() const 
{
  setPars() ;
  return m_das ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Das::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Das::analyticalIntegral
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
  return m_das.integral ( xmin , xmax ) ;
}
// ============================================================================




// ============================================================================
Ostap::Models::CutOffGauss::CutOffGauss 
( const char* name  , 
  const char* title ,
  RooAbsReal& x     , // observable 
  const bool  right , 
  RooAbsReal& x0    , 
  RooAbsReal& sigma ) 
  : RooAbsPdf ( name , title )
  , m_x          ( "x"      , "Observable"            , this , x     ) 
  , m_x0         ( "x0"     , "Threshold parameter"   , this , x0    ) 
  , m_sigma      ( "sigma"  , "Sigma    parameter"    , this , sigma ) 
  , m_cutoff     ( right    ) 
{
  setPars () ;  
}
// ============================================================================
Ostap::Models::CutOffGauss::CutOffGauss 
( const char* name  , 
  const char* title ,
  RooAbsReal& x     , // observable 
  RooAbsReal& x0    , 
  RooAbsReal& sigma , 
  const Ostap::Math::CutOffGauss& cutoff ) 
  : RooAbsPdf ( name , title )
  , m_x          ( "x"      , "Observable"            , this , x     ) 
  , m_x0         ( "x0"     , "Threshold parameter"   , this , x0    ) 
  , m_sigma      ( "sigma"  , "Sigma    parameter"    , this , sigma ) 
  , m_cutoff     ( cutoff ) 
{
  setPars () ;  
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::CutOffGauss::CutOffGauss 
( const Ostap::Models::CutOffGauss& right ,
  const char*                       name  ) 
  : RooAbsPdf    ( right , name ) 
    //
  , m_x          ( "x"     , this , right.m_x     ) 
  , m_x0         ( "x0"    , this , right.m_x0    )  
  , m_sigma      ( "sigma" , this , right.m_sigma ) 
  , m_cutoff     ( right.m_cutoff ) 
{
  setPars () ;  
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Models::CutOffGauss::~CutOffGauss(){} 
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::CutOffGauss*
Ostap::Models::CutOffGauss::clone( const char* name ) const 
{ return new Ostap::Models::CutOffGauss(*this,name) ; }
// ============================================================================
void Ostap::Models::CutOffGauss::setPars () const 
{
  m_cutoff.setX0    ( m_x0    ) ;
  m_cutoff.setSigma ( m_sigma ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::CutOffGauss::evaluate() const 
{
  setPars() ;
  return m_cutoff ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::CutOffGauss::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::CutOffGauss::analyticalIntegral
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
  return m_cutoff.integral ( xmin , xmax ) ;
}
// ============================================================================



// ============================================================================
Ostap::Models::CutOffStudent::CutOffStudent
( const char* name  , 
  const char* title ,
  RooAbsReal& x     , // observable 
  const bool  right , 
  RooAbsReal& x0    , 
  RooAbsReal& nu    , 
  RooAbsReal& sigma ) 
  : RooAbsPdf ( name , title )
  , m_x          ( "x"      , "Observable"            , this , x     ) 
  , m_x0         ( "x0"     , "Threshold parameter"   , this , x0    ) 
  , m_nu         ( "nu"     , "Power     parameter"   , this , nu    ) 
  , m_sigma      ( "sigma"  , "Sigma     parameter"   , this , sigma ) 
  , m_cutoff     ( right    ) 
{
  setPars () ;  
}
// ============================================================================
Ostap::Models::CutOffStudent::CutOffStudent
( const char* name  , 
  const char* title ,
  RooAbsReal& x     , // observable 
  RooAbsReal& x0    , 
  RooAbsReal& nu    , 
  RooAbsReal& sigma , 
  const Ostap::Math::CutOffStudent& cutoff ) 
  : RooAbsPdf ( name , title )
  , m_x          ( "x"      , "Observable"            , this , x     ) 
  , m_x0         ( "x0"     , "Threshold parameter"   , this , x0    ) 
  , m_nu         ( "nu"     , "Power     parameter"   , this , nu    ) 
  , m_sigma      ( "sigma"  , "Sigma     parameter"   , this , sigma ) 
  , m_cutoff     ( cutoff ) 
{
  setPars () ;  
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::CutOffStudent::CutOffStudent
( const Ostap::Models::CutOffStudent& right ,
  const char*                         name  ) 
  : RooAbsPdf    ( right , name ) 
    //
  , m_x          ( "x"     , this , right.m_x     ) 
  , m_x0         ( "x0"    , this , right.m_x0    )  
  , m_nu         ( "nu"    , this , right.m_nu    )  
  , m_sigma      ( "sigma" , this , right.m_sigma ) 
  , m_cutoff     ( right.m_cutoff ) 
{
  setPars () ;  
}
// ============================================================================
// desctructor 
// ============================================================================
Ostap::Models::CutOffStudent::~CutOffStudent(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::CutOffStudent*
Ostap::Models::CutOffStudent::clone( const char* name ) const 
{ return new Ostap::Models::CutOffStudent(*this,name) ; }
// ============================================================================
void Ostap::Models::CutOffStudent::setPars () const 
{
  m_cutoff.setX0    ( m_x0    ) ;
  m_cutoff.setNu    ( m_nu    ) ;
  m_cutoff.setSigma ( m_sigma ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::CutOffStudent::evaluate() const 
{
  setPars() ;
  return m_cutoff ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::CutOffStudent::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::CutOffStudent::analyticalIntegral
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
  return m_cutoff.integral ( xmin , xmax ) ;
}
// ============================================================================


// ============================================================================
// Flat in 1D
// ============================================================================
Ostap::Models::Uniform::Uniform
( const char*          name      , 
  const char*          title     ,
  RooAbsReal&          x         )
  : RooAbsPdf  ( name , title ) 
  , m_dim      ( 1  ) 
  , m_x        ( "!x"  , "Observable" , this , x ) 
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
  , m_dim      ( 2  ) 
  , m_x        ( "!x"  , "Observable" , this , x ) 
  , m_y        ( "!y"  , "Observable" , this , y ) 
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
  , m_dim      ( 3  ) 
  , m_x        ( "!x"  , "x-observable" , this , x ) 
  , m_y        ( "!y"  , "y-observable" , this , y ) 
  , m_z        ( "!z"  , "bservable" , this , z ) 
{}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Uniform::Uniform
( const Ostap::Models::Uniform& right ,
  const char*                   name  ) 
  : RooAbsPdf  ( right , name )
  , m_dim ( right.m_dim  ) 
  , m_x ( "!x" , this , right.m_x ) 
  , m_y ( "!y" , this , right.m_y ) 
  , m_z ( "!z" , this , right.m_z ) 
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
  else if ( 3 == m_dim && matchArgs ( allVars , analVars , m_x       , m_z ) ) { return 2 ; }
  else if ( 3 == m_dim && matchArgs ( allVars , analVars       , m_y , m_z ) ) { return 3 ; }
  else if ( 2 <= m_dim && matchArgs ( allVars , analVars , m_x , m_y       ) ) { return 4 ; }
  else if ( 3 == m_dim && matchArgs ( allVars , analVars             , m_z ) ) { return 5 ; }
  else if ( 2 <= m_dim && matchArgs ( allVars , analVars       , m_y       ) ) { return 6 ; }
  else if ( 1 <= m_dim && matchArgs ( allVars , analVars , m_x             ) ) { return 7 ; }
  //
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Uniform::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  //
  // 3D-integral 
  if      ( 3 == m_dim && 1 == code ) 
  {
    return 
      ( m_x.max ( rangeName ) - m_x.min ( rangeName ) ) *
      ( m_y.max ( rangeName ) - m_y.min ( rangeName ) ) *
      ( m_z.max ( rangeName ) - m_z.min ( rangeName ) ) ;
  }
  // 2D-integral: x,z
  else if ( 3 == m_dim && 2 == code ) 
  {
    return 
      ( m_x.max ( rangeName ) - m_x.min ( rangeName ) ) *
      ( m_z.max ( rangeName ) - m_z.min ( rangeName ) ) ;    
  }
  // 2D-integral: y,z
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
  return 0 ;
}
// ============================================================================




// ============================================================================
// Rice distrbution
// ============================================================================
Ostap::Models::Rice::Rice
( const char*  name      , 
  const char*  title     , 
  RooAbsReal&  x         ,
  RooAbsReal&  nu        ,
  RooAbsReal&  varsigma  ,
  RooAbsReal&  shift     ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"        , "x-observable"       , this , x        )
  , m_nu       ( "nu"       , "nu-parameter"       , this , nu       )
  , m_varsigma ( "varsigma" , "varsigma-parameter" , this , varsigma )
  , m_shift    ( "shift"    , "shift-parameter"    , this , shift    )
  , m_rice () 
{
  setPars() ;
}
// ============================================================================
Ostap::Models::Rice::Rice
( const char*  name      , 
  const char*  title     , 
  RooAbsReal&  x         ,
  RooAbsReal&  nu        ,
  RooAbsReal&  varsigma  ) 
  : Rice ( name , title , x , nu , varsigma , RooFit::RooConst ( 0.0 ) ) 
{}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Rice::Rice
( const Ostap::Models::Rice& right ,
  const char*                   name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"        , this , right.m_x        ) 
  , m_nu       ( "nu"       , this , right.m_nu       ) 
  , m_varsigma ( "varsigma" , this , right.m_varsigma ) 
  , m_shift    ( "shift"    , this , right.m_shift    ) 
  , m_rice     ( right.m_rice ) 
{}  
// ============================================================================
// clone method
// ============================================================================
Ostap::Models::Rice* 
Ostap::Models::Rice::clone ( const char* name ) const 
{ return new Ostap::Models::Rice( *this , name ) ; }
// ============================================================================
Ostap::Models::Rice::~Rice(){}
// ============================================================================
void Ostap::Models::Rice::setPars () const 
{
  m_rice.setNu       ( m_nu       ) ;
  m_rice.setVarsigma ( m_varsigma ) ;
  m_rice.setShift    ( m_shift    ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Rice::evaluate() const 
{
  setPars() ;
  return m_rice( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::Rice::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Rice::analyticalIntegral
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
  return m_rice.integral ( xmin , xmax ) ;
}
// ============================================================================





// ============================================================================
// GIG distrbution
// ============================================================================
Ostap::Models::GenInvGauss::GenInvGauss
( const char*  name      , 
  const char*  title     , 
  RooAbsReal&  x         ,
  RooAbsReal&  theta     ,
  RooAbsReal&  eta       ,
  RooAbsReal&  p         ,
  RooAbsReal&  shift     ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"        , "x-observable"       , this , x        )
  , m_theta    ( "theta"    , "theta-parameter"    , this , theta    )
  , m_eta      ( "eta  "    , "eta-parameter"      , this , eta      )
  , m_p        ( "p"        , "p-parameter"        , this , p        )
  , m_shift    ( "shift"    , "shift-parameter"    , this , shift    )
  , m_gig () 
{
  setPars() ;
}
// ============================================================================
Ostap::Models::GenInvGauss::GenInvGauss
( const char*  name      , 
  const char*  title     , 
  RooAbsReal&  x         ,
  RooAbsReal&  theta     ,
  RooAbsReal&  eta       ,
  RooAbsReal&  p         )
  : GenInvGauss ( name , title , x , theta , eta , p , RooFit::RooConst ( 0.0 ) ) 
{}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::GenInvGauss::GenInvGauss
( const Ostap::Models::GenInvGauss& right ,
  const char*                       name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"        , this , right.m_x        ) 
  , m_theta    ( "theta"    , this , right.m_theta    ) 
  , m_eta      ( "eta"      , this , right.m_eta      ) 
  , m_p        ( "p"        , this , right.m_p        ) 
  , m_shift    ( "shift"    , this , right.m_shift    ) 
  , m_gig      ( right.m_gig ) 
{}  
// ============================================================================
// clone method
// ============================================================================
Ostap::Models::GenInvGauss* 
Ostap::Models::GenInvGauss::clone ( const char* name ) const 
{ return new Ostap::Models::GenInvGauss( *this , name ) ; }
// ============================================================================
Ostap::Models::GenInvGauss::~GenInvGauss(){}
// ============================================================================
void Ostap::Models::GenInvGauss::setPars () const 
{
  m_gig.setTheta ( m_theta ) ;
  m_gig.setEta   ( m_eta   ) ;
  m_gig.setP     ( m_p     ) ;
  m_gig.setShift ( m_shift ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::GenInvGauss::evaluate() const 
{
  setPars() ;
  return m_gig ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::GenInvGauss::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::GenInvGauss::analyticalIntegral
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
  return m_gig.integral ( xmin , xmax ) ;
}
// ============================================================================



// ============================================================================
// SkewGenT 
// ============================================================================
Ostap::Models::SkewGenT::SkewGenT 
( const char*  name  , 
  const char*  title , 
  RooAbsReal&  x     ,
  RooAbsReal&  mu    ,   // location/mean  
  RooAbsReal&  sigma ,   // scale/rms 
  RooAbsReal&  xi    ,   // related to asymmetry 
  RooAbsReal&  r     ,   // shape parameter 
  RooAbsReal&  zeta  )   // shape parameter 
  : RooAbsPdf ( name , title )
    //
  , m_x       ( "!x"     , "x-observable"  , this , x     )
  , m_mu      ( "!mu"    , "location/mean" , this , mu    )
  , m_sigma   ( "!sigma" , "sigma/rms"     , this , sigma )
  , m_xi      ( "!xi"    , "asymmetry"     , this , xi    )
  , m_r       ( "!r"     , "r-shape"       , this , r     )
  , m_zeta    ( "!zeta"  , "zeta-shape"    , this , zeta  )
    //
  , m_sgt     ( 0.1222 , 0.1222 , 0.1222 , 0.1222 , 0.1222 ) 
{
  setPars() ;
}
// ============================================================================
// SkewGenT 
// ============================================================================
Ostap::Models::SkewGenT::SkewGenT 
( const Ostap::Models::SkewGenT& right , 
  const char*                     name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "!x"      , this , right.m_x     ) 
  , m_mu       ( "!mu"     , this , right.m_mu    ) 
  , m_sigma    ( "!mu"     , this , right.m_sigma ) 
  , m_xi       ( "!mu"     , this , right.m_xi    ) 
  , m_r        ( "!mu"     , this , right.m_r     ) 
  , m_zeta     ( "!mu"     , this , right.m_zeta  ) 
    //
  , m_sgt      ( right.m_sgt )
{
  setPars() ;
}
// ============================================================================
// clone method
// ============================================================================
Ostap::Models::SkewGenT* 
Ostap::Models::SkewGenT::clone ( const char* name ) const 
{ return new Ostap::Models::SkewGenT( *this , name ) ; }
// ============================================================================
Ostap::Models::SkewGenT::~SkewGenT(){}
// ============================================================================
void Ostap::Models::SkewGenT::setPars () const 
{
  m_sgt.setMu     ( m_mu    ) ;
  m_sgt.setSigma  ( m_sigma ) ;
  m_sgt.setXi     ( m_xi    ) ;
  m_sgt.setR      ( m_r     ) ;
  m_sgt.setZeta   ( m_zeta  ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::SkewGenT::evaluate() const 
{
  setPars() ;
  return m_sgt ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::SkewGenT::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::SkewGenT::analyticalIntegral
( Int_t       code      , 
  const char* rangeName ) const
{
  assert ( code == 1 ) ;
  if ( 1 != code ){}
  //
  setPars() ;
  //
  return m_sgt.integral ( m_x.min ( rangeName ) , m_x.max ( rangeName ) ) ;
}
// ============================================================================


// ============================================================================
// HORNSdini 
// ============================================================================
Ostap::Models::HORNSdini::HORNSdini
( const char*  name      , 
  const char*  title     , 
  RooAbsReal&  x         ,
  RooAbsReal&  a         ,
  RooAbsReal&  delta     ,
  RooAbsReal&  phi       )
  : RooAbsPdf ( name , title ) 
  , m_x       ( "!x"     , "x-observable" , this , x     )
  , m_a       ( "!a"     , "left horn"    , this , a     )
  , m_delta   ( "!delta" , "width"        , this , delta )
  , m_phi     ( "!phi"   , "modification" , this , phi   )
  , m_horns () 
{
  setPars() ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::HORNSdini::HORNSdini
( const Ostap::Models::HORNSdini& right ,
  const char*                     name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "!x"      , this , right.m_x     ) 
  , m_a        ( "!a"      , this , right.m_a     ) 
  , m_delta    ( "!delta"  , this , right.m_delta ) 
  , m_phi      ( "!phi"    , this , right.m_phi   ) 
  , m_horns      ( right.m_horns ) 
{}  
// ============================================================================
// clone method
// ============================================================================
Ostap::Models::HORNSdini* 
Ostap::Models::HORNSdini::clone ( const char* name ) const 
{ return new Ostap::Models::HORNSdini( *this , name ) ; }
// ============================================================================
Ostap::Models::HORNSdini::~HORNSdini(){}
// ============================================================================
void Ostap::Models::HORNSdini::setPars () const 
{
  m_horns.setA     ( m_a     ) ;
  m_horns.setDelta ( m_delta ) ;
  m_horns.setPhi   ( m_phi   ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::HORNSdini::evaluate() const 
{
  setPars() ;
  return m_horns( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::HORNSdini::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::HORNSdini::analyticalIntegral
( Int_t       code      , 
  const char* rangeName ) const
{
  assert ( code == 1 ) ;
  if ( 1 != code ){}
  //
  setPars() ;
  //
  return m_horns.integral ( m_x.min ( rangeName ) , m_x.max ( rangeName ) ) ;
}
// ============================================================================



// ============================================================================
// HILLdini 
// ============================================================================
Ostap::Models::HILLdini::HILLdini
( const char*  name      , 
  const char*  title     , 
  RooAbsReal&  x         ,
  RooAbsReal&  a         ,
  RooAbsReal&  delta     ,
  RooAbsReal&  phi       )
  : RooAbsPdf ( name , title ) 
  , m_x       ( "!x"     , "x-observable" , this , x     )
  , m_a       ( "!a"     , "left horn"    , this , a     )
  , m_delta   ( "!delta" , "width"        , this , delta )
  , m_phi     ( "!phi"   , "modification" , this , phi   )
  , m_hill () 
{
  setPars() ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::HILLdini::HILLdini
( const Ostap::Models::HILLdini& right ,
  const char*                     name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "!x"      , this , right.m_x     ) 
  , m_a        ( "!a"      , this , right.m_a     ) 
  , m_delta    ( "!delta"  , this , right.m_delta ) 
  , m_phi      ( "!phi"    , this , right.m_phi   ) 
  , m_hill     ( right.m_hill ) 
{}  
// ============================================================================
// clone method
// ============================================================================
Ostap::Models::HILLdini* 
Ostap::Models::HILLdini::clone ( const char* name ) const 
{ return new Ostap::Models::HILLdini( *this , name ) ; }
// ============================================================================
Ostap::Models::HILLdini::~HILLdini(){}
// ============================================================================
void Ostap::Models::HILLdini::setPars () const 
{
  m_hill.setA     ( m_a     ) ;
  m_hill.setDelta ( m_delta ) ;
  m_hill.setPhi   ( m_phi   ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::HILLdini::evaluate() const 
{
  setPars() ;
  return m_hill ( m_x ) ;
}
// ============================================================================
Int_t Ostap::Models::HILLdini::getAnalyticalIntegral
( RooArgSet&  allVars       , 
  RooArgSet&  analVars      ,
  const char* /*rangeName*/ ) const
{
  if ( matchArgs ( allVars , analVars , m_x ) ) { return 1 ; }
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::HILLdini::analyticalIntegral
( Int_t       code      , 
  const char* rangeName ) const
{
  assert ( code == 1 ) ;
  if ( 1 != code ){}
  //
  setPars() ;
  //
  return m_hill.integral ( m_x.min ( rangeName ) , m_x.max ( rangeName ) ) ;
}
// ============================================================================


// ============================================================================
ClassImp(Ostap::Models::Shape1D            ) 
ClassImp(Ostap::Models::Shape2D            ) 
ClassImp(Ostap::Models::Shape3D            ) 
ClassImp(Ostap::Models::Histo1D            ) 
ClassImp(Ostap::Models::Histo2D            ) 
ClassImp(Ostap::Models::Histo3D            ) 
ClassImp(Ostap::Models::Uniform            ) 
ClassImp(Ostap::Models::BreitWigner        ) 
ClassImp(Ostap::Models::BreitWignerMC      ) 
ClassImp(Ostap::Models::BWI                ) 
ClassImp(Ostap::Models::BWPS               ) 
ClassImp(Ostap::Models::BW3L               ) 
ClassImp(Ostap::Models::Flatte             ) 
ClassImp(Ostap::Models::FlatteBugg         ) 
ClassImp(Ostap::Models::LASS               ) 
ClassImp(Ostap::Models::Voigt              ) 
ClassImp(Ostap::Models::PseudoVoigt        ) 
// ClassImp(Ostap::Models::Swanson         ) 
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
ClassImp(Ostap::Models::Novosibirsk        ) 
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
ClassImp(Ostap::Models::GenArgus           ) 
ClassImp(Ostap::Models::Slash              ) 
ClassImp(Ostap::Models::AsymmetricLaplace  ) 
ClassImp(Ostap::Models::Hat                ) 
ClassImp(Ostap::Models::Up                 ) 
ClassImp(Ostap::Models::FupN               ) 
ClassImp(Ostap::Models::Tsallis            ) 
ClassImp(Ostap::Models::QGSM               ) 
ClassImp(Ostap::Models::TwoExpos           ) 
ClassImp(Ostap::Models::DoubleGauss        ) 
ClassImp(Ostap::Models::Gumbel             )
ClassImp(Ostap::Models::Weibull            )
ClassImp(Ostap::Models::RaisingCosine      )
ClassImp(Ostap::Models::QGaussian          )
ClassImp(Ostap::Models::KGaussian          )
ClassImp(Ostap::Models::Hyperbolic         )
ClassImp(Ostap::Models::GenHyperbolic      )
ClassImp(Ostap::Models::Das                )
ClassImp(Ostap::Models::PositiveSpline     ) 
ClassImp(Ostap::Models::MonotonicSpline    ) 
ClassImp(Ostap::Models::ConvexOnlySpline   )
ClassImp(Ostap::Models::ConvexSpline       )
ClassImp(Ostap::Models::CutOffGauss        )
ClassImp(Ostap::Models::CutOffStudent      )
ClassImp(Ostap::Models::Rice               )
ClassImp(Ostap::Models::GenInvGauss        )
ClassImp(Ostap::Models::ExGauss            )
ClassImp(Ostap::Models::NormalLaplace      )
ClassImp(Ostap::Models::PearsonIV          )
ClassImp(Ostap::Models::SkewGenT           )
ClassImp(Ostap::Models::HORNSdini          )
ClassImp(Ostap::Models::HILLdini           )
// ============================================================================
//                                                                      The END 
// ============================================================================

