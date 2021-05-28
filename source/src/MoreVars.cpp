// ============================================================================
// Include files 
// ============================================================================
// ROOT/RooFit
// ============================================================================
#include "RVersion.h"
#include "RooAbsPdf.h"
#include "RooConstVar.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/MoreVars.h"
#include "Ostap/Iterator.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
#include "local_roofit.h"
// ============================================================================
/** @file
 *  Implementation file for objects from the file Ostap/MoreVars.h
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
 *  @date 2020-03-06     
 */
// ============================================================================
namespace 
{
  // ===========================================================================
  const std::string s_EMPTYPARS   = "Vector of coefficients is empty!" ;
  const std::string s_INVALIDPARS = "Invalid parameters!"              ;
  const std::string s_INVALIDPAR  = "Invalid parameter!"               ;
  const std::string s_v1          = "Ostap::MoreRooFit::Bernstein"     ;
  const std::string s_v2          = "Ostap::MoreRooFit::Monotonic"     ;
  const std::string s_v3          = "Ostap::MoreRooFit::Convex"        ;
  const std::string s_v4          = "Ostap::MoreRooFit::ConvexOnly"    ;
  // ===========================================================================
}
// ============================================================================
ClassImp ( Ostap::MoreRooFit::Bernstein  ) ;
ClassImp ( Ostap::MoreRooFit::Monotonic  ) ;
ClassImp ( Ostap::MoreRooFit::Convex     ) ;
ClassImp ( Ostap::MoreRooFit::ConvexOnly ) ;
// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::Bernstein::Bernstein 
( const std::string& name  ,
  const std::string& title ,
  RooAbsReal&        xvar  ,
  const double       xmin  , 
  const double       xmax  ,
  const RooArgList&  pars  ) 
  : RooAbsReal ( name.c_str  () , title.c_str () )
  , m_xvar      ( "x"    , "Dependent"  , this , xvar ) 
  , m_pars      ( "pars" , "Parameters" , this ) 
  , m_bernstein ( ::size ( pars ) - 1 , xmin , xmax ) 
{
  //
  Ostap::Assert ( 1 <= ::size ( pars ) , s_EMPTYPARS  , s_v1 , 510 ) ;
  //
  ::copy_real   ( pars , m_pars , s_INVALIDPAR , s_v1 ) ;
  //
  Ostap::Assert ( m_bernstein.npars() == ::size ( m_pars ), s_INVALIDPARS , s_v1 , 512 ) ;
  //
}
// =============================================================================
// copy constructor 
// =============================================================================
Ostap::MoreRooFit::Bernstein::Bernstein
( const Ostap::MoreRooFit::Bernstein& right , 
  const char*                         name  )
  : RooAbsReal  ( right , name ) 
  , m_xvar      ( "x"    , this , right.m_xvar ) 
  , m_pars      ( "pars" , this , right.m_pars ) 
  , m_bernstein ( right.m_bernstein ) 
{}
// =============================================================================
// default constructor 
// ============================================================================
Ostap::MoreRooFit::Bernstein::Bernstein(){}
// =============================================================================
// destructor  
// ============================================================================
Ostap::MoreRooFit::Bernstein::~Bernstein(){}
// ============================================================================
// clone it!
// ============================================================================
Ostap::MoreRooFit::Bernstein*
Ostap::MoreRooFit::Bernstein::clone ( const char* newname ) const
{ return new Bernstein( *this , newname ) ; }
// ============================================================================
void Ostap::MoreRooFit::Bernstein::setPars() const 
{ ::set_pars ( m_pars , m_bernstein ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Bernstein::evaluate  () const 
{
  setPars () ;
  const double x = m_xvar ;
  return m_bernstein ( x ) ;
}
// ============================================================================
Int_t    Ostap::MoreRooFit::Bernstein::getAnalyticalIntegral
( RooArgSet&  allVars   , 
  RooArgSet&  analVars  , 
  const char* /* rangeName */ ) const 
{
  return matchArgs  ( allVars , analVars , m_xvar ) ? 1 : 0 ;
}
// ===========================================================================
Double_t Ostap::MoreRooFit::Bernstein::analyticalIntegral 
( Int_t code            , 
  const char* rangeName ) const 
{
  setPars() ;
  const double xmin = m_xvar.min ( rangeName ) ;
  const double xmax = m_xvar.max ( rangeName ) ;
  return m_bernstein.integral ( xmin , xmax ) ;
}
// ============================================================================


// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::Monotonic::Monotonic
( const std::string& name       ,
  const std::string& title      ,
  RooAbsReal&        xvar       ,
  const bool         increasing ,
  const double       xmin       , 
  const double       xmax       ,
  RooAbsReal&        a          ,
  RooAbsReal&        b          ,
  const RooArgList&  pars       ) 
  : RooAbsReal ( name.c_str  () , title.c_str () )
  , m_xvar      ( "x"    , "Dependent"  , this , xvar ) 
  , m_a         ( "a"    , "shift/bias" , this , a    ) 
  , m_b         ( "b"    , "scale"      , this , b    ) 
  , m_pars      ( "pars" , "Parameters" , this ) 
  , m_monotonic ( ::size ( pars ) , xmin , xmax , increasing ) 
{
  //
  Ostap::Assert ( 0 <= ::size ( pars ) , s_EMPTYPARS  , s_v2 , 510 ) ;
  //
  ::copy_real   ( pars , m_pars , s_INVALIDPAR , s_v2 ) ;
  //
  Ostap::Assert ( m_monotonic.npars() == ::size ( m_pars ), s_INVALIDPARS , s_v2 , 512 ) ;
  //
}
// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::Monotonic::Monotonic
( const std::string& name       ,
  const std::string& title      ,
  RooAbsReal&        xvar       ,
  const bool         increasing , 
  const double       xmin       , 
  const double       xmax       ,
  const RooArgList&  pars       ) 
  : Monotonic ( name , title , 
                xvar , increasing , xmin , xmax ,
                RooFit::RooConst ( 0.0 ) , 
                RooFit::RooConst ( 1.0 ) ,
                pars ) 
{}
// =============================================================================
// copy constructor 
// =============================================================================
Ostap::MoreRooFit::Monotonic::Monotonic
( const Ostap::MoreRooFit::Monotonic& right , 
  const char*                         name  )
  : RooAbsReal  ( right , name ) 
  , m_xvar      ( "x"    , this , right.m_xvar ) 
  , m_a         ( "a"    , this , right.m_a    ) 
  , m_b         ( "b"    , this , right.m_b    ) 
  , m_pars      ( "pars" , this , right.m_pars ) 
  , m_monotonic ( right.m_monotonic ) 
{}
// =============================================================================
// default constructor 
// ============================================================================
Ostap::MoreRooFit::Monotonic::Monotonic (){}
// =============================================================================
// destructor  
// ============================================================================
Ostap::MoreRooFit::Monotonic::~Monotonic(){}
// ============================================================================
// clone it!
// ============================================================================
Ostap::MoreRooFit::Monotonic*
Ostap::MoreRooFit::Monotonic::clone ( const char* newname ) const
{ return new Monotonic ( *this , newname ) ; }
// ============================================================================
void Ostap::MoreRooFit::Monotonic::setPars() const 
{ ::set_pars (  m_pars , m_monotonic ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Monotonic::evaluate  () const 
{
  setPars () ;
  //
  const double x = m_xvar ;
  if ( x < m_monotonic.xmin () || x > m_monotonic.xmax() ) { return 0 ; }
  const double a = m_a    ;
  const double b = m_b    ;
  //
  return a + b * m_monotonic ( x ) ;
}
// ============================================================================
Int_t    Ostap::MoreRooFit::Monotonic::getAnalyticalIntegral
( RooArgSet&  allVars   , 
  RooArgSet&  analVars  , 
  const char* /* rangeName */ ) const 
{
  return matchArgs  ( allVars , analVars , m_xvar ) ? 1 : 0 ;
}
// ===========================================================================
Double_t Ostap::MoreRooFit::Monotonic::analyticalIntegral 
( Int_t code            , 
  const char* rangeName ) const 
{
  setPars() ;
  //
  const double xmin = std::max ( m_xvar.min ( rangeName ) , m_monotonic.xmin () ) ;
  const double xmax = std::min ( m_xvar.max ( rangeName ) , m_monotonic.xmax () ) ;
  //
  const double a = m_a    ;
  const double b = m_b    ;
  //
  return a * ( xmax - xmin ) + b * m_monotonic.integral ( xmin , xmax ) ;
}
// ============================================================================



// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::Convex::Convex
( const std::string& name       ,
  const std::string& title      ,
  RooAbsReal&        xvar       ,
  const bool         increasing ,
  const bool         convex     ,
  const double       xmin       , 
  const double       xmax       ,
  RooAbsReal&        a          ,
  RooAbsReal&        b          ,
  const RooArgList&  pars       ) 
  : RooAbsReal ( name.c_str  () , title.c_str () )
  , m_xvar      ( "x"    , "Dependent"  , this , xvar ) 
  , m_a         ( "a"    , "shift/bias" , this , a    ) 
  , m_b         ( "b"    , "scale"      , this , b    ) 
  , m_pars      ( "pars" , "Parameters" , this ) 
  , m_convex    ( ::size ( pars ) , xmin , xmax , increasing , convex ) 
{
  //
  Ostap::Assert ( 0 <= ::size ( pars ) , s_EMPTYPARS  , s_v3 , 510 ) ;
  //
  ::copy_real   ( pars , m_pars , s_INVALIDPAR , s_v3 ) ;
  //
  Ostap::Assert ( m_convex.npars() == ::size ( m_pars ), s_INVALIDPARS , s_v3 , 512 ) ;
  //
}
// =============================================================================
// copy constructor 
// =============================================================================
Ostap::MoreRooFit::Convex::Convex
( const Ostap::MoreRooFit::Convex& right , 
  const char*                         name  )
  : RooAbsReal  ( right , name ) 
  , m_xvar      ( "x"    , this , right.m_xvar ) 
  , m_a         ( "a"    , this , right.m_a    ) 
  , m_b         ( "b"    , this , right.m_b    ) 
  , m_pars      ( "pars" , this , right.m_pars ) 
  , m_convex    ( right.m_convex ) 
{}
// =============================================================================
// default constructor 
// ============================================================================
Ostap::MoreRooFit::Convex::Convex (){}
// =============================================================================
// destructor  
// ============================================================================
Ostap::MoreRooFit::Convex::~Convex (){}
// ============================================================================
// clone it!
// ============================================================================
Ostap::MoreRooFit::Convex*
Ostap::MoreRooFit::Convex::clone ( const char* newname ) const
{ return new Convex ( *this , newname ) ; }
// ============================================================================
void Ostap::MoreRooFit::Convex::setPars() const 
{ ::set_pars (  m_pars , m_convex ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Convex::evaluate  () const 
{
  setPars () ;
  //
  const double x = m_xvar ;
  if ( x < m_convex.xmin () || x > m_convex.xmax() ) { return 0 ; }
  const double a = m_a    ;
  const double b = m_b    ;
  //
  return a + b * m_convex ( x ) ;
}
// ============================================================================
Int_t    Ostap::MoreRooFit::Convex::getAnalyticalIntegral
( RooArgSet&  allVars   , 
  RooArgSet&  analVars  , 
  const char* /* rangeName */ ) const 
{
  return matchArgs  ( allVars , analVars , m_xvar ) ? 1 : 0 ;
}
// ===========================================================================
Double_t Ostap::MoreRooFit::Convex::analyticalIntegral 
( Int_t code            , 
  const char* rangeName ) const 
{
  setPars() ;
  //
  const double xmin = std::max ( m_xvar.min ( rangeName ) , m_convex.xmin () ) ;
  const double xmax = std::min ( m_xvar.max ( rangeName ) , m_convex.xmax () ) ;
  //
  const double a = m_a    ;
  const double b = m_b    ;
  //
  return a * ( xmax - xmin ) + b * m_convex.integral ( xmin , xmax ) ;
}
// ============================================================================



// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::ConvexOnly::ConvexOnly
( const std::string& name       ,
  const std::string& title      ,
  RooAbsReal&        xvar       ,
  const bool         convex     ,
  const double       xmin       , 
  const double       xmax       ,
  RooAbsReal&        a          ,
  RooAbsReal&        b          ,
  const RooArgList&  pars       ) 
  : RooAbsReal ( name.c_str  () , title.c_str () )
  , m_xvar      ( "x"    , "Dependent"  , this , xvar ) 
  , m_a         ( "a"    , "shift/bias" , this , a    ) 
  , m_b         ( "b"    , "scale"      , this , b    ) 
  , m_pars      ( "pars" , "Parameters" , this ) 
  , m_convex    ( ::size ( pars ) , xmin , xmax , convex ) 
{
  //
  Ostap::Assert ( 0 <= :: size ( pars ) , s_EMPTYPARS  , s_v4 , 510 ) ;
  //
  ::copy_real   ( pars , m_pars , s_INVALIDPAR , s_v4 ) ;
  //
  Ostap::Assert ( m_convex.npars() == ::size ( m_pars ), s_INVALIDPARS , s_v4 , 512 ) ;
  //
}
// =============================================================================
// copy constructor 
// =============================================================================
Ostap::MoreRooFit::ConvexOnly::ConvexOnly
( const Ostap::MoreRooFit::ConvexOnly& right , 
  const char*                         name  )
  : RooAbsReal  ( right , name ) 
  , m_xvar      ( "x"    , this , right.m_xvar ) 
  , m_a         ( "a"    , this , right.m_a    ) 
  , m_b         ( "b"    , this , right.m_b    ) 
  , m_pars      ( "pars" , this , right.m_pars ) 
  , m_convex    ( right.m_convex ) 
{}
// =============================================================================
// default constructor 
// ============================================================================
Ostap::MoreRooFit::ConvexOnly::ConvexOnly (){}
// =============================================================================
// destructor  
// ============================================================================
Ostap::MoreRooFit::ConvexOnly::~ConvexOnly (){}
// ============================================================================
// clone it!
// ============================================================================
Ostap::MoreRooFit::ConvexOnly*
Ostap::MoreRooFit::ConvexOnly::clone ( const char* newname ) const
{ return new ConvexOnly ( *this , newname ) ; }
// ============================================================================
void Ostap::MoreRooFit::ConvexOnly::setPars() const 
{ ::set_pars (  m_pars , m_convex ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::ConvexOnly::evaluate  () const 
{
  setPars () ;
  //
  const double x = m_xvar ;
  if ( x < m_convex.xmin () || x > m_convex.xmax() ) { return 0 ; }
  const double a = m_a    ;
  const double b = m_b    ;
  //
  return a + b * m_convex ( x ) ;
}
// ============================================================================
Int_t    Ostap::MoreRooFit::ConvexOnly::getAnalyticalIntegral
( RooArgSet&  allVars   , 
  RooArgSet&  analVars  , 
  const char* /* rangeName */ ) const 
{
  return matchArgs  ( allVars , analVars , m_xvar ) ? 1 : 0 ;
}
// ===========================================================================
Double_t Ostap::MoreRooFit::ConvexOnly::analyticalIntegral 
( Int_t code            , 
  const char* rangeName ) const 
{
  setPars() ;
  //
  const double xmin = std::max ( m_xvar.min ( rangeName ) , m_convex.xmin () ) ;
  const double xmax = std::min ( m_xvar.max ( rangeName ) , m_convex.xmax () ) ;
  //
  const double a = m_a    ;
  const double b = m_b    ;
  //
  return a * ( xmax - xmin ) + b * m_convex.integral ( xmin , xmax ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================



