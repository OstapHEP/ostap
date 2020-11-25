// ============================================================================
// Include files 
// ============================================================================
// STD & STL:
// ============================================================================
#include <limits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/PDFs2D.h"
#include "Ostap/Iterator.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
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
  template<class TX, class TY, class FUN>
  void compute_XY ( RooSpan<double> output , FUN& fun , TX x , TY y ) 
  {
    const int n = output.size();
    for ( int i = 0 ; i < n ; ++i ) { output [ i ] = fun ( x [ i ] , y [ i ] ) ; }
  }
  // ==========================================================================
}
#endif
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
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::Poly2DPositive"          ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_positive.npars () , 
                  "Widths/#channels mismatch" , 
                  "Ostap::Models::Poly2DPositive"      ) ;
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
{ ::set_pars ( m_phis , m_positive ) ; }
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
  , m_positive ( n , 
                 std::min ( x.getMin() , y.getMin() ) ,
                 std::max ( x.getMax() , y.getMax() ) )
{
  //
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::Poly2DSymPositive"       ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_positive.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Poly2DSymPositive"       ) ;
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
{ ::set_pars ( m_phis ,  m_positive ) ; }
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
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::PS2DPol          "       ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Ps2DPol"                 ) ;
  //
  setPars () ;
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
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::PS2DPol"                 ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Ps2DPol"                 ) ;
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
{ ::set_pars ( m_phis ,  m_function ) ; }
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


// ============================================================================
//  PS(x)*PS(y)*Polynom 
// ============================================================================
Ostap::Models::PS2DPol2::PS2DPol2
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const Ostap::Math::PhaseSpaceNL& psx , 
  const Ostap::Math::PhaseSpaceNL& psy ,
  const double         mmax      , 
  const unsigned short nX        ,
  const unsigned short nY        ,
  RooArgList&          phis      )
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_function ( psx, psy , mmax , 
                 nX , nY  , 
                 x.getMin() , x.getMax() ,
                 y.getMin() , y.getMax() )
{
  //  
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::PS2DPol2"                ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Ps2DPol2"                ) ;
  //
  setPars () ;
  //
}
// ============================================================================
// generic polinomial
// ============================================================================
Ostap::Models::PS2DPol2::PS2DPol2
( const char*          name       , 
  const char*          title      ,
  RooRealVar&          x          ,
  RooRealVar&          y          ,
  const Ostap::Math::PS2DPol2& ps , 
  RooArgList&          phis       ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_function ( ps ) 
{
  //
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!"   ,
                  "Ostap::Models::PS2DPol2"                  ) ;
  //
  Ostap::Assert ( ::size ( m_phis )   == m_function.npars () , 
                  "Widths/#channels mismatch"                , 
                  "Ostap::Models::Ps2DPol2"                  ) ;
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PS2DPol2::PS2DPol2
( const Ostap::Models::PS2DPol2&  right ,      
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
Ostap::Models::PS2DPol2::~PS2DPol2() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PS2DPol2*
Ostap::Models::PS2DPol2::clone( const char* name ) const 
{ return new Ostap::Models::PS2DPol2 ( *this , name ) ; }
// ============================================================================
void Ostap::Models::PS2DPol2::setPars () const 
{ ::set_pars ( m_phis ,  m_function ) ; }
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PS2DPol2::evaluate() const 
{
  //
  setPars () ;
  //
  return m_function ( m_x , m_y ) ;
}
// ============================================================================
Int_t Ostap::Models::PS2DPol2::getAnalyticalIntegral
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
Double_t Ostap::Models::PS2DPol2::analyticalIntegral 
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
//  PS(x)*PS(y)*Polynom 
// ============================================================================
Ostap::Models::PS2DPol3::PS2DPol3
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const Ostap::Math::PhaseSpaceNL& psx , 
  const Ostap::Math::PhaseSpaceNL& psy ,
  const double         mmax      , 
  const unsigned short nX        ,
  const unsigned short nY        ,
  RooArgList&          phis      )
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_function ( psx, psy , mmax , 
                 nX , nY  ,  
                 x.getMin() , x.getMax() ,
                 y.getMin() , y.getMax() )
{
  //  
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::PS2DPol3"                ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Ps2DPol3"                ) ;
  //
  setPars () ;
  //
}
// ============================================================================
// generic polinomial
// ============================================================================
Ostap::Models::PS2DPol3::PS2DPol3
( const char*          name       , 
  const char*          title      ,
  RooRealVar&          x          ,
  RooRealVar&          y          ,
  const Ostap::Math::PS2DPol3& ps , 
  RooArgList&          phis       ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_function ( ps ) 
{
  //  
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::PS2DPol3"                ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Ps2DPol3"                ) ;
  //
  setPars () ;
  //
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PS2DPol3::PS2DPol3
( const Ostap::Models::PS2DPol3&    right ,      
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
Ostap::Models::PS2DPol3::~PS2DPol3() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PS2DPol3*
Ostap::Models::PS2DPol3::clone( const char* name ) const 
{ return new Ostap::Models::PS2DPol3 ( *this , name ) ; }
// ============================================================================
void Ostap::Models::PS2DPol3::setPars () const 
{ ::set_pars ( m_phis ,  m_function ) ; }
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PS2DPol3::evaluate() const 
{
  //
  setPars () ;
  //
  return m_function ( m_x , m_y ) ;
}
// ============================================================================
Int_t Ostap::Models::PS2DPol3::getAnalyticalIntegral
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
Double_t Ostap::Models::PS2DPol3::analyticalIntegral 
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
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::PS2DPolSym"              ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Ps2DPolSym"              ) ;
  //
  setPars () ;
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
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::PS2DPolSym"              ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Ps2DPolSym"              ) ;
  //
  setPars () ;
  //
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
{ ::set_pars ( m_phis ,  m_function ) ; }
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


// ============================================================================
//  PS(x)*PS(y)*SymPolynom 
// ============================================================================
Ostap::Models::PS2DPol2Sym::PS2DPol2Sym
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const Ostap::Math::PhaseSpaceNL& ps ,
  const double         mmax      ,
  const unsigned short n         ,
  RooArgList&          phis      )
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_function ( ps , mmax , n , x.getMin() , x.getMax() ) 
{
  //
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::PS2DPol2Sym"             ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Ps2DPol2Sym"             ) ;
  //
  setPars () ;
  //
}
// ============================================================================
// generic polinomial
// ============================================================================
Ostap::Models::PS2DPol2Sym::PS2DPol2Sym
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const Ostap::Math::PS2DPol2Sym& ps , 
  RooArgList&          phis      ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_function ( ps ) 
{
  //
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::PS2DPol2Sym"             ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Ps2DPol2Sym"             ) ;
  //
  setPars () ;
  //
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PS2DPol2Sym::PS2DPol2Sym
( const Ostap::Models::PS2DPol2Sym& right ,      
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
Ostap::Models::PS2DPol2Sym::~PS2DPol2Sym() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PS2DPol2Sym*
Ostap::Models::PS2DPol2Sym::clone( const char* name ) const 
{ return new Ostap::Models::PS2DPol2Sym ( *this , name ) ; }
// ============================================================================
void Ostap::Models::PS2DPol2Sym::setPars () const 
{ ::set_pars ( m_phis ,  m_function ) ; }
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PS2DPol2Sym::evaluate() const 
{
  //
  setPars () ;
  //
  return m_function ( m_x , m_y ) ;
}
// ============================================================================
Int_t Ostap::Models::PS2DPol2Sym::getAnalyticalIntegral
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
Double_t Ostap::Models::PS2DPol2Sym::analyticalIntegral 
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
//  PS(x)*PS(y)*SymPolynom 
// ============================================================================
Ostap::Models::PS2DPol3Sym::PS2DPol3Sym
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const Ostap::Math::PhaseSpaceNL& ps ,
  const double         mmax      ,
  const unsigned short n         ,
  RooArgList&          phis      )
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_function ( ps , mmax , n , x.getMin() , x.getMax() ) 
{
  //
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::PS2DPol3Sym"             ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Ps2DPol3Sym"             ) ;
  //
  setPars () ;
  //
}
// ============================================================================
// generic polinomial
// ============================================================================
Ostap::Models::PS2DPol3Sym::PS2DPol3Sym
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  const Ostap::Math::PS2DPol3Sym& ps , 
  RooArgList&          phis      ) 
  : RooAbsPdf ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_function ( ps ) 
{
  //
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::PS2DPol3Sym"             ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Ps2DPol3Sym"             ) ;
  //
  setPars () ;
  //
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PS2DPol3Sym::PS2DPol3Sym
( const Ostap::Models::PS2DPol3Sym& right ,      
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
Ostap::Models::PS2DPol3Sym::~PS2DPol3Sym() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::PS2DPol3Sym*
Ostap::Models::PS2DPol3Sym::clone( const char* name ) const 
{ return new Ostap::Models::PS2DPol3Sym ( *this , name ) ; }
// ============================================================================
void Ostap::Models::PS2DPol3Sym::setPars () const 
{ ::set_pars ( m_phis ,  m_function ) ; }
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::PS2DPol3Sym::evaluate() const 
{
  //
  setPars () ;
  //
  return m_function ( m_x , m_y ) ;
}
// ============================================================================
Int_t Ostap::Models::PS2DPol3Sym::getAnalyticalIntegral
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
Double_t Ostap::Models::PS2DPol3Sym::analyticalIntegral 
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
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::ExpoPS2DPol"             ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::ExpoPS2DPol"             ) ;
  //
  setPars () ;
  //
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
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::ExpoPS2DPol"             ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::ExpoPS2DPol"             ) ;
  //
  setPars () ;
  //
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
  ::set_pars ( m_phis ,  m_function ) ;
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
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::Expo2DPol"               ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Expo2DPol"               ) ;
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
  ::set_pars ( m_phis ,  m_function ) ;
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
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::Expo2DPolSym"            ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_function.npars () , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Expo2DPolSym"            ) ;
  //
  setPars () ;
  //
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
  //
  ::set_pars ( m_phis ,  m_function ) ;
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
( const char*                          name      , 
  const char*                          title     ,
  RooRealVar&                          x         ,
  RooRealVar&                          y         ,
  const Ostap::Math::PositiveSpline2D& spline, 
  RooArgList&          phis      ) 
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_spline   ( spline ) 
{
  //
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::Spline2D"                ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_spline.npars ()   , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Spline2D"                ) ;
  //
  setPars () ;
  //
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
{ ::set_pars ( m_phis ,  m_spline ) ; }
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



// ============================================================================
// generic 2D symmetric spline
// ============================================================================
Ostap::Models::Spline2DSym::Spline2DSym
( const char*                             name      , 
  const char*                             title     ,
  RooRealVar&                             x         ,
  RooRealVar&                             y         ,
  const Ostap::Math::PositiveSpline2DSym& spline, 
  RooArgList&          phis      ) 
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_spline   ( spline ) 
{
  //
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::Spline2DSym"             ) ;
  //
  Ostap::Assert ( ::size ( m_phis ) == m_spline.npars ()   , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Spline2DSym"             ) ;
  //
  setPars () ;
  //
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
{ ::set_pars ( m_phis ,  m_spline ) ; }
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





// ============================================================================
ClassImp(Ostap::Models::Poly2DPositive       ) 
ClassImp(Ostap::Models::Poly2DSymPositive    )
ClassImp(Ostap::Models::PS2DPol              )
ClassImp(Ostap::Models::PS2DPol2             )
ClassImp(Ostap::Models::PS2DPol3             )
ClassImp(Ostap::Models::PS2DPolSym           )
ClassImp(Ostap::Models::PS2DPol2Sym          )
ClassImp(Ostap::Models::PS2DPol3Sym          )
ClassImp(Ostap::Models::ExpoPS2DPol          ) 
ClassImp(Ostap::Models::Expo2DPol            ) 
ClassImp(Ostap::Models::Expo2DPolSym         ) 
ClassImp(Ostap::Models::Spline2D             ) 
ClassImp(Ostap::Models::Spline2DSym          )
// ============================================================================
//                                                                      The END         
// ============================================================================
