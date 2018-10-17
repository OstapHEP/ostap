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
#include "RooArgSet.h"
#include "RooRealVar.h"
// ============================================================================
/** @file 
 *  Implementation file for namespace Ostap::Models
 *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
 *  @date   2011-11-30
 */
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
                            "Ostap::Poly2DPositive"            , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
  { Ostap::throwException ( "Invalid size of parameters vector", 
                            "Ostap::Poly2DSymPositive"         , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
                            "Ostap::PS22DPol"                  , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
                            "Ostap::PS22DPol"                  , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
                            "Ostap::PS2DPol2"                  , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
                            "Ostap::PS2DPol2"                  , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
                            "Ostap::PS2DPol3"                  , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
                            "Ostap::PS2DPol3"                  , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
  //
  setPars () ;
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
                            "Ostap::PS22DPol"                  , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
  { Ostap::throwException ( "Invalid size of parameters vector", 
                            "Ostap::PS22DPol"                  , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
                            "Ostap::PS2DPol2Sym"               , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
                            "Ostap::PS2DPol2Sym"               , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
  //
  setPars () ;
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
                            "Ostap::PS2DPol3Sym"               , 
                            Ostap::StatusCode::FAILURE         ) ; }
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
                            "Ostap::PS2DPol3Sym"               , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
  //
  setPars () ;
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
  { Ostap::throwException ( "Invalid size of parameters vector", 
                            "Ostap::ExpoPS22DPol"              , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
  { Ostap::throwException ( "Invalid size of parameters vector", 
                            "Ostap::ExpoPS22DPol"              , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
  { Ostap::throwException ( "Invalid size of parameters vector", 
                            "Ostap::ExpoPS22DPol"              , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
  { Ostap::throwException ( "Invalid size of parameters vector", 
                            "Ostap::ExpoPS22DPol"              , 
                            Ostap::StatusCode::FAILURE         ) ; }
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
  { Ostap::throwException ( "Invalid size of parameters vector", 
                            "Ostap::Spline2D"                  , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
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
  { Ostap::throwException ( "Invalid size of parameters vector", 
                            "Ostap::Spline2DSym"               , 
                            StatusCode::FAILURE                ) ; }
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
// Rotated models 
// ============================================================================


// ============================================================================
// Rotated Gaussians 
// ============================================================================
Ostap::Models::Rotated2Gaussians::Rotated2Gaussians
( const char*                             name      , 
  const char*                             title     ,
  RooRealVar&                             x         ,
  RooRealVar&                             y         ,
  RooRealVar&                             phi       , // rotation phase 
  // the  first peak 
  RooAbsReal&                             m1         ,
  RooAbsReal&                             sigma1     ,  //
  // the  second peak 
  RooAbsReal&                             m2         ,
  RooAbsReal&                             sigma2     )  //
  : BASE ( name , title , x , y , phi , 
           GAUSS ( m1      .getVal() , sigma1  .getVal() ) , 
           GAUSS ( m2      .getVal() , sigma2  .getVal() ) ) 
    // first signal 
  , m_m01      ( "m01"      , "mass"                       , this , m1      ) 
  , m_sigma1   ( "sigma1"   , "sigma"                      , this , sigma1  )
    // second signal 
  , m_m02      ( "m02"      , "mass"                       , this , m2      ) 
  , m_sigma2   ( "sigma2"   , "sigma"                      , this , sigma2  )
    //
{
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Rotated2Gaussians::Rotated2Gaussians
( const Ostap::Models::Rotated2Gaussians& right ,  
  const char*                             name  ) 
  : BASE ( right , name ) 
    //
  , m_m01      ( "m01"      , this , right.m_m01     ) 
  , m_sigma1   ( "sigma1"   , this , right.m_sigma1  )
    //
  , m_m02      ( "m02"      , this , right.m_m02     ) 
  , m_sigma2   ( "sigma2"   , this , right.m_sigma2  )
    //
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Rotated2Gaussians::~Rotated2Gaussians(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Rotated2Gaussians*
Ostap::Models::Rotated2Gaussians::clone( const char* name ) const 
{ return new Ostap::Models::Rotated2Gaussians ( *this , name ) ; }
// ============================================================================
void Ostap::Models::Rotated2Gaussians::setPars () const 
{
  m_function.setPhase ( m_phi ) ;
  //
  m_function.signal1 () . setM0      ( m_m01     ) ;
  m_function.signal1 () . setSigma   ( m_sigma1  ) ;
  //
  m_function.signal2 () . setM0      ( m_m02     ) ;
  m_function.signal2 () . setSigma   ( m_sigma2  ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Rotated2Gaussians::evaluate() const 
{
  setPars () ;
  return m_function ( m_x , m_y ) ; 
}
//_____________________________________________________________________________
Double_t Ostap::Models::Rotated2Gaussians::analyticalIntegral 
( Int_t       code       , 
  const char* rangeName  ) const 
{
  assert ( 1 == code || 2 == code || 3 == code ) ;
  setPars () ;
  return BASE::analyticalIntegral ( code , rangeName ) ;
}



// ============================================================================
// Rotated Crystal Ball
// ============================================================================
Ostap::Models::Rotated2CrystalBalls::Rotated2CrystalBalls
( const char*                             name      , 
  const char*                             title     ,
  RooRealVar&                             x         ,
  RooRealVar&                             y         ,
  RooRealVar&                             phi       , // rotation phase 
  // the  first peak 
  RooAbsReal&                             m1         ,
  RooAbsReal&                             sigma1     ,  //
  RooAbsReal&                             alphaL1    ,  // alpha_L
  RooAbsReal&                             nL1        ,  //     n_L - 1
  RooAbsReal&                             alphaR1    ,  // alpha_R - 1
  RooAbsReal&                             nR1        , //     n_R
  // the  second peak 
  RooAbsReal&                             m2         ,
  RooAbsReal&                             sigma2     ,  //
  RooAbsReal&                             alphaL2    ,  // alpha_L
  RooAbsReal&                             nL2        ,  //     n_L - 1
  RooAbsReal&                             alphaR2    ,  // alpha_R - 1
  RooAbsReal&                             nR2        ) //     n_R
  : BASE ( name , title , x , y , phi , 
           CB2 ( m1      .getVal() , 
                 sigma1  .getVal() , 
                 alphaL1 .getVal() , 
                 nL1     .getVal() , 
                 alphaR1 .getVal() , 
                 nR1     .getVal() ) , 
           CB2 ( m2      .getVal() , 
                 sigma2  .getVal() , 
                 alphaL2 .getVal() , 
                 nL2     .getVal() , 
                 alphaR2 .getVal() , 
                 nR2     .getVal() ) ) 
    // first signal 
  , m_m01      ( "m01"      , "mass"                       , this , m1      ) 
  , m_sigma1   ( "sigma1"   , "sigma"                      , this , sigma1  )
  , m_alphaL1  ( "alphaL1"  , "(left) alpha = 1 + |alpha|" , this , alphaL1 ) 
  , m_nL1      ( "nL1"      , "(left) n     = 1 + |n|"     , this ,     nL1 ) 
  , m_alphaR1  ( "alphaR1"  , "(left) alpha = 1 + |alpha|" , this , alphaR1 ) 
  , m_nR1      ( "nR1"      , "(left) n     = 1 + |n|"     , this ,     nR1 ) 
    // second signal 
  , m_m02      ( "m02"      , "mass"                       , this , m2      ) 
  , m_sigma2   ( "sigma2"   , "sigma"                      , this , sigma2  )
  , m_alphaL2  ( "alphaL2"  , "(left) alpha = 1 + |alpha|" , this , alphaL2 ) 
  , m_nL2      ( "nL2"      , "(left) n     = 1 + |n|"     , this ,     nL2 ) 
  , m_alphaR2  ( "alphaR2"  , "(left) alpha = 1 + |alpha|" , this , alphaR2 ) 
  , m_nR2      ( "nR2"      , "(left) n     = 1 + |n|"     , this ,     nR2 ) 
    //
{
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Rotated2CrystalBalls::Rotated2CrystalBalls
( const Ostap::Models::Rotated2CrystalBalls& right ,  
  const char*                                name  ) 
  : BASE ( right , name ) 
    //
  , m_m01      ( "m01"      , this , right.m_m01     ) 
  , m_sigma1   ( "sigma1"   , this , right.m_sigma1  )
  , m_alphaL1  ( "alphaL1"  , this , right.m_alphaL1 ) 
  , m_nL1      ( "nL1"      , this , right.m_nL1     ) 
  , m_alphaR1  ( "alphaR1"  , this , right.m_alphaR1 ) 
  , m_nR1      ( "nR1"      , this , right.m_nR1     ) 
    //
  , m_m02      ( "m02"      , this , right.m_m02     ) 
  , m_sigma2   ( "sigma2"   , this , right.m_sigma2  )
  , m_alphaL2  ( "alphaL2"  , this , right.m_alphaL2 ) 
  , m_nL2      ( "nL2"      , this , right.m_nL2     ) 
  , m_alphaR2  ( "alphaR2"  , this , right.m_alphaR2 ) 
  , m_nR2      ( "nR2"      , this , right.m_nR2     )
    //
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Rotated2CrystalBalls::~Rotated2CrystalBalls(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Rotated2CrystalBalls*
Ostap::Models::Rotated2CrystalBalls::clone( const char* name ) const 
{ return new Ostap::Models::Rotated2CrystalBalls ( *this , name ) ; }
// ============================================================================
void Ostap::Models::Rotated2CrystalBalls::setPars () const 
{
  m_function.setPhase ( m_phi ) ;
  //
  m_function.signal1 () . setM0      ( m_m01     ) ;
  m_function.signal1 () . setSigma   ( m_sigma1  ) ;
  m_function.signal1 () . setAlpha_L ( m_alphaL1 ) ;
  m_function.signal1 () . setAlpha_R ( m_alphaR1 ) ;
  m_function.signal1 () . setN_L     ( m_nL1     ) ;
  m_function.signal1 () . setAlpha_R ( m_alphaR1 ) ;
  m_function.signal1 () . setN_R     ( m_nR1     ) ;
  //
  m_function.signal2 () . setM0      ( m_m02     ) ;
  m_function.signal2 () . setSigma   ( m_sigma2  ) ;
  m_function.signal2 () . setAlpha_L ( m_alphaL2 ) ;
  m_function.signal2 () . setAlpha_R ( m_alphaR2 ) ;
  m_function.signal2 () . setN_L     ( m_nL2     ) ;
  m_function.signal2 () . setAlpha_R ( m_alphaR2 ) ;
  m_function.signal2 () . setN_R     ( m_nR2     ) ;
  //
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Rotated2CrystalBalls::evaluate() const 
{
  setPars () ;
  return m_function ( m_x , m_y ) ; 
}
//_____________________________________________________________________________
Double_t Ostap::Models::Rotated2CrystalBalls::analyticalIntegral 
( Int_t       code       , 
  const char* rangeName  ) const 
{
  assert ( 1 == code || 2 == code || 3 == code ) ;
  setPars () ;
  return BASE::analyticalIntegral ( code , rangeName ) ;
}
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
ClassImp(Ostap::Models::Rotated2Gaussians    )
ClassImp(Ostap::Models::Rotated2CrystalBalls )
// ============================================================================
//                                                                      The END         
// ============================================================================
