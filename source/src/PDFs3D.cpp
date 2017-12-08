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
#include "Ostap/PDFs3D.h"
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
 *  @date   2017-11-21
 */
// ============================================================================
// generic polinomial
// ============================================================================
Ostap::Models::Poly3DPositive::Poly3DPositive 
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  RooRealVar&          z         ,
  const unsigned short nX        , 
  const unsigned short nY        ,
  const unsigned short nZ        ,
  RooArgList&           phis      ) 
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_z        ( "z"       , "Observable-Z" , this , z ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_positive ( nX , nY , nZ , 
                 x.getMin() , x.getMax() , 
                 y.getMin() , y.getMax() ,
                 z.getMin() , z.getMax() ) 
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
                            "Ostap::Poly3DPositive"            , 
                            Ostap::StatusCode::FAILURE         ) ; }
  
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Poly3DPositive::Poly3DPositive
( const Ostap::Models::Poly3DPositive&  right ,      
  const char*                              name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_y        ( "y"      , this , right.m_y     ) 
  , m_z        ( "z"      , this , right.m_z     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_positive ( right.m_positive ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Poly3DPositive::~Poly3DPositive(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Poly3DPositive*
Ostap::Models::Poly3DPositive::clone( const char* name ) const 
{ return new Ostap::Models::Poly3DPositive(*this,name) ; }
// ============================================================================
void Ostap::Models::Poly3DPositive::setPars () const 
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
Double_t Ostap::Models::Poly3DPositive::evaluate() const 
{
  setPars () ;
  return m_positive ( m_x , m_y , m_z ) ; 
}
// ============================================================================
Int_t Ostap::Models::Poly3DPositive::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  // ==========================================================================
  if      ( matchArgs ( allVars , analVars , m_x , m_y , m_z ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x , m_y       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars , m_x       , m_z ) ) { return 3 ; }
  else if ( matchArgs ( allVars , analVars ,       m_y , m_z ) ) { return 4 ; }
  else if ( matchArgs ( allVars , analVars , m_x             ) ) { return 5 ; }
  else if ( matchArgs ( allVars , analVars ,       m_y       ) ) { return 6 ; }
  else if ( matchArgs ( allVars , analVars ,             m_z ) ) { return 7 ; }
  // ==========================================================================
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Poly3DPositive::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( 1 <= code && code <=7 ) ;
  //
  return 
    // 3D-integral
    1 == code ? m_positive.integral
    ( m_x.min ( rangeName ) , m_x.max ( rangeName ) ,
      m_y.min ( rangeName ) , m_y.max ( rangeName ) , 
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) : 
    // 2D-integrals
    2 == code ? m_positive.integrateXY 
    ( m_z , 
      m_x.min ( rangeName ) , m_x.max ( rangeName ) ,
      m_y.min ( rangeName ) , m_y.max ( rangeName ) ) :
    3 == code ? m_positive.integrateXZ 
    ( m_y , 
      m_x.min ( rangeName ) , m_x.max ( rangeName ) ,
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) :
    4 == code ? m_positive.integrateYZ 
    ( m_x , 
      m_y.min ( rangeName ) , m_y.max ( rangeName ) ,
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) :
    // 1D-integrals 
    5 == code ? m_positive.integrateX 
    ( m_y , m_z ,  
      m_x.min ( rangeName ) , m_x.max ( rangeName ) ) :
    6 == code ? m_positive.integrateY 
    ( m_x , m_z ,  
      m_y.min ( rangeName ) , m_y.max ( rangeName ) ) :    
    7 == code ? m_positive.integrateZ 
    ( m_x , m_y ,  
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) : 0.0 ;
}




// ============================================================================
// symmetric polinomial
// ============================================================================
Ostap::Models::Poly3DSymPositive::Poly3DSymPositive 
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  RooRealVar&          z         ,
  const unsigned short n         , 
  RooArgList&           phis      ) 
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_z        ( "z"       , "Observable-Z" , this , z ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_positive ( n ,
                 std::min ( x.getMin() , std::min ( y.getMin() , z.getMin() ) ) ,
                 std::max ( x.getMax() , std::max ( y.getMax() , z.getMax() ) ) )
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
                            "Ostap::Poly3DSymPositive"         , 
                            Ostap::StatusCode::FAILURE         ) ; }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Poly3DSymPositive::Poly3DSymPositive
( const Ostap::Models::Poly3DSymPositive&  right ,      
  const char*                              name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_y        ( "y"      , this , right.m_y     ) 
  , m_z        ( "z"      , this , right.m_z     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_positive ( right.m_positive ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Poly3DSymPositive::~Poly3DSymPositive(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Poly3DSymPositive*
Ostap::Models::Poly3DSymPositive::clone( const char* name ) const 
{ return new Ostap::Models::Poly3DSymPositive(*this,name) ; }
// ============================================================================
void Ostap::Models::Poly3DSymPositive::setPars () const 
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
Double_t Ostap::Models::Poly3DSymPositive::evaluate() const 
{
  setPars () ;
  return m_positive ( m_x , m_y , m_z ) ; 
}
// ============================================================================
Int_t Ostap::Models::Poly3DSymPositive::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  // ==========================================================================
  if      ( matchArgs ( allVars , analVars , m_x , m_y , m_z ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x , m_y       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars , m_x       , m_z ) ) { return 3 ; }
  else if ( matchArgs ( allVars , analVars ,       m_y , m_z ) ) { return 4 ; }
  else if ( matchArgs ( allVars , analVars , m_x             ) ) { return 5 ; }
  else if ( matchArgs ( allVars , analVars ,       m_y       ) ) { return 6 ; }
  else if ( matchArgs ( allVars , analVars ,             m_z ) ) { return 7 ; }
  // ==========================================================================
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Poly3DSymPositive::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( 1 <= code && code <=7 ) ;
  //
  return 
    // 3D-integral
    1 == code ? m_positive.integral
    ( m_x.min ( rangeName ) , m_x.max ( rangeName ) ,
      m_y.min ( rangeName ) , m_y.max ( rangeName ) , 
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) : 
    // 2D-integrals
    2 == code ? m_positive.integrateXY 
    ( m_z , 
      m_x.min ( rangeName ) , m_x.max ( rangeName ) ,
      m_y.min ( rangeName ) , m_y.max ( rangeName ) ) :
    3 == code ? m_positive.integrateXZ 
    ( m_y , 
      m_x.min ( rangeName ) , m_x.max ( rangeName ) ,
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) :
    4 == code ? m_positive.integrateYZ 
    ( m_x , 
      m_y.min ( rangeName ) , m_y.max ( rangeName ) ,
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) :
    // 1D-integrals 
    5 == code ? m_positive.integrateX 
    ( m_y , m_z ,  
      m_x.min ( rangeName ) , m_x.max ( rangeName ) ) :
    6 == code ? m_positive.integrateY 
    ( m_x , m_z ,  
      m_y.min ( rangeName ) , m_y.max ( rangeName ) ) :    
    7 == code ? m_positive.integrateZ 
    ( m_x , m_y ,  
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) : 0.0 ;
}
// ============================================================================

// ============================================================================
// mixed symmetry polinomial
// ============================================================================
Ostap::Models::Poly3DMixPositive::Poly3DMixPositive 
( const char*          name      , 
  const char*          title     ,
  RooRealVar&          x         ,
  RooRealVar&          y         ,
  RooRealVar&          z         ,
  const unsigned short n         , 
  const unsigned short nz        , 
  RooArgList&          phis      ) 
  : RooAbsPdf  ( name , title ) 
  , m_x        ( "x"       , "Observable-X" , this , x ) 
  , m_y        ( "y"       , "Observable-Y" , this , y ) 
  , m_z        ( "z"       , "Observable-Z" , this , z ) 
  , m_phis     ( "phis"    , "Coefficients" , this     )
    //
  , m_positive ( n ,  nz ,
                 std::min ( x.getMin() , y.getMin() ) ,
                 std::max ( x.getMax() , y.getMax() ) ,
                 z.getMin() , z.getMax() )
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
                            "Ostap::Poly3DMixPositive"         , 
                            Ostap::StatusCode::FAILURE         ) ; }
  //
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Poly3DMixPositive::Poly3DMixPositive
( const Ostap::Models::Poly3DMixPositive&  right ,      
  const char*                              name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"      , this , right.m_x     ) 
  , m_y        ( "y"      , this , right.m_y     ) 
  , m_z        ( "z"      , this , right.m_z     ) 
  , m_phis     ( "phis"   , this , right.m_phis  ) 
    //
  , m_positive ( right.m_positive ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Poly3DMixPositive::~Poly3DMixPositive(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Poly3DMixPositive*
Ostap::Models::Poly3DMixPositive::clone( const char* name ) const 
{ return new Ostap::Models::Poly3DMixPositive(*this,name) ; }
// ============================================================================
void Ostap::Models::Poly3DMixPositive::setPars () const 
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
Double_t Ostap::Models::Poly3DMixPositive::evaluate() const 
{
  setPars () ;
  return m_positive ( m_x , m_y , m_z ) ; 
}
// ============================================================================
Int_t Ostap::Models::Poly3DMixPositive::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  // ==========================================================================
  if      ( matchArgs ( allVars , analVars , m_x , m_y , m_z ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x , m_y       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars , m_x       , m_z ) ) { return 3 ; }
  else if ( matchArgs ( allVars , analVars ,       m_y , m_z ) ) { return 4 ; }
  else if ( matchArgs ( allVars , analVars , m_x             ) ) { return 5 ; }
  else if ( matchArgs ( allVars , analVars ,       m_y       ) ) { return 6 ; }
  else if ( matchArgs ( allVars , analVars ,             m_z ) ) { return 7 ; }
  // ==========================================================================
  return 0 ;
}
// ============================================================================
Double_t Ostap::Models::Poly3DMixPositive::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  assert ( 1 <= code && code <=7 ) ;
  //
  return 
    // 3D-integral
    1 == code ? m_positive.integral
    ( m_x.min ( rangeName ) , m_x.max ( rangeName ) ,
      m_y.min ( rangeName ) , m_y.max ( rangeName ) , 
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) : 
    // 2D-integrals
    2 == code ? m_positive.integrateXY 
    ( m_z , 
      m_x.min ( rangeName ) , m_x.max ( rangeName ) ,
      m_y.min ( rangeName ) , m_y.max ( rangeName ) ) :
    3 == code ? m_positive.integrateXZ 
    ( m_y , 
      m_x.min ( rangeName ) , m_x.max ( rangeName ) ,
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) :
    4 == code ? m_positive.integrateYZ 
    ( m_x , 
      m_y.min ( rangeName ) , m_y.max ( rangeName ) ,
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) :
    // 1D-integrals 
    5 == code ? m_positive.integrateX 
    ( m_y , m_z ,  
      m_x.min ( rangeName ) , m_x.max ( rangeName ) ) :
    6 == code ? m_positive.integrateY 
    ( m_x , m_z ,  
      m_y.min ( rangeName ) , m_y.max ( rangeName ) ) :    
    7 == code ? m_positive.integrateZ 
    ( m_x , m_y ,  
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) : 0.0 ;
}
// ============================================================================
ClassImp(Ostap::Models::Poly3DPositive       )
ClassImp(Ostap::Models::Poly3DSymPositive    )
ClassImp(Ostap::Models::Poly3DMixPositive    )
// ============================================================================
//                                                                      The END
// ============================================================================

