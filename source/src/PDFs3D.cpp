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
// Local
// ============================================================================
#include "local_roofit.h"
// ============================================================================
/** @file 
 *  Implementation file for namespace Ostap::Models
 *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
 *  @date   2017-11-21
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
  template<class TX, class TY, class TZ , class FUN>
  void compute_XYZ ( RooSpan<double> output , FUN& fun , TX x , TY y , TZ z ) 
  {
    const int n = output.size();
    for ( int i = 0 ; i < n ; ++i ) 
    { output [ i ] = fun ( x [ i ] , y [ i ] , z [ i ] ) ; }
  }
  // ==========================================================================
}
#endif
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
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::Poly3DPositive"          ) ;
  //
  Ostap::Assert ( m_phis.size() == m_positive.npars ()     , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Poly3DPositive"          ) ;
  //
  setPars () ;
  //
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
{ ::set_pars ( m_phis , m_positive ) ; }
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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Poly3DPositive::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x = m_x . getValBatch ( begin , batchSize ) ;
  auto y = m_y . getValBatch ( begin , batchSize ) ;
  auto z = m_z . getValBatch ( begin , batchSize ) ;
  //
  const bool ex = x.empty()  ;
  const bool ey = y.empty()  ;
  const bool ez = z.empty()  ;
  //
  if ( ex && ey && ez ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars() ;
  //
  if      ( !ex &&  ey &&  ez ) 
  { compute_XYZ ( output , m_positive ,        x   , BA ( m_y ) , BA ( m_z ) ) ; }
  else if (  ex && !ey &&  ez ) 
  { compute_XYZ ( output , m_positive , BA ( m_x ) ,        y   , BA ( m_z ) ) ; }
  else if (  ex &&  ey && !ez ) 
  { compute_XYZ ( output , m_positive , BA ( m_x ) , BA ( m_y ) ,        z   ) ; }
  else if (  ex && !ey && !ez ) 
  { compute_XYZ ( output , m_positive , BA ( m_x ) ,        y   ,        z   ) ; }
  else if ( !ex &&  ey && !ez ) 
  { compute_XYZ ( output , m_positive ,        x   , BA ( m_y ) ,        z   ) ; }
  else if ( !ex && !ey &&  ez ) 
  { compute_XYZ ( output , m_positive ,        x   ,        y   , BA ( m_z ) ) ; }
  else 
  { compute_XYZ ( output , m_positive ,        x   ,        y   ,        z   ) ; }
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================

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
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::Poly3DSymPositive"       ) ;
  //
  Ostap::Assert ( m_phis.size() == m_positive.npars ()     , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Poly3DSymPositive"       ) ;
  //
  setPars () ;
  //
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
{ ::set_pars ( m_phis , m_positive ) ; }
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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Poly3DSymPositive::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x = m_x . getValBatch ( begin , batchSize ) ;
  auto y = m_y . getValBatch ( begin , batchSize ) ;
  auto z = m_z . getValBatch ( begin , batchSize ) ;
  //
  const bool ex = x.empty()  ;
  const bool ey = y.empty()  ;
  const bool ez = z.empty()  ;
  //
  if ( ex && ey && ez ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars() ;
  //
  if      ( !ex &&  ey &&  ez ) 
  { compute_XYZ ( output , m_positive ,        x   , BA ( m_y ) , BA ( m_z ) ) ; }
  else if (  ex && !ey &&  ez ) 
  { compute_XYZ ( output , m_positive , BA ( m_x ) ,        y   , BA ( m_z ) ) ; }
  else if (  ex &&  ey && !ez ) 
  { compute_XYZ ( output , m_positive , BA ( m_x ) , BA ( m_y ) ,        z   ) ; }
  else if (  ex && !ey && !ez ) 
  { compute_XYZ ( output , m_positive , BA ( m_x ) ,        y   ,        z   ) ; }
  else if ( !ex &&  ey && !ez ) 
  { compute_XYZ ( output , m_positive ,        x   , BA ( m_y ) ,        z   ) ; }
  else if ( !ex && !ey &&  ez ) 
  { compute_XYZ ( output , m_positive ,        x   ,        y   , BA ( m_z ) ) ; }
  else 
  { compute_XYZ ( output , m_positive ,        x   ,        y   ,        z   ) ; }
  //
  return output ;
}
// ============================================================================
#endif
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
  ::copy_real   ( phis , m_phis , "Invalid phi-parameter!" ,
                  "Ostap::Models::Poly2DMixPositive"       ) ;
  //
  Ostap::Assert ( m_phis.size() == m_positive.npars ()     , 
                  "Widths/#channels mismatch"              , 
                  "Ostap::Models::Poly3DMixPositive"       ) ;
  //
  setPars () ;
  //
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
{ ::set_pars ( m_phis , m_positive ) ; }
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
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
// ============================================================================
RooSpan<double> 
Ostap::Models::Poly3DMixPositive::evaluateBatch 
( std::size_t begin     , 
  std::size_t batchSize ) const 
{ 
  // 
  auto x = m_x . getValBatch ( begin , batchSize ) ;
  auto y = m_y . getValBatch ( begin , batchSize ) ;
  auto z = m_z . getValBatch ( begin , batchSize ) ;
  //
  const bool ex = x.empty()  ;
  const bool ey = y.empty()  ;
  const bool ez = z.empty()  ;
  //
  if ( ex && ey && ez ) { return {} ; }
  //
  auto output = _batchData.makeWritableBatchUnInit ( begin , batchSize ) ;
  //
  setPars() ;
  //
  if      ( !ex &&  ey &&  ez ) 
  { compute_XYZ ( output , m_positive ,        x   , BA ( m_y ) , BA ( m_z ) ) ; }
  else if (  ex && !ey &&  ez ) 
  { compute_XYZ ( output , m_positive , BA ( m_x ) ,        y   , BA ( m_z ) ) ; }
  else if (  ex &&  ey && !ez ) 
  { compute_XYZ ( output , m_positive , BA ( m_x ) , BA ( m_y ) ,        z   ) ; }
  else if (  ex && !ey && !ez ) 
  { compute_XYZ ( output , m_positive , BA ( m_x ) ,        y   ,        z   ) ; }
  else if ( !ex &&  ey && !ez ) 
  { compute_XYZ ( output , m_positive ,        x   , BA ( m_y ) ,        z   ) ; }
  else if ( !ex && !ey &&  ez ) 
  { compute_XYZ ( output , m_positive ,        x   ,        y   , BA ( m_z ) ) ; }
  else 
  { compute_XYZ ( output , m_positive ,        x   ,        y   ,        z   ) ; }
  //
  return output ;
}
// ============================================================================
#endif
// ============================================================================

// ============================================================================
ClassImp(Ostap::Models::Poly3DPositive       )
ClassImp(Ostap::Models::Poly3DSymPositive    )
ClassImp(Ostap::Models::Poly3DMixPositive    )
// ============================================================================
//                                                                      The END
// ============================================================================

