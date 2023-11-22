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
  Ostap::Assert ( ::size ( m_phis ) == m_positive.npars () , 
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
  Ostap::Assert ( ::size ( m_phis ) == m_positive.npars () , 
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
  Ostap::Assert ( ::size ( m_phis ) == m_positive.npars () , 
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




// ============================================================================
// 3D Gaussian
// ============================================================================
Ostap::Models::Gauss3D::Gauss3D
( const char* name    , 
  const char* title   ,
  RooAbsReal& x       ,
  RooAbsReal& y       ,
  RooAbsReal& z       ,
  RooAbsReal& muX     ,
  RooAbsReal& muY     ,
  RooAbsReal& muZ     ,
  RooAbsReal& sigmaX  ,
  RooAbsReal& sigmaY  ,
  RooAbsReal& sigmaZ  ,
  RooAbsReal& phi     , 
  RooAbsReal& theta   , 
  RooAbsReal& psi     ) 
  : RooAbsPdf  ( name    , title ) 
  , m_x        ( "x"     , "Observable-X" , this , x      ) 
  , m_y        ( "y"     , "Observable-Y" , this , y      ) 
  , m_z        ( "z"     , "Observable-Z" , this , z      ) 
  , m_muX      ( "muX"   , "x-locaiton"   , this , muX    ) 
  , m_muY      ( "muY"   , "y-locaiton"   , this , muY    ) 
  , m_muZ      ( "muZ"   , "z-locaiton"   , this , muZ    ) 
  , m_sigmaX   ( "sX"    , "sigma-x"      , this , sigmaX ) 
  , m_sigmaY   ( "sY"    , "sigma-y"      , this , sigmaY ) 
  , m_sigmaZ   ( "sZ"    , "sigma-z"      , this , sigmaZ ) 
  , m_phi      ( "phi"   , "rotation"     , this , phi    ) 
  , m_theta    ( "theta" , "rotation"     , this , theta  ) 
  , m_psi      ( "psi"   , "rotation"     , this , psi    ) 
  , m_gauss3D  ()
{
  setPars () ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::Gauss3D::Gauss3D
( const Ostap::Models::Gauss3D&  right ,      
  const char*                              name  ) 
  : RooAbsPdf  ( right , name ) 
    //
  , m_x        ( "x"     , this , right.m_x      ) 
  , m_y        ( "y"     , this , right.m_y      ) 
  , m_z        ( "z"     , this , right.m_z      ) 
  , m_muX      ( "muX"   , this , right.m_muX    ) 
  , m_muY      ( "muY"   , this , right.m_muY    ) 
  , m_muZ      ( "muZ"   , this , right.m_muZ    ) 
  , m_sigmaX   ( "sX"    , this , right.m_sigmaX ) 
  , m_sigmaY   ( "sY"    , this , right.m_sigmaY ) 
  , m_sigmaZ   ( "sZ"    , this , right.m_sigmaZ ) 
  , m_phi      ( "phi"   , this , right.m_phi    ) 
  , m_theta    ( "theta" , this , right.m_theta  ) 
  , m_psi      ( "psi"   , this , right.m_psi    ) 
    //
  , m_gauss3D  ( right.m_gauss3D ) 
{
  setPars () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Models::Gauss3D::~Gauss3D(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Models::Gauss3D*
Ostap::Models::Gauss3D::clone( const char* name ) const 
{ return new Ostap::Models::Gauss3D(*this,name) ; }
// ============================================================================
void Ostap::Models::Gauss3D::setPars () const 
{ 
  m_gauss3D.setMuX    ( m_muX    ) ;
  m_gauss3D.setMuY    ( m_muY    ) ;
  m_gauss3D.setMuZ    ( m_muZ    ) ;
  m_gauss3D.setSigmaX ( m_sigmaX ) ;
  m_gauss3D.setSigmaY ( m_sigmaY ) ;
  m_gauss3D.setSigmaZ ( m_sigmaZ ) ;
  m_gauss3D.setEuler  ( m_phi , m_theta , m_psi ) ;
}
// ============================================================================
// the actual evaluation of function 
// ============================================================================
Double_t Ostap::Models::Gauss3D::evaluate() const 
{
  //
  setPars () ;
  //
  return m_gauss3D ( m_x , m_y , m_z ) ; 
}
// ============================================================================
Int_t Ostap::Models::Gauss3D::getAnalyticalIntegral
( RooArgSet&     allVars      , 
  RooArgSet&     analVars     ,
  const char* /* rangename */ ) const 
{
  // ==========================================================================
  if      ( matchArgs ( allVars , analVars , m_x , m_y , m_z ) ) { return 1 ; }
  else if ( matchArgs ( allVars , analVars , m_x , m_y       ) ) { return 2 ; }
  else if ( matchArgs ( allVars , analVars , m_x       , m_z ) ) { return 3 ; }
  else if ( matchArgs ( allVars , analVars ,       m_y , m_z ) ) { return 4 ; }
  // else if ( matchArgs ( allVars , analVars , m_x             ) ) { return 5 ; }
  // else if ( matchArgs ( allVars , analVars ,       m_y       ) ) { return 6 ; }
  // else if ( matchArgs ( allVars , analVars ,             m_z ) ) { return 7 ; }
  // ==========================================================================
  return 0 ;
}
//_____________________________________________________________________________
Double_t Ostap::Models::Gauss3D::analyticalIntegral 
( Int_t       code       , 
  const char* rangeName  ) const 
{
  //
  assert ( 1 <= code && code <=7 ) ;
  //
  return 
    // 3D-integral
    1 == code ? m_gauss3D.integral
    ( m_x.min ( rangeName ) , m_x.max ( rangeName ) ,
      m_y.min ( rangeName ) , m_y.max ( rangeName ) , 
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) : 
    // 2D-integrals
    2 == code ? m_gauss3D.integrateXY 
    ( m_z , 
      m_x.min ( rangeName ) , m_x.max ( rangeName ) ,
      m_y.min ( rangeName ) , m_y.max ( rangeName ) ) :
    3 == code ? m_gauss3D.integrateXZ 
    ( m_y , 
      m_x.min ( rangeName ) , m_x.max ( rangeName ) ,
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) :
    4 == code ? m_gauss3D.integrateYZ 
    ( m_x , 
      m_y.min ( rangeName ) , m_y.max ( rangeName ) ,
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) :
    // 1D-integrals 
    5 == code ? m_gauss3D.integrateX 
    ( m_y , m_z ,  
      m_x.min ( rangeName ) , m_x.max ( rangeName ) ) :
    6 == code ? m_gauss3D.integrateY 
    ( m_x , m_z ,  
      m_y.min ( rangeName ) , m_y.max ( rangeName ) ) :    
    7 == code ? m_gauss3D.integrateZ 
    ( m_x , m_y ,  
      m_z.min ( rangeName ) , m_z.max ( rangeName ) ) : 0.0 ;
}
// ============================================================================

// ============================================================================
ClassImp(Ostap::Models::Gauss3D              )
ClassImp(Ostap::Models::Poly3DPositive       )
ClassImp(Ostap::Models::Poly3DSymPositive    )
ClassImp(Ostap::Models::Poly3DMixPositive    )
// ============================================================================
//                                                                      The END
// ============================================================================

