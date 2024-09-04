// ============================================================================
// Include files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <cmath>
#include <array>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/StatEntity.h"
#include "Ostap/WStatEntity.h"
#include "Ostap/Covariance.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::Covariance
 *  @see  Ostap::Math::Covariance
 *  @see  Ostap::Math::WCovariance
 *  @author  Vanya Belyaev Ivan.Belyaev@cern.ch
 *  @date 2024-07-22 
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  const Ostap::Math::Zero    <double> s_zero  {} ;
  const Ostap::Math::Equal_To<double> s_equal {} ; // equality criteria for doubles
  // ==========================================================================
}
// ============================================================================
// constructor from two counters and second moment
// ============================================================================
Ostap::Math::Covariance::Covariance 
( const Ostap::Math::Covariance::Counter& cnt1 ,
  const Ostap::Math::Covariance::Counter& cnt2 ,
  const double                            corr )
  : m_cnt1  ( cnt1 )
  , m_cnt2  ( cnt2 )
  , m_cov2m ( corr ) 
{
  Ostap::Assert ( cnt1.n () == cnt2.n()    , 
		  "Ostap::Math::Covariance: invalid counters!" ,
		  "Ostap::Math::Covariance" ) ;
  //
  const double acorr = std::abs ( corr ) ;
  Ostap::Assert ( acorr <= 1.0 || s_equal ( acorr , 1.0 ) ,
		  "Ostap::Math::Covariance: invalid correlation!" ,
		  "Ostap::Math::Covariance" ) ;
  //
  if      (  1 < corr ) { m_cov2m =  1 ; }
  else if ( -1 > corr ) { m_cov2m = -1 ; }
  //
  if ( empty() ) { m_cov2m = 0 ; }
  else
    {
      const double covx  = m_cnt1.mu2() ;
      const double covy  = m_cnt2.mu2() ;
      if ( s_zero ( covx ) || s_zero ( covy ) ) { m_cov2m = 0 ; }
      Ostap::Assert  ( 0<= covx & 0 <= covy     , 
		       "Ostap::Math::Covariance: invalid variances!" ,
		       "Ostap::Math::Covariance" ) ;
      m_cov2m *= std::sqrt ( covx * covy ) ;
      m_cov2m *= n () ;
    }
}
// ============================================================================
// add two counters 
// ============================================================================
Ostap::Math::Covariance&
Ostap::Math::Covariance::add
( const Ostap::Math::Covariance& right ) 
{
  if      ( right.empty () ) {                 return *this ; }
  else if (       empty () ) { *this = right ; return *this ; }
  //
  const unsigned long long wA =       n () ;
  const unsigned long long wB = right.n () ;
  const unsigned long long wW = wA + wB    ;
  //
  const double xA =       m_cnt1.mean () ;
  const double xB = right.m_cnt1.mean () ;
  const double yA =       m_cnt2.mean () ;
  const double yB = right.m_cnt2.mean () ;
  //
  m_cov2m += right.m_cov2m + ( xB - xA ) * ( yB - yA ) * wA * wB / wW ;
  //
  m_cnt1 += right.m_cnt1 ;
  m_cnt2 += right.m_cnt2 ;
  //
  return *this ;
}
// ==========================================================
double Ostap::Math::Covariance::correlation () const
{
  if ( empty() ) { return 0; }
  //
  const double xv = m_cnt1.mu2() ;
  const double yv = m_cnt2.mu2() ; 
  const double cc = covariance() ;
  //
  if ( s_zero ( cc )    ) { return 0 ; }
  if ( xv < 0 || yv < 0 ) { return 0 ; }
  //
  return cc / std::sqrt ( xv * yv ) ;
}
// ============================================================================
// get the covariance matrix
// ============================================================================
Ostap::Math::Covariance::Matrix
Ostap::Math::covariance
( const Ostap::Math::Covariance& c )
{
  //
  if ( c.empty() ) { return Ostap::Math::Covariance::Matrix() ; }
  //
  static std::array<double,3> buffer; 
  buffer [ 0 ] = c.counter1().mu2() ;
  buffer [ 1 ] = c.covariance () ;
  buffer [ 2 ] = c.counter2().mu2() ;
  return Ostap::Math::Covariance::Matrix ( buffer.begin() , buffer.end() ) ;
}
// ============================================================================
// reset counters
// ============================================================================
void Ostap::Math::Covariance::reset ()
{
  m_cnt1.reset() ;
  m_cnt2.reset() ;
  m_cov2m = 0 ;
}
// ============================================================================
// get the correlation  matrix
// ============================================================================
Ostap::Math::Covariance::Matrix
Ostap::Math::correlation
( const Ostap::Math::Covariance& c )
{
  //
  if ( c.empty() ) { return Ostap::Math::Covariance::Matrix() ; }
  //
  static std::array<double,3> buffer; 
  buffer [ 0 ] = 1 ;
  buffer [ 1 ] = c.correlation () ;
  buffer [ 2 ] = 1 ;
  return Ostap::Math::Covariance::Matrix ( buffer.begin() , buffer.end() ) ;
}
// ============================================================================
// constructor from two counters and second moment
// ============================================================================
Ostap::Math::WCovariance::WCovariance 
( const Ostap::Math::WCovariance::Counter& cnt1 ,
  const Ostap::Math::WCovariance::Counter& cnt2 ,
  const double                             corr )
  : m_cnt1  ( cnt1 )
  , m_cnt2  ( cnt2 )
  , m_cov2m ( corr ) 
{
  Ostap::Assert ( cnt1.weights () == cnt2.weights() , 
		  "Ostap::Math::WCovariance: invalid counters!" ,
		  "Ostap::Math::WCovariance" ) ;
  //
  const double acorr = std::abs ( corr ) ;
  Ostap::Assert ( acorr <= 1.0 || s_equal ( acorr , 1.0 ) ,
		  "Ostap::Math::Covariance: invalid correlation!" ,
		  "Ostap::Math::Covariance" ) ;
  //
  if      (  1 < corr ) { m_cov2m =  1 ; }
  else if ( -1 > corr ) { m_cov2m = -1 ; }
  //
  if ( empty() ) { m_cov2m = 0 ; }
  else
    {
      const double covx  = m_cnt1.mu2() ;
      const double covy  = m_cnt2.mu2() ;
      if ( s_zero ( covx ) || s_zero ( covy ) ) { m_cov2m = 0 ; }
      Ostap::Assert  ( 0<= covx & 0 <= covy     , 
		       "Ostap::Math::Covariance: invalid variances!" ,
		       "Ostap::Math::Covariance" ) ;
      m_cov2m *= std::sqrt ( covx * covy ) ;
      m_cov2m *= w () ;
    }
}
// ============================================================================
// add two counters 
// ============================================================================
Ostap::Math::WCovariance&
Ostap::Math::WCovariance::add
( const Ostap::Math::WCovariance& right ) 
{
  if      ( right.empty () ) {                 return *this ; }
  else if (       empty () ) { *this = right ; return *this ; }
  //
  const long double wA =       w () ;
  const long double wB = right.w () ;
  const long double wW = wA + wB    ;
  //
  m_cnt1 += right.m_cnt1 ;
  m_cnt2 += right.m_cnt2 ;
  //
  const double xA =       m_cnt1.mean () ;
  const double xB = right.m_cnt1.mean () ;
  const double yA =       m_cnt2.mean () ;
  const double yB = right.m_cnt2.mean () ;
  //
  m_cov2m += right.m_cov2m + ( xB - xA ) * ( yB - yA ) * wA * wB / wW ;
  //
  m_cnt1 += right.m_cnt1 ;
  m_cnt2 += right.m_cnt2 ;
  //
  return *this ;
}
// ============================================================================
// reset counters
// ============================================================================
void Ostap::Math::WCovariance::reset ()
{
  m_cnt1.reset() ;
  m_cnt2.reset() ;
  m_cov2m = 0 ;
}
// ==========================================================
double Ostap::Math::WCovariance::correlation () const
{
  if ( empty() ) { return 0; }
  //
  const double xv = m_cnt1.mu2() ;
  const double yv = m_cnt2.mu2() ; 
  const double cc = covariance() ;
  //
  if ( s_zero ( cc )    ) { return 0 ; }
  if ( xv < 0 || yv < 0 ) { return 0 ; }
  //
  return cc / std::sqrt ( xv * yv ) ;
}
// ============================================================================
// get the covariance matrix
// ============================================================================
Ostap::Math::WCovariance::Matrix
Ostap::Math::covariance
( const Ostap::Math::WCovariance& c )
{
  //
  if ( c.empty() ) { return Ostap::Math::WCovariance::Matrix() ; }
  //
  static std::array<double,3> buffer; 
  buffer [ 0 ] = c.counter1().mu2() ;
  buffer [ 1 ] = c.covariance () ;
  buffer [ 2 ] = c.counter2().mu2() ;
  return Ostap::Math::WCovariance::Matrix ( buffer.begin() , buffer.end() ) ;
}
// ============================================================================
// get the correlation  matrix
// ============================================================================
Ostap::Math::WCovariance::Matrix
Ostap::Math::correlation
( const Ostap::Math::WCovariance& c )
{
  //
  if ( c.empty() ) { return Ostap::Math::WCovariance::Matrix() ; }
  //
  static std::array<double,3> buffer; 
  buffer [ 0 ] = 1 ;
  buffer [ 1 ] = c.correlation () ;
  buffer [ 2 ] = 1 ;
  return Ostap::Math::WCovariance::Matrix ( buffer.begin() , buffer.end() ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
