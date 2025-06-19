// ============================================================================
// Include files 
// ============================================================================
//  STD&STL
// ============================================================================
#include <utility>
#include <algorithm>
#include <limits>
#include <cmath> 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Quantile.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
#include "status_codes.h"
#include "local_math.h" 
// ============================================================================
/** @file  
 *  Implementation file for class Ostap::Math::P2Quantile
 *  @date 2020-05-21 
 *  @author Vanya BELYAEV Ivan/Belyaev@itep.ru
 */
// ============================================================================
namespace
{
  const    std::string s_INVALID  = "P2Quantile: invalid quantile!" ;
}
// ============================================================================
// contrcuctor 
// =============================================================================
Ostap::Math::Quantile::Quantile
( const double p )
: m_p    ( p ) 
, m_N    ( 0 )
, m_q    {   }
, m_n    { 0 , 1 , 2 , 3 , 4 }
, m_np   {}
, m_dnp  {}
{
  Ostap::Assert ( s_zero  ( m_p     ) || 
                  s_equal ( m_p , 1 ) || 
                  ( 0 <= m_p  &&  m_p <= 1 ) ,  
                  "Invalid quantile!" ,  
                  "Ostap::Math::QuantileP2" , 
                  INVALID_QUANTILE , __FILE__ , __LINE__ ) ;
  //
  if      ( s_zero  ( m_p     ) ) { m_p = 0 ; }
  else if ( s_equal ( m_p , 1 ) ) { m_p = 1 ; }
  //
  m_np  = std::array<double,5>{ 0 , 2   * m_p , 4  * m_p , 2 * m_p + 2     , 4 * p  } ;
  m_dnp = std::array<double,5>{ 0 , 0.5 * m_p ,      m_p , ( 1 + m_p ) / 2 , 1.0    } ;
  //
}
// =============================================================================
// update the counter 
// =============================================================================
Ostap::Math::Quantile&
Ostap::Math::Quantile::add  
( const  double value ) 
{
  if ( !std::isfinite ( value ) ) { return *this ; }
  //
  if ( m_N  < 5 )
  {
    m_q [ m_N ]  = value ;
    ++m_N ;
    if ( 5 == m_N ) { std::sort ( m_q.begin() , m_q.end () ) ; }
    return *this ; 
  }
  //
  ++m_N ;
  //
  const std::size_t k =     
    value < m_q [ 0 ] ? 0 :
    value < m_q [ 1 ] ? 0 :
    value < m_q [ 2 ] ? 1 :
    value < m_q [ 3 ] ? 2 : 3 ;
  // 
  if      ( 0 == k ) { m_q [ 0 ] = std::min ( m_q [ 0 ] ,  value ) ; }
  else if ( 3 == k ) { m_q [ 4 ] = std::max ( m_q [ 4 ] ,  value ) ; } 
  // 
  for ( std::size_t i = k + 1 ; i < 5 ; ++i ) { ++m_n [ i ]                ; }
  for ( std::size_t i = 0     ; i < 5 ; ++i ) { m_np  [ i ] += m_dnp [ i ] ; }
  //
  
  auto Parabolic = [this] ( const int i , const double d ) -> double 
  {
    const double n2  = this->m_n [ i + 1 ] - this->m_n [ i - 1 ] ;
    const double nm1 = this->m_n [ i     ] - this->m_n [ i - 1 ] ; 
    const double np1 = this->m_n [ i + 1 ] - this->m_n [ i     ] ; 
    const double qp1 = this->m_q [ i + 1 ] - this->m_q [ i     ] ; 
    const double qm1 = this->m_q [ i     ] - this->m_q [ i - 1 ] ; 
    return this->m_q [ i ] +  d / n2 * ( ( nm1 + d ) * qp1 / np1 + 
					 ( np1 + d ) * qm1 / nm1 ) ; 
  }; 
  
  auto Linear = [this] ( const std::size_t i , const int d ) -> double 
  { 
    const double dq = this->m_q [ i + d ] = this->m_q [ i ] ;
    const double dn = this->m_n [ i + d ] + this->m_n [ i ] ;
    return this->m_q [ i ] + d * dq / dn ;
  } ; 
  
  for ( std::size_t i = 1 ; i <= 3 ; ++i )
    {
      const double d = m_np [ i ] - m_n [ i ] ;
      if ( ( (  1 <= d ) && (  1 < m_n [ i + 1 ] - m_n [ i ] ) )  || 
           ( ( -1 >= d ) && ( -1 > m_n [ i - 1 ] - m_n [ i ] ) ) )
        {
          /**  
          const int delta = d < 0 ? -1 : 1 ;
          const double qs = Parabolic ( i , delta );
          if ( m_q [ i - 1] < qs && qs < m_q[i + 1] ) { m_q [ i ] = qs; } 
          else     { m_q [ i ] = Linear ( i , delta );
          m_n[i] += delta ;
          }
          */
        }
    }
  return *this ;
}
// ============================================================================
//                                                                      The END  
// ============================================================================
