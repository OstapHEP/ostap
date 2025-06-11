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
// GSL
// ============================================================================
#include "gsl/gsl_rstat.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/P2Quantile.h"
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
Ostap::Math::QuantileP2::QuantileP2
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
                  ( 1 <= m_p  &&  m_p <= 1 ) ,  
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
Ostap::Math::QuantileP2&
Ostap::Math::QuantileP2::add  
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
       return this->m_q [ i ] +  d / n2 * ( 
             ( nm1 + d ) * qp1 / np1 + 
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
/* Standard constructor for quantile
 *  @param p quatile \f$ 0 < p > 1 \f$
 */
// ============================================================================
Ostap::Math::GSL::P2Quantile::P2Quantile
( const double p )
  : m_ws   ( nullptr )
  , m_p    ( p )  
{
  Ostap::Assert ( 0 < p && p < 1 , s_INVALID ) ;
}
// ======================================================================
// copy contructor 
// ======================================================================
Ostap::Math::GSL::P2Quantile::P2Quantile
( const Ostap::Math::GSL::P2Quantile& right )
  : m_ws   ( nullptr   )
  , m_p    ( right.m_p )  
{
  if ( right.m_ws )
  {
    m_ws = gsl_rstat_quantile_alloc ( m_p ) ;
    //
    std::copy ( right.m_ws->q    , right.m_ws->q    + 5 , m_ws->q    ) ;
    std::copy ( right.m_ws->npos , right.m_ws->npos + 5 , m_ws->npos ) ;
    std::copy ( right.m_ws->np   , right.m_ws->np   + 5 , m_ws->np   ) ;
    std::copy ( right.m_ws->dnp  , right.m_ws->dnp  + 5 , m_ws->dnp  ) ;
    m_ws -> n = right.m_ws-> n ;
  }
}
// ============================================================================
// move constructor
// ============================================================================
Ostap::Math::GSL::P2Quantile::P2Quantile
( Ostap::Math::GSL::P2Quantile&& right )
  : m_ws   ( nullptr   )
  , m_p    ( right.m_p )  
{
  std::swap ( m_ws  , right.m_ws  ) ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::GSL::P2Quantile::~P2Quantile()
{ if  ( m_ws ) { gsl_rstat_quantile_free ( m_ws ) ; m_ws = nullptr ; } }
// ============================================================================
// copy assignement
// ============================================================================
Ostap::Math::GSL::P2Quantile&
Ostap::Math::GSL::P2Quantile::operator=
( const Ostap::Math::GSL::P2Quantile&  right )
{
  if ( &right == this  ) { return *this ; }
  if ( m_ws != nullptr ) { gsl_rstat_quantile_free ( m_ws ) ; m_ws = nullptr ; }
  m_p = right.m_p ;
  //
  if ( right.m_ws )
  {
    m_ws = gsl_rstat_quantile_alloc ( m_p ) ;
    //
    std::copy ( right.m_ws->q    , right.m_ws->q    + 5 , m_ws->q    ) ;
    std::copy ( right.m_ws->npos , right.m_ws->npos + 5 , m_ws->npos ) ;
    std::copy ( right.m_ws->np   , right.m_ws->np   + 5 , m_ws->np   ) ;
    std::copy ( right.m_ws->dnp  , right.m_ws->dnp  + 5 , m_ws->dnp  ) ;
    //
    m_ws -> n = right.m_ws-> n ;
  }
  return *this ;
}
// ============================================================================
// move  assignement
// ============================================================================
Ostap::Math::GSL::P2Quantile&
Ostap::Math::GSL::P2Quantile::operator=
(       Ostap::Math::GSL::P2Quantile&& right )
{
  if ( &right == this  ) { return *this ; }
  std::swap ( m_ws , right.m_ws ) ;
  std::swap ( m_p  , right.m_p  ) ;
  return *this ;
}
// ============================================================================
// get the quantile value
// ============================================================================
double Ostap::Math::GSL::P2Quantile::value() const
{ return m_ws == nullptr ? s_INFINITY : gsl_rstat_quantile_get ( m_ws ) ; }
// ============================================================================
//                                                                      The END  
// ============================================================================
