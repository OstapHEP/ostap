// ============================================================================
#define OSTAP_STATENTITY_CPP 1
// ============================================================================
// include files
// ============================================================================
// STD & STL
// ============================================================================
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>
#include <tuple>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatEntity.h"
// ============================================================================
// Local 
// ============================================================================
#include "format.h"
#include "Exception.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::StatEntity
 *  @date 26/06/2001
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
/* The full contructor from all important values
 * @see StatEntity::format
 * @param entries number of entries
 * @param mu   the mean value 
 * @param mu2  the second central moment/variance/dispersion 
 * @param minv the minimum value 
 * @param maxv the maximum value 
 */  
// ============================================================================
Ostap::StatEntity::StatEntity
( const unsigned long entries ,
  const double        mu      ,
  const double        mu2     ,
  const double        minv    ,
  const double        maxv    )
  : m_n   ( entries )
  , m_mu  ( mu      )
  , m_mu2 ( mu2     )
  , m_min ( minv    )
  , m_max ( maxv    )
{
  // reset empty counter 
  if ( 0 == m_n ) { reset() ; }
  else 
  {
    Ostap::Assert ( m_min <= mu && mu <= m_max , 
                    "Ostap::StatEntity: invalid minv/mu/maxv" ,
                    "Ostap::StatEntity" ) ; 
    Ostap::Assert ( 0 <= m_mu2 , 
                    "Ostap::StatEntity: invalid second moment",
                    "Ostap::StatEntity" ) ;
  }
}
// ============================================================================
/* add a value : the main method 
 * @see https://arxiv.org/abs/1510.04923v1
 * @param value (INPUT) value to be added
 * @return self reference 
 */
// ============================================================================
Ostap::StatEntity& 
Ostap::StatEntity::add ( const double value ) 
{
  if ( 0 == m_n ) 
  {
    m_n   = 1 ;
    m_mu  = value ;
    m_mu2 = 0 ;
    m_min = value ;
    m_max = value ;
    //
    return *this ;      // RETURN
  }
  // update the regular case 
  m_n   += 1       ;                                                   // UPDATE 
  const long double delta   = value - m_mu ;             
  const long double delta_n = delta / m_n  ;
  m_mu  += delta_n ;                                                   // UPDATE
  m_mu2  = m_mu2 * ( m_n - 1 ) / m_n + delta_n * ( delta - delta_n ) ; // UPDATE
  //
  m_min  = std::min ( m_min , value ) ;                                // UPDATE       
  m_max  = std::max ( m_max , value ) ;                                // UPDATE 
  //
  return *this ;
}
// ============================================================================
// r.m.s of flag
// ============================================================================
double Ostap::StatEntity::rms    () const
{ return 2 > m_n || 0 >= m_mu2 ? 0.0L : std::sqrt ( m_mu2 * 1.0L ) ; }
// ============================================================================
// sum of values 
// ============================================================================
double Ostap::StatEntity::sum    () const 
{ return 0 < m_n ? m_mu * m_n : 0.0L ; }
// ============================================================================
// sum of value^2
// ============================================================================
double Ostap::StatEntity::sum2   () const
{ return 0 < m_n ? ( m_mu2 + 1.0L * m_mu * m_mu ) * m_n : 0.0L ; }
// ============================================================================
// error in mean value of flag
// ============================================================================
double Ostap::StatEntity::meanErr() const
{ return 0 < m_n && 0 < m_mu2 ? std::sqrt ( 1.0L * m_mu2 / m_n ) : 0 ; }
// ============================================================================
// interprete the content as efficiency
// ============================================================================
double Ostap::StatEntity::efficiency    () const
{
  //
  if ( 1 > m_n || 0 > m_mu || 1 < m_mu ) { return -1 ; }
  if ( 0 != m_min && 1 != m_min        ) { return -1 ; }
  if ( 0 != m_max && 1 != m_max        ) { return -1 ; }
  //
  return m_mu ;
}
// ============================================================================
// evaluate "the binomial" error in efficiency (note quotes) 
// ============================================================================
double Ostap::StatEntity::efficiencyErr () const
{
  if ( 0 > efficiency() ) { return -1 ; }
  //
  long double nA = 1.0L * m_n *          m_mu   ; // "accepted"
  long double nR = 1.0L * m_n * ( 1.0L - m_mu ) ; // "rejected" 
  // treat properly the bins with eff=0 and  eff=100% 
  if ( 1 > std::abs ( nA ) ) { nA = 1 ; } 
  if ( 1 > std::abs ( nR ) ) { nR = 1 ; } 
  //
  return std::sqrt ( nA * nR / m_n ) / m_n ;
}
// ============================================================================
/* increment with other counter (useful for Monitoring/Adder )
 * @code
 * const StatEntity second = ... ;
 * StatEntity first = ... ;
 * first += second ;
 * @endcode
 * @param other counter to be added
 * @return self-reference
 * @see Pebay, P., Terriberry, T.B., Kolla, H. et al. Comput Stat (2016) 31: 1305. 
 * @see https://doi.org/10.1007/s00180-015-0637-z
 */
// ============================================================================
Ostap::StatEntity& 
Ostap::StatEntity::add ( const Ostap::StatEntity& other )
{
  /// trivial updates:
  if      ( 0 == other.m_n ) { return *this ; }
  else if ( 0 ==       m_n ) 
  {
    m_n   = other.m_n   ;
    m_mu  = other.m_mu  ;
    m_mu2 = other.m_mu2 ;
    m_min = other.m_min ;
    m_max = other.m_max ;
    //
    return *this ;
  }
  //
  const unsigned long long N     = m_n + other.m_n   ;
  const long double        fA    = m_n * 1.0L / N    ;
  const long double        fB    = 1.0L - fA         ;
  const long double        delta = 1.0L * m_mu - other.m_mu ;
  //
  m_n   = N                                 ;                       // UPDATE 
  m_mu  = m_mu  * fA + other.m_mu  * fB     ;                       // UPDATE 
  m_mu2 = m_mu2 * fA + other.m_mu2 * fB + fA * fB * delta * delta ; // UPDATE
  //
  m_min = std::min ( m_min , other.m_min )  ;                       // UPDATE           
  m_max = std::max ( m_max , other.m_max )  ;                       // UPDATE 
  //
  return *this ;
}
// ============================================================================
// comparison
// ============================================================================
bool Ostap::StatEntity::operator< ( const Ostap::StatEntity& right ) const
{
  return  &right == this ? false :
    std::tie (       m_n ,       m_mu ,       m_mu2 ,       m_min ,       m_max ) < 
    std::tie ( right.m_n , right.m_mu , right.m_mu2 , right.m_min , right.m_max ) ;
}
// ============================================================================
// comparison
// ============================================================================
bool Ostap::StatEntity::operator==( const Ostap::StatEntity& right ) const
{
  return  &right == this ? true :
    std::tie (       m_n ,       m_mu ,       m_mu2 ,       m_min ,       m_max ) ==
    std::tie ( right.m_n , right.m_mu , right.m_mu2 , right.m_min , right.m_max ) ;
}
// ============================================================================
// reset all quantities
// ============================================================================
void Ostap::StatEntity::reset()
{
  m_n   = 0 ;
  m_mu  = 0 ;
  m_mu2 = 0 ;
  m_min =   std::numeric_limits<double>::max() ;
  m_max = - std::numeric_limits<double>::max() ;
}
// ============================================================================
// representation as string
// ============================================================================
std::string Ostap::StatEntity::toString () const
{
  std::ostringstream ost ;
  fillStream ( ost )  ;
  return ost.str () ;
}
// ============================================================================
// printout to std::ostream
// ============================================================================
std::ostream& Ostap::StatEntity::fillStream ( std::ostream& o ) const
{
  return
    o << Ostap::format ( "#=%-14.8g sum=%-14.8g" , n () , sum() )
      << Ostap::format ( " mean=%10.4g +- %-10.5g min/max=%10.4g/%-10.4g" ,
                         mean() , rms () , min() , max() ) ;
}
// ============================================================================

// ============================================================================
// The END
// ============================================================================
