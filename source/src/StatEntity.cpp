// ============================================================================
#define OSTAP_STATENTITY_CPP 1
// ============================================================================
// include files
// ============================================================================
// STD & STL
// ============================================================================
#include <sstream>
#include <string>
#include <cmath>
#include <limits>
#include <tuple>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/StatEntity.h"
// ============================================================================
// Local 
// ============================================================================
#include "format.h"
#include "Exception.h"
#include "status_codes.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::StatEntity
 *  @date 26/06/2001
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  const Ostap::Math::Zero    <double> s_zero  {} ;
  // ==========================================================================
}
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
  Ostap::Assert ( std::isfinite ( m_mu  ) &&
                  std::isfinite ( m_mu2 ) &&
                  std::isfinite ( m_min ) &&
                  std::isfinite ( m_max )            ,                   
                  "Invalid parameters!"              , 
                  "Ostap::Math::StatEntity"          ,
                  INVALID_PARS , __FILE__ , __LINE__ ) ;
  // empty counter  (ignore min/max)
  if ( empty () )
    {
      Ostap::Assert ( s_zero ( m_mu ) && s_zero ( m_mu2 )  , 
                      "Invalid mu/mu2 for empty counter!"  ,
                      "Ostap::StatEntity"                  ,
                      INVALID_PARS  , __FILE__ , __LINE__  ) ;
      m_mu  = 0 ;
      m_mu2 = 0 ;
      /// redefine/ignore min/max 
      m_min =   std::numeric_limits<double>::max() ;
      m_max = - std::numeric_limits<double>::max() ;
    }
  else
    {
      Ostap::Assert ( m_min <= mu && mu <= m_max         , 
                      "Invalid minv/mu/maxv"             ,
                      "Ostap::StatEntity"                ,
                      INVALID_PARS , __FILE__ , __LINE__ ) ; 
    }
  //
  if ( s_zero ( m_mu2 ) ) { m_mu2 = 0 ; }
  //
  Ostap::Assert ( ( !empty () ) || ( empty() && !m_mu2 ) ,		  
                  "Inconsistent mu2/empty!"              ,
                  "Ostap::StatEntity"                    ,
                  INVALID_PARS , __FILE__ , __LINE__     ) ;  
  //
  Ostap::Assert ( 0 <= m_mu2                         , 
                  "Invalid second moment"            ,
                  "Ostap::StatEntity"                , 
                  INVALID_PARS , __FILE__ , __LINE__ ) ;
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
  /// ignore non-finite values 
  if ( !std::isfinite ( value ) ) { return *this ; }
  //
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
  //
  // update the regular case 
  //
  const unsigned long long N     = m_n + 1             ;  
  const long double        fA    = m_n * 1.0L / N      ;
  const long double        fB    = 1.0L - fA           ;
  const long double        delta = 1.0L * value - m_mu ;
  //
  m_n  += 1                                     ; // UPDATE 
  m_mu  = fA * m_mu  + fB * value               ; // UPDATE 
  m_mu2  = fA * m_mu2 + fA * fB * delta * delta ; // UPDATE
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
{ return 0 < m_n ? 1.0L * m_mu * m_n : 0.0L ; }
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
  const unsigned long long N     = m_n + other.m_n          ;
  const long double        fA    = m_n * 1.0L / N           ;
  const long double        fB    = 1.0L - fA                ;
  const long double        delta = 1.0L * other.m_mu - m_mu ;
  //
  m_n   = m_n + other.m_n                   ;                       // UPDATE 
  m_mu  = fA * m_mu  + fB * other.m_mu      ;                       // UPDATE 
  m_mu2 = fA * m_mu2 + fB * other.m_mu2 + fA * fB * delta * delta ; // UPDATE
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
// swap two counters
// ============================================================================
void Ostap::StatEntity::swap ( Ostap::StatEntity& right )
{
  std::swap ( m_n   , right.m_n   ) ;
  std::swap ( m_mu  , right.m_mu  ) ;
  std::swap ( m_mu2 , right.m_mu2 ) ;
  std::swap ( m_min , right.m_min ) ;
  std::swap ( m_max , right.m_max ) ;
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
    o << Ostap::format ( "#=%-10.5g sum=%+-10.5g" , n () , sum() )
      << Ostap::format ( " mean/rms=%+10.5g/%-10.5g min/max=%+10.5g/%+-10.5g" ,
                         mean () , rms () , min () , max () ) ;
}
// ============================================================================

// ============================================================================
//                                                                      The END
// ============================================================================
