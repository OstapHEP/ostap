// ============================================================================
// Include files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <cmath>
#include <sstream>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/WStatEntity.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// local
// ============================================================================
#include "format.h"
#include "status_codes.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::WStatEntity
 *  @see StatEntity
 *  @author  Vanya Belyaev Ivan.Belyaev@itep.ru
 *  @date 2014-04-07 
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  const Ostap::Math::Zero    <double> s_zero  {} ;
  // ==========================================================================
}
// ============================================================================
// constructor from StatEntity of values 
// ============================================================================
Ostap::WStatEntity::WStatEntity ( const Ostap::StatEntity& values ) 
  : m_mu      ( values.mu      () ) 
  , m_mu2     ( values.mu2     () ) 
  , m_values  ( values            ) 
  , m_weights ( values.n () , 1 , 0 , 1 , 1 )  // weights are trivial!
{}
// ============================================================================
// full constructor
// ============================================================================
Ostap::WStatEntity::WStatEntity
( const double             mu      ,
  const double             mu2     ,
  const Ostap::StatEntity& values  ,
  const Ostap::StatEntity& weights )
  : m_mu      ( mu      )
  , m_mu2     ( mu2     )
  , m_values  ( values  ) 
  , m_weights ( weights ) 
{
  Ostap::Assert ( std::isfinite ( m_mu   ) &&
                  std::isfinite ( m_mu2  )           , 
                  "Invalid parameters!"              , 
                  "Ostap::Math::WStatEntity"         ,
                  INVALID_PARS , __FILE__ , __LINE__ ) ;
  //
  if ( empty () )
    {
      Ostap::Assert ( s_zero ( m_mu ) && s_zero  ( m_mu2 ) ,		      
                      "Ostap::WStatWEntity: invalid mu/mu2 for empty counter!" ,
                      "Ostap::WStatEntity"                 ,
                      INVALID_PARS , __FILE__ , __LINE__   ) ;
      m_mu  = 0 ;
      m_mu2 = 0 ;
    }
  //
  Ostap::Assert ( m_values.n () == m_weights.n () , 
                  "Ostap::WStatWEntity: inconsistent values/weights counters!" ,
                  "Ostap::WStatEntity" ,
                  INVALID_PARS , __FILE__ , __LINE__ ) ;
  //
  if ( s_zero ( m_mu2 ) ) { m_mu2 = 0 ; }
  //
}
// ============================================================================
// finite values? 
// ============================================================================
bool Ostap::WStatEntity::isfinite () const
{
  return
    std::isfinite ( m_mu  ) &&
    std::isfinite ( m_mu2 ) &&
    m_values .isfinite ()   &&
    m_weights.isfinite ()   ;
}
// ============================================================================
// Is it OK ?
// ============================================================================
bool Ostap::WStatEntity::ok () const
{
  if  ( empty() ) { return !m_mu && !m_mu2 && m_values.empty() && m_weights.empty() ; }
  return isfinite () && m_values.ok() && m_weights.ok () ;
}   
// ============================================================================
// update statistics 
// ============================================================================
Ostap::WStatEntity&
Ostap::WStatEntity::add   
( const double value  ,  
  const double weight )
{
  //
  /// ignore non-finite values and zero weights 
  if ( !weight || !std::isfinite ( value  ) || !std::isfinite ( weight ) ) { return *this ; }
  //
  if ( 0 == n() ) 
    {
      m_mu       = value  ;
      m_mu2      = 0      ;
      m_values  += value  ; 
      m_weights += weight ;    
      //
      return *this ;
    }
  //
  const long double wA    = n () * m_weights.mu() ;
  const long double wB    = weight                ; 
  //
  const long double W     = wA + wB               ;
  //
  // if      ( !W ) { /* soem action is needed ? */ } 
  if ( !wA )
    {
      m_mu  = value ;
      m_mu2 = 0     ;
    }
  else 
    {
      const long double fA    = wA / W                ;
      const long double fB    = 1.0L - fA             ;
      const long double delta = 1.0L * value - m_mu   ;
      //
      m_mu  = fA * m_mu  + fB * value              ; // UPDATE 
      m_mu2 = fA * m_mu2 + fA * fB * delta * delta ; // UPDATE
    }
  //
  m_values  += value  ;                          // UPDATE 
  m_weights += weight ;                          // UPDATE 
  //
  return *this ;
}
// ============================================================================
/* add another counter 
 * @see Pebay, P., Terriberry, T.B., Kolla, H. et al. Comput Stat (2016) 31: 1305. 
 * @see https://doi.org/10.1007/s00180-015-0637-z
 */
// ============================================================================
Ostap::WStatEntity& 
Ostap::WStatEntity::add
( const Ostap::WStatEntity& other ) 
{
  // treat the trivial cases  
  if      ( 0 == other.n () ) {                  return *this ; }
  else if ( 0 ==       n () ) { *this == other ; return *this ; }
  //
  const long double wA    = 1.0L *       n () *       m_weights.mu () ; //       sumw () 
  const long double wB    = 1.0L * other.n () * other.m_weights.mu () ; // other.sumw ()     ;
  //
  const long double W     = wA + wB                  ;
  const long double fA    = wA / W                   ;
  const long double fB    = 1.0L - fA                ; // wB / W            ;
  const long double delta = 1.0L * other.m_mu - m_mu ;
  //
  m_mu  = fA * m_mu  + fB * other.m_mu   ;                           // UPDATE 
  m_mu2 = fA * m_mu2 + fB * other.m_mu2  + fA * fB * delta * delta ; // UPDATE 
  //
  m_values  += other.m_values  ; // UPDATE 
  m_weights += other.m_weights ; // UPDATE 
  //
  return *this ;
}
// ============================================================================
// sum_i weight_i*value_i
// ============================================================================
double Ostap::WStatEntity::sum   () const  // sum_i weight_i * value_i
{ return 0 < n () ? m_mu * sumw () : 0.0 ; }
// ============================================================================
// sum_i weight_i*(value_i**2)
// ============================================================================
double Ostap::WStatEntity::sum2  () const  // sum_i weight_i * (value_i**2)
{ return 1 < n () ? ( dispersion () + std::pow ( mean () , 2 ) ) * sumw () : 0.0 ; }
// ============================================================================
// get the sample mean
// ============================================================================
double Ostap::WStatEntity::meanErr () const
{
  const double neff = nEff() ;
  if ( s_zero ( neff ) ) { return 0 ; }
  const double v = dispersion() / neff ;
  return v <= 0 ? 0.0 : std::sqrt ( v ) ;
}
// ============================================================================
// calculate rms 
// ============================================================================
double Ostap::WStatEntity::rms () const
{
  if ( 2 > n() ) { return 0 ; }
  const long double d = dispersion () ;
  if ( 0 >= d || s_zero ( d ) ) { return 0 ; }
  return std::sqrt ( d ) ;
}
// ============================================================================
// calculate effective number of entries 
// ============================================================================
double Ostap::WStatEntity::nEff () const
{
  if ( 0 == n ()      ) { return 0 ; }
  const double sw2 =  sumw2() ;
  return s_zero ( sw2 ) ? 0.0 : std::pow ( sumw () , 2 ) / sw2  ;
}
// ============================================================================
// reset statistic 
// ============================================================================
void Ostap::WStatEntity::reset () 
{
  m_mu  = 0 ;
  m_mu2 = 0 ;
  m_values .reset() ;
  m_weights.reset() ;
}
// ============================================================================
// printout 
// ============================================================================
std::ostream& Ostap::WStatEntity::fillStream ( std::ostream& o ) const 
{
  return 
    o << Ostap::format ( "#=%+-10.5g sum=%+-10.5g " , nEff() , sum () ) 
      << Ostap::format ( " mean/rms=%+10.5g/%-10.5g min/max=%+10.5g/%+-10.5g" ,
                         mean () , rms () , m_values.min() , m_values.max () ) ;
}
// ============================================================================
// convert to string 
// ============================================================================
std::string Ostap::WStatEntity::toString () const
{
  std::ostringstream ost ;
  fillStream ( ost )  ;
  return ost.str () ;
}
// ============================================================================
// comparison
// ============================================================================
bool Ostap::WStatEntity::operator< ( const Ostap::WStatEntity& right ) const
{
  return  &right == this ? false :
    std::tie (       m_mu ,       m_mu2 ,       m_values,       m_weights ) < 
    std::tie ( right.m_mu , right.m_mu2 , right.m_values, right.m_weights ) ;
}
// ============================================================================
// comparison
// ============================================================================
bool Ostap::WStatEntity::operator==( const Ostap::WStatEntity& right ) const
{
  return  &right == this ? true :
    std::tie (       m_mu ,       m_mu2 ,       m_values,       m_weights ) == 
    std::tie ( right.m_mu , right.m_mu2 , right.m_values, right.m_weights ) ;
}
// ============================================================================
// swap two counters
// ============================================================================
void Ostap::WStatEntity::swap ( Ostap::WStatEntity& right )
{
  std::swap   ( m_mu      , right.m_mu      ) ;
  std::swap   ( m_mu2     , right.m_mu2     ) ;
  Ostap::swap ( m_values  , right.m_values  ) ;
  Ostap::swap ( m_weights , right.m_weights ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
