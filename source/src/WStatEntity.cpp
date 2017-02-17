// $ID:$ 
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
// ============================================================================
// local
// ============================================================================
#include "format.h"
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
// empty constructor 
// ============================================================================
Ostap::WStatEntity::WStatEntity()
  : m_sum     ( 0 ) 
  , m_sum2    ( 0 )
  , m_values  (   )
  , m_weights (   )
{}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::WStatEntity::WStatEntity ( const Ostap::WStatEntity& right ) 
  : m_sum     ( right.m_sum     ) 
  , m_sum2    ( right.m_sum2    ) 
  , m_values  ( right.m_values  ) 
  , m_weights ( right.m_weights )  
{}
// ============================================================================
// constructor from StatEntity of values 
// ============================================================================
Ostap::WStatEntity::WStatEntity ( const Ostap::StatEntity& values ) 
  : m_sum     ( values.sum      () ) 
  , m_sum2    ( values.sum2     () ) 
  , m_values  ( values             ) 
  , m_weights ( values.nEntries () , 
                values.nEntries () , 
                values.nEntries () , 1 , 1 ) 
{}
// ============================================================================
// update statistics 
// ============================================================================
Ostap::WStatEntity&
Ostap::WStatEntity::add   
( const double value  ,  
  const double weight )
{
  m_sum     += weight * value         ;
  m_sum2    += weight * value * value ;  
  //
  if ( !s_zero ( weight ) ) { m_values += value ; }
  //
  m_weights += weight ;
  //
  return *this ;
}
// ============================================================================
// get the sample mean
// ============================================================================
double Ostap::WStatEntity::mean () const
{
  return ( 0 == nEntries () 
           || s_zero ( m_sum            ) 
           || s_zero ( m_weights.sum () ) ) ? 0.0 
    : m_sum / m_weights.sum () ;
}
// ============================================================================
// get the sample mean
// ============================================================================
double Ostap::WStatEntity::meanErr () const
{
  const double neff = nEff() ;
  if ( s_zero ( neff ) ) { return 0 ; }
  //
  const double v = dispersion() / neff ;
  //
  return v <= 0 ? 0.0 : std::sqrt ( v ) ;
}
// ============================================================================
// calculate dispersion 
// ============================================================================
double Ostap::WStatEntity::dispersion () const
{ 
  //
  if ( 1 >= nEntries() || s_zero ( m_weights.sum() ) ) { return 0 ; }
  //
  const long double a1 = m_sum2   ;
  const long double a2 = mean  () ;
  //
  return a1 / m_weights.sum () - std::pow ( a2 , 2 ) ;
}
// ============================================================================
// calculate rms 
// ============================================================================
double Ostap::WStatEntity::rms () const
{
  const long double d = dispersion () ;
  //
  if ( 0 >= d || s_zero ( d ) ) { return 0 ; }
  return std::sqrt ( d ) ;
}
// ============================================================================
// calculate effective number of entries 
// ============================================================================
double Ostap::WStatEntity::nEff () const
{
  //
  if ( 0 == nEntries() || s_zero ( m_weights.sum2 () ) ) { return 0 ; }
  //
  return std::pow ( m_weights.sum() , 2 ) /  m_weights.sum2 () ;
}
// ============================================================================
// reset statistic 
// ============================================================================
void Ostap::WStatEntity::reset () 
{
  m_sum  = 0 ;
  m_sum2 = 0 ;
  m_values .reset() ;
  m_weights.reset() ;
}
// ============================================================================
// printout 
// ============================================================================
std::ostream& Ostap::WStatEntity::fillStream ( std::ostream& o ) const 
{
  return 
    o << Ostap::format ( "#=%-14.8g Sum=%-14.8g " , nEff() , m_sum ) 
      << Ostap::format ( " Mean=%10.4g +- %-10.5g Min/Max=%10.4g/%-10.4g" ,
                         mean() ,  rms() ,  m_values.min() ,  m_values.max() ) ;
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
// The END 
// ============================================================================
