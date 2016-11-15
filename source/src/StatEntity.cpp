// ============================================================================
#define GAUDIKERNEL_STATENTITY_CPP 1
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
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatEntity.h"
// ============================================================================
// Local 
// ============================================================================
#include "format.h"
// ============================================================================
/** @file
 *  Implementation file for class StatEntity
 *  @date 26/06/2001
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
// The full contructor from all important values
// ============================================================================
Ostap::StatEntity::StatEntity
( const unsigned long entries ,
  const double        flag    ,
  const double        flag2   ,
  const double        minFlag ,
  const double        maxFlag )
  : m_se_nEntries            ( entries )
  , m_se_accumulatedFlag     ( flag    )
  , m_se_accumulatedFlag2    ( flag2   )
  , m_se_minimalFlag         ( minFlag )
  , m_se_maximalFlag         ( maxFlag )
{}
// ============================================================================
// mean value of flag
// ============================================================================
double Ostap::StatEntity::mean   () const
{
  if ( 0 >= nEntries() ) { return 0 ;}
  const long double f1 = m_se_accumulatedFlag ;
  const long double f2 = m_se_nEntries ;
  return f1 / f2 ;
}
// ============================================================================
// r.m.s of flag
// ============================================================================
double Ostap::StatEntity::rms    () const
{
  if ( 0 >= nEntries() ) { return 0 ; }
  const long double f1  = m_se_accumulatedFlag  ;
  const long double f2  = f1 / nEntries () ;
  const long double f3  = m_se_accumulatedFlag2 ;
  const long double f4  = f3 / nEntries () ;
  const long double result = f4 - f2 * f2  ;
  return ( 0 > result ) ? 0 : std::sqrt ( result ) ;
}
// ============================================================================
// error in mean value of flag
// ============================================================================
double Ostap::StatEntity::meanErr() const
{
  if ( 0 >= nEntries () ) { return 0 ; }
  const long double f1  = m_se_accumulatedFlag  ;
  const long double f2  = f1 / nEntries () ;
  const long double f3  = m_se_accumulatedFlag2 ;
  const long double f4  = f3 / nEntries () ;
  const long double result = f4 - f2 * f2  ;
  if ( 0 > result      ) { return 0 ; }
  //
  return std::sqrt ( result / nEntries () ) ;
}
// ============================================================================
// interprete the content as efficiency
// ============================================================================
double Ostap::StatEntity::efficiency    () const
{
  if ( 1 > nEntries () || 0 > sum() || sum() > nEntries() ) { return -1 ; }
  const long double fMin = min () ;
  if ( 0 != fMin && 1 != fMin ) { return -1 ; }
  const long double fMax = max () ;
  if ( 0 != fMax && 1 != fMax ) { return -1 ; }
  return mean() ;
}
// ============================================================================
// evaluate the binomial error in efficiency
// ============================================================================
double Ostap::StatEntity::efficiencyErr () const
{
  if ( 0 > efficiency() ) { return -1 ; }
  //
  long double n1 = sum () ;
  // treat properly the bins with eff=0
  if ( 0 == n1 ) { n1 = 1 ; } ///< @attention treat properly the bins with eff=0
  const long double n3 = nEntries  () ;
  long double       n2 = n3 - sum () ;
  // treat properly the bins with eff=100%
  if ( 1 > std::abs ( n2 ) ) { n2 = 1 ; } ///< @attention treat properly the bins with eff=100%
  //
  return std::sqrt ( n1 * n2 / n3 ) / n3 ;
}
// ============================================================================
// increment with other entity
// ============================================================================
Ostap::StatEntity& 
Ostap::StatEntity::operator+=( const Ostap::StatEntity& other )
{
  m_se_nEntries         += other.m_se_nEntries ;
  m_se_accumulatedFlag  += other.m_se_accumulatedFlag  ;
  m_se_accumulatedFlag2 += other.m_se_accumulatedFlag2 ;
  m_se_minimalFlag = std::min ( m_se_minimalFlag , other.m_se_minimalFlag ) ;
  m_se_maximalFlag = std::max ( m_se_maximalFlag , other.m_se_maximalFlag ) ;
  //
  return *this ;
}
// ============================================================================
// comparison
// ============================================================================
bool Ostap::StatEntity::operator< ( const Ostap::StatEntity& se ) const
{
  if      ( &se == this                     ) { return false ; }
  else if ( nEntries () <   se.nEntries ()  ) { return true  ; }
  else if ( nEntries () ==  se.nEntries () &&
            sum      () <   se.sum      ()  ) { return true  ; }
  else if ( nEntries () ==  se.nEntries () &&
            sum      () ==  se.sum      () &&
            min      () <   se.min      ()  ) { return true  ; }
  else if ( nEntries () ==  se.nEntries () &&
            sum      () ==  se.sum      () &&
            min      () ==  se.min      () &&
            max      () <   se.max      ()  ) { return true  ; }
  else if ( nEntries () ==  se.nEntries () &&
            sum      () ==  se.sum      () &&
            min      () ==  se.min      () &&
            max      () ==  se.max      () &&
            sum2     () <   se.sum2     ()  ) { return true  ; }
  ///
  return false;
}
// ============================================================================
// comparison
// ============================================================================
bool Ostap::StatEntity::operator==( const Ostap::StatEntity& se ) const
{
  return ( &se == this ) || 
    ( m_se_nEntries         == se.m_se_nEntries         && 
      m_se_accumulatedFlag  == se.m_se_accumulatedFlag  && 
      m_se_accumulatedFlag2 == se.m_se_accumulatedFlag2 && 
      m_se_minimalFlag      == se.m_se_minimalFlag      && 
      m_se_maximalFlag      == se.m_se_maximalFlag       ) ;
}
// ============================================================================
// increment a flag
// ============================================================================
unsigned long long Ostap::StatEntity::add ( const double Flag )
{
  /// accumulate the flag
  m_se_accumulatedFlag   += Flag        ; // accumulate the flag
  /// evaluate min/max
  m_se_minimalFlag = std::min ( m_se_minimalFlag , Flag ) ; // evaluate min/max
  m_se_maximalFlag = std::max ( m_se_maximalFlag , Flag ) ; // evaluate min/max
  // accumulate statistics, but avoid FPE for very small flags...
  static const double s_min1 = 2 * std::sqrt ( std::numeric_limits<double>::min() ) ;
  if ( s_min1 < Flag || -s_min1 > Flag )
  { m_se_accumulatedFlag2  += Flag * Flag ; }// accumulate statistics:
  //
  return ++m_se_nEntries ;
}
// ============================================================================
// reset all quantities
// ============================================================================
void Ostap::StatEntity::reset()
{
  //
  m_se_nEntries            =        0 ;
  m_se_accumulatedFlag     =        0 ;
  m_se_minimalFlag         =       std::numeric_limits<double>::max() ;
  m_se_maximalFlag         =  -1 * std::numeric_limits<double>::max() ;
  m_se_accumulatedFlag2    =        0 ;
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
    o << Ostap::format ( "#=%-14.8g Sum=%-14.8g" , nEntries() , sum() )
      << Ostap::format ( " Mean=%10.4g +- %-10.5g Min/Max=%10.4g/%-10.4g" ,
                         mean() , rms() , min() , max() ) ;
}
// ============================================================================

// ============================================================================
// The END
// ============================================================================
