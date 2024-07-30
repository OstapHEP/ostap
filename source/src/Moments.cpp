// ============================================================================
// Include files 
// ============================================================================
//  STD&STL
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Moments.h"
// ============================================================================
// Local
// ============================================================================
#include  "Exception.h"
#include  "local_math.h"
// ============================================================================
/** @file 
 *  Implementation file for classes from the file Ostap/Moments.h
 *  @date 2020-06-07 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
/*  @var s_INVALID_CENTRAL_MOMENT ; 
 *  the invalid value of the central moment 
 */ 
// ============================================================================
const double Ostap::Math::Moments::s_INVALID_MOMENT = -1 * s_INFINITY ;
// ===========================================================================
// virtual destructor
// ===========================================================================
Ostap::Math::Statistic::~Statistic(){}
// ===========================================================================
// virtual destructor
// ===========================================================================
Ostap::Math::WStatistic::~WStatistic(){}
// ===========================================================================
// virtual destructor
// ===========================================================================
Ostap::Math::Moment::~Moment(){}
// ===========================================================================
// virtual destructor
// ===========================================================================
Ostap::Math::WMoment::~WMoment(){}
// ===========================================================================


// ===========================================================================
// HarmoincMean:: accumulate only non-zero entries 
// ===========================================================================
Ostap::Math::HarmonicMean&
Ostap::Math::HarmonicMean::add ( const double x )
{
  if ( !s_zero ( x ) ) { m_inv.add ( 1/x ) ; }
  return *this ;
}
// ===========================================================================
// WHarmoincMean:: accumulate only non-zero entries 
// ===========================================================================
Ostap::Math::WHarmonicMean&
Ostap::Math::WHarmonicMean::add
( const double x ,
  const double w ) 
{
  if ( !s_zero ( x ) ) { m_inv.add ( 1/x , w ) ; }
  return *this ;
}


// ===========================================================================
Ostap::Math::PowerMean::PowerMean
( const double p )
  : m_p   ( p )
  , m_pow (   )
{}
// ===========================================================================
// add two counters togather if p is common 
// ===========================================================================
Ostap::Math::PowerMean&
Ostap::Math::PowerMean::add ( const Ostap::Math::PowerMean& x )
{
  Ostap::Assert ( s_equal ( m_p , x.m_p ) ,
		  "Cannot add counters with non-eual values of 'p'" , 
		  "Ostap::Math::PowerMean"  ) ;
  m_pow.add ( x.m_pow ) ;
  return *this ;
}

// ===========================================================================
Ostap::Math::WPowerMean::WPowerMean
( const double p )
  : m_p   ( p )
  , m_pow (   )
{}
// ===========================================================================
// add two counters togather if p is common 
// ===========================================================================
Ostap::Math::WPowerMean&
Ostap::Math::WPowerMean::add ( const Ostap::Math::WPowerMean& x )
{
  Ostap::Assert ( s_equal ( m_p , x.m_p ) ,
		  "Cannot add counters with non-eual values of 'p'" , 
		  "Ostap::Math::WPowerMean"  ) ;
  m_pow.add ( x.m_pow ) ;
  return *this ;
}
// ===========================================================================
Ostap::Math::WPowerMean&
Ostap::Math::WPowerMean::add
( const double x  ,
  const double w  )
{
  if ( 0 < x && !s_zero ( w ) ) { m_pow.add ( std::pow ( x , m_p ) , w ) ; }
  return *this ;
}

// ===========================================================================
Ostap::Math::LehmerMean::LehmerMean
( const double p )
  : m_p ( p )
  , m_lp   ()
  , m_lpm1 ()
{}
// ===========================================================================
Ostap::Math::WLehmerMean::WLehmerMean
( const double p )
  : m_p    ( p )
  , m_lp   ()
  , m_lpm1 ()
{}
// ===========================================================================
// add two counters togather if p is common 
// ===========================================================================
Ostap::Math::LehmerMean&
Ostap::Math::LehmerMean::add ( const Ostap::Math::LehmerMean& x )
{
  Ostap::Assert ( s_equal ( m_p , x.m_p ) ,
		  "Cannot add counters with non-eual values of 'p'" , 
		  "Ostap::Math::LehmerMean"  ) ;
  m_lp   .add ( x.m_lp   ) ;
  m_lpm1 .add ( x.m_lpm1 ) ;
  return *this ;
}
// ===========================================================================
Ostap::Math::WLehmerMean&
Ostap::Math::WLehmerMean::add
( const Ostap::Math::WLehmerMean& x )
{
  Ostap::Assert ( s_equal ( m_p , x.m_p ) ,
		  "Cannot add counters with non-eual values of 'p'" , 
		  "Ostap::Math::WLehmerMean"  ) ;
  m_lp   .add ( x.m_lp   ) ;
  m_lpm1 .add ( x.m_lpm1 ) ;
  return *this ;
}
// ===========================================================================
// accumulate only positive entries with non-zero weight 
// ===========================================================================
Ostap::Math::WLehmerMean&
Ostap::Math::WLehmerMean::add
( const double x ,
  const double w ) 
{
  if ( ( 0 < x ) && !s_zero ( x ) && !s_zero ( w ) )
    {
      m_lp  .add ( std::pow ( x , m_p     ) , w ) ;
      m_lpm1.add ( std::pow ( x , m_p - 1 ) , w ) ;	    
    }
  return *this ;
}


// ===========================================================================
// default constructor
// ===========================================================================
Ostap::Math::MinValue::MinValue ()
  : m_min ( std::numeric_limits<double>::max() ) 
  , m_cnt ()
{}
// ===========================================================================
// default constructor
// ===========================================================================
Ostap::Math::MaxValue::MaxValue ()
  : m_max ( - std::numeric_limits<double>::max() ) 
  , m_cnt ()
{}
// ===========================================================================
// default constructor
// ===========================================================================
Ostap::Math::MinMaxValue::MinMaxValue ()
  : m_min (   std::numeric_limits<double>::max() ) 
  , m_max ( - std::numeric_limits<double>::max() ) 
  , m_cnt ()
{}
// ===========================================================================
// default constructor
// ===========================================================================
Ostap::Math::WMinMaxValue::WMinMaxValue ()
  : m_min (   std::numeric_limits<double>::max() ) 
  , m_max ( - std::numeric_limits<double>::max() ) 
  , m_cnt ()
{}
// ======================================================================


// int test()
// {
//   Ostap::Math::Moment_<5> m1 ;
//   for ( unsigned i = 0 ; i < 100 ; ++i ) { m1 += i ; }
//   //
//   m1.M(0) ;
//   m1.M(1) ;
//   m1.M(2) ;
//   Ostap::Math::Moment_<5> m2 ;
//   for ( unsigned i = 0 ; i < 100 ; ++i ) { m2 += i ; }
//   //
//   m2.M(0) ;
//   m2.M(1) ;
//   m2.M(2) ;
//   Ostap::Math::Moment_<2> m4 ;
//   for ( unsigned i = 0 ; i < 100 ; ++i ) { m4 += i ; }
//   const auto m3  = m1 + m2 ;
//   const auto a1 = Ostap::Math::Moments::central_moment<2> ( m3 ) ;
//   const auto a2 = Ostap::Math::Moments::central_moment<2> ( m4 ) ;
//   return 0 ;
// }


// ============================================================================
//                                                                      The END 
// ============================================================================
