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
// constructor 
// ===========================================================================
/* full constructor
 *  @param size number of entries 
 *  @param sumw sum of weights 
 *  @param sumw sum of squared weights 
 */
// ===========================================================================
Ostap::Math::WMoment_<0>::WMoment_ 
( const unsigned long long size  ,
  const double             sumw  ,
  const double             sumw2 )
  : m_size ( size  )
  , m_w    ( sumw  )
  , m_w2   ( sumw2 ) 
{
  if ( 0 == m_size )
    {
      Ostap::Assert ( s_zero ( m_w ) && s_zero ( m_w2 ) ,  
		      "Non-zero sum of (squared) weights for empty counter!" , 
		      "Ostap::Math::WMoment_<0>"  ) ;
      m_w  = 0 ;
      m_w2 = 0 ;       
    }
  //
  if ( s_zero ( m_w2 ) )
    {
      Assert ( s_zero ( m_w )  ,
	       "Non zero sum of weigth for zero sumw2 !" , 
	       "Ostap::Math::WMoment_<0>"  ) ;      
      m_w  = 0 ;
      m_w2 = 0 ;
    }
  //
  if ( s_zero ( m_w ) ) { m_w = 0 ; }
  //
  Ostap::Assert ( 0 <= m_w2  ,      
		  "Negatoive sum of squared weights!" , 
		  "Ostap::Math::WMoment_<0>"  ) ;
  
}
// ===========================================================================



// ===========================================================================
Ostap::Math::GeometricMean::GeometricMean
( const Ostap::Math::GeometricMean::Counter& cnt )
  : m_log ( cnt )
{}
// ===========================================================================
Ostap::Math::HarmonicMean::HarmonicMean
( const Ostap::Math::HarmonicMean::Counter& cnt )
  : m_inv ( cnt )
{}
// ===========================================================================
Ostap::Math::ArithmeticMean::ArithmeticMean
( const Ostap::Math::ArithmeticMean::Counter& cnt )
  : Counter ( cnt ) 
{}
// ===========================================================================
Ostap::Math::PowerMean::PowerMean
( const double p )
  : m_p   ( p )
  , m_pow (   )
{}
// ===========================================================================
Ostap::Math::PowerMean::PowerMean
( const double                           p   ,
  const Ostap::Math::PowerMean::Counter& cnt ) 
  : m_p   ( p   )
  , m_pow ( cnt )
{}
// ===========================================================================
Ostap::Math::LehmerMean::LehmerMean
( const double p )
  : m_p    ( p )
  , m_lp   (   )
  , m_lpm1 (   )
{}
// ===========================================================================
Ostap::Math::LehmerMean::LehmerMean
( const double                           p     ,
  const Ostap::Math::LehmerMean::Counter& cnt1 , 
  const Ostap::Math::LehmerMean::Counter& cnt2 ) 
  : m_p    ( p    )
  , m_lp   ( cnt1 )
  , m_lpm1 ( cnt2 )
{
  Ostap::Assert  ( m_lp.size () == m_lpm1.size () ,
		   "Inconsistent structure of two conuters!" ,
		   "Ostap::Math::LehmerMean" ) ;
  
}
// ===========================================================================
Ostap::Math::WGeometricMean::WGeometricMean
( const Ostap::Math::WGeometricMean::Counter& cnt )
  : m_log ( cnt )
{}
// ===========================================================================
Ostap::Math::WHarmonicMean::WHarmonicMean
( const Ostap::Math::WHarmonicMean::Counter& cnt )
  : m_inv ( cnt )
{}
// ===========================================================================
Ostap::Math::WArithmeticMean::WArithmeticMean
( const Ostap::Math::WArithmeticMean::Counter& cnt )
  : Counter ( cnt ) 
{}
// ===========================================================================
Ostap::Math::WPowerMean::WPowerMean
( const double p )
  : m_p   ( p )
  , m_pow (   )
{}
// ===========================================================================
Ostap::Math::WPowerMean::WPowerMean
( const double                            p   ,
  const Ostap::Math::WPowerMean::Counter& cnt ) 
  : m_p   ( p   )
  , m_pow ( cnt )
{}
// ===========================================================================
Ostap::Math::WLehmerMean::WLehmerMean
( const double p )
  : m_p    ( p )
  , m_lp   (   )
  , m_lpm1 (   )
{}
// ===========================================================================
Ostap::Math::WLehmerMean::WLehmerMean
( const double                             p    ,
  const Ostap::Math::WLehmerMean::Counter& cnt1 , 
  const Ostap::Math::WLehmerMean::Counter& cnt2 ) 
  : m_p    ( p    )
  , m_lp   ( cnt1 )
  , m_lpm1 ( cnt2 )
{
  Ostap::Assert  ( m_lp.size () == m_lpm1.size () ,
		   "Inconsistent structure of two conuters!" ,
		   "Ostap::Math::WLehmerMean" ) ;
  Ostap::Assert  ( s_equal ( m_lp.w  () , m_lpm1.w  () ) , 
		   "Inconsistent structure of two conuters!" ,
		   "Ostap::Math::WLehmerMean" ) ;
  Ostap::Assert  ( s_equal ( m_lp.w2 () , m_lpm1.w2 () ) , 
		   "Inconsistent structure of two conuters!" ,
		   "Ostap::Math::WLehmerMean" ) ;
}


// ===========================================================================
// GeometricMean::accumulate only positive entries 
// ===========================================================================
Ostap::Math::GeometricMean&
Ostap::Math::GeometricMean::add
( const double         x )
{
  if ( 0 < x && !s_zero ( x ) ) { m_log.add ( std::log2 ( x ) ) ; }
  return *this ;
}
// ===========================================================================
// HarmonicMean:: accumulate only non-zero entries 
// ===========================================================================
Ostap::Math::HarmonicMean&
Ostap::Math::HarmonicMean::add ( const double x )
{
  if ( !s_zero ( x ) ) { m_inv.add ( 1/x ) ; }
  return *this ;
}
// ===========================================================================
// PowerMean::accumulate only positive entries 
// ===========================================================================
Ostap::Math::PowerMean&
Ostap::Math::PowerMean::add
( const double         x )
{
  if ( ( 0 < x ) && !s_zero ( x ) ) { m_pow.add ( std::pow ( x , m_p ) ) ; }
  return *this ;
}
// ===========================================================================
// LehmerMean: accumulate only positive entries
// ===========================================================================
Ostap::Math::LehmerMean&
Ostap::Math::LehmerMean::add ( const double x )
{
  if ( ( 0 < x ) && !s_zero( x ) )
    {
      m_lp  .add ( std::pow ( x , m_p     ) ) ;
      m_lpm1.add ( std::pow ( x , m_p - 1 ) ) ;	    
    }
  return *this ;
}
// ===========================================================================
// WGeometricMean:: accumulate only positive entries 
// ===========================================================================
Ostap::Math::WGeometricMean&
Ostap::Math::WGeometricMean::add
( const double x ,
  const double w ) 
{
  if ( ( 0 < x ) && !s_zero ( x ) ) { m_log.add ( std::log2 ( x ) , w ) ; }
  return *this ;
}
// ===========================================================================
// WHarmonicMean:: accumulate only non-zero entries 
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
// WPowerMean:: accumulate only positive entries 
// ===========================================================================
Ostap::Math::WPowerMean&
Ostap::Math::WPowerMean::add
( const double x  ,
  const double w  )
{
  if ( ( 0 < x ) && !s_zero ( x ) ) { m_pow.add ( std::pow ( x , m_p ) , w ) ; }
  return *this ;
}
// ===========================================================================
// LehmerMean: accumulate only positive entries
// ===========================================================================
Ostap::Math::WLehmerMean&
Ostap::Math::WLehmerMean::add
( const double x ,
  const double w ) 
{
  if ( ( 0 < x ) && !s_zero ( x ) )
    {
      m_lp  .add ( std::pow ( x , m_p     ) , w ) ;
      m_lpm1.add ( std::pow ( x , m_p - 1 ) , w ) ;	    
    }
  return *this ;
}
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

// ===========================================================================
// default constructor
// ===========================================================================
Ostap::Math::MinMaxValue::MinMaxValue ()
  : m_min (   std::numeric_limits<double>::max() ) 
  , m_max ( - std::numeric_limits<double>::max() ) 
  , m_cnt ()
{}
// ===========================================================================
// constructor
// ===========================================================================
Ostap::Math::MinMaxValue::MinMaxValue
( const double                             minv ,
  const double                             maxv ,
  const Ostap::Math::MinMaxValue::Counter& cnt  )
  : m_min ( minv )  
  , m_max ( maxv ) 
  , m_cnt ()
{
  Ostap::Assert ( ( empty()
		    && m_min ==   std::numeric_limits<double>::max()
		    && m_max == - std::numeric_limits<double>::max() )
		  || ( !empty() && m_min <= m_max ) ,
		  "Invalid min/max/empty structure!" ,
		  "Ostap::Math::MinMaxValue!" ) ;
}
// ===========================================================================
// default constructor
// ===========================================================================
Ostap::Math::WMinMaxValue::WMinMaxValue ()
  : m_min (   std::numeric_limits<double>::max() ) 
  , m_max ( - std::numeric_limits<double>::max() ) 
  , m_cnt ()
{}
// ======================================================================
// constructor
// ===========================================================================
Ostap::Math::WMinMaxValue::WMinMaxValue
( const double                              minv ,
  const double                              maxv ,
  const Ostap::Math::WMinMaxValue::Counter& cnt  )
  : m_min ( minv )  
  , m_max ( maxv ) 
  , m_cnt ()
{
  Ostap::Assert ( ( empty()
		    && m_min ==   std::numeric_limits<double>::max()
		    && m_max == - std::numeric_limits<double>::max() ) 
		  || ( !empty() && m_min <= m_max ) ,		      
		  "Invalid min/max/empty structure!" ,
		  "Ostap::Math::WMinMaxValue!" ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
