// ============================================================================
// Include files 
// ============================================================================
//  STD&STL
// ============================================================================
#include <limits> 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Moments.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// Local
// ============================================================================
#include  "local_math.h"
#include  "status_codes.h"
// ============================================================================
/** @file 
 *  Implementation file for classes from the file Ostap/Moments.h
 *  @date 2020-06-07 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ===========================================================================
namespace
{
  // =========================================================================
  static_assert ( std::numeric_limits<Ostap::Math::Moment_<0>::size_type>::is_specialized  ,
		  "Ostap::Math::Moment_<0>::size_type is not specialized!" ) ;
  // =========================================================================
  static_assert ( std::numeric_limits<double>::is_specialized   ,
		  "std""numeric_limits<double> is not specialized!" ) ;
  // =========================================================================
  /** @var S_INVALID_MOMENT 
   *  invalid vaolue for moments
   */
  constexpr double s_INVALID_MOMENT { std::numeric_limits<double>::quiet_NaN() } ;
  // =========================================================================  
} // =========================================================================
// ===========================================================================
// get the invalid moment 
// ===========================================================================
double Ostap::Math::Moments::invalid_moment ()
{ return ::s_INVALID_MOMENT ; }  
// ===========================================================================
// virtual destructor
// ===========================================================================
Ostap::Math::Moment::~Moment(){}
// ===========================================================================
// virtual destructor
// ===========================================================================
Ostap::Math::WMoment::~WMoment(){}
// ===========================================================================
// return the value for invalid moment 
// ===========================================================================
double Ostap::Math::Moment::invalid_moment () const
{ return Ostap::Math::Moments::invalid_moment () ; }
// ===========================================================================  
// return the value for invalid moment 
// ===========================================================================
double Ostap::Math::WMoment::invalid_moment () const
{ return Ostap::Math::Moments::invalid_moment () ; } 
// ===========================================================================
// constructor 
// ===========================================================================
/*  full constructor
 *  @param size number of entries 
 *  @param sumw sum of weights 
 *  @param sumw sum of squared weights 
 */
// ===========================================================================
// (default) constructor 
// ===========================================================================
Ostap::Math::Moment_<0>::Moment_
( const Ostap::Math::Moment_<0>::size_type size )
  : m_size ( size )
{}
// ===========================================================================
// constructor from mu and previous 
// ===========================================================================
Ostap::Math::Moment_<1>::Moment_
( const Ostap::Math::Moment_<0>& prev , 
  const double                   mu   ,
  const double                   xmin ,
  const double                   xmax )
  : m_prev ( prev ) 
  , m_mu   ( mu   )
  , m_min  { std::min ( xmin ,   std::numeric_limits<double>::max () ) } 
  , m_max  { std::max ( xmax , - std::numeric_limits<double>::max () ) } 
{}
// ===========================================================================
Ostap::Math::WMoment_<0>::WMoment_
( const Ostap::Math::WMoment_<0>::size_type size  ,  
  const double                              sumw  ,
  const double                              sumw2 , 
  const double                              wmin  , 
  const double                              wmax  ) 
  : m_size ( size  )
  , m_w    ( sumw  )
  , m_w2   ( sumw2 )
  , m_wmin { std::min ( wmin ,   std::numeric_limits<double>::max () ) } 
  , m_wmax { std::max ( wmax , - std::numeric_limits<double>::max () ) } 
{
  Ostap::Assert ( std::isfinite ( m_w ) && std::isfinite ( m_w2 ) ,
		  "Invalid sumw/sumw2!"              , 
		  "Ostap::Math::WMoment_<0>"         ,
		  INVALID_PARS , __FILE__ , __LINE__ ) ;
  
  if ( 0 == m_size )
    {
      Ostap::Assert ( s_zero ( m_w ) && s_zero ( m_w2 )  ,  
		      "Non-zero sum of (squared) weights for empty counter!" , 
		      "Ostap::Math::WMoment_<0>"         ,
		      INVALID_PARS , __FILE__ , __LINE__ ) ;
      m_w  = 0 ;
      m_w2 = 0 ;       
    }
  //
  if ( s_zero ( m_w2 ) )
    {
      Ostap::Assert ( s_zero ( m_w )                     ,
		      "Non zero sum of weigth for zero sumw2 !" , 
		      "Ostap::Math::WMoment_<0>"         ,
		      INVALID_PARS , __FILE__ , __LINE__ ) ;      
      m_w  = 0 ;
      m_w2 = 0 ;
    }
  //
  if ( s_zero ( m_w ) ) { m_w = 0 ; }
  //
  Ostap::Assert ( 0 <= m_w2                          ,      
		  "Negative sum of squared weights!" , 
		  "Ostap::Math::WMoment_<0>"         , 
		  INVALID_PARS , __FILE__ , __LINE__ ) ;      
}
 // ===========================================================================
 /// constructor from mu and previous moment 
 // ===========================================================================
 Ostap::Math::WMoment_<1>::WMoment_
 ( const Ostap::Math::WMoment_<0>& prev ,
   const double                    mu   ,
  const double                    xmin ,
  const double                    xmax ) 
  : m_prev ( prev ) 
  , m_mu   ( mu   )
  , m_min  { std::min ( xmin ,   std::numeric_limits<double>::max () ) } 
  , m_max  { std::max ( xmax , - std::numeric_limits<double>::max () ) } 
{}   
// ============================================================================
 
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
  : m_cnt ( cnt ) 
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
		   "Inconsistent structure of two counters!" ,
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
  : m_cnt ( cnt ) 
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
		   "Inconsistent structure of two counters!" ,
		   "Ostap::Math::WLehmerMean" ) ;
  Ostap::Assert  ( s_equal ( m_lp.w  () , m_lpm1.w  () ) , 
		   "Inconsistent structure of two counters!" ,
		   "Ostap::Math::WLehmerMean" ) ;
  Ostap::Assert  ( s_equal ( m_lp.w2 () , m_lpm1.w2 () ) , 
		   "Inconsistent structure of two counters!" ,
		   "Ostap::Math::WLehmerMean" ) ;
}


// ===========================================================================
// GeometricMean::accumulate only positive entries 
// ===========================================================================
Ostap::Math::GeometricMean&
Ostap::Math::GeometricMean::add
( const double         x )
{
  if ( !std::isfinite ( x ) ) { return *this ; }
  if ( 0 < x && !s_zero ( x ) ) { m_log.add ( std::log2 ( x ) ) ; }
  return *this ;
}
// ===========================================================================
// HarmonicMean:: accumulate only non-zero entries 
// ===========================================================================
Ostap::Math::HarmonicMean&
Ostap::Math::HarmonicMean::add ( const double x )
{
  if ( !std::isfinite ( x ) ) { return *this ; }
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
  if ( !std::isfinite ( x ) ) { return *this ; }
  if ( ( 0 < x ) && !s_zero ( x ) ) { m_pow.add ( std::pow ( x , m_p ) ) ; }
  return *this ;
}
// ===========================================================================
// LehmerMean: accumulate only positive entries
// ===========================================================================
Ostap::Math::LehmerMean&
Ostap::Math::LehmerMean::add ( const double x )
{
  if ( !std::isfinite ( x ) ) { return *this ; }
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
  if ( !w || !std::isfinite ( x ) || !std::isfinite ( w ) ) { return *this ; }
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
  if ( !w || !std::isfinite ( x ) || !std::isfinite ( w ) ) { return *this ; }
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
  if ( !w || !std::isfinite ( x ) || !std::isfinite ( w ) ) { return *this ; }
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
  if ( !w || !std::isfinite ( x ) || !std::isfinite ( w ) ) { return *this ; }
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
                  "Ostap::Math::PowerMean"  ,
                  INVALID_ORDER , __FILE__  , __LINE__ ) ;
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
                  "Ostap::Math::LehmerMean" ,
                  INVALID_ORDER , __FILE__  , __LINE__ ) ;
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
                  "Ostap::Math::WLehmerMean" ,
                  INVALID_ORDER , __FILE__ , __LINE__ ) ;
  m_lp   .add ( x.m_lp   ) ;
  m_lpm1 .add ( x.m_lpm1 ) ;
  return *this ;
}
// ===========================================================================

// ============================================================================
// Finiteness 
// ============================================================================
bool Ostap::Math::GeometricMean::isfinite () const
{
  return m_log.isfinite ()
    && m_log.mean  () < std::numeric_limits<double>::max_exponent ;
}
// ============================================================================
// Finiteness
// ============================================================================
bool Ostap::Math::HarmonicMean::isfinite () const
{ return m_inv.isfinite () && m_inv.mean() ; }
// ============================================================================
// Finiteness
// ============================================================================
bool Ostap::Math::PowerMean::isfinite () const
{
  if ( !m_pow.isfinite() || !m_p ) { return false ; }
  const double x   = m_pow.mean()    ;
  if ( x <= 0                  ) { return false ; }
  const double l2x = std::log2 ( x ) ;
  return l2x / m_p < std::numeric_limits<double>::max_exponent ;
}
// ============================================================================
// Finiteness
// ============================================================================
bool Ostap::Math::LehmerMean::isfinite () const
{ return m_lp.isfinite () && m_lpm1.isfinite() && m_lpm1.mean () ; }
// ============================================================================
// Finiteness 
// ============================================================================
bool Ostap::Math::WGeometricMean::isfinite () const
{
  return m_log.isfinite ()
    && m_log.mean  () < std::numeric_limits<double>::max_exponent ;
}
// ============================================================================
// Finiteness
// ============================================================================
bool Ostap::Math::WHarmonicMean::isfinite () const
{ return m_inv.isfinite () && m_inv.mean() ; }
// ============================================================================
// Finiteness
// ============================================================================
bool Ostap::Math::WPowerMean::isfinite () const
{
  if ( !m_pow.isfinite() || !m_p ) { return false ; }
  const double x   = m_pow.mean()    ;
  if ( x <= 0                    ) { return false ; }
  const double l2x = std::log2 ( x ) ;
  return l2x / m_p < std::numeric_limits<double>::max_exponent ;
}
// ============================================================================
// Finiteness
// ============================================================================
bool Ostap::Math::WLehmerMean::isfinite () const
{ return m_lp.isfinite () && m_lpm1.isfinite() && m_lpm1.mean () ; }



// ============================================================================
//                                                                      The END 
// ============================================================================
