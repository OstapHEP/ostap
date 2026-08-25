// ============================================================================
// Incldue files 
// ============================================================================
// STD&STL
// ============================================================================
#include <memory>
#include <map> 
#include <cstring>
#include <numeric>
#include <algorithm>
#include <tuple>
#include <array>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/ECDF.h"
#include "Ostap/EPDF.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// Local
// ============================================================================
#include "status_codes.h"
#include "local_math.h"
#include "local_hash.h"
#include "syncedcache.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::ECDF
 *  @see Ostap::Math::ECDF
 *  @date 2024-09-16 
 *  @author Vanya BELYAEV 
 */
// ============================================================================
namespace
{
  // ==========================================================================
  /// actual kernel function type 
  typedef double (*FUNPTR)(double) ; 
  typedef std::tuple<
    Ostap::Math::DensityEstimator::Kernel,
    FUNPTR,
    double, 
    double, 
    double, 
    double> RECORD ;
  //
  typedef std::array<RECORD,int(Ostap::Math::DensityEstimator::Kernel::Last)> TABLE ;  
  const TABLE s_TABLE { 
    RECORD ( Ostap::Math::DensityEstimator::Kernel::Uniform,     &Ostap::Math::k_uniform     ,1,1./3   ,1./2                   ,0.929) ,
    RECORD ( Ostap::Math::DensityEstimator::Kernel::Triangular,  &Ostap::Math::k_triangular  ,1,1./6   ,2./3                   ,0.986) ,
    RECORD ( Ostap::Math::DensityEstimator::Kernel::Epanechnikov,&Ostap::Math::k_epanechnikov,1,1./5   ,3./5                   ,1.000) ,
    RECORD ( Ostap::Math::DensityEstimator::Kernel::Quartic,     &Ostap::Math::k_quartic     ,1,1./7   ,5./7                   ,0.994) ,
    RECORD ( Ostap::Math::DensityEstimator::Kernel::Triweight,   &Ostap::Math::k_triweight   ,1,1./9   ,350./429               ,0.987) ,
    RECORD ( Ostap::Math::DensityEstimator::Kernel::Tricube,     &Ostap::Math::k_tricube     ,1,35./243,175./447               ,0.998) ,
    RECORD ( Ostap::Math::DensityEstimator::Kernel::Gaussian,    &Ostap::Math::k_gaussian    ,0,1      ,0.5/std::sqrt(s_pi)    ,0.951) ,
    RECORD ( Ostap::Math::DensityEstimator::Kernel::Cosine,      &Ostap::Math::k_cosine      ,1,1-8/(s_pi*s_pi) , s_pi*s_pi/16 ,0.999) ,
    RECORD ( Ostap::Math::DensityEstimator::Kernel::Logistic,    &Ostap::Math::k_logistic    ,0,s_pi*s_pi/3 , 1.6              ,0.887) ,
    RECORD ( Ostap::Math::DensityEstimator::Kernel::Sigmoid,     &Ostap::Math::k_sigmoid     ,0,s_pi*s_pi/4 , 2/(s_pi*s_pi)    ,0.843) ,
  } ;
  // ==========================================================================
  typedef std::map<Ostap::Math::DensityEstimator::Kernel,std::string> KMAP ;
  const KMAP s_KMAP = {
    { Ostap::Math::DensityEstimator::Kernel::Uniform      , "Uniform"      } , 
    { Ostap::Math::DensityEstimator::Kernel::Rectangular  , "Rectangualar" } , 
    { Ostap::Math::DensityEstimator::Kernel::Boxcar       , "Boxcar"       } ,
    { Ostap::Math::DensityEstimator::Kernel::Triangular   , "Triangular"   } ,
    { Ostap::Math::DensityEstimator::Kernel::Epanechnikov , "Epanechnikov" } , 
    { Ostap::Math::DensityEstimator::Kernel::Parabolic    , "Parabolic"    } , 
    { Ostap::Math::DensityEstimator::Kernel::Quartic      , "Quartic"      } , 
    { Ostap::Math::DensityEstimator::Kernel::Biweight     , "Biweight"     } , 
    { Ostap::Math::DensityEstimator::Kernel::Triweight    , "Triweight"    } , 
    { Ostap::Math::DensityEstimator::Kernel::Tricube      , "Tricube"      } , 
    { Ostap::Math::DensityEstimator::Kernel::Gaussian     , "Gaussian"     } , 
    { Ostap::Math::DensityEstimator::Kernel::Cosine       , "Cosine"       } , 
    { Ostap::Math::DensityEstimator::Kernel::Logistic     , "Logistic"     } ,
    { Ostap::Math::DensityEstimator::Kernel::Sigmoid      , "Sigmoid"      } } ;
  // ==========================================================================
} //                                             The end of anonymous namespace 
// ============================================================================
double Ostap::Math::DensityEstimator::kernel 
( const double                                u , 
  const Ostap::Math::DensityEstimator::Kernel k )
{
  //
  Ostap::Assert ( Kernel::First <= k && k <  Kernel::Last , 
                  "Invalid Kernel!"       , 
                  "Ostap::Math::DensityEstimator::kernel"  ,
                  INVALID_KERNEL , __FILE__, __LINE__ ) ;
  //
  switch ( k )
  {
  case Kernel::Uniform      : return Ostap::Math::k_uniform      ( u ) ;
    // case Kernel::Rectangular  : return Ostap::Math::k_uniform      ( u ) ;
    // case Kernel::Boxcar       : return Ostap::Math::k_uniform      ( u ) ;
  case Kernel::Triangular   : return Ostap::Math::k_triangular   ( u ) ;
  case Kernel::Epanechnikov : return Ostap::Math::k_epanechnikov ( u ) ;
    // case Kernel::Parabolic    : return Ostap::Math::k_epanechnikov ( u ) ;
  case Kernel::Quartic      : return Ostap::Math::k_quartic      ( u ) ;
    // case Kernel::Biweight     : return Ostap::Math::k_quartic      ( u ) ;
  case Kernel::Triweight    : return Ostap::Math::k_triweight    ( u ) ;
  case Kernel::Tricube      : return Ostap::Math::k_tricube      ( u ) ;
  case Kernel::Gaussian     : return Ostap::Math::k_gaussian     ( u ) ;
  case Kernel::Cosine       : return Ostap::Math::k_cosine       ( u ) ;
  case Kernel::Logistic     : return Ostap::Math::k_logistic     ( u ) ;
  case Kernel::Sigmoid      : return Ostap::Math::k_sigmoid      ( u ) ;
  default                   : return Ostap::Math::k_epanechnikov ( u ) ; // DEFAULT!   
  } 
  return Ostap::Math::k_epanechnikov ( u ) ; // DEFAULT!   
}
// ============================================================================

// ============================================================================
// Cosine kernel 
// ============================================================================
double Ostap::Math::k_cosine       ( const double u ) 
{ return std::abs ( u ) <= 1 ? 0.25 * s_pi * std::cos ( 0.5 * s_pi * u ) : 0.0 ;  } 
// ============================================================================
// Logistic Kernel
// ============================================================================
double Ostap::Math::k_logistic     ( const double u )
{ return std::abs ( u ) < s_EXP_OVERFLOW ? 1 / ( 2.0 * std::cosh ( u ) + 2 ) :  0.0 ;  }
// ============================================================================
// sigmoid kernel 
// ============================================================================
double Ostap::Math::k_sigmoid      ( const double u ) 
{ return std::abs ( u ) < s_EXP_OVERFLOW ? s_1_pi / std::cosh ( u ) : 0.0 ; }
// ============================================================================
// Gaussian kernel 
// ============================================================================
double Ostap::Math::k_gaussian    ( const double u ) { return gauss_pdf ( u ) ; }
// ============================================================================
// get the "optimal" value for the smoothing parameter 
// ============================================================================
double Ostap::Math::DensityEstimator::hopt
( const  Ostap::Math::ECDF& data )
{
  if ( 2 > data.N() ) { return -1 ; }
  //
  Ostap::StatEntity cnt {} ;
  for ( const auto v : data ) { cnt += v ; } 
  //
  double value = cnt.rms() ;
  //
  if ( 4 <= data.N () )
    {
      // the first * third  quartiles
      const std::size_t i1 = static_cast<std::size_t> ( 0.25 * data.N () ) ;
      const std::size_t i3 = static_cast<std::size_t> ( 0.75 * data.N () ) ;
      const double x1 = data [ i1 ] ;
      const double x3 = data [ i3 ] ;
      if ( x1 < x3 ) { value = std::min ( value , ( x3 - x1 ) / 1.34 ) ;} 
    }
  //
  return 0.9 * value * std::pow ( data.N() , -0.2 );
}
// ============================================================================
// get the "optimal" value for the smoothing parameter 
// ============================================================================
double Ostap::Math::DensityEstimator::hopt
( const  Ostap::Math::WECDF& data )
{
  if ( 2 > data.N() || 2 > data.nEff() ) { return -1 ; }
  //
  Ostap::WStatEntity cnt {} ;
  for ( const auto v : data ) { cnt.add ( v.first , v.second ) ; } 
  //
  double value = cnt.rms() ;
  //
  if ( 4 <= data.N () && 4 <= data.nEff() )
    {      
      // the first * third  quartiles, well... not really... @todo FIX ME!
      const std::size_t i1 = static_cast<std::size_t> ( 0.25 * data.N () ) ;
      const std::size_t i3 = static_cast<std::size_t> ( 0.75 * data.N () ) ;
      const double x1 = data [ i1 ].first ;
      const double x3 = data [ i3 ].first ;
      if ( x1 < x3 ) { value = std::min ( value , ( x3 - x1 ) / 1.34 ) ;} 
    }
  //
  return 0.9 * value * std::pow ( data.nEff () , -0.2 );
}
// ============================================================================
// get the name for the given Kernel 
// ============================================================================
const std::string& Ostap::Math::DensityEstimator::name
( const Ostap::Math::DensityEstimator::Kernel k )
{
  KMAP::const_iterator found = s_KMAP.find ( k ) ;
  /// unknown kernel 
  static const std::string s_UNKNOWN { "<Unknown>" } ;
  return s_KMAP.end() == found ? found->second : s_UNKNOWN ; 
}  
// ============================================================================
/* create the empirical PDF from empirical CDF 
 *  @attention data are not copied!
 */
// ============================================================================
Ostap::Math::EPDF::EPDF
( const Ostap::Math::ECDF&                    cdf ,
  const Ostap::Math::DensityEstimator::Kernel k   ,
  const double                                h   ) 
  : m_cdf   ( cdf )
  , m_k     ( k   )
  , m_h     ( h   )
{
  Ostap::Assert ( Ostap::Math::DensityEstimator::Kernel::First <= k &&
                  Ostap::Math::DensityEstimator::Kernel::Last  >  k &&
                  static_cast<std::size_t> ( k ) < s_TABLE.size () , 
                  "Invalid Kernel!"   , 
                  "Ostap::Math::EPDF" ,
                  INVALID_KERNEL , __FILE__, __LINE__ ) ;
  // check the smoothing parameter 
  if ( m_h <= 0 ) { m_h = Ostap::Math::DensityEstimator::hopt ( m_cdf ) ; }
  Ostap::Assert ( 0 < m_h             ,
                  "Invalid smoothing parameter" ,
                  "Ostap::Math::EPDF" ,
                  INVALID_SMOOTH , __FILE__ , __LINE__ ) ;  
}
// =============================================================================
// get the PDF
// =============================================================================
double Ostap::Math::EPDF::evaluate ( const double x ) const
{
  FUNPTR k = std::get<1> ( s_TABLE [ static_cast<std::size_t> ( m_k ) ] ) ;
  double s = std::get<2> ( s_TABLE [ static_cast<std::size_t> ( m_k ) ] ) ;
  if  ( s <= 0 ) { s = 5 ; } ;
  //
  const double xmn = x - s * m_h ;  
  const double xmx = x + s * m_h ;
  //
  if ( m_cdf.xmax () < xmn || m_cdf.xmin () > xmx ) { return 0 ; }
  Ostap::Math::ECDF::iterator imin = std::lower_bound ( m_cdf.begin () , m_cdf.end () , xmn ) ;
  if ( imin == m_cdf.end() ) { return 0 ; }
  Ostap::Math::ECDF::iterator imax = std::upper_bound ( imin + 1 , m_cdf.end () , xmx ) ;
  double       value = 0 ;
  const double ih    = 1 / m_h ;
  for ( Ostap::Math::ECDF::iterator i = imin ; i < imax ; ++i )
    {
      const double u = ih * ( x - *i ) ;
      value += (*k) ( u ) ;
    }
  return ih * value / m_cdf.N() ;
}
// ============================================================================
/* create the emppirical PDF from empirical CDF 
 *  @attention data are not copied!
 */
// ============================================================================
Ostap::Math::WEPDF::WEPDF
( const Ostap::Math::WECDF&                   cdf ,
  const Ostap::Math::DensityEstimator::Kernel k   ,
  const double                                h   ) 
  : m_cdf ( cdf )
  , m_k   ( k   )
  , m_h   ( h   )
{
  Ostap::Assert ( Ostap::Math::DensityEstimator::Kernel::First <= k &&
                  Ostap::Math::DensityEstimator::Kernel::Last  >  k &&
                  static_cast<std::size_t> ( k ) < s_TABLE.size () , 
                  "Invalid Kernel!"    , 
                  "Ostap::Math::WEPDF" ,
                  INVALID_KERNEL , __FILE__, __LINE__ ) ;
  // check the smoothing parameter 
  if ( m_h <= 0 ) { m_h = Ostap::Math::DensityEstimator::hopt ( m_cdf ) ; }
  Ostap::Assert ( 0 < m_h             ,
                  "Invalid smoothing parameter" ,
                  "Ostap::Math::WEPDF" ,
                  INVALID_SMOOTH , __FILE__ , __LINE__ ) ;  
}
// =============================================================================
// get the PDF
// =============================================================================
double Ostap::Math::WEPDF::evaluate ( const double x ) const
{
  FUNPTR k = std::get<1> ( s_TABLE [ static_cast<std::size_t> ( m_k ) ] ) ;
  double s = std::get<2> ( s_TABLE [ static_cast<std::size_t> ( m_k ) ] ) ;
  if  ( s <= 0 ) { s = 5 ; } ;
  //
  const double xmn = x - s * m_h ;  
  const double xmx = x + s * m_h ;
  //
  if ( m_cdf.xmax () < xmn || m_cdf.xmin () > xmx ) { return 0 ; }
  //
  const Ostap::Math::WECDF::COMPARE cmp {} ;
  //
  Ostap::Math::WECDF::iterator imin = std::lower_bound ( m_cdf.begin () , m_cdf.end () , xmn , cmp ) ;
  if ( imin == m_cdf.end() ) { return 0 ; }
  Ostap::Math::WECDF::iterator imax = std::upper_bound ( imin + 1 , m_cdf.end () , xmx , cmp ) ;
  double       value = 0 ;
  const double ih    = 1 / m_h ;
  for ( Ostap::Math::WECDF::iterator i = imin ; i < imax ; ++i )
    {
      const double w = i->second ;
      if ( !w ) { continue ; } 
      const double u = ih * ( x - i->first ) ;
      value += w * (*k) ( u ) ;
    }
  return ih * value / m_cdf.sumw ();
}
// ============================================================================
// update smoothing parameters
// ============================================================================
bool Ostap::Math::EPDF::setH      ( const double h )
{
  if ( ( 0 < h ) && s_equal ( m_h , h ) ) { return false ; }
  m_h = h ;
  // check the smoothing parameter 
  if ( m_h <= 0 ) { m_h = Ostap::Math::DensityEstimator::hopt ( m_cdf ) ; }
  return true ;
}
// ============================================================================
// update kernel kernel
// ============================================================================
bool Ostap::Math::EPDF::setKernel
( const Ostap::Math::DensityEstimator::Kernel k )
{
  Ostap::Assert ( Ostap::Math::DensityEstimator::Kernel::First <= k &&
                  Ostap::Math::DensityEstimator::Kernel::Last  >  k &&
                  static_cast<std::size_t> ( k ) < s_TABLE.size () , 
                  "Invalid Kernel!"    , 
                  "Ostap::Math::EPDF::setKernel" ,
                  INVALID_KERNEL , __FILE__, __LINE__ ) ;
  if ( k == m_k ) { return false ; }
  m_k = k ;
  return true ;
}
// ============================================================================
// update smoothing parameters
// ============================================================================
bool Ostap::Math::WEPDF::setH      ( const double h )
{
  if ( ( 0 < h ) && s_equal ( m_h , h ) ) { return false ; }
  m_h = h ;
  // check the smoothing parameter 
  if ( m_h <= 0 ) { m_h = Ostap::Math::DensityEstimator::hopt ( m_cdf ) ; }
  return true ;
}
// ============================================================================
// update kernel kernel
// ============================================================================
bool Ostap::Math::WEPDF::setKernel
( const Ostap::Math::DensityEstimator::Kernel k )
{
  Ostap::Assert ( Ostap::Math::DensityEstimator::Kernel::First <= k &&
                  Ostap::Math::DensityEstimator::Kernel::Last  >  k &&
                  static_cast<std::size_t> ( k ) < s_TABLE.size () , 
                  "Invalid Kernel!"    , 
                  "Ostap::Math::WEPDF::setKernel" ,
                  INVALID_KERNEL , __FILE__, __LINE__ ) ;
  if ( k == m_k ) { return false ; }
  m_k = k ;
  return true ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================

