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
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
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
  typedef std::array<RECORD,Ostap::Math::DensityEstimator::Last+1> TABLE ;  
  const TABLE s_TABLE { 
    RECORD ( Ostap::Math::DensityEstimator::Uniform,     &Ostap::Math::k_uniform     ,1,1./3   ,1./2                   ,0.929) ,
    RECORD ( Ostap::Math::DensityEstimator::Triangular,  &Ostap::Math::k_triangular  ,1,1./6   ,2./3                   ,0.986) ,
    RECORD ( Ostap::Math::DensityEstimator::Epanechnikov,&Ostap::Math::k_epanechnikov,1,1./5   ,3./5                   ,1.000) ,
    RECORD ( Ostap::Math::DensityEstimator::Quartic,     &Ostap::Math::k_quartic     ,1,1./7   ,5./7                   ,0.994) ,
    RECORD ( Ostap::Math::DensityEstimator::Triweight,   &Ostap::Math::k_triweight   ,1,1./9   ,350./429               ,0.987) ,
    RECORD ( Ostap::Math::DensityEstimator::Tricube,     &Ostap::Math::k_tricube     ,1,35./243,175./447               ,0.998) ,
    RECORD ( Ostap::Math::DensityEstimator::Gaussian,    &Ostap::Math::k_gaussian    ,0,1      ,0.5/std::sqrt(M_PI)    ,0.951) ,
    RECORD ( Ostap::Math::DensityEstimator::Cosine,      &Ostap::Math::k_cosine      ,1,1-8/(M_PI*M_PI) , M_PI*M_PI/16 ,0.999) ,
    RECORD ( Ostap::Math::DensityEstimator::Logistic,    &Ostap::Math::k_logistic    ,0,M_PI*M_PI/3 , 1.6              ,0.887) ,
    RECORD ( Ostap::Math::DensityEstimator::Sigmoid,     &Ostap::Math::k_sigmoid     ,0,M_PI*M_PI/4 , 2/(M_PI*M_PI)    ,0.843) ,
  } ;
  // ==========================================================================
}
// ============================================================================
double Ostap::Math::DensityEstimator::kernel 
( const double                                u , 
  const Ostap::Math::DensityEstimator::Kernel k )
{
  //
  Ostap::Assert ( Uniform <= k && k <=  Last , 
		  "Invalid Kernel!" , 
		  "Ostap::Math::DensityEstimator::kernel"  ,
		  INVALID_KERNEL , __FILE__, __LINE__ ) ;
  //
  return 
    Gaussian     == k ? Ostap::Math::k_gaussian     ( u ) :
    Epanechnikov == k ? Ostap::Math::k_epanechnikov ( u ) :
    Triangular   == k ? Ostap::Math::k_triangular   ( u ) :
    Quartic      == k ? Ostap::Math::k_quartic      ( u ) : 
    Triweight    == k ? Ostap::Math::k_triweight    ( u ) :
    Tricube      == k ? Ostap::Math::k_tricube      ( u ) : 
    Cosine       == k ? Ostap::Math::k_cosine       ( u ) : 
    Logistic     == k ? Ostap::Math::k_logistic     ( u ) : 
    Sigmoid      == k ? Ostap::Math::k_sigmoid      ( u ) : Ostap::Math::k_uniform ( u ) ;  
}
// ============================================================================

// ============================================================================
// Cosine kernel 
// ============================================================================
double Ostap::Math::k_cosine       ( const double u ) 
{ return std::abs ( u ) <= 1 ? 0.25 * M_PI * std::cos ( 0.5 * M_PI * u ) : 0.0 ;  } 
// ============================================================================
// Logistic Kernel
// ============================================================================
double Ostap::Math::k_logistic     ( const double u )
{ return std::abs ( u ) < s_EXP_OVERFLOW ? 1 / ( 2.0 * std::cosh ( u ) + 2 ) :  0.0 ;  }
// ============================================================================
// sigmoid kernel 
// ============================================================================
double Ostap::Math::k_sigmoid      ( const double u ) 
{ return std::abs ( u ) < s_EXP_OVERFLOW ? s_PIi / std::cosh ( u ) : 0.0 ; }
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
  Ostap::Assert ( Ostap::Math::DensityEstimator::Uniform <= k &&
		  Ostap::Math::DensityEstimator::Last    >= k && k < s_TABLE.size () , 
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
  FUNPTR k = std::get<1> ( s_TABLE[m_k] ) ;
  double s = std::get<2> ( s_TABLE[m_k] ) ;
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
  Ostap::Assert ( Ostap::Math::DensityEstimator::Uniform <= k &&
		  Ostap::Math::DensityEstimator::Last    >= k && k < s_TABLE.size () , 
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
  FUNPTR k = std::get<1> ( s_TABLE[m_k] ) ;
  double s = std::get<2> ( s_TABLE[m_k] ) ;
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
  Ostap::Assert ( Ostap::Math::DensityEstimator::Uniform <= k &&
		  Ostap::Math::DensityEstimator::Last    >= k && k < s_TABLE.size () , 
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
  Ostap::Assert ( Ostap::Math::DensityEstimator::Uniform <= k &&
		  Ostap::Math::DensityEstimator::Last    >= k && k < s_TABLE.size () , 
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

