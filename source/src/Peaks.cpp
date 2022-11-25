// ============================================================================
// Include files 
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <array>
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_sf_exp.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/Peaks.h"
#include "Ostap/MoreMath.h"
#include "Ostap/qMath.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/Polynomials.h"
#include "Ostap/ToStream.h"
// ============================================================================
//  Local
// ============================================================================
#include "Exception.h"
#include "local_math.h"
#include "local_gsl.h"
#include "local_hash.h"
#include "gauss.h"      
#include "Integrator1D.h"      
#include "syncedcache.h"  // the cache 
// ============================================================================
/** @file 
 *  implmentation file for classes from the file Ostap/Peaks.h
 *  @author 
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-04-19
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  /** evaluate the helper function \f[ f = \frac{\sinh (x) }{x} \f]
  *  it allows to calculate Novosibirsk's function in efficient and regular way
  *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
  *  @date 2010-05-23
  */
  double x_sinh ( const double x , double precision = s_APRECISION )
  {
    //
    if      ( s_equal   ( x , 0 )    ) { return 1 ; }
    else if ( std::fabs ( x ) < 0.1  )
    {
      double result = 1.0  ;
      double delta  = x    ;
      //
      precision = std::fabs ( precision ) ;
      precision = std::min  ( precision , std::fabs ( s_APRECISION_TAIL ) ) ;
      unsigned int n = 1 ;
      //
      do
      {
        delta  *= x * x / ( ( n + 1 )  * ( n + 2 ) ) ;
        result += delta ;
        n      += 2 ;
      }
      while ( std::fabs ( delta ) > 0.1 * precision && n < 10000 ) ;
      //
      return result ;
    }
    //
    if ( 100 < std::fabs ( x )  ) { return  s_INFINITY ; }
    //
    // the generic evaluation
    //
    return std::sinh ( x ) / x ;
  }
  // ==========================================================================
  // Crystal Ball & Co 
  // ==========================================================================
  /** @var s_TRUNC
   *  trunkating parameter for CrystalBall-functions
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double      s_TRUNC = 15.0 ;
  // ==========================================================================
  /** evaluate very simple power-law integral
   *  
   *  \f[ I = \int_{x_{low}}^{x_{high}} \left( \frac{A}{B+Cx}\right)^{N} dx \f]
   * 
   *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
   */
  double tail_integral 
  ( const double A    , 
    const double B    , 
    const double C    , 
    const double N    , 
    const double low  , 
    const double high )
  {
    //
    // few really very simple cases:
    if      ( s_equal ( N , 0 ) ) { return high - low ; }
    else if ( s_equal ( A , 0 ) ) { return 0 ; }
    else if ( s_equal ( C , 0 ) ) { return std::pow ( A / B , N ) * ( high - low ) ; }
    //
    // again the trivial cases 
    if ( s_equal ( low , high ) ) { return 0 ; }
    else if      ( low > high   ) { return -1 * tail_integral ( A , B , C , N , high , low ) ; }
    //
    //  y = (B+C*x)/A
    //
    const double y_low  = ( B + C * low  ) / A ;
    const double y_high = ( B + C * high ) / A ;
    //
    // the special case 
    if ( s_equal ( N , 1 ) ) { return A / C * my_log ( y_high / y_low ) ; } // RETURN 
    //
    // the regular case 
    return A / C * ( std::pow ( y_high , 1 - N ) - 
                     std::pow ( y_low  , 1 - N ) ) / ( 1 - N ) ;
  }
  // ========================================================================== 
  // Atlas/Zeus& Co 
  // ==========================================================================
  /** @var S_ATLAL 
   *  magic constant - integral for Atlas function 
   *  @see Ostap::Math::Atlas 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date 2015-08-21
   */   
  const double s_ATLAS = 3.052369876253939 ;
  // ==========================================================================
  // Sinh-asinh 
  // ==========================================================================
  inline double shash
  ( const double x   , 
    const double eps , 
    const double dlt ) 
  {
    const double y = eps + dlt * std::asinh ( x ) ;
    return 
      (     GSL_LOG_DBL_MAX < y ) ?    s_INFINITY :
      ( -1* GSL_LOG_DBL_MAX > y ) ? -1*s_INFINITY : std::sinh ( y ) ;
  }
  // // ==========================================================================
  // /en-t'
  // // ==========================================================================
  // inline double student_cdf (  const double t , const double nu ) 
  // {
  //   const double xt    = nu / ( t * t + nu ) ;
  //   const double value = 0.5 * gsl_sf_beta_inc ( 0.5 * nu , 0.5 , xt ) ;
  //   return t >= 0 ? 1 - value : value ;
  // }
  // ==========================================================================
  /** product of Gaussian PDF and Mill's ratio 
   *  \f$ f(x; a ) = \phi ( x ) R( a + x ) \f$ 
   */
  // inline double gauss_mills 
  // ( const double x , 
  //  const double a ) 
  // { 
  //   return Ostap::Math::gauss_pdf ( x ) * Ostap::Math::mills_normal ( a + x ) ; 
  // } ;
  // ==========================================================================
}
// ============================================================================

// ============================================================================
// BifurcatedGauss
// ============================================================================
/*  constructor from all parameters
 *  @param peak    the peak posiion
 *  @param sigma_L (left sigma)
 *  @param sigma_R (right sigma)
 */
// ============================================================================
  Ostap::Math::BifurcatedGauss::BifurcatedGauss
  ( const double peak   ,
    const double sigmaL ,
    const double sigmaR )
  : m_peak    ( peak )
  , m_sigmaL  ( std::fabs ( sigmaL ) )
  , m_sigmaR  ( std::fabs ( sigmaR ) )
//
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::BifurcatedGauss::~BifurcatedGauss (){}
// ============================================================================
// evaluate Bifurcated Gaussian
// ============================================================================
double Ostap::Math::BifurcatedGauss::evaluate ( const double x ) const
{
  const double dx   = x - m_peak ;
  const double norm = s_SQRTPIHALF * ( sigmaL() + sigmaR() ) ;
  //
  return
    dx < 0 ?
    my_exp ( -0.5 * dx * dx / ( sigmaL () * sigmaL () ) ) / norm :
    my_exp ( -0.5 * dx * dx / ( sigmaR () * sigmaR () ) ) / norm ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::BifurcatedGauss::integral () const { return 1 ; }
// ============================================================================
//  get CDF 
// ============================================================================
double Ostap::Math::BifurcatedGauss::cdf ( const double x ) const 
{
  // left half-gaussian
  if ( x <= m_peak )
  {
  const double sigma = sigmaL () ;
  const double sf    = s_SQRT2i  / sigma  ;
  const double nf    = sigma    / ( sigmaL() + sigmaR() ) ;
  const double b     = ( x - m_peak ) * sf ;
    return std::erfc ( -b ) * nf ; // RETURN
  }
  //
  const double bias = sigmaL() / ( sigmaL() + sigmaR() ) ;
  return bias + integral ( m_peak , x ) ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::BifurcatedGauss::integral ( const double low  ,
                                                const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }                         // RETURN
  if (           low > high   ) { return - integral ( high , low ) ; } // RETURN
  //
  // left half-gaussian
  if       ( high <= m_peak )
  {
    const double sigma = sigmaL () ;
    const double sf    = s_SQRT2i  / sigma  ;
    const double nf    = sigma    / ( sigmaL() + sigmaR() ) ;
    const double a     = ( low  - m_peak ) * sf ;
    const double b     = ( high - m_peak ) * sf ;
    return  ( std::erf ( b ) -  std::erf ( a ) ) * nf ; // RETURN
  }
  // right half-gaussian
  else if ( low >= m_peak )
  {
    const double sigma = sigmaR () ;
    const double sf    = s_SQRT2i  / sigma  ;
    const double nf    = sigma    / ( sigmaL() + sigmaR() ) ;
    const double a     = ( low  - m_peak ) * sf ;
    const double b     = ( high - m_peak ) * sf ;
    return  ( std::erf ( b ) -  std::erf ( a ) ) * nf ; // RETURN
  }
  // split into two intervals 
  return
    integral ( low    , m_peak ) +
    integral ( m_peak , high   ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::BifurcatedGauss::tag () const 
{
  static const std::string s_name = "BiFurcatedGauss" ;
  return Ostap::Utils::hash_combiner ( s_name , m_peak , m_sigmaL , m_sigmaR ) ; 
}
// ============================================================================

// ============================================================================
bool Ostap::Math::BifurcatedGauss::setSigmaL ( const double value )
{
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( m_sigmaL , value_ ) ) { return false ; }
  m_sigmaL = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::BifurcatedGauss::setSigmaR ( const double value )
{
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( m_sigmaR , value_ ) ) { return false ; }
  m_sigmaR = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::BifurcatedGauss::setPeak( const double value )
{
  if ( s_equal ( m_peak , value ) ) { return false ; }
  m_peak   = value  ;
  //
  return true ;
}
// ============================================================================
/* constructor from all parameters
 *  @param peak     the peak posiion
 *  @param sigmaL   the sigma for first component
 *  @param fraction the fraction of first component 
 *  @param scale    the ratio of sigmas for second and first components
 */
// ============================================================================
Ostap::Math::DoubleGauss::DoubleGauss
( const double peak     ,
  const double sigma    , 
  const double fraction , 
  const double scale    ) 
  : m_peak     ( peak               ) 
  , m_sigma    ( std::abs ( sigma ) )
  , m_fraction ( std::min ( std::max ( fraction , 0.0 ), 1.0 ) ) 
  , m_scale    ( std::abs ( scale ) )
{}
// ============================================================================
bool Ostap::Math::DoubleGauss::setPeak   ( const double value ) 
{
  if ( s_equal ( value , m_peak ) ) { return false ; }
  m_peak = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::DoubleGauss::setSigma ( const double value ) 
{
  const double value_ =  std::abs (value ) ;
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  m_sigma = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::DoubleGauss::setFraction ( const double value ) 
{
  const double value_ = std::min ( std::max ( value , 0.0 ), 1.0 ) ;
  if ( s_equal ( value_ , m_fraction ) ) { return false ; }
  m_fraction = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::DoubleGauss::setScale ( const double value ) 
{
  const double value_ =  std::abs (value ) ;
  if ( s_equal ( value_ , m_scale ) ) { return false ; }
  m_scale = value_ ;
  return true ;
}
// ============================================================================
// evaluate  double Gaussian
// ============================================================================
double Ostap::Math::DoubleGauss::pdf ( const double x ) const 
{
  const double mu       = m_peak     ;
  const double sigma    = m_sigma    ;
  const double scale    = m_scale    ;
  const double fraction = m_fraction ;
  //
  const double sigma2   =  scale *  sigma ;
  //
  const double dx1      = ( x - mu ) / sigma  ;
  const double dx2      = ( x - mu ) / sigma2 ;
  //
  const double  f1 = fraction ;
  const double  f2 = 1 - f1   ;
  //
  static const double s_norm = 1.0 / std::sqrt ( 2.0 * M_PI  ) ;
  //
  return 
    s_norm * ( f1 * std::exp ( -0.5 * dx1 * dx1 ) / sigma  +
               f2 * std::exp ( -0.5 * dx2 * dx2 ) / sigma2 ) ;  
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::DoubleGauss::integral 
( const double xmin , 
  const double xmax ) const 
{
  const double mu       = m_peak     ;
  const double sigma    = m_sigma    ;
  const double scale    = m_scale    ;
  const double fraction = m_fraction ;
  //
  const double sigma2   =  scale *  sigma ;
  //
  const double  f1 = fraction ;
  const double  f2 = 1 - f1 ;
  //
  static const double s_isqrt2 = 1.0 / std::sqrt ( 2.0 ) ;
  //
  const double ixscale1 = s_isqrt2 / sigma  ;
  const double ixscale2 = s_isqrt2 / sigma2 ;
  //
  const double r1 = 
    std::erf ( ( xmax - mu ) * ixscale1 ) - 
    std::erf ( ( xmin - mu ) * ixscale1 ) ;
  //
  const double r2 = 
    std::erf ( ( xmax - mu ) * ixscale2 ) - 
    std::erf ( ( xmin - mu ) * ixscale2 ) ;
  //
  return 0.5 * ( f1 * r1 + f2 * r2  ) ;
} 
// ============================================================================
// get cdf 
// ============================================================================
double Ostap::Math::DoubleGauss::cdf ( const double x )  const 
{
  const double mu       = m_peak     ;
  const double sigma    = m_sigma    ;
  const double scale    = m_scale    ;
  const double fraction = m_fraction ;
  //
  const double sigma2   =  scale *  sigma ;
  //
  const double  f1 = fraction ;
  const double  f2 =  1 - f1  ;
  //
  static const double s_isqrt2 = 1.0 / std::sqrt ( 2.0 ) ;
  //
  const double ixscale1 = s_isqrt2 / sigma  ;
  const double ixscale2 = s_isqrt2 / sigma2 ;
  //
  const double r1 = std::erf ( ( x - mu ) * ixscale1 ) ;
  const double r2 = std::erf ( ( x - mu ) * ixscale2 ) ;
  //
  return  0.5 * ( f1 * ( r1 + 1  ) + f2 * ( r2  + 1 ) ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::DoubleGauss::tag () const 
{
  static const std::string s_name = "DoubleGauss" ;
  return Ostap::Utils::hash_combiner ( s_name , m_peak , m_sigma , m_fraction , m_scale ) ; 
}
// ============================================================================





// ============================================================================
// Gauss
// ============================================================================
/*  constructor from all parameters
 *  @param peak    the peak posiion
 *  @param sigma_L (left sigma)
 *  @param sigma_R (right sigma)
 */
// ============================================================================
Ostap::Math::Gauss::Gauss
( const double peak   ,
  const double sigma  )
  : m_peak  ( peak )
  , m_sigma ( std::fabs ( sigma ) )
    //
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Gauss::~Gauss (){}
// ============================================================================
// evaluate Bifurcated Gaussian
// ============================================================================
double Ostap::Math::Gauss::evaluate ( const double x ) const
{
  const double dx   = ( x - m_peak ) / m_sigma ;
  const double norm = s_SQRTPIHALF * m_sigma  ;
  //
  return my_exp ( -0.5 * dx * dx ) / norm ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::Gauss::integral () const { return 1 ; }
// ============================================================================
//  get CDF 
// ============================================================================
double Ostap::Math::Gauss::cdf ( const double x ) const 
{
  const double dx =  s_SQRT2i * ( x - m_peak ) / m_sigma ;
  return 0.5 * ( 1 + std::erf ( dx ) ) ;
}
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::Gauss::integral ( const double low  ,
                                      const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }                         // RETURN
  //
  const double c = s_SQRT2i / m_sigma ;
  const double l = c * ( low  - m_peak ) ;
  const double h = c * ( high - m_peak ) ;
  //
  return 0.5 * ( std::erf ( h ) - std::erf ( l ) ) ;
}
// ============================================================================
bool Ostap::Math::Gauss::setSigma ( const double value )
{
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( m_sigma , value_ ) ) { return false ; }
  m_sigma = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Gauss::setPeak ( const double value )
{
  if ( s_equal ( m_peak , value ) ) { return false ; }
  m_peak   = value  ;
  return true ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Gauss::tag () const 
{ 
  static const std::string s_name = "Gauss" ;
  return Ostap::Utils::hash_combiner ( s_name , m_peak , m_sigma ) ; 
}
// ============================================================================





// ============================================================================
/*  constructor from all agruments 
 *  @param mu     location/peak posiiton 
 *  @param alpha  "scale" parameter 
 *  @param beta   "shape" parameter 
 */
// ============================================================================
Ostap::Math::GenGaussV1::GenGaussV1 
( const double mu    ,
  const double alpha , 
  const double beta  ) 
  : m_mu     ( mu                 ) 
  , m_alpha  ( std::abs ( alpha ) ) 
  , m_beta   ( std::abs ( beta  ) ) 
  , m_gbeta1 ( 0 )
  , m_gbeta2 ( 0 )
{
  //
  setBeta ( beta ) ;
  //
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::GenGaussV1::~GenGaussV1(){}
// ============================================================================
bool Ostap::Math::GenGaussV1::setMu        ( const double value ) 
{
  if ( s_equal ( value , m_mu) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGaussV1::setAlpha     ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_alpha ) ) { return false ; }
  m_alpha = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGaussV1::setBeta     ( const double value ) 
{
  //
  const double value_ = std::max ( std::abs ( value ) , 1.5/GSL_SF_GAMMA_XMAX ) ;
  //
  if ( s_equal ( value_ , m_beta ) ) { return false ; }
  //
  m_beta = value_ ;
  //
  if   ( beta() * GSL_SF_GAMMA_XMAX < 6 ) 
  { 
    m_gbeta1  = 0 ;
    m_gbeta2  = std::lgamma ( 3 / beta() ) ;    
    m_gbeta2 -= std::lgamma ( 1 / beta() ) ;
    m_gbeta2  = my_exp      ( m_gbeta2   ) ;
  }
  else 
  { 
    m_gbeta1 = 1. / std::tgamma ( 1 / beta() )            ;
    m_gbeta2 =      std::tgamma ( 3 / beta() ) * m_gbeta1 ;
  }
  //
  return true ;
}
// ============================================================================
double Ostap::Math::GenGaussV1::pdf ( const double x ) const 
{
  //
  const double delta  = std::abs ( x - m_mu )         ;
  const double delta1 =            delta  / m_alpha   ;
  const double delta2 = std::pow ( delta1 , m_beta  ) ;
  //
  if ( delta2 > 60 || 0 == m_gbeta1 || beta() * GSL_SF_GAMMA_XMAX < 4 ) 
  {
    double result  = std::log    ( 0.5 * beta() / alpha() ) ;
    result        -= delta2 ;
    result        -= std::lgamma ( 1            / beta()  ) ;
    return my_exp ( result ) ;
  }
  //
  double result   = 0.5 * beta() / alpha() ;
  result         *= my_exp   ( -delta2 ) ;
  result         *= m_gbeta1  ;
  //
  return result ;
}
// ============================================================================
double Ostap::Math::GenGaussV1::cdf ( const double x ) const 
{
  //
  const double delta  = std::abs ( x - m_mu )         ;
  const double delta1 =            delta  / m_alpha   ;
  const double delta2 = std::pow ( delta1 , m_beta  ) ;
  //
  const double c = 0.5 * gsl_sf_gamma_inc_P ( 1/beta() , delta2 ) ;
  //
  return x < m_mu ?  0.5 - c : 0.5 + c ;
}
// ============================================================================ 
double Ostap::Math::GenGaussV1::integral ( const double low  , 
                                           const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================ 
double  Ostap::Math::GenGaussV1::variance () const 
{ return alpha() * alpha() *             m_gbeta2   ; }
// ============================================================================
double  Ostap::Math::GenGaussV1::sigma    () const 
{ return alpha()           * std::sqrt ( m_gbeta2 ) ; }
// ============================================================================
double  Ostap::Math::GenGaussV1::kurtosis () const 
{
  double result  =   std::lgamma ( 5 / beta() ) ;
  result        +=   std::lgamma ( 1 / beta() ) ;
  result        -= 2*std::lgamma ( 3 / beta() ) ;
  //
  return my_exp ( result ) - 3 ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::GenGaussV1::tag () const 
{ 
  static const std::string s_name = "GenGaussV1" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_alpha , m_beta ) ; 
}
// ============================================================================



// ============================================================================
/* constructor from all agruments 
 *  @param xi     location/peak posiiton 
 *  @param alpha  "scale" parameter 
 *  @param kappa  "shape" parameter 
 */
// ============================================================================
Ostap::Math::GenGaussV2::GenGaussV2 
( const double xi    ,
  const double alpha , 
  const double kappa )  // kappa=0 correponds to gaussian  
  : m_xi      ( xi                 ) 
  , m_alpha   ( std::abs ( alpha ) ) 
  , m_kappa   (            kappa   ) 
{
  //
  setKappa ( kappa ) ;
  //
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::GenGaussV2::~GenGaussV2(){}
// ============================================================================
bool Ostap::Math::GenGaussV2::setXi        ( const double value ) 
{
  if ( s_equal ( value , m_xi) ) { return false ; }
  m_xi = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGaussV2::setAlpha     ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_alpha ) ) { return false ; }
  m_alpha = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenGaussV2::setKappa  ( const double value ) 
{
  //
  double value_ = value ;
  //
  if ( s_equal ( value_ , 0       ) ) { value_ = 0 ; }
  if ( s_equal ( value_ , m_kappa ) ) { return false ; }
  //
  m_kappa = value_;
  //
  return true ;
}
// ============================================================================
double  Ostap::Math::GenGaussV2::y ( const double x ) const
{
  if  ( s_equal( m_kappa , 0 ) ) { return ( x - xi() ) / alpha() ; }
  //
  const double delta =  - ( x - xi () ) * kappa() / alpha () ;
  //
  return 
    delta > 1 ?
    -std::log   ( 1 + delta ) / kappa() :
    -std::log1p (     delta ) / kappa() ;  
}
// ============================================================================
double Ostap::Math::GenGaussV2::pdf ( const double x ) const 
{
  //
  if      ( s_equal ( m_kappa , 0  ) ) {}
  else if ( m_kappa * x >= m_kappa * m_xi + m_alpha ) { return 0 ; } // cover both cases(?)
  //
  const double y_   = y ( x ) ;
  //
  const double gau  = my_exp ( -0.5 * y_ * y_ ) / s_SQRT2PI ;
  //
  return gau / ( alpha() - kappa() * ( x - xi () ) ) ;
}
// ============================================================================
double Ostap::Math::GenGaussV2::cdf ( const double x ) const 
{
  //
  if      ( s_equal ( m_kappa , 0 ) ) {}
  else if ( kappa () > 0 && ( m_kappa * x >= m_kappa * m_xi + m_alpha ) ) { return 1 ; }
  else if ( kappa () < 0 && ( m_kappa * x >= m_kappa * m_xi + m_alpha ) ) { return 0 ; }
  //
  const double y_   = y ( x ) ;
  //
  const double e_   = std::erf ( y_ * s_SQRT2i ) ;
  //
  return 0.5 * ( 1 + e_ ) ;
}
// ============================================================================
double Ostap::Math::GenGaussV2::integral ( const double low  , 
                                           const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
double Ostap::Math::GenGaussV2::mean () const 
{
  if ( s_equal ( kappa () , 0 ) ) { return xi () ; }
  //
  const double k2 = 0.5 * kappa() * kappa() ;
  //
  return xi() - 0.5 * alpha() * kappa() * exprel ( k2 ) ;
}
// ============================================================================
double Ostap::Math::GenGaussV2::variance () const 
{
  if ( s_equal ( kappa() , 0 ) ) { return alpha () * alpha() ; }
  //
  const double k2 = kappa() * kappa() ;
  //
  return alpha () * alpha() * std::exp ( k2 ) * exprel ( k2 ) ;
}
// ============================================================================
double Ostap::Math::GenGaussV2::sigma () const 
{ return std::sqrt ( variance() ) ; }
// ============================================================================
double Ostap::Math::GenGaussV2::skewness () const 
{
  const double k2     = kappa () * kappa() ;
  //
  const double a1     = exprel (     k2 ) ;
  const double a3     = exprel ( 3 * k2 ) ;
  //
  const double a      = std::pow ( a1 , 1.5 ) ;
  //
  const double result = 3 *  ( a1 - a3 ) / a ;
  //
  return kappa() * result ;
  //  
}
// ============================================================================
double Ostap::Math::GenGaussV2::kurtosis () const 
{
  //
  const double ek2 = my_exp ( kappa() * kappa() ) ;
  //
  return  Ostap::Math::POW ( ek2 , 4 )  
    + 2 * Ostap::Math::POW ( ek2 , 3 ) 
    + 3 * Ostap::Math::POW ( ek2 , 2 ) - 6 ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::GenGaussV2::tag () const 
{ 
  static const std::string s_name = "GenGaussV2" ;
  return Ostap::Utils::hash_combiner ( s_name , m_xi , m_alpha , m_kappa ) ; 
}
// ============================================================================



// ============================================================================
/*  constructor from all agruments 
 *  @param xi     location/peak posiiton 
 *  @param omega  "scale" parameter 
 *  @param alpha  "shape" parameter 
 */
// ============================================================================
Ostap::Math::SkewGauss::SkewGauss
( const double xi    ,
  const double omega , 
  const double alpha )  // alpha=0 correponds to gaussian  
  : m_xi      ( xi                 ) 
  , m_omega   ( std::abs ( omega ) ) 
  , m_alpha   (            alpha   ) 
{}
// ============================================================================
// desctructor
// ============================================================================
Ostap::Math::SkewGauss::~SkewGauss(){}
// ============================================================================
bool Ostap::Math::SkewGauss::setXi        ( const double value ) 
{
  if ( s_equal ( value , m_xi) ) { return false ; }
  m_xi = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::SkewGauss::setOmega     ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_omega ) ) { return false ; }
  m_omega = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::SkewGauss::setAlpha  ( const double value ) 
{
  //
  if ( s_equal ( value , m_alpha ) ) { return false ; }
  m_alpha = value ;
  if ( s_equal ( 0     , m_alpha ) ) { m_alpha = 0  ; }
  //
  return true ;
}
// ============================================================================
double Ostap::Math::SkewGauss::pdf ( const double x ) const 
{
  const double y = ( x - m_xi ) / m_omega ;
  return 2* gauss_pdf ( y ) * gauss_cdf ( m_alpha * y ) / m_omega ;
}
// ============================================================================
double Ostap::Math::SkewGauss::cdf ( const double x ) const 
{
  const double y = ( x - m_xi ) / m_omega ;
  return gauss_cdf ( y ) - 2 * owen ( y , m_alpha ) ;
}
// ============================================================================
double Ostap::Math::SkewGauss::integral ( const double low  , 
                                          const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
double Ostap::Math::SkewGauss::mean () const 
{
  static const double s_c = std::sqrt ( 2.0L / M_PI ) ;
  const double delta = m_alpha / std::sqrt ( 1 + m_alpha * m_alpha ) ;
  return m_xi + m_omega * delta * s_c ;
}
// ============================================================================
double Ostap::Math::SkewGauss::variance () const 
{
  const double delta = m_alpha / std::sqrt( 1 + m_alpha * m_alpha ) ;
  return m_omega * m_omega * ( 1 - 2 * delta * delta , M_PI ) ;
}
// ============================================================================
double Ostap::Math::SkewGauss::skewness () const 
{
  static const double s_c1 =  ( 4.0- M_PI ) / 2.0 ;
  static const double s_c2 =  std::sqrt ( 2.0 / M_PI ) ;
  //
  const double delta = m_alpha / std::sqrt( 1 + m_alpha * m_alpha ) ;
  return s_c1 * std::pow ( delta * s_c2 , 3 ) /
    std::pow ( 1 - 2 * delta * delta / M_PI , 1.5  ) ;
}
// ============================================================================
double Ostap::Math::SkewGauss::kurtosis () const 
{
  static const double s_c1 =  2.0 * ( M_PI - 3 ) ;
  static const double s_c2 =  std::sqrt ( 2.0 / M_PI ) ;
  //
  const double delta = m_alpha / std::sqrt ( 1 + m_alpha * m_alpha ) ;
  return s_c1 * std::pow ( delta * s_c2 , 4 ) / std::pow ( 1 - 2 * delta * delta / M_PI , 2 ) ;
}
// ============================================================================
double Ostap::Math::SkewGauss::sigma  () const 
{ return std::sqrt ( variance () ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::SkewGauss::tag () const 
{ 
  static const std::string s_name = "SkewGauss" ;
  return Ostap::Utils::hash_combiner ( s_name , m_xi , m_omega , m_alpha ) ; 
}
// ============================================================================


// ============================================================================
// constructor from all parameters 
// ============================================================================
Ostap::Math::ExGauss::ExGauss
( const double mu       , 
  const double varsigma , 
  const double k        ) 
  : m_mu ( mu ) 
  , m_varsigma ( std::abs ( varsigma ) )
  , m_k  ( k  )
{}
// ============================================================================
double Ostap::Math::ExGauss::evaluate          ( const double x ) const 
{
  //
  const double z     = ( x - m_mu ) / m_varsigma ;
  const bool k_zero  = s_zero ( m_k ) ;
  //
  const double kk    = std::abs ( m_k ) ;
  //
  return 
    k_zero  ? Ostap::Math::gauss_pdf   ( z ) / m_varsigma : 
    m_k > 0 ? Ostap::Math::gauss_mills ( z , 1.0/kk - z ) / ( kk * m_varsigma ) :
    m_k < 0 ? Ostap::Math::gauss_mills ( z , 1.0/kk + z ) / ( kk * m_varsigma ) :
    Ostap::Math::gauss_pdf   ( z ) / m_varsigma  ;
}
// ============================================================================
bool Ostap::Math::ExGauss::setMu        ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::ExGauss::setVarsigma ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_varsigma , avalue ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ============================================================================
bool Ostap::Math::ExGauss::setK ( const double value ) 
{
  if ( s_equal ( m_k , value ) ) { return false ; }
  m_k = s_zero ( value ) ? 0.0 : value ;
  return true ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::ExGauss::integral   () const { return 1 ; }
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::ExGauss::integral   
( const double low  ,
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low )        { return -integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get CDF
// ============================================================================
double Ostap::Math::ExGauss::cdf ( const double x ) const 
{
  const double z     = ( x - m_mu ) / m_varsigma ;
  const bool k_zero  = s_zero ( m_k ) ;
  //
  const double gauss = Ostap::Math::gauss_cdf ( z ) ;
  const double kk    = std::abs ( m_k ) ;
  //
  return 
    k_zero  ? gauss  :
    m_k > 0 ? gauss - Ostap::Math::gauss_mills ( z , 1 / kk - z ) :
    m_k < 0 ? gauss + Ostap::Math::gauss_mills ( z , 1 / kk + z ) :
    gauss ;
}
// ============================================================================
// mean value 
// ============================================================================
double Ostap::Math::ExGauss::mean        () const 
{ return m_mu + m_k * m_varsigma ; }
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::ExGauss::variance    () const 
{ return m_varsigma * m_varsigma * ( 1 + m_k * m_k ) ; }
// ============================================================================
// RMS value 
// ============================================================================
double Ostap::Math::ExGauss::rms        () const
{ return std::sqrt ( variance () )  ; }
// ============================================================================
// skewness 
// ============================================================================
double Ostap::Math::ExGauss::skewness    () const 
{ return cumulant ( 3 ) / std::pow ( cumulant ( 2 ) , 1.5 ) ; }
// ============================================================================
// (excess) kurtosis  
// ============================================================================
double Ostap::Math::ExGauss::kurtosis    () const 
{
  const double k4 = cumulant ( 4 ) ;
  const double k2 = cumulant ( 2 ) ;
  const double s2 = variance ()  ;
  //
  return ( k4 + 3 * k2 * k2 ) / ( s2 * s2 ) - 3 ;
}
// ============================================================================
// get cumulant 
// ============================================================================
double Ostap::Math::ExGauss::cumulant    ( const unsigned short r ) const 
{
  return 
    0 == r ? 0.0         :
    1 == r ? mean     () :
    2 == r ? variance () :
    s_zero ( m_k ) ? 0.0 :
    std::tgamma ( r ) * std::pow ( m_k * m_varsigma , r ) ;  
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::ExGauss::tag () const 
{ 
  static const std::string s_name = "ExGauss" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_varsigma , m_k  ) ; 
}
// ============================================================================

// ============================================================================
// constructor 
// ============================================================================
Ostap::Math::NormalLaplace::NormalLaplace 
( const double mu       ,
  const double varsigma ,
  const double kL       , 
  const double kR       ) 
  : m_mu ( mu ) 
  , m_varsigma ( std::abs ( varsigma ) ) 
  , m_kL       ( std::abs ( kL       ) ) 
  , m_kR       ( std::abs ( kR       ) ) 
{}
// ============================================================================
bool Ostap::Math::NormalLaplace::setMu        ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::NormalLaplace::setVarsigma ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_varsigma , avalue ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ============================================================================
bool Ostap::Math::NormalLaplace::setKL ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_kL , avalue ) ) { return false ; }
  m_kL = s_zero ( avalue ) ? 0.0 : avalue ;
  return true ;
}
// ============================================================================
bool Ostap::Math::NormalLaplace::setKR ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_kR , avalue ) ) { return false ; }
  m_kR = s_zero ( avalue ) ? 0.0 : avalue ;
  return true ;
}
// ============================================================================
double Ostap::Math::NormalLaplace::evaluate ( const double x ) const 
{
  //
  const double z = ( x - m_mu ) / m_varsigma ;
  const bool l_zero = s_zero ( m_kL ) ;
  const bool r_zero = s_zero ( m_kR ) ;
  //
  return 
    l_zero && r_zero ? Ostap::Math::gauss_pdf   ( z ) / m_varsigma :
    l_zero           ? Ostap::Math::gauss_mills ( z , 1/m_kR - z ) / ( m_kR * m_varsigma ) :
    r_zero           ? Ostap::Math::gauss_mills ( z , 1/m_kL + z ) / ( m_kL * m_varsigma ) :
    ( Ostap::Math::gauss_mills ( z , 1/m_kR - z ) +
      Ostap::Math::gauss_mills ( z , 1/m_kL + z ) ) / ( ( m_kL + m_kR ) * m_varsigma ) ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::NormalLaplace::integral   () const { return 1.0 ; } 
// ============================================================================
// get the integral between low and high limits
// ============================================================================
double Ostap::Math::NormalLaplace::integral  
( const double low  ,
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low )        { return -integral ( high , low ) ; }
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get CDF
// ============================================================================
double Ostap::Math::NormalLaplace::cdf ( const double x ) const 
{
  const double z = ( x - m_mu ) / m_varsigma ;
  const bool l_zero = s_zero ( m_kL ) ;
  const bool r_zero = s_zero ( m_kR ) ;
  //
  const double gauss = Ostap::Math::gauss_cdf ( z ) ;
  //
  return 
    l_zero && r_zero ? gauss : 
    l_zero           ? gauss - Ostap::Math::gauss_mills ( z , 1 / m_kR - z ) :
    r_zero           ? gauss + Ostap::Math::gauss_mills ( z , 1 / m_kL + z ) :
    gauss  - ( Ostap::Math::gauss_mills ( z , 1 / m_kR - z ) * m_kR - 
               Ostap::Math::gauss_mills ( z , 1 / m_kL + z ) * m_kL ) / ( m_kL + m_kR ) ;
}
// ============================================================================
// get cumulant 
// ============================================================================
double Ostap::Math::NormalLaplace::cumulant    ( const unsigned short r ) const 
{
  return 
    0 == r ? 0.0         :
    1 == r ? mean     () :
    2 == r ? variance () :
    std::tgamma ( r ) * ( std::pow ( m_kR * m_varsigma , r )  +
                          std::pow ( m_kL * m_varsigma , r ) ) ; 
}
// ============================================================================
// mean value 
// ============================================================================
double Ostap::Math::NormalLaplace::mean        () const 
{  return m_mu + m_varsigma * ( m_kR - m_kL ) ; }
// ============================================================================
// variance  
// ============================================================================
double Ostap::Math::NormalLaplace::variance    () const 
{  return m_varsigma * m_varsigma * ( 1 + m_kR * m_kR  + m_kL * m_kL ) ; }
// ============================================================================
// RMS
// ============================================================================
double Ostap::Math::NormalLaplace::rms        () const 
{ return std::sqrt ( variance () )  ; }
// ============================================================================
// skewness 
// ============================================================================
double Ostap::Math::NormalLaplace::skewness    () const 
{ return cumulant ( 3 ) / std::pow ( cumulant ( 2 ) , 1.5 ) ; }
// ============================================================================
// (excess) kurtosis  
// ============================================================================
double Ostap::Math::NormalLaplace::kurtosis    () const 
{
  const double k4 = cumulant ( 4 ) ;
  const double k2 = cumulant ( 2 ) ;
  const double s2 = variance ()  ;
  //
  return ( k4 + 3 * k2 * k2 ) / ( s2 * s2 ) - 3 ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::NormalLaplace::tag () const 
{ 
  static const std::string s_name = "NormalLaplace" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_varsigma , m_kL , m_kR ) ; 
}
// ============================================================================

  
  
 
  
  
    


// ============================================================================
// Bukin
// ============================================================================
/*  constructor from all parameters
 *  @param peak  the peak position
 *  @param sigma the effective sigma, defined as FWHM/2.35
 *  @param xi    the asymmetry parameter
 *  @param rhoL  the left  tail paramter
 *  @param rhoR  the right tail paramter
 */
// ============================================================================
Ostap::Math::Bukin::Bukin
( const double peak   ,
  const double sigma  ,
  const double xi     ,
  const double rho_L  ,
  const double rho_R  )
//
  : m_peak      ( M_PI + peak  )
  , m_sigma     ( M_PI + sigma )
  , m_xi        ( M_PI + xi    )
  , m_rho_L     ( M_PI + rho_L )
  , m_rho_R     ( M_PI + rho_R )
//
  , m_x1        ( M_PI )
  , m_x2        ( M_PI )
//
  , m_workspace ()
{
  //
  setXi    ( xi    ) ; // must be the first
  //
  setPeak  ( peak  ) ;
  //
  setSigma ( sigma ) ;
  //
  setRhoL  ( rho_L ) ;
  //
  setRhoR  ( rho_R ) ;
  //
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Bukin::~Bukin () {}
// ============================================================================
bool Ostap::Math::Bukin::setPeak  ( const double value )
{
  if ( s_equal ( value , m_peak ) ) { return false ; }
  //
  m_peak   = value ;
  //
  const double xi_ = m_xi / std::sqrt ( 1 + m_xi * m_xi ) ;
  m_x1     = m_peak + m_sigma * s_Bukin * ( xi_ - 1 ) ;
  m_x2     = m_peak + m_sigma * s_Bukin * ( xi_ + 1 ) ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Bukin::setSigma ( const double value )
{
  //
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  //
  m_sigma  = value_ ;
  //
  const double xi_ = m_xi/std::sqrt ( 1 + m_xi * m_xi ) ;
  m_x1 = m_peak + m_sigma * s_Bukin * ( xi_ - 1 ) ;
  m_x2 = m_peak + m_sigma * s_Bukin * ( xi_ + 1 ) ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Bukin::setXi    ( const double value )
{
  // no need for update
  if ( s_equal ( value , m_xi ) ) { return false ; } // no need for update
  //
  m_xi     = value ;
  //
  const double xi      = m_xi    ;
  const double xi2     =   xi*xi ;
  const double xi2sqrt = std::sqrt ( 1 + xi2 ) ;
  //
  const double alpha = 2 * xi  * xi2sqrt / s_Bukin   ;
  const double beta  = 2 * xi  * ( xi - xi2sqrt )    ;
  // well, it is actually alpha/beta:
  const double ab    = xi2sqrt / ( xi - xi2sqrt ) / s_Bukin ;
  //
  m_A     = alpha             ;
  //
  m_B2    = 1/Ostap::Math::log1p_x ( beta )  ;
  m_B2   *= m_B2              ;
  m_B2   *= ab*ab             ;
  //
  //
  const double delta = xi + xi2sqrt - 1 ;
  const double tail  =
    0.5 * s_Bukin * xi2sqrt * ( 1 + xi + xi2sqrt ) / ( xi + xi2sqrt ) / Ostap::Math::log1p_x  ( delta ) ;
  //
  // left  tail parameter
  //
  m_L  = tail ;
  m_L /= ( xi2sqrt - xi ) ;
  m_L /= ( xi2sqrt - xi ) ;
  //
  // right tail parameter
  //
  m_R  = tail ;
  m_R /= ( xi2sqrt + xi ) ;
  m_R /= ( xi2sqrt + xi ) ;
  //
  // central region
  //
  const double xi_ = m_xi / xi2sqrt ;
  m_x1 = m_peak + m_sigma * s_Bukin * ( xi_ - 1 ) ;
  m_x2 = m_peak + m_sigma * s_Bukin * ( xi_ + 1 ) ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Bukin::setRhoL  ( const double value )
{
  //
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( value_ , m_rho_L  ) ) { return false ; }
  //
  m_rho_L    = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Bukin::setRhoR  ( const double value )
{
  //
  const double value_ = std::fabs ( value ) ;
  if ( s_equal ( value_ , m_rho_R  ) ) { return false ; }
  //
  m_rho_R    = value_ ;
  //
  return true ;
}
// ============================================================================
// evaluate Bukin's function
// ============================================================================
double Ostap::Math::Bukin::pdf ( const double x ) const
{
  //
  //  left tail :
  //
  if       ( m_x1 >= x )
  {
    const double dx  = x - m_x1               ;
    const double dx2 = dx / ( m_peak - m_x1 ) ;
    return  0.5 * my_exp (   m_L * dx / m_sigma  - m_rho_L * m_rho_L * dx2 * dx2 ) ;
  }
  //
  // right tail :
  //
  else if ( m_x2 <=  x )
  {
    const double dx  = x - m_x2               ;
    const double dx2 = dx / ( m_peak - m_x2 ) ;
    return 0.5 * my_exp ( - m_R * dx / m_sigma  - m_rho_R * m_rho_R * dx2 * dx2 ) ;
  }
  //
  // central region
  //
  const double dx   = ( x - m_peak ) / m_sigma ;
  const double A    = Ostap::Math::log1p_x ( m_A * dx ) ;
  //
  return my_exp ( - s_ln2 * dx * dx * A * A * m_B2 ) ;
}
// =========================================================================
// get the integral between low and high limits
// =========================================================================
double Ostap::Math::Bukin::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0        ; } // RETURN
  else if (           low > high   ) { return - integral ( high , low  ) ; } // RETURN
  //
  // split into reasonable sub-intervals
  //
  if ( low < m_x1    && m_x1   < high )
  { return integral (  low , m_x1   ) + integral ( m_x1   , high ) ; }
  if ( low < m_x2    && m_x2   < high )
  { return integral (  low , m_x2   ) + integral ( m_x2   , high ) ; }
  if ( low < m_peak  && m_peak < high )
  { return integral (  low , m_peak ) + integral ( m_peak , high ) ; }
  //
  const bool in_tail = 
    ( high < m_x1 - 5 * std::abs ( m_x2 - m_x1 ) )  || 
    ( low  > m_x2 + 5 * std::abs ( m_x2 - m_x1 ) ) ; 
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Bukin> s_integrator {} ;
  static char s_message[] = "Integral(Bukin)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () , 
      &F     ,  
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size () ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
}
// =========================================================================
// get the integral
// =========================================================================
double Ostap::Math::Bukin::integral () const
{
  //
  // Tails
  //
  static const Ostap::Math::GSL::Integrator1D<Bukin> s_integrator {} ;
  static char s_message1[] = "Integral(Bukin/left)"  ;
  static char s_message2[] = "Integral(Bukin/right)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  //
  int    ierror1  =  0 ;
  double result1  =  1 ;
  double error1   = -1 ;
  std::tie ( ierror1 , result1 , error1 ) = s_integrator.gaqil_integrate
    ( tag () , 
      &F     , 
      m_x1   ,                       // low edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION_TAIL    ,          // relative precision
      m_workspace.size()   ,          // size of workspace
      s_message1          , 
      __FILE__ , __LINE__ ) ;
  //
  int    ierror2  =  0 ;
  double result2  =  1 ;
  double error2   = -1 ;
  std::tie ( ierror2 , result2 , error2 ) = s_integrator.gaqiu_integrate
    ( tag () , 
      &F     , 
      m_x2   ,                       // high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION_TAIL    ,          // relative precision
      m_workspace.size()   ,          // size of workspace
      s_message2          , 
      __FILE__ , __LINE__ ) ;
  //
  return result1 + result2 + integral ( m_x1 , m_x2 ) ;
  //
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Bukin::tag () const 
{ 
  static const std::string s_name = "Bukin" ;
  return Ostap::Utils::hash_combiner ( s_name , m_peak , m_sigma , m_xi , m_rho_L , m_rho_R ) ; 
}
// ============================================================================




// ============================================================================
// Novosibirsk function
// ============================================================================
/*  constructor from all parameters
 *  @param m0    the peak posiion
 *  @param sigma the effective sigma
 *  @param tau   the tail paramter
 */
// ============================================================================
Ostap::Math::Novosibirsk::Novosibirsk
( const double m0    ,
  const double sigma ,
  const double tau   )
  : m_m0        ( m0                   )
  , m_sigma     ( std::fabs ( sigma )  )
  , m_tau       (             tau      )
    //
  , m_workspace ()
{
  //
  m_lambda = x_sinh ( m_tau * s_Novosibirsk ) ;
  //
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Novosibirsk::~Novosibirsk() {}
// ============================================================================
// set parameter m0
// ============================================================================
bool Ostap::Math::Novosibirsk::setM0    ( const double value )
{
  if ( s_equal ( m_m0 ,  value ) ) { return false ; }
  m_m0 = value ;
  return true ;
}
// ============================================================================
// set parameter sigma
// ============================================================================
bool Ostap::Math::Novosibirsk::setSigma ( const double value )
{
  const double avalue = std::fabs ( value ) ;
  if ( s_equal ( avalue , m_sigma ) ) { return false ; }
  m_sigma    = value ;
  return true ;
}
// ============================================================================
// set parameter tau
// ============================================================================
bool Ostap::Math::Novosibirsk::setTau ( const double value )
{
  if ( s_equal ( value , m_tau ) ) { return false ; }
  m_tau      = value ;
  m_lambda   = x_sinh ( m_tau * s_Novosibirsk ) ;
  return true ;
}
// ============================================================================
// evaluate Novosibirsk's function
// ============================================================================
double Ostap::Math::Novosibirsk::pdf  ( const double x ) const
{
  const double dx     = ( x - m_m0 ) / m_sigma ;
  const double arg    = m_lambda * dx * m_tau ;
  if ( arg <= -1 || s_equal ( arg , -1 ) ) { return 0 ; } // RETURN
  const double l      = Ostap::Math::log1p_x  ( arg ) * m_lambda * dx ;
  const double result = l * l ; // + m_tau * m_tau ;
  //
  return  my_exp ( -0.5 * result ) * s_SQRT2PIi / m_sigma ;
}
// =========================================================================
// get the integral between low and high limits
// =========================================================================
double Ostap::Math::Novosibirsk::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                  0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high ,
                                                           low   ) ; } // RETURN
  //
  // split into reasonable sub intervals
  //
  if ( low <  m_m0  && m_m0 < high ) { return integral ( low , m_m0 ) + integral ( m_m0 , high ) ; }
  //
  {
    const double x1 = m_m0 +  3 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_m0 -  3 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  {
    const double x1 = m_m0 +  5 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_m0 -  5 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }  
  //
  {
    const double x1 = m_m0 + 10 * m_sigma ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_m0 - 10 * m_sigma ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  {
    const double x1 = m_m0 + 15 * m_sigma ;
    if ( 0 < m_tau && low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_m0 - 15 * m_sigma ;
    if ( 0 > m_tau && low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
  }
  //
  const double x1     = m_m0 - 15 * m_sigma  ;
  const double x2     = m_m0 + 15 * m_sigma  ;
  const double x_low  = std::min ( x1 , x2 ) ;
  const double x_high = std::max ( x1 , x2 ) ;
  //
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Novosibirsk> s_integrator {} ;
  static char s_message[] = "Integral(Novosibirsk)" ;
  //
  const bool in_tail = high <= x_low || x_high <= low ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () , 
      &F     , 
      low    , high             ,    // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()              ,           // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// =========================================================================
// get the integral
// =========================================================================
double Ostap::Math::Novosibirsk::integral () const
{
  //
  if ( s_zero ( m_tau ) ) { return 1 ; }
  //
  const double tau1 = std::max ( 1.0 , std::abs ( m_tau ) ) ;
  const double tau2 = 1 ;
  //
  const double x_low  = m_m0 - ( 0 <= m_tau ?  5 * tau2 : 15 * tau1 ) * m_sigma ;
  const double x_high = m_m0 + ( 0 <= m_tau ? 15 * tau1 :  5 * tau2 ) * m_sigma ;
  //
  // use GSL to evaluate the tails:
  //
  static const Ostap::Math::GSL::Integrator1D<Novosibirsk> s_integrator {} ;
  static char s_message1[] = "Integral(Novosibirsk/left)"  ;
  static char s_message2[] = "Integral(Novosibirs/right)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  //
  int    ierror1  =  0 ;
  double result1  =  1 ;
  double error1   = -1 ;
  std::tie ( ierror1 , result1 , error1 ) = s_integrator.gaqil_integrate
    ( tag () , 
      &F     , 
      x_low  ,                       // low edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION_TAIL   ,          // absolute precision
      s_APRECISION_TAIL   ,          // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message1          , 
      __FILE__ , __LINE__ ) ;
  //
  int    ierror2  =  0 ;
  double result2  =  1 ;
  double error2   = -1 ;
  std::tie ( ierror2 , result2 , error2 ) = s_integrator.gaqiu_integrate
    ( tag () , 
      &F     , 
      x_high ,                       // high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION_TAIL    ,          // absolute precision
      s_RPRECISION_TAIL    ,          // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message2          , 
      __FILE__ , __LINE__ ) ;
  //
  return result1 + result2 + integral ( x_low ,  x_high ) ;
  //  
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Novosibirsk::tag () const 
{ 
  static const std::string s_name = "Novosibirsk" ;
  return Ostap::Utils::hash_combiner ( s_name , m_m0 , m_sigma , m_tau ) ; 
}
// ============================================================================


// ============================================================================
// Crystal Ball & Co
// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 */
// ============================================================================
Ostap::Math::CrystalBall::CrystalBall
( const double m0     ,
  const double sigma  ,
  const double alpha  ,
  const double n      )
  : m_m0         ( m0 )
  , m_sigma      (  1 )
  , m_alpha      (  2 )
  , m_n          (  2 )
//
  , m_A  ( -1000 ) 
  , m_B  ( -1000 ) 
  , m_C  ( -1000 ) 
{
  //
  setM0     ( m0     ) ;
  setAlpha  ( alpha  ) ;
  setSigma  ( sigma  ) ;
  setN      ( n      ) ;
  //
  m_A = my_exp ( -0.5 * m_alpha * m_alpha ) ;
  m_B =  0.5 * ( 1 + std::erf ( - m_alpha * s_SQRT2i ) ) ;
  if   ( !s_equal ( m_n , 0 ) && !s_equal ( m_alpha , 0 )  ) 
  { m_C  = ( m_n + 1 )  / aa ()  / m_n  * s_SQRT2PIi ; }
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::CrystalBall::~CrystalBall (){}
// ============================================================================
bool  Ostap::Math::CrystalBall::setM0 ( const double value )
{
  //
  if ( s_equal ( value , m_m0 ) ) { return false ; }
  //
  m_m0       = value ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBall::setSigma ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  //
  m_sigma    = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBall::setAlpha  ( const double value )
{
  //
  if ( s_equal ( value , m_alpha ) ) { return false ; }
  //
  m_alpha    = value  ;
  //
  m_A        = my_exp ( -0.5 * alpha () * alpha ()) ;
  //
  // 
  if   ( s_equal ( n () , 0 ) || s_equal ( m_alpha , 0 )  ) { m_C = -1000 ; }
  else { m_C  = np1 () / aa () / n()  * s_SQRT2PIi ; }
  //
  m_B  = 0.5 * ( 1 + std::erf ( - m_alpha * s_SQRT2i ) ) ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBall::setN      ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_n     ) ) { return false ; }
  //
  m_n        = value_ ;
  if ( s_equal ( m_n , 0 ) ) { m_n = 0 ; }
  //
  if   ( s_equal ( n () , 0 ) || s_equal ( m_alpha , 0 )  ) { m_C = -1000 ; }
  else { m_C  = np1()  / aa () / n() * s_SQRT2PIi ; }
  //
  return true ;
}
// ============================================================================
//  evaluate CrystalBall's function
// ============================================================================
double Ostap::Math::CrystalBall::pdf ( const double x ) const
{
  //
  const double dx    = ( x - m_m0 ) / m_sigma ;
  //
  // the tail
  //
  if  ( dx < -m_alpha )
  {
    const double frac = np1 () / ( np1 () - aa() * ( m_alpha + dx ) )  ;
    return std::pow ( frac , np1 () ) * m_A * s_SQRT2PIi / sigma() ;
  }
  //
  // the peak
  //
  return my_exp ( -0.5 * dx * dx ) * s_SQRT2PIi / sigma() ;
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CrystalBall::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high ,
                                                           low  ) ; } // RETURN
  //
  const double x0 = m_m0 - m_alpha * m_sigma ;
  //
  // split into proper subintervals
  //
  if      ( low < x0 && x0 < high )
  { return integral ( low , x0 ) + integral ( x0 , high ) ; }
  //
  // Z = (x-x0)/sigma 
  //
  const double zlow  = ( low  - m_m0 ) / sigma() ;
  const double zhigh = ( high - m_m0 ) / sigma() ;
  //
  // peak
  //
  if ( x0 <= low  )
  { return s_SQRT2PIi * details::gaussian_int ( 0.5   , 0 , zlow  , zhigh ) ; }
  //
  // tail
  //
  const double A =   np1 () ;
  const double B =   np1 () ;
  const double C = - aa  () ;
  //
  const double result = s_SQRT2PIi * m_A * 
    tail_integral ( A , B , C , np1 () , zlow + alpha() , zhigh + alpha() ) ;
  //
  return result ;
}
// =========================================================================
// get the integral
// =========================================================================
double Ostap::Math::CrystalBall::integral () const
{
  /// the regular case 
  if ( 0 < m_C ) { return m_C + m_B ; }
  //
  /// trunkate it! 
  const double left = ( 0 < m_alpha ) ?  (-m_alpha-s_TRUNC) : -s_TRUNC ;
  // 
  return m_B + integral ( m0 () + left     * sigma() , 
                          m0 () - alpha () * sigma() ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::CrystalBall::tag () const 
{
  static const std::string s_name = "CrystalBall" ;
  return Ostap::Utils::hash_combiner ( s_name , m_m0 , m_sigma , m_alpha , m_n ) ; 
}
// ============================================================================


// ============================================================================
// Needham function
// ============================================================================
/* constructor from all parameters
 *  @param m0     m0       parameter
 *  @param sigma  sigma    parameter
 *  @param a0     a0       parameter
 *  @param a1     a1       parameter
 *  @param a2     a2       parameter
 */
// ============================================================================
Ostap::Math::Needham::Needham
( const double m0    ,
  const double sigma ,
  const double a0    ,
  const double a1    ,
  const double a2    )
/// @see Ostap::Math:CrystalBall
  : m_cb  ( m0 , sigma , 1 , 0 ) // Ostap::Math:CrystalBall
  , m_a0  ( std::abs ( a0 )  )
  , m_a1  (            a1    )
  , m_a2  (            a2    )
{
  m_cb.setAlpha ( alpha () ) ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Needham::~Needham(){}
// ============================================================================
bool Ostap::Math::Needham::setA0 ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_a0 ) ) { return false ; }
  m_a0      = value_ ;
  return m_cb.setAlpha ( alpha () ) ;
}
// ============================================================================
bool Ostap::Math::Needham::setA1 ( const double value )
{
  if ( s_equal ( value , m_a1 ) ) { return false ; }
  m_a1 = value ;
  return m_cb.setAlpha ( alpha () ) ;
}
// ============================================================================
bool Ostap::Math::Needham::setA2 ( const double value )
{
  if ( s_equal ( value , m_a2 ) ) { return false ; }
  m_a2 = value ;
  return m_cb.setAlpha ( alpha () ) ;
}
// ===========================================================================
// evaluate Needham's function
// ===========================================================================
double Ostap::Math::Needham::pdf ( const double x ) const
{ return m_cb ( x ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Needham::tag () const 
{
  static const std::string s_name = "Needham" ;
  return Ostap::Utils::hash_combiner ( s_name , m_cb.tag() ,  m_a0 , m_a1 , m_a2 ) ; 
}
// ============================================================================


// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 */
// ============================================================================
Ostap::Math::CrystalBallRightSide::CrystalBallRightSide
( const double m0    ,
  const double sigma ,
  const double alpha ,
  const double n     )
  : m_cb         ( m0 , sigma , alpha , n ) 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::CrystalBallRightSide::~CrystalBallRightSide (){}
// ============================================================================
//  evaluate CrystalBall's function
// ============================================================================
double Ostap::Math::CrystalBallRightSide::pdf ( const double x ) const
{
  const double y = 2 * m0 ()  - x ;
  //
  return  m_cb.pdf ( y ) ;  
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CrystalBallRightSide::integral
( const double low ,
  const double high ) const
{ return m_cb.integral ( 2 * m0 () - high  , 2 * m0 () - low ) ; }
// =========================================================================
// get the integral
// =========================================================================
double Ostap::Math::CrystalBallRightSide::integral () const 
{ return m_cb.integral () ; }
// =========================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::CrystalBallRightSide::tag () const 
{  
  static const std::string s_name = "CrystalBallRightSide" ;
  return Ostap::Utils::hash_combiner ( s_name , m_cb.tag() ,  -1 ) ;
}
// ============================================================================


// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 */
// ============================================================================
Ostap::Math::CrystalBallDoubleSided::CrystalBallDoubleSided
( const double m0      ,
  const double sigma   ,
  const double alpha_L ,
  const double n_L     ,
  const double alpha_R ,
  const double n_R     )
  : m_m0         (  m0 )
  , m_sigma      (   1 )
  , m_alpha_L    (   2 )
  , m_n_L        (   2 )
  , m_alpha_R    (   2 )
  , m_n_R        (   2 )
//
  , m_AL         ( -1000 ) 
  , m_AR         ( -1000 )
  , m_B          ( -1000 ) 
  , m_TL         ( -1000 ) 
  , m_TR         ( -1000 ) 
{
  //
  setM0       ( m0      ) ;
  setSigma    ( sigma   ) ;
  setAlpha_L  ( alpha_L ) ;
  setAlpha_R  ( alpha_R ) ;
  setN_L      ( n_L     ) ;
  setN_R      ( n_R     ) ;
  //
  m_AL = my_exp ( -0.5 * m_alpha_L * m_alpha_L ) ;
  m_AR = my_exp ( -0.5 * m_alpha_R * m_alpha_R ) ;
  m_B  = 0.5 *  ( std::erf (  m_alpha_R * s_SQRT2i ) - 
                  std::erf ( -m_alpha_L * s_SQRT2i ) ) ;
  //
  if   ( !s_equal ( m_n_L , 0 ) && !s_equal ( m_alpha_L , 0 )  ) 
  { m_TL  = ( m_n_L + 1 )  / std::abs ( m_alpha_L )  / m_n_L  * s_SQRT2PIi ; }
  if   ( !s_equal ( m_n_R , 0 ) && !s_equal ( m_alpha_R , 0 )  ) 
  { m_TR  = ( m_n_R + 1 )  / std::abs ( m_alpha_R )  / m_n_R  * s_SQRT2PIi ; }
  //
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::CrystalBallDoubleSided::~CrystalBallDoubleSided(){}
// ============================================================================
bool Ostap::Math::CrystalBallDoubleSided::setM0 ( const double value )
{
  if ( s_equal ( value , m_m0 ) ) { return false ; }
  //
  m_m0       = value ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBallDoubleSided::setSigma ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  //
  m_sigma    = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBallDoubleSided::setAlpha_L ( const double value )
{
  if ( s_equal ( value , m_alpha_L ) ) { return false ; }
  //
  m_alpha_L  = value  ;
  m_AL       = my_exp ( -0.5 * m_alpha_L * m_alpha_L ) ;
  m_B        = 0.5 *  ( std::erf (  m_alpha_R * s_SQRT2i ) - 
                        std::erf ( -m_alpha_L * s_SQRT2i ) ) ;
  //
  if   ( s_equal ( m_n_L , 0 ) || s_equal  ( m_alpha_L , 0 )  ) {  m_TL = -1000 ; }
  else { m_TL  = ( m_n_L + 1 )  / std::abs ( m_alpha_L )  / m_n_L  * s_SQRT2PIi ; }
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBallDoubleSided::setAlpha_R ( const double value )
{
  if ( s_equal ( value , m_alpha_R ) ) { return false ; }
  //
  m_alpha_R  = value  ;
  m_AR       = my_exp ( -0.5 * m_alpha_R * m_alpha_R ) ;
  m_B        = 0.5 *  ( std::erf (  m_alpha_R * s_SQRT2i ) - 
                        std::erf ( -m_alpha_L * s_SQRT2i ) ) ;
  //
  if   ( s_equal ( m_n_R , 0 ) || s_equal ( m_alpha_R , 0 )  ) { m_TR = -1000 ; }
  else { m_TR  = ( m_n_R + 1 )  / std::abs ( m_alpha_R )  / m_n_R  * s_SQRT2PIi ; }
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBallDoubleSided::setN_L     ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_n_L    ) ) { return false ; }
  //
  m_n_L      = value_ ;
  if ( s_equal ( m_n_L , 0 ) ) { m_n_L = 0 ; }
  //
  if   ( s_equal ( m_n_L , 0 ) || s_equal ( m_alpha_L , 0 )  ) {  m_TL = -1000 ; }
  else { m_TL  = ( m_n_L + 1 )  / std::abs ( m_alpha_L )  / m_n_L  * s_SQRT2PIi ; }
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::CrystalBallDoubleSided::setN_R     ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_n_R    ) ) { return false ; }
  //
  m_n_R      = value_ ;
  if ( s_equal ( m_n_R , 0 ) ) { m_n_R = 1 ; }
  //
  if   ( s_equal ( m_n_R , 0 ) || s_equal ( m_alpha_R , 0 )  ) { m_TR = -1000 ; }
  else { m_TR  = ( m_n_R + 1 )  / std::abs ( m_alpha_R )  / m_n_R  * s_SQRT2PIi ; }
  //
  return true ;
}
// ============================================================================
//  evaluate CrystalBall's function
// ============================================================================
double Ostap::Math::CrystalBallDoubleSided::pdf ( const double x ) const
{
  //
  const double dx   = ( x - m_m0 ) / m_sigma ;
  //
  // the left tail
  if      ( dx  < -m_alpha_L )  // left tail
  {
    const double np1  = n_L() + 1 ;
    const double frac = np1 / ( np1 - std::abs ( m_alpha_L ) * ( m_alpha_L + dx ) )  ;
    return std::pow ( frac , np1 ) * m_AL * s_SQRT2PIi / sigma() ;
  }
  // the right tail
  else if  ( dx >  m_alpha_R )  // right tail
  {
    const double np1  = n_R () + 1 ;
    const double frac = np1 / ( np1 - std::abs ( m_alpha_R ) * ( m_alpha_R - dx ) )  ;
    return std::pow ( frac , np1 ) * m_AR * s_SQRT2PIi / sigma() ;
  }
  //
  // the peak
  //
  return my_exp ( -0.5 * dx * dx ) * s_SQRT2PIi / sigma() ; 
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CrystalBallDoubleSided::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return 0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high , low ) ; }
  //
  const double x_low  = m_m0 - m_alpha_L * m_sigma ;
  const double x_high = m_m0 + m_alpha_R * m_sigma ;
  //
  // split into proper subintervals
  //
  if ( low < x_low  && x_low  < high )
  { return integral ( low , x_low  ) + integral ( x_low  , high ) ; }
  if ( low < x_high && x_high < high )
  { return integral ( low , x_high ) + integral ( x_high , high ) ; }
  //
  // Z = (x-x0)/sigma 
  //
  const double zlow  = ( low  - m_m0 ) / sigma() ;
  const double zhigh = ( high - m_m0 ) / sigma() ;
  //
  // the peak
  //
  if ( x_low <= low && high <= x_high )
  { return  s_SQRT2PIi * details::gaussian_int ( 0.5   , 0 , zlow  , zhigh ) ; }
  //
  // left tail 
  //
  if ( high <= x_low ) 
  {
    const double np1 = n_L () + 1 ;
    //
    const double A   = np1 ;
    const double B   = np1 ;
    const double C   = - std::abs ( alpha_L () ) ;
    //
    return s_SQRT2PIi * m_AL * 
      tail_integral ( A , B , C , np1 , zlow + alpha_L() , zhigh + alpha_L() ) ;
  }
  //
  // right tail 
  // 
  if ( low  >= x_high ) 
  {
    //
    const double np1 = n_R () + 1 ;
    //
    const double A   = np1 ;
    const double B   = np1 ;
    const double C   = std::abs ( alpha_R () ) ;
    //
    return s_SQRT2PIi * m_AR * 
      tail_integral ( A , B , C , np1 , zlow - alpha_R() , zhigh - alpha_R() ) ;
  }
  //
  return 0 ;
}
// ============================================================================
// get the (truncated)  integral
// ============================================================================
double Ostap::Math::CrystalBallDoubleSided::integral () const 
{
  //
  if      ( 0 < m_TL && 0 <= m_TR ) { return m_TL + m_TR + m_B ; }
  else if ( 0 < m_TR ) 
  {
    /// truncate it! 
    const double left = ( 0 < alpha_L() ) ?  (-alpha_L()-s_TRUNC) : -s_TRUNC ;
    // 
    return m_TR + m_B + integral ( m0 () + left       * sigma() , 
                                   m0 () - alpha_L () * sigma() ) ;
    
  }
  else if ( 0 < m_TL ) 
  {
    /// truncate it! 
    const double right = ( 0 < alpha_R() ) ?  ( alpha_R () + s_TRUNC) : + s_TRUNC ;
    // 
    return m_TL + m_B + integral ( m0 () + alpha_R () * sigma() , 
                                   m0 () + right      * sigma() ) ;
    
  }
  //
  /// truncate both
  const double left  = ( 0 < alpha_L() ) ?  (-alpha_L () - s_TRUNC ) : -s_TRUNC ;
  /// truncate it! 
  const double right = ( 0 < alpha_R() ) ?  ( alpha_R () + s_TRUNC ) : + s_TRUNC ;
  //
  return integral ( m0 () - left  * sigma () , m0 () + right * sigma () ) ;

}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::CrystalBallDoubleSided::tag () const 
{ 
  static const std::string s_name = "CrystalBallDoubleSide" ;
  return Ostap::Utils::hash_combiner ( s_name    , 
                             m_m0      , m_sigma , 
                             m_alpha_L , m_n_L   , 
                             m_alpha_R , m_n_R   ) ; 
}
// ============================================================================




// ============================================================================
// apolonios 
// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 *  @param b     b-parameter 
 */
// ============================================================================
Ostap::Math::Apollonios::Apollonios
( const double m0    ,
  const double sigma ,
  const double alpha ,
  const double n     ,
  const double bp    )
  : m_m0         ( m0 )
  , m_sigma      (  1 )
  , m_alpha      (  2 )
  , m_n          (  2 )
  , m_b          (  2 )
//
  , m_A  ( -1000 ) 
//
  , m_workspace () 
{
  //
  setM0     ( m0    ) ;
  setAlpha  ( alpha ) ;
  setSigma  ( sigma ) ;
  setN      ( n     ) ;
  setB      ( bp    ) ;
  //
  m_A = my_exp ( -b () * a1 () ) ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Apollonios::~Apollonios(){}
// ============================================================================
bool  Ostap::Math::Apollonios::setM0 ( const double value )
{
  //
  if ( s_equal ( value , m_m0 ) ) { return false ; }
  //
  m_m0       = value ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apollonios::setSigma ( const double value )
{
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  //
  m_sigma    = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apollonios::setAlpha  ( const double value )
{
  //
  if ( s_equal ( value , m_alpha ) ) { return false ; }
  //
  m_alpha    = value  ;
  //
  m_A = my_exp ( -b() * a1 () ) ; 
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apollonios::setN      ( const double value )
{
  const double value_ = std::fabs ( value );
  if ( s_equal ( value_ , m_n     ) ) { return false ; }
  //
  m_n        = value_ ;
  if ( s_equal ( m_n , 0 ) ) { m_n = 0 ; }
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apollonios::setB  ( const double value )
{
  //
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_b ) ) { return false ; }
  //
  m_b    = value_  ;
  //
  if ( s_equal ( m_b , 0 ) ) { m_b = 0 ; }
  if ( s_equal ( m_b , 1 ) ) { m_b = 1 ; }
  //
  m_A = my_exp ( -b () * a1 () ) ;
  //
  return true ;
}
// ============================================================================
//  evaluate Apollonios' function
// ============================================================================
double Ostap::Math::Apollonios::pdf ( const double x ) const
{
  //
  const double dx    = ( x - m_m0 ) / m_sigma ;
  //
  // the tail
  //
  if  ( dx < -m_alpha )
  {
    const double frac = np1 () / ( np1 () - ( m_alpha + dx ) * aa () )  ;
    return std::pow ( frac , np1 () ) * m_A * s_SQRT2PIi / sigma() ;
  }
  //
  // the peak
  //
  return my_exp ( -b() * std::sqrt ( 1 + dx*dx ) ) * s_SQRT2PIi / sigma() ;  
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::Apollonios::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high ,
                                                           low  ) ; } // RETURN
  //
  const double x0 = m_m0 - m_alpha * m_sigma ;
  //
  // split into proper subintervals
  //
  if      ( low < x0 && x0 < high )
  { return integral ( low , x0 ) + integral ( x0 , high ) ; }
  //
  // Z = (x-x0)/sigma 
  //
  const double zlow  = ( low  - m_m0 ) / sigma() ;
  const double zhigh = ( high - m_m0 ) / sigma() ;
  //
  // peak
  //
  if ( x0 <= low  )
  {
    //
    // use GSL to evaluate the integral
    //
    static const Ostap::Math::GSL::Integrator1D<Apollonios> s_integrator {} ;
    static char s_message[] = "Integral(Apollonios)" ;
    //
    const auto F = s_integrator.make_function ( this ) ;
    int    ierror   =  0 ;
    double result   =  1 ;
    double error    = -1 ;
    std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
      ( tag () , 
        &F     , 
        low    , high  ,               // low & high edges
        workspace ( m_workspace ) ,    // workspace
        s_APRECISION         ,          // absolute precision
        s_RPRECISION         ,          // relative precision
        m_workspace.size()              ,          // size of workspace
        s_message           , 
        __FILE__ , __LINE__ ) ;
    //  
    return result ;
  }
  //
  // tail
  //
  const double A = np1 () ;
  const double B = np1 () ;
  const double C = - std::abs ( alpha () * b () ) / a1 ()  ;
  //
  const double result = s_SQRT2PIi * m_A * 
    tail_integral ( A , B , C , np1 () , zlow + alpha() , zhigh + alpha() ) ;
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Apollonios::tag () const 
{ 
  static const std::string s_name = "Apollonios" ;
  return Ostap::Utils::hash_combiner ( s_name , m_m0 , m_sigma , m_alpha , m_n , m_b ) ; 
}
// ============================================================================

// ============================================================================
// apolonios2 
// ============================================================================
/*  constructor from all parameters
 *  @param m0 m0 parameter
 *  @param alpha alpha parameter
 *  @param n     n-parameter
 *  @param b     b-parameter 
 */
// ============================================================================
Ostap::Math::Apollonios2::Apollonios2
( const double m0     ,
  const double sigmaL ,
  const double sigmaR ,
  const double beta   )
  : m_m0         (  0 )
  , m_sigmaL     (  1 )
  , m_sigmaR     (  1 )
  , m_beta       (  1 )
  , m_workspace () 
{
  //
  setM0     ( m0     ) ;
  setSigmaL ( sigmaL ) ;
  setSigmaR ( sigmaR ) ;
  setBeta   ( beta   ) ;
  //
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Apollonios2::~Apollonios2(){}
// ============================================================================
bool  Ostap::Math::Apollonios2::setM0 ( const double value )
{
  //
  if ( s_equal ( value , m_m0 ) ) { return false ; }
  //
  m_m0       = value ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apollonios2::setSigmaL ( const double value )
{
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_sigmaL ) ) { return false ; }
  //
  m_sigmaL    = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apollonios2::setSigmaR ( const double value )
{
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_sigmaR ) ) { return false ; }
  //
  m_sigmaR    = value_ ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Apollonios2::setBeta ( const double value )
{
  //
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_beta ) ) { return false ; }
  //
  m_beta    = value_  ;
  //
  if ( s_equal ( m_beta , 0 ) ) { m_beta = 0 ; }
  if ( s_equal ( m_beta , 1 ) ) { m_beta = 1 ; }
  //
  return true ;
}
// ============================================================================
//  evaluate Apollonios' function
// ============================================================================
double Ostap::Math::Apollonios2::pdf ( const double x ) const
{
  //
  const double dx = 
    ( x < m_m0 ) ?
    ( x - m_m0 ) / m_sigmaL :
    ( x - m_m0 ) / m_sigmaR ;
  //
  // the peak
  //
  return my_exp ( beta() * ( beta()  - std::sqrt ( b2 () + dx * dx ) ) ) * s_SQRT2PIi / sigma()  ;  
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::Apollonios2::integral
( const double low ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  else if (           low > high   ) { return - integral ( high ,
                                                           low  ) ; } // RETURN
  //
  const double xR = m_m0 + 4.0 * m_sigmaR ;
  if ( low < xR && xR < high ) 
  { return integral ( low , xR ) + integral ( xR , high ) ; }
  //
  const double xL = m_m0 - 4.0 * m_sigmaL ;
  if ( low < xL && xL < high ) 
  { return integral ( low , xL ) + integral ( xL , high ) ; }
  //
  const double in_tail = ( low >= xR || high <= xL ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Apollonios2> s_integrator {} ;
  static char s_message[] = "Integral(Apollonios2)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () , 
      &F     , 
      low , high  ,                  // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Apollonios2::tag () const 
{ 
  static const std::string s_name = "Apollonios2" ;
  return Ostap::Utils::hash_combiner ( s_name , m_m0 , m_sigmaL , m_sigmaR , m_beta ) ; 
}
// ============================================================================


// ============================================================================
/*  constructor with all parameters
 *  @param mean  \f$\mu\f$-parameter 
 *  @param sigma \f$\sigma\f$-parameter 
 */
// ============================================================================
Ostap::Math::Atlas::Atlas
( const double mean  ,
  const double sigma ) 
  : m_mean      ( mean               ) 
  , m_sigma     ( std::abs ( sigma ) )    
  , m_workspace ()
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Atlas::~Atlas(){}
// ============================================================================
// get variance:  very good numerical approximation d
// ============================================================================
double Ostap::Math::Atlas::variance () const { return 3 * m_sigma * m_sigma ; }
// ============================================================================
// get rms :  very good numerical approximation 
// ============================================================================
double Ostap::Math::Atlas::rms      () const { return s_SQRT3     * m_sigma ; }
// ============================================================================
bool Ostap::Math::Atlas::setMean ( const double value ) 
{
  if ( s_equal ( value , m_mean ) ) { return false ; }
  m_mean = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Atlas::setSigma ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_sigma ) ) { return false ; }
  m_sigma = value_ ;
  return true ;
}
// ============================================================================
// evaluate atlas function 
// ============================================================================
double Ostap::Math::Atlas::pdf        ( const double x ) const 
{
  const double dx = std::abs  ( x - m_mean ) / m_sigma ;
  if ( s_zero ( dx ) ) { return 1 ; }                        // return 1 
  const double x2 = std::pow ( dx , 1.0 + 1 / ( 1 + 0.5 * dx ) ) ;
  return std::exp ( -0.5 * x2 ) / ( s_ATLAS * m_sigma )  ;
}
// ============================================================================
double Ostap::Math::Atlas::integral
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low ,high ) ) { return 0 ; }
  else if ( low > high            ) { return -integral ( high , low ) ; }
  //
  // split 
  if ( low < m_mean && m_mean < high ) 
  { return integral ( low , m_mean ) + integral ( m_mean , high ) ; }
  //
  const double left  = m_mean - 5 * m_sigma ;  
  if ( low < left   &&  left < high ) 
  { return integral ( low , left   ) + integral ( left   , high ) ; }
  //
  const double right = m_mean + 5 * m_sigma ;  
  if ( low < right  && right  < high ) 
  { return integral ( low , right  ) + integral ( right  , high ) ; }
  //
  const bool in_tail = ( high <= left || low >= right ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Atlas> s_integrator {} ;
  static char s_message[] = "Integral(Atlas)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () , 
      &F     , 
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// overall integral, not exact but precise enough...
// ============================================================================
double Ostap::Math::Atlas::integral () const { return 1 ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Atlas::tag () const 
{ 
  static const std::string s_name = "Atlas" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mean , m_sigma ) ; 
}
// ============================================================================



// ============================================================================
// Sech 
// ============================================================================
/*  constructor with all parameters
 *  @param mean  \f$\mu\f$-parameter 
 *  @param sigma \f$\sigma\f$-parameter
 */
// ============================================================================
Ostap::Math::Sech::Sech  
( const double mean  ,
  const double sigma ) 
  : m_mean  (             mean    ) 
  , m_sigma (  std::abs ( sigma ) )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Sech::~Sech(){}
// ============================================================================
// evaluate sech function 
// ============================================================================
double Ostap::Math::Sech::pdf ( const double x ) const 
{
  const long double y = ( x - m_mean ) * M_PI_2 / m_sigma ;
  return 
    GSL_LOG_DBL_MAX < std::abs ( y )  ? 0 : 
    0.5 / ( m_sigma * std::cosh ( y ) ) ;
}
// ============================================================================
bool Ostap::Math::Sech::setMean ( const double value ) 
{
  if ( s_equal ( value , m_mean  ) ) { return false ; }
  m_mean  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Sech::setSigma ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_sigma  ) ) { return false ; }
  m_sigma  = value_ ;
  return true ;
}
// ============================================================================
// get integral from low to high 
// ============================================================================
double Ostap::Math::Sech::integral 
( const double low  ,
  const double high ) const 
{ return s_equal ( low , high ) ? 0.0 : cdf ( high ) - cdf ( low ) ; }
// ============================================================================
// get integral from -infinity to + infinity 
// ============================================================================
double Ostap::Math::Sech::integral () const { return 1 ; } 
// ============================================================================
// evaluate CDF function 
// ============================================================================
double Ostap::Math::Sech::cdf ( const double x ) const 
{
  const long double y = ( x - m_mean ) * M_PI_2 / m_sigma ;
  return
    ( GSL_LOG_DBL_MAX < y ) ? 1 :
    ( GSL_LOG_DBL_MIN > y ) ? 0 : 
    std::atan (  std::exp ( y ) ) / M_PI_2 ;
}
// ============================================================================
// get quantile (0<p<1)
// ============================================================================
double Ostap::Math::Sech::quantile ( const double p ) const 
{ return 
    0 >= p || s_zero  ( p     ) ? -s_INFINITY :
    1 <= p || s_equal ( p , 1 ) ? +s_INFINITY : 
    m_mean + m_sigma * 2 / M_PI * std::log( std::tan ( M_PI * p /  2 ) ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Sech::tag () const 
{ 
  static const std::string s_name = "Sech" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mean , m_sigma ) ; 
}
// ============================================================================

// ============================================================================
/*  constructor with all parameters
 *  @param mean  \f$\mu\f$-parameter 
 *  @param alpha \f$\alpha\f$-parameter 
 *  @param beta  \f$\beta\f$-parameter 
 */
// ============================================================================
Ostap::Math::Losev::Losev
( const double mu    , 
  const double alpha , 
  const double beta  ) 
  : m_mu        ( mu                 ) 
  , m_alpha     ( std::abs ( alpha ) ) 
  , m_beta      ( std::abs ( beta  ) ) 
  , m_norm      ( -1 )
  , m_workspace () 
{}
// ============================================================================
bool Ostap::Math::Losev::setMu ( const double value ) 
{
  if ( s_equal ( value , m_mu  ) ) { return false ; }
  m_mu  = value ;
  return true ;
}
// =============================================================================
bool Ostap::Math::Losev::setAlpha ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_alpha  ) ) { return false ; }
  m_alpha  = v  ;
  m_norm   = -1 ;
  return true ;
}
// =============================================================================
bool Ostap::Math::Losev::setBeta ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_beta  ) ) { return false ; }
  m_beta   = v  ;
  m_norm   = -1 ;
  return true ;
}
// =============================================================================
// the mode of the distribution 
// =============================================================================
double Ostap::Math::Losev::mode () const 
{ return m_mu + std::log ( m_alpha / m_beta ) / ( m_alpha + m_beta ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Losev::tag () const 
{ 
  static const std::string s_name = "Losev" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_alpha , m_beta ) ; 
}
// =============================================================================
// evaluate the function 
// =============================================================================
double Ostap::Math::Losev::pdf ( const double x ) const 
{
  if ( m_norm <= 0 ) 
  {
    const double sumab = m_alpha + m_beta ;
    m_norm = sumab * std::sin ( M_PI * m_beta / sumab ) / M_PI ;
  }
  //
  const double dx = x - m_mu ;
  return 0 <= dx  ? 
    m_norm * std::exp ( -m_beta  * dx ) / ( 1 + std::exp  ( - ( m_alpha + m_beta ) * dx ) ) :
    m_norm * std::exp (  m_alpha * dx ) / ( 1 + std::exp  (   ( m_alpha + m_beta ) * dx ) ) ;  
}
// ============================================================================
/*  get the integral between low and high values 
 *  \f$ \int_{low}^{high}f(x) dx\f$
 */
// ============================================================================
double Ostap::Math::Losev::integral 
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low ,high ) ) { return 0 ; }
  else if ( low > high            ) { return -integral ( high , low ) ; }
  //  
  // split 
  const double left  = m_mu - 6 * m_alpha;  
  if ( low < left && left < high ) 
  { return integral ( low , left ) + integral ( left   , high ) ; }
  //
  const double right = m_mu + 6 * m_beta ;  
  if ( low < right && right < high ) 
  { return integral ( low , right ) + integral ( right , high ) ; }
  //
  const bool in_tail = ( high <= left || low >= right ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Losev> s_integrator {} ;
  static char s_message[] = "Integral(Losev)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () , 
      &F     , 
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// Logistic
// ============================================================================
/*  constructor with all parameters
 *  @param mean  \f$\mu\f$-parameter 
 *  @param sigma \f$\sigma\f$-parameter
 */
// ============================================================================
Ostap::Math::Logistic::Logistic
( const double mean  ,
  const double sigma ) 
  : m_mean  (             mean    ) 
  , m_sigma (  std::abs ( sigma ) )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Logistic::~Logistic(){}
// ============================================================================
// evaluate sech function 
// ============================================================================
double Ostap::Math::Logistic::pdf ( const double x ) const 
{
  const double s = m_sigma * s_SQRT3overPI ;
  const long double y = ( x - m_mean ) / ( 2 * s ) ;
  if  ( GSL_LOG_DBL_MAX < std::abs( y ) ) { return 0 ; }
  const long double c = std::cosh ( y ) ;
  return 0.25 / c / c / s  ;
}
// ============================================================================
bool Ostap::Math::Logistic::setMean ( const double value ) 
{
  if ( s_equal ( value , m_mean  ) ) { return false ; }
  m_mean  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Logistic::setSigma ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_sigma  ) ) { return false ; }
  m_sigma  = value_ ;
  return true ;
}
// ============================================================================
// get integral from low to high 
// ============================================================================
double Ostap::Math::Logistic::integral 
( const double low  ,
  const double high ) const 
{ return s_equal ( low , high ) ? 0.0 : cdf ( high ) - cdf ( low ) ; }
// ============================================================================
// get integral from -infinity to + infinity 
// ============================================================================
double Ostap::Math::Logistic::integral () const { return 1 ; } 
// ============================================================================
// evaluate CDF function 
// ============================================================================
double Ostap::Math::Logistic::cdf ( const double x ) const 
{
  const double s = m_sigma * s_SQRT3overPI ;
  const long double y = ( x - m_mean ) / ( 2 * s ) ;
  return 0.5 * ( 1 + std::tanh ( y ) ) ;
}
// ============================================================================
// get parameter s 
// ============================================================================
double Ostap::Math::Logistic::s() const 
{ return m_sigma * s_SQRT3overPI ; }
// ============================================================================
// quantile function  (0<p<1)
// ============================================================================
double Ostap::Math::Logistic::quantile ( const double p ) const 
{ return
    0 >= p || s_zero  ( p     ) ? -s_INFINITY :
    1 <= p || s_equal ( p , 1 ) ? +s_INFINITY : 
    m_mean + m_sigma * s_SQRT3overPI * std::log ( p / ( 1 - p ) ) ; }
// ===========================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Logistic::tag () const 
{ 
  static const std::string s_name = "Logistic" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mean , m_sigma ) ;
}
// ============================================================================








// ============================================================================
// Student-T 
// ============================================================================
/*  constructor from mass, resolution and "n"-parameter 
 *  @param M     mass 
 *  @param sigma width parameter
 *  @param N     n-parameter  ( actually  n=1+|N| ) 
 */
// ============================================================================
Ostap::Math::StudentT::StudentT 
( const double mass  , 
  const double sigma ,
  const double n     ) 
//
  : m_M    (      std::abs ( mass  ) )
  , m_s    (      std::abs ( sigma ) )
  , m_n    ( -1 )
  , m_norm ( -1 ) 
{
  setN ( n ) ;  
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::StudentT::~StudentT (){}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::StudentT::setM ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_M ) ) { return false ; }
  //
  m_M = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::StudentT::setSigma ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_s ) ) { return false ; }
  //
  m_s = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::StudentT::setN ( const double x )
{
  //
  const double v = 1 + std::abs ( x ) ;
  //
  if ( m_norm < 0 ) 
  {
    m_norm  = gsl_sf_gamma ( 0.5 * ( v + 1 ) ) / gsl_sf_gamma ( 0.5 * v ) ;  
    m_norm /= std::sqrt    ( M_PI * v ) ;
  }
  //
  if ( s_equal ( v , m_n ) ) { return false ; }
  //
  m_n = v ;
  //
  m_norm  = gsl_sf_gamma ( 0.5 * ( v + 1 ) ) / gsl_sf_gamma ( 0.5 * v ) ;  
  m_norm /= std::sqrt    ( M_PI * v ) ;
  //
  return true ;
}
// ==========================================================================
double Ostap::Math::StudentT::pdf ( const double x ) const
{
  //
  const double y = ( x - M () ) / sigma () ;
  //
  const double f = std::pow (  1 + y * y / nu () ,  -0.5 * ( nu () + 1 ) ) ;
  //
  return m_norm * f / sigma () ; // sigma comes from dx = dy * sigma 
}
// ============================================================================
double Ostap::Math::StudentT::cdf ( const double y ) const
{
  //
  const double  t    = ( y - M () ) / sigma () ;
  return Ostap::Math::student_cdf ( t , nu() ) ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::StudentT::integral() const { return 1 ; }
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::StudentT::integral
( const double low  , 
  const double high ) const 
{
  //
  if ( s_equal ( low , high ) ) { return                 0.0 ; } // RETURN
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::StudentT::tag () const 
{ 
  static const std::string s_name = "StudentT" ;
  return Ostap::Utils::hash_combiner ( s_name , m_M , m_s , m_n ) ; 
}
// ============================================================================

// ============================================================================
// Bifurcated Student-T 
// ============================================================================
/*  constructor from mass, resolution and "n"-parameter 
 *  @param M     mass 
 *  @param sigma width parameter
 *  @param N     n-parameter  ( actually  n=1+|N| ) 
 */
// ============================================================================
Ostap::Math::BifurcatedStudentT::BifurcatedStudentT 
( const double mass   , 
  const double sigmaL ,
  const double sigmaR ,
  const double nL     , 
  const double nR     ) 
//
  : m_M     (      std::abs ( mass   ) )
  , m_sL    (      std::abs ( sigmaL ) )
  , m_sR    (      std::abs ( sigmaR ) )
  , m_nL    ( -1 )
  , m_nR    ( -1 )
  , m_normL ( -1 ) 
  , m_normR ( -1 ) 
{
  setNL ( nL ) ;  
  setNR ( nR ) ;  
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::BifurcatedStudentT::~BifurcatedStudentT (){}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setM ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_M ) ) { return false ; }
  //
  m_M = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setSigmaL ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_sL ) ) { return false ; }
  //
  m_sL = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setSigmaR ( const double x )
{
  //
  const double v = std::abs ( x ) ;
  if ( s_equal ( v , m_sR ) ) { return false ; }
  //
  m_sR = v ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setNL ( const double x )
{
  //
  const double v = 1 + std::abs ( x ) ;
  //
  if ( m_normL < 0 ) 
  {
    m_normL = gsl_sf_gamma  ( 0.5 * ( v + 1 ) ) / gsl_sf_gamma ( 0.5 * v ) ;  
    m_normL /= std::sqrt    ( M_PI * v ) ;
  }
  //
  if ( s_equal ( v , m_nL ) ) { return false ; }
  //
  m_nL =  v ;
  //
  m_normL = gsl_sf_gamma  ( 0.5 * ( v + 1 ) ) / gsl_sf_gamma ( 0.5 * v ) ;  
  m_normL /= std::sqrt    ( M_PI * v ) ;
  //
  return true ;
}
// ============================================================================
// set the proper parameters
// ============================================================================
bool Ostap::Math::BifurcatedStudentT::setNR ( const double x )
{
  //
  const double v = 1 + std::abs ( x ) ;
  //
  if ( m_normR < 0 ) 
  {
    m_normR  = std::tgamma ( 0.5 * ( v + 1 ) ) / std::tgamma ( 0.5 * v ) ;  
    m_normR /= std::sqrt   ( M_PI * v ) ;
  }
  //
  if ( s_equal ( v , m_nR ) ) { return false ; }
  //
  m_nR = v ;
  //
  m_normR  = std::tgamma ( 0.5 * ( v + 1 ) ) / std::tgamma ( 0.5 * v ) ;  
  m_normR /= std::sqrt   ( M_PI * v ) ;
  //
  return true ;
}
// ==========================================================================
double Ostap::Math::BifurcatedStudentT::pdf ( const double x ) const
{
  //
  const double y = ( x <= M() ) ? 
    ( x - M () ) / sigmaL () : ( x - M () ) / sigmaR () ;
  //
  const double f = ( x <= M() ) ? 
    std::pow (  1 + y * y / nuL () ,  -0.5 * ( nuL () + 1 ) ) :
    std::pow (  1 + y * y / nuR () ,  -0.5 * ( nuR () + 1 ) ) ;
  //
  const double n_1 = m_normL       / sigmaL () ;
  const double n_2 = m_normR       / sigmaR () ;
  const double n_t = 2 * n_1 * n_2 / ( n_1 + n_2 ) ;
  //
  return n_t * f ; 
}
// ============================================================================
double Ostap::Math::BifurcatedStudentT::cdf ( const double y ) const
{
  //
  const double n_1 = m_normL / sigmaL () ;
  const double n_2 = m_normR / sigmaR () ;
  //
  if ( y <= M() ) 
  {
    const double  t    = ( y - M () ) / sigmaL () ;
    return     2 * n_2 / ( n_1 + n_2 ) * Ostap::Math::student_cdf (  t , nuL () ) ;  
  }
  //
  const   double  t    = ( y - M () ) / sigmaR () ;
  return   1 - 2 * n_1 / ( n_1 + n_2 ) * Ostap::Math::student_cdf ( -t , nuR () ) ;  
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::BifurcatedStudentT::integral() const { return 1 ; }
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::BifurcatedStudentT::integral
( const double low  , 
  const double high ) const 
{
  //
  if ( s_equal ( low , high ) ) { return 0.0 ; } // RETURN
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::BifurcatedStudentT::tag () const 
{  
  static const std::string s_name = "BiFurcatedStudentT" ;
  return Ostap::Utils::hash_combiner ( s_name , m_M , m_sL , m_sR  , m_nL , m_nR ) ; 
}
// ============================================================================



// ============================================================================
/*  constructor from all parameters 
 *  @param mu    location parameter 
 *  @param sigma width/scale parameter 
 *  @param n     n-parameter 
 *  @param kappa asymmetry parameter 
 */
// ============================================================================
Ostap::Math::PearsonIV::PearsonIV
( const double mu       , 
  const double varsigma , 
  const double n        , 
  const double kappa     )
  : m_mu       ( mu ) 
  , m_varsigma ( std::abs ( varsigma ) )
  , m_n        ( std::abs ( n        ) ) 
  , m_kappa    ( kappa )
  , m_C        ( -1 ) 
{
  setN ( n ) ;
}
// ===========================================================================
// get value of the function 
// ===========================================================================
double Ostap::Math::PearsonIV::evaluate ( const double x ) const 
{
  const double y = ( x - m_mu ) / m_varsigma ;
  const double s = m_C * std::pow ( 1 + y * y , - m() ) / m_varsigma ;
  return s_zero ( m_kappa ) ? s : s * std::exp ( -m_kappa * std::atan ( y ) ) ;
}
// ===========================================================================
// set location parameter
// ===========================================================================
bool Ostap::Math::PearsonIV::setMu ( const double value )
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ===========================================================================
// set width/scale parameter
// ===========================================================================
bool Ostap::Math::PearsonIV::setVarsigma ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_varsigma , avalue ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ===========================================================================
// set n-parameter
// ===========================================================================
bool Ostap::Math::PearsonIV::setN ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_n , avalue ) && 0 < m_C ) { return false ; }
  m_n = avalue ;
  //
  // m_C = std::norm ( Ostap::Math::gamma ( std::complex<double> ( m() , 0.5 * nu () ) ) /
  //                   Ostap::Math::gamma ( m() ) ) / std::beta  ( m() - 0.5 , 0.5 ) ;
  //
  m_C = Ostap::Math::pearsonIV_g2 ( m() , 0.5 * nu() ) / std::beta  ( m() - 0.5 , 0.5 ) ;
  //
  return true ;
}   
// ===========================================================================
// set asymmetry parameter
// ===========================================================================
bool Ostap::Math::PearsonIV::setKappa ( const double value )
{
  if ( s_equal ( m_kappa , value ) && 0 < m_C ) { return false ; }
  m_kappa = value ;
  //
  // m_C     = std::norm ( Ostap::Math::gamma ( std::complex<double> ( m() , 0.5 * nu () ) ) /
  //                       Ostap::Math::gamma ( m() ) ) / std::beta  ( m() - 0.5 , 0.5 ) ;
  //
  m_C = Ostap::Math::pearsonIV_g2 ( m() , 0.5 * nu() ) / std::beta  ( m() - 0.5 , 0.5 ) ;
  //
  return true ;
}  
// ===========================================================================
// get the integral
// ===========================================================================
double Ostap::Math::PearsonIV::integral () const { return 1 ; }
// ===========================================================================
// get the integral between low and high limits
// ===========================================================================
double Ostap::Math::PearsonIV::integral
( const double low  ,
  const double high ) const 
{
  //
  if ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low ) { return - integral ( high , low ) ; }
  //
  const bool symmetric = s_zero ( nu () ) ;
  //
  const double m0      = ( 1 < m() ) ? 0.5 * ( mode () + mean () ) : mode () ;
  const double width   = ( 2 * m() <= 3 ) ? 0.5 * 
    std::max ( m_varsigma , 0.5 * infection_width () ) : rms () ;
  //
  { // split at mode 
    const double x1 = mode () ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = x1 + 2 * width ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
    const double x3 = x1 - 2 * width ;
    if ( low < x2 && x2 < high ) { return integral ( low , x3 ) + integral ( x3 , high ) ; }
  }
  //
  if ( !symmetric && ( 1 < m () ) ) 
  {
    // split at mean 
    const double x1 = mean ()  ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = x1 + 2 * width ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
    const double x3 = x1 - 2 * width ;
    if ( low < x2 && x2 < high ) { return integral ( low , x3 ) + integral ( x3 , high ) ; }
  }
  //
  // more splits 
  { 
    const double x2 = m0 - 6 * width  ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
    const double x3 = m0 + 6 * width  ;
    if ( low < x3 && x3 < high ) { return integral ( low , x3 ) + integral ( x3 , high ) ; }
  }
  //
  if ( !symmetric && ( 0 < nu() ) )
  {
    const double xx = m0 - 12 * width  ;
    if ( low < xx && xx < high ) { return integral ( low , xx ) + integral ( xx , high ) ; }
  }
  //
  if ( !symmetric && ( 0 > nu() ) )
  {
    const double xx = m0 + 12 * width  ;
    if ( low < xx && xx < high ) { return integral ( low , xx ) + integral ( xx , high ) ; }
  }
  //
  //
  const bool in_tail = ( high <= m0 - 10 * width ) || ( low >=  m0 + 10 * width ) ;             
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PearsonIV> s_integrator {} ;
  static char s_message[] = "Integral(PEarsonIV)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
             ( tag () , 
               &F     , 
               low    , high  ,               // low & high edges
               workspace ( m_workspace ) ,    // workspace
               in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
               in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
               m_workspace.size()   ,          // size of workspace
               s_message           , 
               __FILE__ , __LINE__ ) ;
  //
  return result ; 
}
// ============================================================================
// mode 
// ============================================================================
double Ostap::Math::PearsonIV::mode () const // mode of the distribution 
{ return m_mu - 0.5 * nu() * a() / m() ; }
// ============================================================================
// mean
// ============================================================================
double Ostap::Math::PearsonIV::mean () const // mean of the distribution 
{
  return 
    s_zero ( nu() ) ? m_mu                                  :
    1 < m ()        ? m_mu - 0.5 * a() * nu() / ( m() - 1 ) :
    std::copysign ( std::numeric_limits<double>::infinity() , -nu() ) ;
}
// ============================================================================
// (central) moment 
// ============================================================================
double Ostap::Math::PearsonIV::moment ( const unsigned short k ) const 
{
  if      ( 0 == k ) { return 1 ; }
  else if ( 1 == k ) { return 0 ; }
  //
  const bool odd = ( 1 == k % 2 )  ;
  if ( odd && s_zero ( nu () ) ) { return 0 ; }
  //
  if ( r () + 1 <= k ) 
  {
    return odd ?
      std::copysign ( std::numeric_limits<double>::infinity () , -nu() ) :
      std::numeric_limits<double>::infinity () ;
  }
  //
  const double r2  = std::pow ( r  () , 2 ) ;
  const double nu2 = std::pow ( nu () , 2 ) ;
  //
  double m2 = 1  ;
  double m1 = 0  ;
  double m  = m1 ;
  for  ( unsigned short kk = 2 ; kk <= k ; ++kk ) 
  {
    const double c = a () * ( kk - 1 ) / ( r2 * ( r () - ( kk - 1 ) ) ) ;
    m  = -2 * nu() * r () * m1 + a () * ( r2 + nu2 ) * m2 ;
    m *= c  ;
    m2 = m1 ;
    m1 = m  ;
  }
  return m ;
}
// ============================================================================
// variance                        (for m>3.5)
// ============================================================================
double Ostap::Math::PearsonIV::variance () const  // variance of the distribution 
{ return 2 * m () <= 3 ? std::numeric_limits<double>::infinity () : moment ( 2 ) ; }
// ============================================================================
// rms                              (for m>3.5)
// ============================================================================
double Ostap::Math::PearsonIV::rms     () const  // rms of the distribution 
{ return 2 * m () <= 3 ? std::numeric_limits<double>::infinity () : std::sqrt ( moment ( 2 ) ) ; }
// ============================================================================
// skewness ( for m>2) 
// ============================================================================
double Ostap::Math::PearsonIV::skewness () const  // skewness 
{
  return 
    s_zero ( nu () ) ? 0.0 :
    m () <= 2 ? std::copysign ( std::numeric_limits<double>::infinity () , -nu() ) :
    moment ( 3 ) / std::pow ( moment ( 2 ) , 1.5 ) ;
}
// ============================================================================
// (excessive) kurtosis ( for m>5/2) 
// ============================================================================
double Ostap::Math::PearsonIV::kurtosis () const  // (excessive) kurtosis 
{
  return 
    2 * m () <= 5 ? std::numeric_limits<double>::infinity () :
    moment ( 4 ) / std::pow ( moment ( 2 ) , 2) - 3 ;
}
// ============================================================================
// beta1 parameter of Pearson family (m>2) 
// ============================================================================
double Ostap::Math::PearsonIV::beta1 () const  // beta1 parameter of Pearson family 
{ 
  return 
    s_zero ( nu () ) ? 0.0 :
    m () <= 2 ? std::numeric_limits<double>::infinity () :
    std::pow ( moment ( 3 ) , 2 ) / std::pow ( moment ( 2 ) , 3 ) ;
}
// ============================================================================
// beta2 parameter of Pearson family (m>2) 
// ============================================================================
double Ostap::Math::PearsonIV::beta2 () const  // beta2 parameter of Pearson family 
{ 
  return 
    2 * m () <= 5 ? std::numeric_limits<double>::infinity () :
    moment ( 4 ) / std::pow ( moment ( 2 ) , 2) ;
}
// ============================================================================
/* distance between two infection points:
 *  distance between two points with \f$ f^{\prime\prime}=0\f$.
 *  the twp points are equidstance fro mthe mode 
 */
// ============================================================================
double Ostap::Math::PearsonIV::infection_width () const 
{
  return a () / m() * std::sqrt ( ( 4 * std::pow ( m() , 2 ) + 
                                    std::pow ( nu () , 2 ) ) / ( 2 * m() + 1 ) ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::PearsonIV::tag () const 
{  
  static const std::string s_name = "PearsonIV" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_varsigma , m_n , m_kappa ) ; 
}
// ============================================================================


// ============================================================================
/*  constructor with all parameters
 *  @param location \f$\mu\f$-parameter       \f$-\inf<\mu<+\inf\f$
 *  @param scale    \f$\sigma\f$-parameter    \f$0<\sigma\f$
 *  @param epsilon  \f$\epsilon\f$-parameter  \f$-\inf<\epsilon<+\inf\f$
 *  @param delta    \f$\delta\f$-parameter    \f$0<\epsilon<+\inf\f$
 */
// ============================================================================
Ostap::Math::SinhAsinh::SinhAsinh 
( const double location  ,
  const double scale     , 
  const double epsilon   , 
  const double delta     )
  : m_mu       (            location   ) 
  , m_sigma    ( std::abs ( scale    ) ) 
  , m_epsilon  (            epsilon    ) 
  , m_delta    ( std::abs ( delta    ) ) 
{}
// ============================================================================
Ostap::Math::SinhAsinh::~SinhAsinh(){}
// ============================================================================
bool Ostap::Math::SinhAsinh::setMu      ( const double value ) 
{
  if ( s_equal ( value , m_mu  ) ) { return false ; }
  m_mu  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::SinhAsinh::setSigma      ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_sigma  ) ) { return false ; }
  m_sigma  = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::SinhAsinh::setEpsilon ( const double value ) 
{
  if ( s_equal ( value , m_epsilon  ) ) { return false ; }
  m_epsilon  = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::SinhAsinh::setDelta ( const double value ) 
{
  const double value_ = std::abs ( value ) ;
  if ( s_equal ( value_ , m_delta  ) ) { return false ; }
  m_delta  = value_ ;
  return true ;
}
// ============================================================================
// evaluate sinhasinh-distributions
// ============================================================================
double Ostap::Math::SinhAsinh::pdf ( const double x ) const 
{
  //
  const double y = ( x - mu () ) / sigma()  ;
  const double z = shash ( y , epsilon() , delta() )  ;
  //
  const double r = s_SQRT2PIi * delta() 
    // * std::sqrt  ( ( 1.0 + z * z ) / ( 1.0 + y * y )  ) 
    *    std::hypot ( 1 , z )          / std::hypot ( 1 , y )  
    *    my_exp ( -0.5 * z * z ) ;
  //
  return  r / sigma() ;
}
// ============================================================================
// evaluate sinhasinh cimulative distribution
// ============================================================================
double Ostap::Math::SinhAsinh::cdf ( const double x ) const 
{
  //
  const double y = ( x - mu() ) / sigma()  ;
  const double z = shash ( y , epsilon () , delta () )  ;
  //
  return gsl_cdf_ugaussian_P ( z ) ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::SinhAsinh::integral ( const double low  ,
                                          const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  return cdf ( high ) - cdf ( low ) ;
}
// ======================================================================
// get median
// ======================================================================
double Ostap::Math::SinhAsinh::median () const 
{ return m_mu - m_sigma * std::sinh ( m_epsilon / m_delta ) ; }
// ============================================================================
// get the mean for the distgribution
// ============================================================================
double Ostap::Math::SinhAsinh::mean   () const 
{
  const double d1 = 0.5 * ( 1 - m_delta ) / m_delta ;
  const double d2 = 0.5 * ( 1 - m_delta ) / m_delta ;
  //
  static const double s_const1 = std::pow ( std::exp ( 1 ) , 0.25 ) / std::sqrt ( 8 *M_PI ) ;
  //
  const double a = std::sinh ( m_epsilon / m_delta ) * s_const1 *
    ( Ostap::Math::bessel_Knu ( d1 , 0.25 ) + 
      Ostap::Math::bessel_Knu ( d2 , 0.25 ) ) ;
  //
  return m_mu - m_sigma * a ;
}
// ============================================================================
// get the variance for the distribution
// ============================================================================
double Ostap::Math::SinhAsinh::variance () const
{
  //
  const double d1 = 0.5 * ( 1 + m_delta ) / m_delta ;
  const double d2 = 0.5 * ( 1 - m_delta ) / m_delta ;
  //
  static const double s_const1 = std::pow ( std::exp ( 1 ) , 0.25 ) / std::sqrt ( 8 *M_PI ) ;
  //
  const double a = std::sinh (     m_epsilon / m_delta ) * s_const1 *
    ( Ostap::Math::bessel_Knu ( d1 , 0.25 ) + 
      Ostap::Math::bessel_Knu ( d2 , 0.25 ) ) ;
  //
  const double p1 = 0.5 * ( 2 + m_delta ) / m_delta ;
  const double p2 = 0.5 * ( 2 - m_delta ) / m_delta ;
  //
  static const double s_const2 = s_const1 / 2 ;
  //
  const double b = std::cosh ( 2 * m_epsilon / m_delta ) * s_const2 *
    ( Ostap::Math::bessel_Knu ( p1 , 0.25 ) + 
      Ostap::Math::bessel_Knu ( p2 , 0.25 ) ) ;
  //
  return m_sigma * m_sigma * ( b - a * a - 0.5 ) ;
}
// ============================================================================
// get the RMS for the distribution
// ============================================================================
double Ostap::Math::SinhAsinh::rms () const
{ return std::sqrt ( variance () ) ; }

// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::SinhAsinh::tag () const 
{  
  static const std::string s_name = "SinhAsinh" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_sigma , m_epsilon , m_delta ) ; 
}
// ============================================================================

// ============================================================================
// Johnson-SU
// ============================================================================
/*  constructor with all parameters
 *  @param xi     \f$\xi\f$-parameter       \f$-\inf<\xi<+\inf\f$
 *  @param lambda \f$\lambda\f$-parameter   \f$   0<\lambda<+\inf\f$
 *  @param delta  \f$\delta\f$-parameter    \f$   0<\delta<+\inf\f$
 *  @param gamma  \f$\gamma\f$-parameter    \f$-\inf<\epsilon<+\inf\f$
 */
// ============================================================================
Ostap::Math::JohnsonSU::JohnsonSU  
( const double xi      ,   // related to location 
  const double lambda  ,   // related to variance
  const double delta   ,   // shape 
  const double gamma   )  // shape 
  : m_xi      (            xi) 
  , m_lambda  ( std::abs ( lambda ) ) 
  , m_delta   ( std::abs ( delta  ) ) 
  , m_gamma   (            gamma    ) 
{}
// ============================================================================
// Destructor
// ============================================================================
Ostap::Math::JohnsonSU::~JohnsonSU(){}
// ============================================================================
// get the mean value
// ============================================================================
double Ostap::Math::JohnsonSU::mean() const 
{
  const double d = 
    std::exp   ( 0.5 / ( m_delta * m_delta ) ) * 
    std::sinh  ( m_gamma / m_delta           ) ;
  //
  return m_xi - m_lambda * d ;
}
// ============================================================================
// get the variance
// ============================================================================
double Ostap::Math::JohnsonSU::variance () const 
{
  const double d1 = std::exp ( 1.0 / ( m_delta * m_delta ) );
  //
  const double d2 = ( d1 - 1 ) * ( d1 * std::cosh ( 2  *m_gamma / m_delta ) + 1 ) ;
  //
  return 0.5 * m_lambda * m_lambda * d2 ;
}
// ============================================================================
bool Ostap::Math::JohnsonSU::setXi( const double value ) 
{
  if ( s_equal ( value , m_xi ) ) { return false ; }
  m_xi = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::JohnsonSU::setGamma( const double value ) 
{
  if ( s_equal ( value , m_gamma ) ) { return false ; }
  m_gamma = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::JohnsonSU::setLambda ( const double value ) 
{
  const double value_ = std::abs ( value ) ;  
  if ( s_equal ( value_ , m_lambda ) ) { return false ; }
  m_lambda = value_ ;
  return true ;
}
// ============================================================================
bool Ostap::Math::JohnsonSU::setDelta ( const double value ) 
{
  const double value_ = std::abs ( value ) ;  
  if ( s_equal ( value_ , m_delta ) ) { return false ; }
  m_delta = value_ ;
  return true ;
}
// ============================================================================
// evaluate JohnsonSU-distributions
// ============================================================================
double Ostap::Math::JohnsonSU::pdf        ( const double x ) const 
{
  // get z 
  const long double dx  = ( x - m_xi ) / m_lambda ;
  const long double z   = m_gamma + m_delta * std::asinh ( dx ) ;
  //
  const long double res = std::exp ( -0.5 * z * z ) / std::sqrt ( 1 + dx * dx ) ;
  //
  return res * m_delta / ( m_lambda * s_SQRT2PI ) ;
}
// ============================================================================
// evaluate JohnsonSU-distributions
// ============================================================================
double Ostap::Math::JohnsonSU::cdf        ( const double x ) const 
{
  // get z 
  const long double dx  = ( x - m_xi ) / m_lambda ;
  const long double z   = m_gamma + m_delta * std::asinh ( dx ) ;
  //
  return gsl_cdf_ugaussian_P ( z ) ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::JohnsonSU::integral 
( const double low  ,
  const double high ) const 
{ return  s_equal ( low , high ) ? 0.0 : ( cdf ( high ) - cdf ( low ) ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::JohnsonSU::tag () const 
{ 
  static const std::string s_name = "JohnsonSU" ;
  return Ostap::Utils::hash_combiner ( s_name , m_xi , m_lambda , m_delta , m_gamma ) ; 
}
// ============================================================================





// ============================================================================
/*  Constructor from location and mean 
 *  @param mu location 
 *  @param scale the scale, scale>0
 */
// ============================================================================
Ostap::Math::Slash::Slash
( const double mu    , // location 
  const double scale ) // scale 
  : m_mu    ( mu ) 
  , m_scale ( std::abs ( scale ) ) 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::Slash::~Slash(){}
// ============================================================================
bool Ostap::Math::Slash::setMu    ( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Slash::setScale ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_scale ) ) { return false ; }
  m_scale = v;
  return true ;
}
// ============================================================================
namespace 
{
  // ==========================================================================
  const long double s_slash = 0.5L / std::sqrt( 2.0L * M_PI );
  // ==========================================================================
  //  (phi(0)-phi(x))/x^2
  inline long double _slash_pdf_ ( const long double x ) 
  {
    if        ( s_zero ( x ) ) { return s_slash ; }
    else if   ( 0.1 < std::abs ( x ) ) 
    { return  ( 2 * s_slash -  Ostap::Math::gauss_pdf ( x ) ) / ( x * x ) ; }
    //
    const long double z = - 0.5L * x * x ;
    return  s_slash * ( std::expm1 ( z ) / z ) ;
  }
  // ==========================================================================
  //  Phi(x) - (phi(0)-phi(x))/x
  inline long double _slash_cdf_ ( const long double x ) 
  {
    return 
      s_equal ( x , 0 ) ? 0.5 : 
      Ostap::Math::gauss_cdf ( x ) - x * _slash_pdf_ ( x ) ;
  }
}
// ============================================================================
// evaluate slash function
// ============================================================================
double Ostap::Math::Slash::pdf ( const double x ) const
{
  const double y = ( x - m_mu ) / m_scale ;
  return _slash_pdf_ ( y ) / m_scale ;
}
// ============================================================================
// evaluate slash CDF 
// ============================================================================
double Ostap::Math::Slash::cdf ( const double x ) const
{
  const double y = ( x - m_mu ) / m_scale ;
  return _slash_cdf_ ( y ) ;
}
// ============================================================================
// get integral from low to high
// ============================================================================
double Ostap::Math::Slash::integral
( const double low  ,
  const double high ) const 
{ return s_equal ( low ,  high ) ? 0.0 : cdf ( high ) - cdf ( low ) ; }
// ===========================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Slash::tag () const 
{ 
  static const std::string s_name = "Slash" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_scale ) ; 
}
// ============================================================================




// ============================================================================
/*  constructor with all arguments 
 *  @param mu  the mean/mode/median of the distribution
 *  @param s   the width-parameteter
 */
// ============================================================================
Ostap::Math::RaisingCosine::RaisingCosine
( const double mu , 
  const double s  ) 
  : m_mu ( mu ) 
  , m_s  ( std::abs ( s ) ) 
{}
// ============================================================================
bool Ostap::Math::RaisingCosine::setS ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_s ) ) { return false ; }
  m_s = v;
  return true ;
}
// ============================================================================
bool Ostap::Math::RaisingCosine::setMu ( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
// evaluate raising cosine distribution
// ============================================================================
double Ostap::Math::RaisingCosine::pdf ( const double x ) const 
{
  return 
    x <= m_mu - m_s ? 0.0 : 
    x >= m_mu + m_s ? 0.0 : 
    ( 1  + std::cos ( M_PI * ( x  - m_mu ) / m_s ) ) / ( 2 * m_s )  ;
}
// ============================================================================
// variance  
// ============================================================================
double Ostap::Math::RaisingCosine::variance () const 
{
  static const double s_c1 = ( 1./3 - 2 / ( M_PI * M_PI ) ) ;
  return m_s * m_s * s_c1 ;
}
// ============================================================================
// rms
// ============================================================================
double Ostap::Math::RaisingCosine::rms () const 
{
  static const double s_c2 = std::sqrt ( 1./3 - 2 / ( M_PI * M_PI ) ) ;
  return m_s * s_c2 ;
}
// ============================================================================
// kurtosis
// ============================================================================
double Ostap::Math::RaisingCosine::kurtosis () const 
{
  static const double s_k = 
    1.2 * ( 90. - std::pow ( M_PI , 4 ) ) / std::pow ( M_PI*M_PI  - 6. , 2 ) ;
  return  s_k ;
}
// ============================================================================
// get CDF 
// ============================================================================
double Ostap::Math::RaisingCosine::cdf      ( const double x ) const 
{
  if      ( x <= m_mu - m_s ) { return 0 ; }
  else if ( x >= m_mu - m_s ) { return 1 ; }
  //
  const double y = ( x - m_mu ) / m_s ;
  return 0.5 * ( 1 + y + std::sin ( y  * M_PI ) / M_PI ) ;
}
// ============================================================================
// evaluate the integral
// ============================================================================
double Ostap::Math::RaisingCosine::integral
( const double low  , 
  const double high ) const 
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( low > high             ) { return -integral ( high , low ) ; }
  else if ( high <  m_mu - m_s     ) { return 0 ; }
  else if ( low  >  m_mu + m_s     ) { return 0 ; }
  //                                            
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::RaisingCosine::tag () const 
{ 
  static const std::string s_name = "RasisingCosine" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_s ) ; 
}
// ============================================================================


// ============================================================================
/*  constructor from all parameters 
 *  @param mu  peak location
 *  @param lambdaL ``left''  exponential slope  (lambdaL>0)
 *  @param lambdaR ``right'' exponential slope  (lambdaR>0)
 */
// ============================================================================
Ostap::Math::AsymmetricLaplace::AsymmetricLaplace 
( const double mu      , // location 
  const double lambdaL , // left  exponential slope 
  const double lambdaR ) // right exponential slope 
  : m_mu ( mu ) 
  , m_lambdaL ( std::abs ( lambdaL ) ) 
  , m_lambdaR ( std::abs ( lambdaR ) ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::AsymmetricLaplace::~AsymmetricLaplace(){}
// ============================================================================
bool Ostap::Math::AsymmetricLaplace::setMu    ( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::AsymmetricLaplace::setLambdaL ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_lambdaL ) ) { return false ; }
  m_lambdaL = v;
  return true ;
}
// ============================================================================
bool Ostap::Math::AsymmetricLaplace::setLambdaR ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_lambdaR ) ) { return false ; }
  m_lambdaR = v;
  return true ;
}
// ============================================================================
// evaluate asymmetic laplace function
// ============================================================================
double  Ostap::Math::AsymmetricLaplace::pdf ( const double x ) const 
{
  const long double l = 1.0L / ( m_lambdaL + m_lambdaR ) ;
  return 
    x < m_mu ? 
    l * std::exp (   ( x - m_mu ) / m_lambdaL ) :
    l * std::exp ( - ( x - m_mu ) / m_lambdaR ) ;
}
// ============================================================================
// evaluate CDF
// ============================================================================
double  Ostap::Math::AsymmetricLaplace::cdf ( const double x ) const 
{
  const long double l = 1.0L / ( m_lambdaL + m_lambdaR ) ;
  return 
    x < m_mu ? 
    m_lambdaR     * l * std::exp (   ( x - m_mu ) / m_lambdaL ) :
    1 - m_lambdaL * l * std::exp ( - ( x - m_mu ) / m_lambdaR ) ;
}
// ============================================================================
// get integral from low to high
// ============================================================================
double Ostap::Math::AsymmetricLaplace::integral
( const double low  ,
  const double high ) const 
{ return s_equal ( low ,  high ) ? 0.0 : cdf ( high ) - cdf ( low ) ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::AsymmetricLaplace::tag () const 
{ 
  static const std::string s_name = "AsymmetricLaplace" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_lambdaL , m_lambdaR ) ; 
}
// ============================================================================


// ============================================================================
/*  constructor from all arguments            
 *  @param mean  the mean/mode/location of the peak 
 *  @param q     q-value   (q<3, for q>3, it will be reflected)
 *  @param scale 
 */
// ============================================================================
Ostap::Math::QGaussian::QGaussian 
( const double mean  ,  // mean/mode/location 
  const double scale ,  // scale/sigma
  const double q     )  // q-parameter (shape)  
  : m_mean  ( mean               )
  , m_scale ( std::abs ( scale ) )
  , m_q     ( 1                  )
  , m_cq    ( s_SQRTPI           ) 
{
  setQ ( q ) ; 
}
// ============================================================================
// evaluate  pdf  for q-Gaussian distribution
// ============================================================================
// evaluate PDF of q-Gaussian distribution 
// ============================================================================
double Ostap::Math::QGaussian::pdf ( const  double x ) const 
{
  //
  if ( 1 == m_q || s_equal ( m_q , 1 ) ) { return gauss_pdf ( x , m_mean , m_scale ) ; }
  //
  const double dx  = ( x - m_mean ) / m_scale ;
  //
  static const double s_sq2 = std::sqrt ( 2.0 ) ;
  //
  return Ostap::Math::tsallis_qexp ( -0.5 * dx * dx , m_q ) / 
    ( s_sq2 * m_scale * m_cq ) ;  
}
// ============================================================================
// set mean 
// ============================================================================
bool Ostap::Math::QGaussian::setMean ( const double value ) 
{
  if ( s_equal ( value , m_mean ) ) { return false ; }
  m_mean = value ;
  return true ;
}
// ============================================================================
// set q
// ============================================================================
bool Ostap::Math::QGaussian::setQ ( const double value ) 
{
  if ( value > 3 ) { return setQ ( 6  - value ) ; } // ATTENTION! 
  //
  if ( s_equal ( value , m_q ) ) { return false ; }
  //
  m_q  = value ;
  //
  m_cq = s_SQRTPI ; 
  //
  if      ( 1 > m_q ) 
  {
    const long double q  = m_q ;
    const long double g1 = std::lgamma ( 1.0             / ( 1.0L - q ) ) ;
    const long double g2 = std::lgamma ( 0.5 * ( 3 - q ) / ( 1.0L - q ) ) ;
    m_cq *= 2 * std::exp ( g1  - std::log ( 3.0L - q ) - 
                           0.5 * std::log ( 1.0L - q ) - g2 ) ;
  }
  else if ( 1 < m_q )
  {
    const long double q  = m_q ;
    const long double g1 = std::lgamma ( 1.0                / ( q - 1.0L ) ) ;
    const long double g2 = std::lgamma ( 0.5 * ( 3.0L - q ) / ( q - 1.0L ) ) ;
    m_cq *=     std::exp ( g2 - 0.5 * std::log ( q - 1.0L ) - g1 ) ;
    //
  }
  //
  return true ;
}
// ============================================================================
// set scale
// ============================================================================
bool Ostap::Math::QGaussian::setScale ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_scale ) ) { return false ; }
  m_scale = v ;
  return true ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::QGaussian::integral
( const double low  , 
  const double high ) const 
{
  ///
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if (           low > high   ) {  return -integral ( high , low ) ; }
  ///
  if ( 1 == m_q || s_equal ( m_q , 1 ) ) 
  {
    return 
      gauss_cdf  ( high , m_mean , m_scale ) - 
      gauss_cdf  (  low , m_mean , m_scale ) ;
  }
  //
  if ( m_q  > 1 ) 
  {
    if ( low < m_mean && m_mean < high ) 
    {
      const double dx1 = m_mean -  low    ;
      const double dx2 = high   -  m_mean ;
      return 
        dx1 < dx2 ? 
        2 * integral ( low    , m_mean ) + integral ( 2 * m_mean - low ,              high ) :
        2 * integral ( m_mean , high   ) + integral (              low , 2 * m_mean - high ) ;
    }
    else if ( high - low > 3 * m_scale ) 
    {
      const double mid = 0.5 * ( low  + high ) ;
      return integral ( low , mid ) + integral ( mid , high ) ;
    } 
    //
  }
  
  double xlow  = low  ;
  double xhigh = high ;
  
  if ( m_q < 1 ) 
  {
    static const double s_sq2 = std::sqrt ( 2.0 ) ;
    const double win  = s_sq2 * m_scale / std::sqrt ( 1.0L - m_q ) ;
    const double xmin = m_mean - win ;
    const double xmax = m_mean + win ;
    if ( high <= xmin || low >= xmax ) { return 0 ; } // RETURN
    xlow  = std::max ( xmin , xlow  ) ;
    xhigh = std::min ( xmax , xhigh ) ;
  }
  //  are we already in the tail? 
  const bool in_tail = 
    std::min ( std::abs ( xhigh - m_mean ) , std::abs ( m_mean - xlow ) )  > 8 * m_scale ;   
  //
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<QGaussian> s_integrator {} ;
  static char s_message[] = "Integral(QGaussian)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () , 
      &F     , 
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()   ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::QGaussian::tag () const 
{ 
  static const std::string s_name = "QGaussian" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mean , m_q , m_scale ) ; 
}
// ============================================================================





// ============================================================================
/*  constructor from all arguments            
 *  @param mean  the mean/mode/location of the peak 
 *  @param kappa  kappa-value 
 *  @param scale 
 */
// ============================================================================
Ostap::Math::KGaussian::KGaussian 
( const double mean  ,  // mean/mode/location 
  const double scale ,  // scale/sigma
  const double kappa )  // kappa-parameter (shape)  
  : m_mean  ( mean               )
  , m_scale ( std::abs ( scale ) )
  , m_k     ( 100   )
  , m_kappa ( kappa )
{
  setKappa ( kappa ) ; 
}
// ============================================================================
// evaluate  pdf  for k-Gaussian distribution
// ============================================================================
double Ostap::Math::KGaussian::pdf ( const  double x ) const 
{
  //
  if ( 0 == m_k || s_zero ( m_k ) ) { return gauss_pdf ( x , m_mean , m_scale ) ; }
  //
  const double dx  = ( x - m_mean ) / m_scale ;
  //
  return m_Zk / m_scale * Ostap::Math::kaniadakis_kexp ( -0.5 * dx * dx , m_k ) ;
}
// ============================================================================
// set mean 
// ============================================================================
bool Ostap::Math::KGaussian::setMean ( const double value ) 
{
  if ( s_equal ( value , m_mean ) ) { return false ; }
  m_mean = value ;
  return true ;
}
// ============================================================================
// set kappa 
// ============================================================================
bool Ostap::Math::KGaussian::setKappa ( const double value ) 
{
  //
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_kappa ) && std::abs ( m_k ) <= 1 ) { return false ; }
  //
  m_kappa = avalue                ;
  m_k     = std::tanh ( m_kappa ) ;
  //
  if ( s_zero ( m_k ) ) { m_Zk = s_SQRT2PIi ; }
  else 
  {
    m_Zk = std::sqrt ( m_k / M_PI ) 
      * ( 1 + 0.5 * m_k ) 
      * std::exp ( std::lgamma ( 0.5/m_k + 0.25 ) - std::lgamma ( 0.5/m_k - 0.25 ) ) ;  
  }
  //
  return true ;
}
// ============================================================================
// set scale
// ============================================================================
bool Ostap::Math::KGaussian::setScale ( const double value ) 
{
  const double v = std::abs ( value ) ;
  if ( s_equal ( v , m_scale ) ) { return false ; }
  m_scale = v ;
  return true ;
}
// ============================================================================
// get variance 
// ============================================================================
double Ostap::Math::KGaussian::variance   () const 
{
  if ( 0 == m_k || s_zero ( m_k ) ) { return m_scale * m_scale ; }
  
  const double f1 = std::exp ( std::lgamma ( 0.5/m_k + 0.25 ) - 
                               std::lgamma ( 0.5/m_k - 0.25 ) ) ;
  const double f2 = 4 * m_k * ( 2 + m_k ) / ( ( 2 - m_k ) * ( 4 - 9 * m_k * m_k ) ) ;
  return 2 * m_scale * m_scale * f2 * f1 * f1 ;  
}
// ============================================================================
// get RMS
// ============================================================================
double Ostap::Math::KGaussian::rms () const 
{ return std::sqrt ( variance () ) ; }  
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::KGaussian::integral
( const double low  , 
  const double high ) const 
{
  ///
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if (           low > high   ) {  return -integral ( high , low ) ; }
  ///
  if ( 0 == m_k || s_zero ( m_k ) ) 
  {
    return 
      gauss_cdf  ( high , m_mean , m_scale ) - 
      gauss_cdf  (  low , m_mean , m_scale ) ;
  }
  //
  // split into some "reasonable" intervals 
  if ( low <  m_mean  && m_mean < high ) 
  { return integral ( low , m_mean ) + integral ( m_mean , high ) ; }
  //
  { 
    const double x1 = m_mean + 3 * m_scale ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_mean - 3 * m_scale ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }    
  }
  //
  { 
    const double x1 = m_mean + 5 * m_scale ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_mean - 5 * m_scale ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }    
  }
  //
  { 
    const double x1 = m_mean + 10 * m_scale ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_mean - 10 * m_scale ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }    
  }
  //
  { 
    const double x1 = m_mean + 15 * m_scale ;
    if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
    const double x2 = m_mean - 15 * m_scale ;
    if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
    
  }
  //
  const double x_low  = m_mean - 15 * m_scale  ;
  const double x_high = m_mean + 15 * m_scale  ;
  //
  const bool in_tail = high <= x_low || x_high <= low ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<KGaussian> s_integrator {} ;
  static char s_message[] = "Integral(KGaussian)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () , 
      &F     , 
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()   ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::KGaussian::tag () const 
{ 
  static const std::string s_name = "KGaussian" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mean , m_kappa , m_scale ) ; 
}
// ============================================================================



namespace 
{
  // ==========================================================================
  const double z_SMALL = 1.e-6 ;
  // ==========================================================================
  inline double _knu_ ( const double z , const double nu ) 
  {
    const double zh  = 0.5 * z ;
    const double zh2 = zh  * zh ;
    const double gn  = std::tgamma ( nu ) ;
    return gn * std::pow ( zh , -nu ) * ( 1 + zh2 /(1-nu) + 0.5 * zh2 * zh2 / ( ( 1 - nu ) * ( 2 - nu ) ) ) ;
  }
  // ==========================================================================
  ///  evaluate \f$ K_{\nu}(z) \f$ for small values of z 
  inline double knu  ( const double z , const double nu ) 
  { 
    return 
      z < z_SMALL && !s_zero ( nu ) ?
          0.5 * ( _knu_( z , nu ) + _knu_  ( z , -nu ) ) 
          : Ostap::Math::bessel_Knu ( nu , z ) ;
  }
  // ===========================================================================
  /** evaluate \f$ z^{\nu} K^{*}_{\nu}(z)}\f$ for small values of z,
   *  where \f$ K^*_{\nu}(z)\f$ is a scaled modified Bessel function 
   */
  inline double z_knu_scaled ( const double z , 
                               const double nu ) 
  {
    //
    if ( s_zero ( z ) ) 
    { return ( nu <= 0 ) ? 0.0 : std::pow ( 2 , nu - 1 ) * std::tgamma ( nu  ) ; }
    //
    if      ( z > z_SMALL ) { return std::pow ( z , nu ) * Ostap::Math::bessel_Knu_scaled ( nu , z ) ; }
    //
    if      ( nu >  0.2 ) { return 0.5 * std::pow ( 2       ,  nu ) * std::tgamma (  nu ) ; }
    else if ( nu < -0.2 ) { return 0.5 * std::pow ( 2/(z*z) , -nu ) * std::tgamma ( -nu ) ; }
    //
    if ( s_zero ( nu ) ) { return z * ( -M_EULER - std::log ( 0.5 * z ) ) ; }
    //
    const double zh  = 0.5 * z ;
    const double zh2 = zh  * zh ;
    const double gn1 = std::tgamma (  nu ) ;
    const double gn2 = std::tgamma ( -nu ) ;
    //
    // const double g1  = gn1 * std::pow ( zh , -nu ) * ( 1 + zh2 /(1-nu) + 0.5 * zh2 * zh2 / ( ( 1 - nu ) * ( 2 - nu ) ) ) ;
    // const double g2  = gn2 * std::pow ( zh ,  nu ) * ( 1 + zh2 /(1+nu) + 0.5 * zh2 * zh2 / ( ( 1 + nu ) * ( 2 + nu ) ) ) ;
    //
    const double g1  = gn1 * ( 1 + zh2 /(1-nu) + 0.5 * zh2 * zh2 / ( ( 1 - nu ) * ( 2 - nu ) ) ) ;
    const double g2  = gn2 * ( 1 + zh2 /(1+nu) + 0.5 * zh2 * zh2 / ( ( 1 + nu ) * ( 2 + nu ) ) ) ;
    ///
    return 0 <= nu ? 
      0.5 * ( g1 + std::pow ( zh,  2 * nu ) * g2 ) :
      0.5 * ( g2 + std::pow ( zh, -2 * nu ) * g1 ) ;
  }
  // ==========================================================================
  /// calculate \f$ z K_{nu+1}(z)/K_{nu}(z) \f$ 
  inline double _AL2_ ( const  double nu , const double z ) 
  {
    // 
    // return z * 
    //   Ostap::Math::bessel_Knu_scaled ( nu + 1 , z ) / 
    //   Ostap::Math::bessel_Knu_scaled ( nu     , z ) ;
    //
    if ( z_SMALL  <=  z ) 
    {
      return z * 
        Ostap::Math::bessel_Knu_scaled ( nu + 1 , z ) / 
        Ostap::Math::bessel_Knu_scaled ( nu     , z ) ;
    }
    //
    if ( s_equal ( nu , -1 ) )
    {
      const double zh = 0.5 * z ;
      const double zlog = std::log ( zh ) ;
      return z * z * ( -M_EULER - zlog ) / ( 1 + z * zh * zlog ) ;
    }
    //
    else if ( s_equal ( nu , 0 ) )
    {
      const double zh = 0.5 * z ;
      const double zlog = std::log ( zh ) ;
      return ( 1 + zh * zh * ( 1 + 2 * zlog ) ) / ( -M_EULER + ( 1 - M_EULER ) * zh * zh - zlog ) ;
    }
    //
    else if ( nu < -1.15 ) 
    { 
      return 0.5 * z * z / abs ( nu ) ;   
    }
    else if ( nu < -1 ) 
    {
      return z * 
        Ostap::Math::bessel_Knu_scaled ( nu + 1 , z ) / 
        Ostap::Math::bessel_Knu_scaled ( nu     , z ) ;
    }
    else if ( nu < -0.2 ) 
    {
      const double d  = std::abs  ( nu ) ;
      const double xh = 0.5 * z ;
      return 2 * std::tgamma ( 1 - d ) / std::tgamma ( d ) * std::pow ( xh , 2*d ) ;
    }
    else if ( std::abs ( nu ) <= 0.2 ) 
    {
      return z * knu ( z , nu + 1 ) / knu ( z , nu ) ;
    }
    //
    return 2 * nu ;
  }
  // ==========================================================================
}
// ============================================================================



// ============================================================================
/* constructor from mu, sigma, zeta and kappa 
 *  @param mu    related to location 
 *  @param beta  related to asymmetry
 *  @param sigma related to width 
 *  @param zeta  related to what ?
 */
// ============================================================================
Ostap::Math::Hyperbolic::Hyperbolic
( const double mu     ,   // related to location 
  const double sigma  ,   // related to withs  
  const double zeta   ,   // shape parameter
  const double kappa  )   // related to asymmetry 
  : m_mu    ( mu    ) 
  , m_sigma ( -1    )
  , m_zeta  ( -1    )
  , m_kappa ( kappa )
  , m_AL    ( -1    )  
  , m_N     ( -1    )
{
  setSigma ( sigma ) ;
  setZeta  ( zeta  ) ;
}
// ============================================================================
bool Ostap::Math::Hyperbolic::setMu    ( const double value ) 
{
  if ( s_equal ( value , m_mu  ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Hyperbolic::setSigma ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_sigma ) ) { return false ; }
  m_sigma = avalue ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Hyperbolic::setZeta ( const double value ) 
{
  //
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_zeta ) &&  ( 0 < m_AL ) && ( 0 < m_N ) ) { return false ; }
  m_zeta = avalue ;
  //
  m_AL = std::sqrt ( _AL2_ ( 1 , m_zeta ) ) ;
  m_N  = 1 / ( z_knu_scaled ( m_zeta , 1 ) ) ; 
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::Hyperbolic::setKappa ( const double value ) 
{
  if ( s_equal ( value , m_kappa ) ) { return false ; }
  m_kappa = value ;
  return true ;
}
// ============================================================================
/*  set "standard" parameters 
 *  @param mu     mu-parameter, location 
 *  @param beta   beta-parameter, asymmetry 
 *  @param gamma  alpha-parameter, shape 
 *  @param delta  delta-parameter, scale 
 *
 *  \f$ \alpha = \sqrt{\beta^2+\gamma^2} \f$ 
 */ 
// ============================================================================
bool Ostap::Math::Hyperbolic::setStandard
( const double mu     ,
  const double beta   , 
  const double gamma  , 
  const double delta  )
{
  bool modified = !s_equal ( m_mu , mu )  ;
  //
  m_mu     = mu    ;
  //
  const double _zeta   = std::abs ( delta ) * std::abs ( gamma ) ;
  if ( !s_equal ( m_zeta , _zeta ) ) { modified = true ; }
  m_zeta = _zeta ;
  //
  if ( modified ) { m_AL = std::sqrt ( _AL2_ ( 1  , m_zeta ) ) ; }
  //
  const double _sigma = m_AL / std::abs ( gamma ) ;
  if ( s_equal ( m_sigma , _sigma ) ) { modified = true ; }
  m_sigma = _sigma ;
  //
  if ( modified ) { m_N = 1 / ( s_SQRT2PI * z_knu_scaled ( m_zeta , 1 ) ) ; }
  //
  const double _kappa = beta / m_sigma ;
  if ( s_equal ( m_kappa , _kappa ) ) { modified = true ; }
  m_kappa = _kappa ;
  //
  return modified ;
}
// ============================================================================
// calculate the mean of the distribution  
// ============================================================================
double Ostap::Math::Hyperbolic::mean () const 
{ return m_mu + m_kappa * m_sigma ; }
// ============================================================================
// get the actual mode of the distribution
// ============================================================================
double Ostap::Math::Hyperbolic::mode () const 
{ return m_mu + m_kappa * m_sigma * m_zeta / ( m_AL * m_AL ) ; }
// ============================================================================
// get the variance/dispersion 
// ============================================================================
double Ostap::Math::Hyperbolic::variance () const 
{
  //
  const double s2 = sigma2 () ;
  const double k2 = kappa2 () ;
  const double z2 = zeta2  () ;
  //
  return s2 + k2 * s2 * ( _AL2_ ( 1 + 1 , m_zeta ) /  ( m_AL * m_AL ) - 1.0 ) ;
}
// ============================================================================
// evaluate  pdf  for the Hyperbolic distribution
// ============================================================================
double Ostap::Math::Hyperbolic::pdf ( const double x ) const 
{
  //
  const double dx =  ( x - m_mu ) / m_sigma ;
  //
  const double a2 = m_AL * m_AL            ;
  const double ka = m_kappa * m_kappa + a2 ;
  //
  const double q  = 
    - std::sqrt ( ka * ( m_zeta * m_zeta / a2 + dx * dx ) )  
    + m_kappa * dx    
    + m_zeta          ;  // from normalization 
  //
  const double aa = 0.5 * a2 / ( m_sigma * std::sqrt ( ka ) ) ;
  //
  return m_N * std::exp ( q ) * aa ;
}
// ============================================================================
// get the integral between low and high limits
// =========================================================================
double Ostap::Math::Hyperbolic::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0        ; } // RETURN
  else if (           low > high   ) { return - integral ( high , low  ) ; } // RETURN
  //
  const double m1    = mode () ;
  const double m2    = mean () ;
  const double mmin  = std::min ( m1 , m2 ) ;
  const double mmax  = std::max ( m1 , m2 ) ;
  const double mlow  = mmin - 5 * m_sigma ;
  const double mhigh = mmax + 5 * m_sigma ;
  //
  const double mc [] = { mmin - 3.0 * m_sigma , 
                         mmax + 3.0 * m_sigma , 
                         mlow , mhigh         } ;
  //
  for ( const double c : mc ) 
  { if ( low < c  && c < high ) { return integral ( low , c ) + integral ( c , high ) ; } }
  //
  // in tails 
  const bool in_tail = ( high <= mlow ) || ( low >= mhigh ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Hyperbolic> s_integrator {} ;
  static char s_message[] = "Integral(Hyperbolic)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () , 
      &F     ,  
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size () ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Hyperbolic::tag () const 
{ 
  static const std::string s_name = "Hyperbolic" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_sigma , m_zeta , m_kappa  ) ; 
}
// ============================================================================




// ============================================================================
/*  constructor from mu, sigma, zeta and kappa 
 *  @param mu     related to location 
 *  @param sigma  related to width 
 *  @param zeta   related to kurtosis 
 *  @param kappa  related to asymmetry
 *  @param lambda shape-related    
 */
// ============================================================================
Ostap::Math::GenHyperbolic::GenHyperbolic
( const double mu     ,  // related to location 
  const double sigma  ,  // related to width 
  const double zeta   ,  // related to kurtosis
  const double kappa  ,  // related to asymmetry 
  const double lambda ) 
  : m_mu     ( mu    ) 
  , m_sigma  ( std::abs ( sigma ) )
  , m_zeta   ( zeta  ) 
  , m_kappa  ( kappa ) 
  , m_lambda ( lambda ) 
  , m_AL     ( -1 ) 
  , m_N      ( -1 ) 
{
  setSigma  ( sigma  ) ;  
  setLambda ( lambda ) ;  
  setZeta   ( zeta   ) ;  
}
// ============================================================================
bool Ostap::Math::GenHyperbolic::setMu    ( const double value ) 
{
  if ( s_equal ( value , m_mu  ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenHyperbolic::setSigma ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_sigma ) ) { return false ; }
  m_sigma = avalue ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenHyperbolic::setKappa ( const double value ) 
{
  if ( s_equal ( value , m_kappa ) ) { return false ; }
  m_kappa = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::GenHyperbolic::setZeta ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_zeta ) &&  ( 0 < m_AL ) && ( 0 < m_N ) ) { return false ; }
  m_zeta = avalue ;
  //
  m_AL = std::sqrt ( _AL2_ ( m_lambda , m_zeta ) ) ;
  m_N  = 1 / ( s_SQRT2PI * z_knu_scaled ( m_zeta , m_lambda ) ) ;
  //
  return true ;
}
// ============================================================================
bool Ostap::Math::GenHyperbolic::setLambda ( const double value ) 
{
  if ( s_equal ( value , m_lambda ) && ( 0 < m_AL ) && ( 0 < m_N ) ) { return false ; }
  m_lambda = value ;
  //
  m_AL = std::sqrt ( _AL2_ ( m_lambda , m_zeta ) ) ;
  m_N  = 1 / ( s_SQRT2PI * z_knu_scaled ( m_zeta , m_lambda ) ) ;
  //
  return true ;
}
// ============================================================================
/*  set "standard" parameters 
 *  @param mu     mu-parameter, location 
 *  @param beta   beta-parameter, asymmetry 
 *  @param gamma  alpha-parameter, shape 
 *  @param delta  delta-parameter, scale 
 *  @param lambda lambda-parameter, shape   
 *
 *  \f$ \alpha = \sqrt{\beta^2+\gamma^2} \f$ 
 *  \f[ \begin{array}{ll} 
 *    \delta \ge 0, left| \beta \right|<  \alpha &  if~\lambda > 0 \\ 
 *    \delta >   0, left| \beta \right|<  \alpha &  if~\lambda = 0 \\ 
 *    \delta >   0, left| \beta \right|\le\alpha &  if~\lambda < 0 
 *  \end{array}\f] 
 */ 
// ============================================================================
bool Ostap::Math::GenHyperbolic::setStandard
( const double mu     ,
  const double beta   , 
  const double gamma  , 
  const double delta  ,
  const double lambda ) 
{
  bool modified = !s_equal ( m_mu , mu ) || !s_equal ( m_lambda , lambda ) ;
  //
  m_mu     = mu     ;
  m_lambda = lambda ;
  //
  const double _zeta   = std::abs ( delta ) * std::abs ( gamma ) ;
  if ( !s_equal ( m_zeta , _zeta ) ) { modified = true ; }
  m_zeta = _zeta ;
  //
  if ( modified ) { m_AL = std::sqrt ( _AL2_ ( m_lambda , m_zeta ) ) ; }
  //
  const double _sigma = m_AL / std::abs ( gamma ) ;
  if ( s_equal ( m_sigma , _sigma ) ) { modified = true ; }
  m_sigma = _sigma ;
  //
  if ( modified ) { m_N = 1 / ( s_SQRT2PI * z_knu_scaled ( m_zeta , m_lambda ) ) ; }
  //
  const double _kappa = beta / m_sigma ;
  if ( s_equal ( m_kappa , _kappa ) ) { modified = true ; }
  m_kappa = _kappa ;
  //
  return modified ;
}
// ============================================================================
// evaluate  pdf  for Generalised Hyperbolic distribution
// ============================================================================
double Ostap::Math::GenHyperbolic::pdf ( const  double x ) const 
{  
  const double dx   =  ( x - m_mu ) / m_sigma ;
  //
  const double k2   = m_kappa * m_kappa ;
  const double k2pA = k2 + m_AL * m_AL  ;
  const double z_A  = m_zeta    / m_AL  ;
  //
  const double arg2 =  k2pA * ( dx * dx + z_A * z_A ) ;
  const double arg  =  std::sqrt ( arg2 ) ;
  //
  // NB: we use scaled bessel function here!
  const double kfun = Ostap::Math::bessel_Knu_scaled ( m_lambda - 0.5 , arg ) ;
  //
  const double f   = 
    + std::log ( kfun )             // scaled bessel function 
    - arg                           // "unscale" it 
    + m_zeta                        // from normalzation 
    + m_kappa * dx                  // asymmetry factor  
    + ( m_lambda - 0.5 ) * std::log ( arg * m_sigma * m_sigma / k2pA ) ;
  //
  return m_N * std::exp ( f ) * std::pow ( gamma2() , m_lambda ) ;
}
// ============================================================================
// get the integral between low and high limits
// =========================================================================
double Ostap::Math::GenHyperbolic::integral
( const double low  ,
  const double high ) const
{
  //
  if      ( s_equal ( low , high ) ) { return                 0.0        ; } // RETURN
  else if (           low > high   ) { return - integral ( high , low  ) ; } // RETURN
  //
  const double m1    = mean () ;
  const double mlow  = m1 - 5 * m_sigma ;
  const double mhigh = m1 + 5 * m_sigma ;
  //
  const double mc [] = { m1 - 3.0 * m_sigma , 
                         m1 + 3.0 * m_sigma , 
                         mlow , mhigh       } ;
  //
  for ( const double c : mc ) 
  { if ( low < c  && c < high ) { return integral ( low , c ) + integral ( c , high ) ; } }
  //
  // in tails 
  const bool in_tail = ( high <= mlow ) || ( low >= mhigh ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<GenHyperbolic> s_integrator {} ;
  static char s_message[] = "Integral(GenHyperbolic)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () , 
      &F     ,  
      low    , high  ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size () ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
}
// ==============================================================================
// get mean value 
// ==============================================================================
double Ostap::Math::GenHyperbolic::mean        () const 
{ return m_mu + m_kappa * m_sigma ; }
// ==============================================================================
// get variance 
// ==============================================================================
double Ostap::Math::GenHyperbolic::variance    () const 
{
  //
  const double s2 = sigma2 () ;
  const double k2 = kappa2 () ;
  const double z2 = zeta2  () ;
  //
  return s2 + k2 * s2 * ( _AL2_ ( m_lambda + 1 , m_zeta ) /  ( m_AL * m_AL ) - 1.0 ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::GenHyperbolic::tag () const 
{ 
  static const std::string s_name = "GHD" ;
  return Ostap::Utils::hash_combiner ( s_name   , 
                             m_mu     , 
                             m_sigma  , 
                             m_kappa  , 
                             m_zeta   , 
                             m_lambda ) ; 
}
// ============================================================================







// ============================================================================
/*  constructor with full parameters 
 *  @param mu peak location 
 *  @param sigma sigma for Gaussian Core 
 *  @param kL    left tail parameter 
 *  @param kR    right tail parameter 
 */
// ============================================================================
Ostap::Math::Das::Das 
( const double mu     , // location parameter 
  const double sigma  , // width parameter 
  const double kL     , // left tails 
  const double kR     ) // right tail 
  : m_mu    ( mu ) 
  , m_sigma ( std::abs ( sigma ) ) 
  , m_kL    ( std::abs ( kL    ) ) 
  , m_kR    ( std::abs ( kR    ) ) 
{}
// ============================================================================
bool Ostap::Math::Das::setMu    ( const double value ) 
{
  if ( s_equal ( value , m_mu  ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Das::setSigma ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_sigma ) ) { return false ; }
  m_sigma = avalue ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Das::setKL ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_kL ) ) { return false ; }
  m_kL = avalue ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Das::setKR ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_kR ) ) { return false ; }
  m_kR = avalue ;
  return true ;
}
// ============================================================================
// evaluate  pdf 
// ============================================================================
double Ostap::Math::Das::pdf ( const  double x ) const 
{
  const double dx = ( x - m_mu ) / m_sigma ;
  ///
  static const double s_N = 1.0 / std::sqrt ( 2.0 * M_PI ) ;
  //
  return 
    ( dx  <= - m_kL ) ? 
    s_N * std::exp ( m_kL * ( 0.5 * m_kL + dx ) ) / m_sigma : 
    ( dx  >=   m_kR ) ? 
    s_N * std::exp ( m_kR * ( 0.5 * m_kR - dx ) ) / m_sigma : 
    s_N * std::exp ( -0.5 * dx * dx             ) / m_sigma ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::Das::integral () const 
{
  //
  static const double s_N = 1.0 / std::sqrt ( 2.0 * M_PI ) ;
  //
  return 
    // gaussian core 
    Ostap::Math::gauss_int ( -m_kL , m_kR ) 
    // left  tail
    + s_N * std::exp ( -0.5 * m_kL * m_kL ) / m_kL 
    // right tail
    + s_N * std::exp ( -0.5 * m_kR * m_kR ) / m_kR ;  
}
// =============================================================================
// get the integral 
// =============================================================================
double Ostap::Math::Das::integral 
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if (           low > high   ) { return - integral ( high , low ) ; }
  //
  const double sL = m_mu - m_kL * m_sigma ;
  if ( low < sL && sL < high ) { return integral ( low , sL ) + integral ( sL , high ) ; }
  const double sR = m_mu + m_kR * m_sigma ;
  if ( low < sR && sR < high ) { return integral ( low , sR ) + integral ( sR , high ) ; }
  //
  static const double s_N = 1.0 / std::sqrt ( 2.0 * M_PI ) ;
  //
  // left tail 
  if ( std::max ( low , high ) <= sL )
  {
    const double k2h = 0.5* m_kL * m_kL ;
    const double kS  = m_kL / m_sigma   ;
    return s_N * ( std::exp ( k2h + ( high - m_mu ) * kS ) - 
                   std::exp ( k2h + ( low  - m_mu ) * kS ) ) / m_kL ;
  }
  /// right tail
  if ( std::min ( low , high ) >= sR )
  {
    const double k2h = 0.5* m_kR * m_kR ;
    const double kS  = m_kR / m_sigma   ;
    return s_N * ( std::exp ( k2h - ( low  - m_mu ) * kS ) - 
                   std::exp ( k2h - ( high - m_mu ) * kS ) ) / m_kR ;
  }
  // gaussian core 
  return Ostap::Math::gauss_int ( low ,  high , m_mu , m_sigma );
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Das::tag () const 
{ 
  static const std::string s_name = "Das" ;
  return Ostap::Utils::hash_combiner ( s_name   , 
                             m_mu     , 
                             m_sigma  , 
                             m_kL     , 
                             m_kR     ) ;
}
// ============================================================================





// // ============================================================================
// /*  constructor with full parameters 
//  *  @param mu    related to location 
//  *  @param sigma related to RSM/scale/width 
//  *  @param xi    related to asymmetry/skewness
//  *  @param r     shape parameter 
//  *  @param alpha shape parameter     
//  */
// // ============================================================================
// Ostap::Math::SkewGenT::SkewGenT  
// ( const double mu     ,    // location parameter 
//   const double sigma  ,    // width parameter 
//   const double xi     ,    // asymmetry/skewness parameter 
//   const double r      ,    // shape parameter 
//   const double alpha  )    // shape parameter 
//   : m_mu    ( mu    ) 
//   , m_sigma ( sigma ) 
//   , m_xi    ( xi ) 
//   , m_r     ( std::abs ( r     ) ) 
//   , m_alpha ( std::abs ( alpha ) )   
// {
//   setMu    ( mu    ) ;
//   setSigma ( sigma ) ;
//   setXi    ( xi    ) ;
//   setR     ( r     ) ;
//   setAlpha ( alpha ) ;
// }
// // ======================================================================
// // set mu-parameter
// // ======================================================================
// bool Ostap::Math::SkewGenT::setMu
// ( const double value ) 
// {
//   if ( s_equal ( value , m_mu ) ) { return false ; }
//   m_mu = value ;
//   return true ;
// }
// // ======================================================================
// // set sigma-parameter
// // ======================================================================
// bool Ostap::Math::SkewGenT::setSigma
// ( const double value ) 
// {
//   const double avalue = std::abs ( value ) ;
//   if ( s_equal ( avalue , m_sigma ) ) { return false ; }
//   m_sigma = avalue ;
//   return true ;
// }
// // ======================================================================
// // set xi-parameter
// // ======================================================================
// bool Ostap::Math::SkewGenT::setXi
// ( const double value ) 
// {
//   if ( s_equal ( value , m_xi ) && -1<= m_lambda && m_lambda <= 1 ) { return false ; }
//   m_xi     = value ;
//   m_lambda = std::tanh ( value ) ; 
//   return true ;
// }
// // ======================================================================
// // set r-parameter
// // ======================================================================
// bool Ostap::Math::SkewGenT::setR
// ( const double value ) 
// {
//   const double avalue = std::abs ( value ) ;
//   if ( s_equal ( avalue , m_r ) 
//        && 0    < m_p 
//        && 0    < m_norm 
//        && 0    < m_qip 
//        && -100 != m_b2  
//        && -100 != m_b3 ) { return false ; }
//   //
//   m_r = avalue  ;
//   m_p = 1.0/m_r ;
//   //
//   const double q_= ( m_alpha + 2 ) * m_r ;
//   const double p_= m_p ;
//   //
//   m_qip          = std::pow  ( q_ ,  m_r ) ;
//   const double lnb1  = Ostap::Math::lnbeta ( m_r , q_ ) ;
//   const double lnqp  = m_r * std::log ( q_ ) ;
//   m_norm = std::exp ( - lnqp - lnb1 ) ;
//   m_b2   = std::exp ( Ostap::Math::lnbeta ( 2 * m_r , q_ -     m_r ) - lnb1 ) ;
//   m_b3   = std::exp ( Ostap::Math::lnbeta ( 3 * m_r , q_ - 2 * m_r ) - lnb1 ) ;
//   //
//   m_mc   = m_qip * m_b2 ;
//   //
//   return true ;
// }
// // ======================================================================
// // set alpha-parameter
// // ======================================================================
// bool Ostap::Math::SkewGenT::setAlpha
// ( const double value ) 
// {
//   const double avalue = std::abs ( value ) ;
//   if ( s_equal ( avalue , m_alpha ) 
//        && 0    < m_q 
//        && 0    < m_norm 
//        && 0    < m_qip
//        && -100 != m_b2  
//        && -100 != m_b3   ) { return false ; }
//   //
//   m_alpha = avalue  ;
//   m_q     = ( m_alpha + 2 ) * m_r ;
//   //
//   m_qip          = std::pow  ( m_q , m_r ) ;
//   //
//   const double lnb1  = Ostap::Math::lnbeta ( m_r , m_q ) ;
//   const double lnqp  = m_r * std::log ( m_q ) ;
//   //
//   // m_norm = std::exp ( -lbqp - lnb1 ) ;
//   m_norm = std::exp ( - lnb1 ) ;
//   m_b2   = std::exp ( Ostap::Math::lnbeta ( 2 * m_r , m_q -     m_r ) - lnb1 ) ;
//   m_b3   = std::exp ( Ostap::Math::lnbeta ( 3 * m_r , m_q - 2 * m_r ) - lnb1 ) ;
//   //
//   return true ;
// }
// // ============================================================================
// // helper bias parameter for location 
// // ============================================================================
// double Ostap:: Math::SkewGenT::m ( const double vv ) const
// { return s_zero ( m_lambda ) ? 0 : 2 * vv * m_sigma * m_lambda * m_qip * m_b2  ; }
// // ============================================================================
// // helper bias parameter for location 
// // ============================================================================
// double Ostap:: Math::SkewGenT::m () const
// { return s_zero ( m_lambda ) ? 0 : m ( v () ) ; }
// // ============================================================================
// // helper scale parameter 
// // ============================================================================
// double Ostap:: Math::SkewGenT::v () const 
// {  
//   return 1 / std::sqrt ( ( 3 * m_lambda * m_lambda + 1 )  * m_b3 
//                          - 4 * m_lambda * m_lambda * m_b2 * m_b2 )  ;
// }
// // ============================================================================
// /// evaluate the pdf 
// // ============================================================================
// double Ostap::Math::SkewGenT::pdf ( const double x ) const 
// { 
//   const double v_ = v (    ) ;
//   const double m_ = m ( v_ ) ; 
//   //
//   const double dx = ( x - m_mu + m_ ) / ( v_ * m_sigma ) ;
//   const double t  = dx / ( m_lambda * std::copysign ( 1.0 , dx ) + 1 ) ;
//   //
//   return 
//     m_norm * m_p * std::pow ( std::pow ( std::abs ( t ) , m_p ) + 1 , -m_r - m_q ) / ( 2 * v_ * m_sigma ) ;
// }
// // ============================================================================
// // integral 
// // ============================================================================
// double Ostap::Math::SkewGenT::integral () const { return 1 ; }
// // ============================================================================
// // integral
// // ============================================================================
// double Ostap::Math::SkewGenT::integral
// ( const double low  , 
//   const double high ) const 
// {
//   //
//   if      ( s_equal ( low , high ) ) { return 0 ; }
//   else if ( high < low             ) { return - integral ( high , low ) ; }
//   //
//   // split into reasonable sub intervals
//   //
//   const double mm = m_mu - m ()  ;
//   //
//   if ( low <  mm && mm < high ) { return integral ( low , mm ) + integral ( mm , high ) ; }
//   //
//   {
//     const double x1 = mm + 3 * m_sigma ;
//     if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
//     const double x2 = mm - 3 * m_sigma ;
//     if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
//   }
//   //
//   {
//     const double x1 = mm  + 5 * m_sigma ;
//     if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
//     const double x2 = mm  - 5 * m_sigma ;
//     if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
//   }  
//   //
//   {
//     const double x1 = mm + 10 * m_sigma ;
//     if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
//     const double x2 = mm - 10 * m_sigma ;
//     if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
//   }
//   //
//   {
//     const double x1 = mm + 15 * m_sigma ;
//     if ( low < x1 && x1 < high ) { return integral ( low , x1 ) + integral ( x1 , high ) ; }
//     const double x2 = mm -  15 * m_sigma ;
//     if ( low < x2 && x2 < high ) { return integral ( low , x2 ) + integral ( x2 , high ) ; }
//   }
//   //
//   const double x1     = mm - 15 * m_sigma  ;
//   const double x2     = mm + 15 * m_sigma  ;
//   const double x_low  = std::min ( x1 , x2 ) ;
//   const double x_high = std::max ( x1 , x2 ) ;
//   //
//   // use GSL to evaluate the integral
//   //
//   static const Ostap::Math::GSL::Integrator1D<SkewGenT> s_integrator {} ;
//   static char s_message[] = "Integral(SkewGenT)" ;
//   //
//   const bool in_tail = high <= x_low || x_high <= low ;
//   //
//   const auto F = s_integrator.make_function ( this ) ;
//   int    ierror   =  0 ;
//   double result   =  1 ;
//   double error    = -1 ;
//   std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
//     ( tag () , 
//       &F     , 
//       low    , high             ,    // low & high edges
//       workspace ( m_workspace ) ,    // workspace
//       in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
//       in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
//       m_workspace.size()              ,           // size of workspace
//       s_message           , 
//       __FILE__ , __LINE__ ) ;
//   //
//   return result ;
// }
// // ============================================================================
// // get the tag 
// // ============================================================================
// std::size_t Ostap::Math::SkewGenT::tag () const 
// { 
//   static const std::string s_name = "SkewGenT" ;
//   return Ostap::Utils::hash_combiner ( s_name   , 
//                              m_mu     , 
//                              m_sigma  , 
//                              m_xi     , 
//                              m_r      , 
//                              m_alpha  ) ;
// }
// // ============================================================================





// ===========================================================================
// constructor with location and scale parmaeters 
// ===========================================================================
Ostap::Math::Hat::Hat
( const double mu       , 
  const double varsigma ) 
  : m_mu       ( mu )
  , m_varsigma ( std::abs ( varsigma ) ) 
{}
// ============================================================================
// set mu
// ============================================================================
bool Ostap::Math::Hat::setMu
( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
// set varsigma
// ============================================================================
bool Ostap::Math::Hat::setVarsigma
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_varsigma ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::Hat::evaluate 
( const double x ) const
{
  static const double s_norm = 1/0.443993816168079313833061405603 ;
  const double z = ( x - m_mu ) / m_varsigma ;
  return 1 <= std::abs ( z ) ? 0.0 : Ostap::Math::hat ( z )  * s_norm / m_varsigma ;
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::Hat::integral () const { return 1 ; }
// ============================================================================
// integral
// ============================================================================
double Ostap::Math::Hat::integral
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  //
  const double mn = ( low  - m_mu ) / m_varsigma ;
  const double mx = ( high - m_mu ) / m_varsigma ;
  //
  if      ( mx <= -1            ) { return 0 ; }
  else if ( mn >=  1            ) { return 0 ; }
  else if ( mn <= -1 && mx >= 1 ) { return 1 ; }
  //
  const double xmn = std::max ( low  , m_mu - m_varsigma ) ;
  const double xmx = std::min ( high , m_mu + m_varsigma ) ;
  //
   // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Hat> s_integrator {} ;
  static char s_message[] = "Integral(Hat)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () , 
      &F     ,  
      xmn    , xmx              , // low & high edges
      workspace ( m_workspace ) , // workspace
      s_APRECISION              , // absolute precision
      s_RPRECISION              , // relative precision
      m_workspace.size ()       , // size of workspace
      s_message                 , 
      __FILE__ , __LINE__       ) ;
  //
  return result ;
}
// ============================================================================
// get the variance of the distribution 
// ============================================================================
double Ostap::Math::Hat::variance () const 
{ return m_varsigma * m_varsigma * 0.15811363626379668 ; }
// ============================================================================
// get the RMS of the distribution 
// ============================================================================
double Ostap::Math::Hat::rms      () const 
{ return m_varsigma * 0.3976350541184676 ; }
// ============================================================================
// get the (excess) kurtosis 
// ============================================================================
double Ostap::Math::Hat::kurtosis () const 
{ return -0.8807206646393597 ; }
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Hat::tag () const 
{ 
  static const std::string s_name = "Hat" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_varsigma ) ;
}
// ============================================================================


// ===========================================================================
// constructor with location and scale parmaeters 
// ===========================================================================
Ostap::Math::Up::Up
( const double mu       , 
  const double varsigma ) 
  : m_mu       ( mu )
  , m_varsigma ( std::abs ( varsigma ) ) 
{}
// ============================================================================
// set mu
// ============================================================================
bool Ostap::Math::Up::setMu
( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
// set varsigma
// ============================================================================
bool Ostap::Math::Up::setVarsigma
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_varsigma ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::Up::evaluate 
( const double x ) const
{
  const double z = ( x - m_mu ) / m_varsigma ;
  return 1 <= std::abs ( z ) ? 0.0 : eval ( z ) / m_varsigma ;
}
// ===========================================================================
// evaluate the "standard" up function 
// ===========================================================================
double Ostap::Math::Up::eval ( const double z )  const 
{
  //
  auto fourrier = []( unsigned int k ) -> double  
    { return k == 0 || 1 == k % 2 ?  Ostap::Math::up_F ( M_PI * k ) : 0.0 ; } ;
  static const std::array<double,120> s_fourrier = 
    detail::make_array( fourrier , std::make_index_sequence<120>() ) ;
  //
  return 1 <= std::abs ( z ) ? 0.0 : Ostap::Math::Clenshaw::cosine_sum 
    ( s_fourrier.begin () , s_fourrier.end () , z * M_PI ) ; 
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::Up::integral () const { return 1 ; }
// ============================================================================
// integral
// ============================================================================
double Ostap::Math::Up::integral
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  //
  const double mn = ( low  - m_mu ) / m_varsigma ;
  const double mx = ( high - m_mu ) / m_varsigma ;
  //
  if      ( mx <= -1            ) { return 0 ; }
  else if ( mn >=  1            ) { return 0 ; }
  else if ( mn <= -1 && mx >= 1 ) { return 1 ; }
  //
  const double xmn = std::max ( low  , xmin () ) ;
  const double xmx = std::min ( high , xmax () ) ;
  //
   // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Up> s_integrator {} ;
  static char s_message[] = "Integral(Up)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () , 
      &F     ,  
      xmn    , xmx              , // low & high edges
      workspace ( m_workspace ) , // workspace
      s_APRECISION              , // absolute precision
      s_RPRECISION              , // relative precision
      m_workspace.size ()       , // size of workspace
      s_message                 , 
      __FILE__ , __LINE__       ) ;
  //
  return result ;
}
// ============================================================================
// get the variance of the distribution 
// ============================================================================
double Ostap::Math::Up::variance () const 
{ return m_varsigma * m_varsigma / 9 ;}
// ============================================================================
// get the RMS of the distribution 
// ============================================================================
double Ostap::Math::Up::rms      () const 
{ return m_varsigma / 3 ; }
// ============================================================================
// get the (excess) kurtosis 
// ============================================================================
double Ostap::Math::Up::kurtosis () const 
{
  static const double s_kurtosis = 19.0 * 9 * 9 / ( std::pow ( 3 , 3 ) * 5 * 5 ) - 3 ;
  return s_kurtosis ;
}
// ============================================================================
// get the value of the derivative 
// ============================================================================
double Ostap::Math::Up::derivative ( const double x ) const 
{
  const double z = ( x - m_mu ) / m_varsigma ;
  return 1 <= std::abs ( z )  ? 0.0 : 
    2 * ( eval ( 2 * z + 1 )  - eval ( 2 * z - 1 ) ) / m_varsigma ;
}  
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Up::tag () const 
{ 
  static const std::string s_name = "Up" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_varsigma ) ;
}
// ============================================================================


// ============================================================================
namespace 
{  
  // ==========================================================================
  typedef std::array<double,120>           RESULT ;
  typedef std::map<unsigned short,RESULT>  MAP    ;
  typedef SyncedCache<MAP>                 CACHE  ;
  // =========================================================================
  CACHE s_FupN_cache {} ;
  // =========================================================================
}
// ===========================================================================
// constructor with location and scale parmaeters 
// ===========================================================================
Ostap::Math::FupN::FupN
( const unsigned short N        , 
  const double         mu       , 
  const double         varsigma ) 
  : m_N        ( N  )
  , m_mu       ( mu )
  , m_varsigma ( std::abs ( varsigma ) ) 
{
  CACHE::Lock lock { s_FupN_cache.mutex() } ;
  auto it = s_FupN_cache->find  ( m_N  ) ;
  if ( s_FupN_cache->end() == it ) 
  {
    auto fourrier = [this]( unsigned int k ) -> double
      { return Ostap::Math::fupN_F ( this->m_N , M_PI * k / ( this->m_N + 1 ) ) ; } ;
    //
    const RESULT res = detail::make_array( fourrier , std::make_index_sequence<120>() ) ;
    s_FupN_cache->insert ( std::make_pair ( m_N , res ) ) ;  
  }
}
// ============================================================================
// set mu
// ============================================================================
bool Ostap::Math::FupN::setMu
( const double value ) 
{
  if ( s_equal ( value , m_mu ) ) { return false ; }
  m_mu = value ;
  return true ;
}
// ============================================================================
// set varsigma
// ============================================================================
bool Ostap::Math::FupN::setVarsigma
( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_varsigma ) ) { return false ; }
  m_varsigma = avalue ;
  return true ;
}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::FupN::evaluate 
( const double x ) const
{
  const double z = ( x - m_mu ) / m_varsigma ;
  return 0.5 * ( m_N + 2 )  <= std::abs ( z ) ? 0.0 : eval ( z ) / m_varsigma ;
}
// ===========================================================================
// evaluate the "standard" fupN function 
// ===========================================================================
double Ostap::Math::FupN::eval ( const double z )  const 
{
  //
  auto it = s_FupN_cache->find  ( m_N  ) ;
  Ostap::Assert ( s_FupN_cache->end() != it , 
                  "Cache does not exist!"   , 
                  "Ostap::Math::FupN"       ) ;
  //
  return 0.5 * ( m_N + 2 ) <= std::abs ( z ) ? 0.0 : Ostap::Math::Clenshaw::cosine_sum 
    ( it->second.begin () , 
      it->second.end   () ,  M_PI * z / ( m_N + 1 ) ) / ( m_N + 1 ) ; 
}
// ============================================================================
// integral 
// ============================================================================
double Ostap::Math::FupN::integral () const { return 1 ; }
// ============================================================================
// integral
// ============================================================================
double Ostap::Math::FupN::integral
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  //
  const double mn = ( low  - m_mu ) / m_varsigma ;
  const double mx = ( high - m_mu ) / m_varsigma ;
  //
  const double nn = 0.5 * ( m_N + 2 ) ;
  //
  if      ( mx <= -1               ) { return 0 ; }
  else if ( mn >=  1               ) { return 0 ; }
  else if ( mn <= -nn  && mx >= nn ) { return 1 ; }
  //
  const double xmn = std::max ( low  , xmin () ) ;
  const double xmx = std::min ( high , xmax () ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<FupN> s_integrator {} ;
  static char s_message[] = "Integral(FupN)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   =  0 ;
  double result   =  1 ;
  double error    = -1 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () , 
      &F     ,  
      xmn    , xmx              , // low & high edges
      workspace ( m_workspace ) , // workspace
      s_APRECISION              , // absolute precision
      s_RPRECISION              , // relative precision
      m_workspace.size ()       , // size of workspace
      s_message                 , 
      __FILE__ , __LINE__       ) ;
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::FupN::tag () const 
{ 
  static const std::string s_name = "FupN" ;
  return Ostap::Utils::hash_combiner ( s_name , m_N , m_mu , m_varsigma ) ;
}
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
